
/*!
 ***********************************************************************
 * \file
 *    configfile.c
 * \brief
 *    Configuration handling.
 * \author
 *  Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Stephan Wenger           <stewe@cs.tu-berlin.de>
 * \note
 *    In the future this module should hide the Parameters and offer only
 *    Functions for their access.  Modules which make frequent use of some parameters
 *    (e.g. picture size in macroblocks) are free to buffer them on local variables.
 *    This will not only avoid global variable and make the code more readable, but also
 *    speed it up.  It will also greatly facilitate future enhancements such as the
 *    handling of different picture sizes in the same sequence.                         \n
 *                                                                                      \n
 *    For now, everything is just copied to the inp_par structure (gulp)
 *
 **************************************************************************************
 * \par Configuration File Format
 **************************************************************************************
 * Format is line oriented, maximum of one parameter per line                           \n
 *                                                                                      \n
 * Lines have the following format:                                                     \n
 * \<ParameterName\> = \<ParameterValue\> # Comments \\n                                    \n
 * Whitespace is space and \\t
 * \par
 * \<ParameterName\> are the predefined names for Parameters and are case sensitive.
 *   See configfile.h for the definition of those names and their mapping to
 *   cfgparams->values.
 * \par
 * \<ParameterValue\> are either integers [0..9]* or strings.
 *   Integers must fit into the wordlengths, signed values are generally assumed.
 *   Strings containing no whitespace characters can be used directly.  Strings containing
 *   whitespace characters are to be inclosed in double quotes ("string with whitespace")
 *   The double quote character is forbidden (may want to implement something smarter here).
 * \par
 * Any Parameters whose ParameterName is undefined lead to the termination of the program
 * with an error message.
 *
 * \par Known bug/Shortcoming:
 *    zero-length strings (i.e. to signal an non-existing file
 *    have to be coded as "".
 *
 * \par Rules for using command files
 *                                                                                      \n
 * All Parameters are initially taken from DEFAULTCONFIGFILENAME, defined in configfile.h.
 * If an -f \<config\> parameter is present in the command line then this file is used to
 * update the defaults of DEFAULTCONFIGFILENAME.  There can be more than one -f parameters
 * present.  If -p <ParameterName = ParameterValue> parameters are present then these
 * override the default and the additional config file's settings, and are themselves
 * overridden by future -p parameters.  There must be whitespace between -f and -p commands
 * and their respective parameters
 ***********************************************************************
 */

#define INCLUDED_BY_CONFIGFILE_C

#include <sys/stat.h>

#include "global.h"
#include "configfile.h"
#include "fmo.h"
#include "conformance.h"
#include "mc_prediction.h"
#include "mv_search.h"
#include "img_io.h"
#include "ratectl.h"

char *GetConfigFileContent (char *Filename);
static void ParseContent            (InputParameters *p_Inp, Mapping *Map, char *buf, int bufsize);
static void PatchInp                (VideoParameters *p_Vid, InputParameters *p_Inp);
static int  ParameterNameToMapIndex (Mapping *Map, char *s);
static int  InitEncoderParams       (Mapping *Map);
static int  TestEncoderParams       (Mapping *Map, int bitdepth_qp_scale[3]);
static int  DisplayEncoderParams    (Mapping *Map);

static const int mb_width_cr[4] = {0,8, 8,16};
static const int mb_height_cr[4]= {0,8,16,16};

#define MAX_ITEMS_TO_PARSE  10000


/*!
 ***********************************************************************
 * \brief
 *   print help message and exit
 ***********************************************************************
 */
void JMHelpExit (void)
{
  fprintf( stderr, "\n   lencod [-h] [-d defenc.cfg] {[-f curenc1.cfg]...[-f curencN.cfg]}"
    " {[-p EncParam1=EncValue1]..[-p EncParamM=EncValueM]}\n\n"
    "## Parameters\n\n"

    "## Options\n"
    "   -h :  prints function usage\n"
    "   -d :  use <defenc.cfg> as default file for parameter initializations.\n"
    "         If not used then file defaults to encoder.cfg in local directory.\n"
    "   -f :  read <curencM.cfg> for reseting selected encoder parameters.\n"
    "         Multiple files could be used that set different parameters\n"
    "   -p :  Set parameter <EncParamM> to <EncValueM>.\n"
    "         See default encoder.cfg file for description of all parameters.\n\n"

    "## Supported video file formats\n"
    "   RAW:  .yuv -> YUV 4:2:0\n\n"

    "## Examples of usage:\n"
    "   lencod\n"
    "   lencod  -h\n"
    "   lencod  -d default.cfg\n"
    "   lencod  -f curenc1.cfg\n"
    "   lencod  -f curenc1.cfg -p InputFile=\"e:\\data\\container_qcif_30.yuv\" -p SourceWidth=176 -p SourceHeight=144\n"
    "   lencod  -f curenc1.cfg -p FramesToBeEncoded=30 -p QPISlice=28 -p QPPSlice=28 -p QPBSlice=30\n");

  exit(-1);
}

/*!
 ************************************************************************
 * \brief
 *    Reads Input File Size 
 *
 ************************************************************************
 */
int64 getVideoFileSize(int video_file)
{
   int64 fsize;   

   lseek(video_file, 0, SEEK_END); 
   fsize = tell((int) video_file); 
   lseek(video_file, 0, SEEK_SET); 

   return fsize;
}

/*!
 ************************************************************************
 * \brief
 *    Updates the number of frames to encode based on the file size
 *
 ************************************************************************
 */
void getNumberOfFrames (InputParameters *p_Inp, VideoDataFile *input_file)
{
  int64 fsize = getVideoFileSize(input_file->f_num);
  int64 isize = (int64) p_Inp->source.size;
  int maxBitDepth = imax(p_Inp->source.bit_depth[0], p_Inp->source.bit_depth[1]);

  isize <<= (maxBitDepth > 8)? 1: 0;
  p_Inp->no_frames   = (int) (((fsize - p_Inp->infile_header)/ isize) - p_Inp->start_frame);
}

/*!
 ************************************************************************
 * \brief
 *    Updates images max values
 *
 ************************************************************************
 */
static void updateMaxValue(FrameFormat *format)
{
  format->max_value[0] = (1 << format->bit_depth[0]) - 1;
  format->max_value_sq[0] = format->max_value[0] * format->max_value[0];
  format->max_value[1] = (1 << format->bit_depth[1]) - 1;
  format->max_value_sq[1] = format->max_value[1] * format->max_value[1];
  format->max_value[2] = (1 << format->bit_depth[2]) - 1;
  format->max_value_sq[2] = format->max_value[2] * format->max_value[2];
}

/*!
 ************************************************************************
 * \brief
 *    Update output format parameters (resolution & bit-depth) given input
 *
 ************************************************************************
 */
static void updateOutFormat(InputParameters *p_Inp)
{
  FrameFormat *output = &p_Inp->output;
  FrameFormat *source = &p_Inp->source;
  output->yuv_format  = (ColorFormat) p_Inp->yuv_format;
  source->yuv_format  = (ColorFormat) p_Inp->yuv_format;

  if (p_Inp->src_resize == 0)
  {
    output->width  = source->width;
    output->height = source->height;
  }

  if (p_Inp->yuv_format == YUV400) // reset bitdepth of chroma for 400 content
  {
    source->bit_depth[1] = 8;
    output->bit_depth[1] = 8;
    source->width_cr  = 0;
    source->height_cr = 0;
    output->width_cr  = 0;
    output->height_cr = 0;
  }
  else
  {
    source->width_cr  = (source->width  * mb_width_cr [output->yuv_format]) >> 4;
    source->height_cr = (source->height * mb_height_cr[output->yuv_format]) >> 4;
    output->width_cr  = (output->width  * mb_width_cr [output->yuv_format]) >> 4;
    output->height_cr = (output->height * mb_height_cr[output->yuv_format]) >> 4;
  }

  // source size
  source->size_cmp[0] = source->width * source->height;
  source->size_cmp[1] = source->width_cr * source->height_cr;
  source->size_cmp[2] = source->size_cmp[1];
  source->size        = source->size_cmp[0] + source->size_cmp[1] + source->size_cmp[2];
  source->mb_width    = source->width  / MB_BLOCK_SIZE;
  source->mb_height   = source->height / MB_BLOCK_SIZE;
  source->pic_unit_size_on_disk = (imax(source->bit_depth[0], source->bit_depth[1]) > 8) ? 16 : 8;
  source->pic_unit_size_shift3 = source->pic_unit_size_on_disk >> 3;


  // output size (excluding padding)
  output->size_cmp[0] = output->width * output->height;
  output->size_cmp[1] = output->width_cr * output->height_cr;
  output->size_cmp[2] = output->size_cmp[1];
  output->size        = output->size_cmp[0] + output->size_cmp[1] + output->size_cmp[2];
  output->mb_width    = output->width  / MB_BLOCK_SIZE;
  output->mb_height   = output->height / MB_BLOCK_SIZE;


  // both chroma components have the same bitdepth
  source->bit_depth[2] = source->bit_depth[1];
  output->bit_depth[2] = output->bit_depth[1];
  
  // if no bitdepth rescale ensure bitdepth is same
  if (p_Inp->src_BitDepthRescale == 0) 
  {    
    output->bit_depth[0] = source->bit_depth[0];
    output->bit_depth[1] = source->bit_depth[1];
    output->bit_depth[2] = source->bit_depth[2];
  }
  output->pic_unit_size_on_disk = (imax(output->bit_depth[0], output->bit_depth[1]) > 8) ? 16 : 8;
  output->pic_unit_size_shift3 = output->pic_unit_size_on_disk >> 3;
  
  output->frame_rate = source->frame_rate / (p_Inp->frame_skip + 1);
  output->color_model = source->color_model;

  updateMaxValue(source);
  updateMaxValue(output);
}


/*!
 ***********************************************************************
 * \brief
 *    Parse the command line parameters and read the config files.
 * \param p_Vid
 *    VideoParameters structure for encoding
 * \param p_Inp
 *    InputParameters structure as input configuration
 * \param ac
 *    number of command line parameters
 * \param av
 *    command line parameters
 ***********************************************************************
 */
void Configure (VideoParameters *p_Vid, InputParameters *p_Inp, int ac, char *av[])
{
  char *content = NULL;
  int CLcount, ContentLen, NumberParams;
  char *filename=DEFAULTCONFIGFILENAME;

  if (ac==2)
  {
    if (0 == strncmp (av[1], "-v", 2))
    {
      printf("JM-" VERSION "\n");
      exit(0);
    }
    if (0 == strncmp (av[1], "-V", 2))
    {
      printf("JM " JM ": compiled " __DATE__ " " __TIME__ "\n");
#if ( IMGTYPE == 0 )
      printf("support for more than 8 bits/pel disabled\n");
#endif
#if ( ENABLE_FIELD_CTX == 0 )
      printf("CABAC field coding disabled\n");
#endif
#if ( ENABLE_HIGH444_CTX == 0 )
      printf("CABAC High 4:4:4 profile coding disabled\n");
#endif
      exit(0);
    }

    if (0 == strncmp (av[1], "-h", 2))
    {
      JMHelpExit();
    }
  }

  memset (&cfgparams, 0, sizeof (InputParameters));
  //Set default parameters.
  printf ("Setting Default Parameters...\n");
  InitEncoderParams(Map);

  // Process default config file
  CLcount = 1;

  if (ac>=3)
  {
    if (0 == strncmp (av[1], "-d", 2))
    {
      filename=av[2];
      CLcount = 3;
    }
    if (0 == strncmp (av[1], "-h", 2))
    {
      JMHelpExit();
    }
  }
  printf ("Parsing Configfile %s", filename);
  content = GetConfigFileContent (filename);
  if (NULL==content)
    error (errortext, 300);
  ParseContent (p_Inp, Map, content, strlen(content));
  printf ("\n");
  free (content);

  // Parse the command line

  while (CLcount < ac)
  {
    if (0 == strncmp (av[CLcount], "-h", 2))
    {
      JMHelpExit();
    }

    if (0 == strncmp (av[CLcount], "-f", 2) || 0 == strncmp (av[CLcount], "-F", 2))  // A file parameter?
    {
      content = GetConfigFileContent (av[CLcount+1]);
      if (NULL==content)
        error (errortext, 300);
      printf ("Parsing Configfile %s", av[CLcount+1]);
      ParseContent (p_Inp, Map, content, strlen (content));
      printf ("\n");
      free (content);
      CLcount += 2;
    } 
    else
    {
      if (0 == strncmp (av[CLcount], "-p", 2) || 0 == strncmp (av[CLcount], "-P", 2))  // A config change?
      {
        // Collect all data until next parameter (starting with -<x> (x is any character)),
        // put it into content, and parse content.

        ++CLcount;
        ContentLen = 0;
        NumberParams = CLcount;

        // determine the necessary size for content
        while (NumberParams < ac && av[NumberParams][0] != '-')
          ContentLen += strlen (av[NumberParams++]);        // Space for all the strings
        ContentLen += 1000;                     // Additional 1000 bytes for spaces and \0s


        if ((content = malloc (ContentLen))==NULL) no_mem_exit("Configure: content");;
        content[0] = '\0';

        // concatenate all parameters identified before

        while (CLcount < NumberParams)
        {
          char *source = &av[CLcount][0];
          char *destin = &content[strlen (content)];

          while (*source != '\0')
          {
            if (*source == '=')  // The Parser expects whitespace before and after '='
            {
              *destin++=' '; *destin++='='; *destin++=' ';  // Hence make sure we add it
            } 
            else
              *destin++=*source;
            source++;
          }
          *destin = '\0';
          CLcount++;
        }
        printf ("Parsing command line string '%s'", content);
        ParseContent (p_Inp, Map, content, strlen(content));
        free (content);
        printf ("\n");
      }
      else
      {
        snprintf (errortext, ET_SIZE, "Error in command line, ac %d, around string '%s', missing -f or -p parameters?", CLcount, av[CLcount]);
        error (errortext, 300);
      }
    }
  }
  printf ("\n");
  PatchInp(p_Vid, p_Inp);

  memcpy (&cfgparams, p_Inp, sizeof (InputParameters));

  if (p_Inp->DisplayEncParams)
    DisplayEncoderParams(Map);
}

/*!
 ***********************************************************************
 * \brief
 *    allocates memory buf, opens file Filename in f, reads contents into
 *    buf and returns buf
 * \param Filename
 *    name of config file
 * \return
 *    if successfull, content of config file
 *    NULL in case of error. Error message will be set in errortext
 ***********************************************************************
 */
char *GetConfigFileContent (char *Filename)
{
  long FileSize;
  FILE *f;
  char *buf;

  if (NULL == (f = fopen (Filename, "r")))
  {
      snprintf (errortext, ET_SIZE, "Cannot open configuration file %s.", Filename);
      return NULL;
  }

  if (0 != fseek (f, 0, SEEK_END))
  {
    snprintf (errortext, ET_SIZE, "Cannot fseek in configuration file %s.", Filename);
    return NULL;
  }

  FileSize = ftell (f);
  if (FileSize < 0 || FileSize > 100000)
  {
    snprintf (errortext, ET_SIZE, "Unreasonable Filesize %ld reported by ftell for configuration file %s.", FileSize, Filename);
    return NULL;
  }
  if (0 != fseek (f, 0, SEEK_SET))
  {
    snprintf (errortext, ET_SIZE, "Cannot fseek in configuration file %s.", Filename);
    return NULL;
  }

  if ((buf = malloc (FileSize + 1))==NULL) no_mem_exit("GetConfigFileContent: buf");

  // Note that ftell() gives us the file size as the file system sees it.  The actual file size,
  // as reported by fread() below will be often smaller due to CR/LF to CR conversion and/or
  // control characters after the dos EOF marker in the file.

  FileSize = fread (buf, 1, FileSize, f);
  buf[FileSize] = '\0';


  fclose (f);
  return buf;
}


/*!
 ***********************************************************************
 * \brief
 *    Parses the character array buf and writes global variable input, which is defined in
 *    configfile.h.  This hack will continue to be necessary to facilitate the addition of
 *    new parameters through the Map[] mechanism (Need compiler-generated addresses in map[]).
 * \param p_Inp
 *    InputParameters of configuration
 * \param Map
 *    Mapping structure to specify the name and value mapping relation
 * \param buf
 *    buffer to be parsed
 * \param bufsize
 *    buffer size of buffer
 ***********************************************************************
 */
void ParseContent (InputParameters *p_Inp, Mapping *Map, char *buf, int bufsize)
{
  char *items[MAX_ITEMS_TO_PARSE] = {NULL};
  int MapIdx;
  int item = 0;
  int InString = 0, InItem = 0;
  char *p = buf;
  char *bufend = &buf[bufsize];
  int IntContent;
  double DoubleContent;
  int i;

  // Stage one: Generate an argc/argv-type list in items[], without comments and whitespace.
  // This is context insensitive and could be done most easily with lex(1).

  while (p < bufend)
  {
    switch (*p)
    {
    case 13:
      ++p;
      break;
    case '#':                 // Found comment
      *p = '\0';              // Replace '#' with '\0' in case of comment immediately following integer or string
      while (*p != '\n' && p < bufend)  // Skip till EOL or EOF, whichever comes first
        ++p;
      InString = 0;
      InItem = 0;
      break;
    case '\n':
      InItem = 0;
      InString = 0;
      *p++='\0';
      break;
    case ' ':
    case '\t':              // Skip whitespace, leave state unchanged
      if (InString)
        p++;
      else
      {                     // Terminate non-strings once whitespace is found
        *p++ = '\0';
        InItem = 0;
      }
      break;

    case '"':               // Begin/End of String
      *p++ = '\0';
      if (!InString)
      {
        items[item++] = p;
        InItem = ~InItem;
      }
      else
        InItem = 0;
      InString = ~InString; // Toggle
      break;

    default:
      if (!InItem)
      {
        items[item++] = p;
        InItem = ~InItem;
      }
      p++;
    }
  }

  item--;

  for (i=0; i<item; i+= 3)
  {
    if (0 > (MapIdx = ParameterNameToMapIndex (Map, items[i])))
    {
      //snprintf (errortext, ET_SIZE, " Parsing error in config file: Parameter Name '%s' not recognized.", items[i]);
      //error (errortext, 300);
      printf ("\n\tParsing error in config file: Parameter Name '%s' not recognized.", items[i]);
      continue;
    }
    if (strcasecmp ("=", items[i+1]))
    {
      snprintf (errortext, ET_SIZE, " Parsing error in config file: '=' expected as the second token in each line.");
      error (errortext, 300);
    }

    // Now interpret the Value, context sensitive...

    switch (Map[MapIdx].Type)
    {
    case 0:           // Numerical
      if (1 != sscanf (items[i+2], "%d", &IntContent))
      {
        snprintf (errortext, ET_SIZE, " Parsing error: Expected numerical value for Parameter of %s, found '%s'.", items[i], items[i+2]);
        error (errortext, 300);
      }
      * (int *) (Map[MapIdx].Place) = IntContent;
      printf (".");
      break;
    case 1:
      if (items[i + 2] == NULL)
        memset((char *) Map[MapIdx].Place, 0, Map[MapIdx].char_size);
      else
        strncpy ((char *) Map[MapIdx].Place, items [i+2], Map[MapIdx].char_size);
      printf (".");
      break;
    case 2:           // Numerical double
      if (1 != sscanf (items[i+2], "%lf", &DoubleContent))
      {
        snprintf (errortext, ET_SIZE, " Parsing error: Expected numerical value for Parameter of %s, found '%s'.", items[i], items[i+2]);
        error (errortext, 300);
      }
      * (double *) (Map[MapIdx].Place) = DoubleContent;
      printf (".");
      break;
    default:
      error ("Unknown value type in the map definition of configfile.h",-1);
    }
  }
  memcpy (p_Inp, &cfgparams, sizeof (InputParameters));
}

/*!
 ***********************************************************************
 * \brief
 *    Returns the index number from Map[] for a given parameter name.
 * \param Map
 *    Mapping structure
 * \param s
 *    parameter name string
 * \return
 *    the index number if the string is a valid parameter name,         \n
 *    -1 for error
 ***********************************************************************
 */
static int ParameterNameToMapIndex (Mapping *Map, char *s)
{
  int i = 0;

  while (Map[i].TokenName != NULL)
    if (0==strcasecmp (Map[i].TokenName, s))
      return i;
    else
      i++;
  return -1;
}

/*!
 ***********************************************************************
 * \brief
 *    Sets initial values for encoding parameters.
 * \return
 *    -1 for error
 ***********************************************************************
 */
static int InitEncoderParams(Mapping *Map)
{
  int i = 0;

  while (Map[i].TokenName != NULL)
  {
    if (Map[i].Type == 0)
        * (int *) (Map[i].Place) = (int) Map[i].Default;
    else if (Map[i].Type == 2)
    * (double *) (Map[i].Place) = Map[i].Default;
      i++;
  }
  return -1;
}

/*!
 ***********************************************************************
 * \brief
 *    Validates encoding parameters.
 * \return
 *    -1 for error
 ***********************************************************************
 */
static int TestEncoderParams(Mapping *Map, int bitdepth_qp_scale[3])
{
  int i = 0;

  while (Map[i].TokenName != NULL)
  {
    if (Map[i].param_limits == 1)
    {
      if (Map[i].Type == 0)
      {
        if ( * (int *) (Map[i].Place) < (int) Map[i].min_limit || * (int *) (Map[i].Place) > (int) Map[i].max_limit )
        {
          snprintf(errortext, ET_SIZE, "Error in input parameter %s. Check configuration file. Value should be in [%d, %d] range.", Map[i].TokenName, (int) Map[i].min_limit,(int)Map[i].max_limit );
          error (errortext, 400);
        }

      }
      else if (Map[i].Type == 2)
      {
        if ( * (double *) (Map[i].Place) < Map[i].min_limit || * (double *) (Map[i].Place) > Map[i].max_limit )
        {
          snprintf(errortext, ET_SIZE, "Error in input parameter %s. Check configuration file. Value should be in [%.2f, %.2f] range.", Map[i].TokenName,Map[i].min_limit ,Map[i].max_limit );
          error (errortext, 400);
        }
      }
    }
    else if (Map[i].param_limits == 2)
    {
      if (Map[i].Type == 0)
      {
        if ( * (int *) (Map[i].Place) < (int) Map[i].min_limit )
        {
          snprintf(errortext, ET_SIZE, "Error in input parameter %s. Check configuration file. Value should not be smaller than %d.", Map[i].TokenName, (int) Map[i].min_limit);
          error (errortext, 400);
        }
      }
      else if (Map[i].Type == 2)
      {
        if ( * (double *) (Map[i].Place) < Map[i].min_limit )
        {
          snprintf(errortext, ET_SIZE, "Error in input parameter %s. Check configuration file. Value should not be smaller than %2.f.", Map[i].TokenName,Map[i].min_limit);
          error (errortext, 400);
        }
      }
    }
    else if (Map[i].param_limits == 3) // Only used for QPs
    {
      
      if (Map[i].Type == 0)
      {
        int cur_qp = * (int *) (Map[i].Place);
        int min_qp = (int) (Map[i].min_limit - bitdepth_qp_scale[0]);
        int max_qp = (int) Map[i].max_limit;
        
        if (( cur_qp < min_qp ) || ( cur_qp > max_qp ))
        {
          snprintf(errortext, ET_SIZE, "Error in input parameter %s. Check configuration file. Value should be in [%d, %d] range.", Map[i].TokenName, min_qp, max_qp );
          error (errortext, 400);
        }
      }
    }

    i++;
  }
  return -1;
}



/*!
 ***********************************************************************
 * \brief
 *    Outputs encoding parameters.
 * \return
 *    -1 for error
 ***********************************************************************
 */
static int DisplayEncoderParams(Mapping *Map)
{
  int i = 0;

  printf("******************************************************\n");
  printf("*               Encoder Parameters                   *\n");
  printf("******************************************************\n");
  while (Map[i].TokenName != NULL)
  {
    if (Map[i].Type == 0)
      printf("Parameter %s = %d\n",Map[i].TokenName,* (int *) (Map[i].Place));
    else if (Map[i].Type == 1)
      printf("Parameter %s = ""%s""\n",Map[i].TokenName,(char *)  (Map[i].Place));
    else if (Map[i].Type == 2)
      printf("Parameter %s = %.2f\n",Map[i].TokenName,* (double *) (Map[i].Place));
      i++;
  }
  printf("******************************************************\n");
  return -1;
}

/*!
 ************************************************************************
 * \brief
 *    read the slice group configuration file. Returns without action
 *    if type is not 0, 2 or 6
 ************************************************************************
 */
void read_slice_group_info(VideoParameters *p_Vid, InputParameters *p_Inp)
{
  FILE * sgfile=NULL;
  int i;
  int ret;
  unsigned PicSizeInMapUnits;

  if ((p_Inp->slice_group_map_type != 0) && (p_Inp->slice_group_map_type != 2) && (p_Inp->slice_group_map_type != 6))
  {
    // nothing to do
    return;
  }

  // do we have a file name (not only NULL character)
  if (strlen (p_Inp->SliceGroupConfigFileName) <= 1)
    error ("No slice group config file name specified", 500);
    
  // open file
  sgfile = fopen(p_Inp->SliceGroupConfigFileName,"r");

  if ( NULL==sgfile )
  {
    snprintf(errortext, ET_SIZE, "Error opening slice group file %s", p_Inp->SliceGroupConfigFileName);
    error (errortext, 500);
  }

  switch (p_Inp->slice_group_map_type)
  {
  case 0:
    p_Inp->run_length_minus1=(int *)malloc(sizeof(int)*(p_Inp->num_slice_groups_minus1+1));
    if ( NULL==p_Inp->run_length_minus1 )
    {
      fclose(sgfile);
      no_mem_exit("read_slice_group_info: p_Inp->run_length_minus1");
    }

    // each line contains one 'run_length_minus1' value
    for(i=0; i <= p_Inp->num_slice_groups_minus1;i++)
    {
      ret  = fscanf(sgfile,"%d",(p_Inp->run_length_minus1+i));
      if ( 1!=ret )
      {
        fclose(sgfile);
        snprintf(errortext, ET_SIZE, "Error while reading slice group config file (line %d)", i+1);
        error (errortext, 500);
      }
      // scan remaining line
      ret = fscanf(sgfile,"%*[^\n]");
    }
    break;

  case 2:
    // determine expected frame size in map units
    PicSizeInMapUnits = (p_Inp->output.width >> 4) * (p_Inp->output.height >> 4);
    if (p_Inp->MbInterlace||p_Inp->PicInterlace) 
      PicSizeInMapUnits >>= 1;

    p_Inp->top_left=(unsigned *)malloc(sizeof(unsigned)*p_Inp->num_slice_groups_minus1);
    p_Inp->bottom_right=(unsigned *)malloc(sizeof(unsigned)*p_Inp->num_slice_groups_minus1);
    
    if (NULL==p_Inp->top_left)
    {
      fclose(sgfile);
      no_mem_exit("PatchInp: p_Inp->top_left");
    }
    
    if (NULL==p_Inp->bottom_right)
    {
      fclose(sgfile);
      no_mem_exit("PatchInp: p_Inp->bottom_right");
    }

    // every two lines contain 'top_left' and 'bottom_right' value
    for(i=0;i<p_Inp->num_slice_groups_minus1;i++)
    {
      ret = fscanf(sgfile,"%ud",(p_Inp->top_left+i));
      if ( 1!=ret )
      {
        fclose(sgfile);
        snprintf(errortext, ET_SIZE, "Error while reading slice group config file (line %d)", 2*i +1);
        error (errortext, 500);
      }
      if (p_Inp->top_left[i] > PicSizeInMapUnits)
      {
        fprintf(stderr, "Warning: slice group # %d top_left exceeds picture size (will be clipped)\n", i);
      }
      // scan remaining line
      ret = fscanf(sgfile,"%*[^\n]");
      ret = fscanf(sgfile,"%ud",(p_Inp->bottom_right+i));
      if ( 1!=ret )
      {
        fclose(sgfile);
        snprintf(errortext, ET_SIZE, "Error while reading slice group config file (line %d)", 2*i + 2);
        error (errortext, 500);
      }
      if (p_Inp->bottom_right[i] > PicSizeInMapUnits)
      {
        fprintf(stderr, "Warning: slice group # %d bottom_right exceeds picture size (will be clipped)\n", i);
      }
      // scan remaining line
      ret = fscanf(sgfile,"%*[^\n]");
    }
    break;
  case 6:
    {
      int tmp;
      int frame_mb_only;
      int mb_width, mb_height, mapunit_height;

      frame_mb_only = !(p_Inp->PicInterlace || p_Inp->MbInterlace);
      mb_width  = (p_Inp->output.width + p_Vid->auto_crop_right)>>4;
      mb_height = (p_Inp->output.height + p_Vid->auto_crop_bottom)>>4;
      mapunit_height = mb_height / (2-frame_mb_only);

      p_Inp->slice_group_id=(byte * ) malloc(sizeof(byte)*mapunit_height*mb_width);
      if (NULL==p_Inp->slice_group_id)
      {
        fclose(sgfile);
        no_mem_exit("PatchInp: p_Inp->slice_group_id");
      }

      // each line contains slice_group_id for one Macroblock
      for (i=0;i<mapunit_height*mb_width;i++)
      {
        ret = fscanf(sgfile,"%d", &tmp);
        p_Inp->slice_group_id[i]= (byte) tmp;
        if ( 1!=ret )
        {
          fclose(sgfile);
          snprintf(errortext, ET_SIZE, "Error while reading slice group config file (line %d)", i + 1);
          error (errortext, 500);
        }
        if ( *(p_Inp->slice_group_id+i) > p_Inp->num_slice_groups_minus1 )
        {
          fclose(sgfile);
          snprintf(errortext, ET_SIZE, "Error while reading slice group config file: slice_group_id not allowed (line %d)", i + 1);
          error (errortext, 500);
        }
        // scan remaining line
        ret = fscanf(sgfile,"%*[^\n]");
      }
    }
    break;
  default:
    // we should not get here
    error ("Wrong slice group type while reading config file", 500);
    break;
  }

  // close file again
  fclose(sgfile);
}

/*!
 ***********************************************************************
 * \brief
 *    Checks the input parameters for consistency.
 ***********************************************************************
 */
static void PatchInp (VideoParameters *p_Vid, InputParameters *p_Inp)
{
  int i;
  int storedBplus1;
  int bitdepth_qp_scale[3];

  if (p_Inp->src_BitDepthRescale)
  {
    bitdepth_qp_scale [0] = 6*(p_Inp->output.bit_depth[0] - 8);
    bitdepth_qp_scale [1] = 6*(p_Inp->output.bit_depth[1] - 8);
    bitdepth_qp_scale [2] = 6*(p_Inp->output.bit_depth[2] - 8);
  }
  else
  {
    bitdepth_qp_scale [0] = 6*(p_Inp->source.bit_depth[0] - 8);
    bitdepth_qp_scale [1] = 6*(p_Inp->source.bit_depth[1] - 8);
    bitdepth_qp_scale [2] = 6*(p_Inp->source.bit_depth[2] - 8);
  }

  TestEncoderParams(Map, bitdepth_qp_scale);

  if (p_Inp->source.frame_rate == 0.0)
    p_Inp->source.frame_rate = (double) INIT_FRAME_RATE;

  ParseVideoType(&p_Inp->input_file1);
  ParseFrameNoFormatFromString (&p_Inp->input_file1);

  // Read resolution from file name
  if (p_Inp->source.width == 0 || p_Inp->source.height == 0)
  {
    if (ParseSizeFromString (&p_Inp->input_file1, &(p_Inp->source.width), &(p_Inp->source.height), &(p_Inp->source.frame_rate)) == 0)
    {
      snprintf(errortext, ET_SIZE, "File name does not contain resolution information.");    
      error (errortext, 500);
    }
  }

#if (!ENABLE_FIELD_CTX)
  if ( (p_Inp->PicInterlace || p_Inp->MbInterlace) && p_Inp->symbol_mode )
  {
    snprintf(errortext, ET_SIZE, "Recompile with ENABLE_FIELD_CTX set to one to enable interlaced coding with CABAC.");    
    error (errortext, 500);
  }
#endif
#if (!ENABLE_HIGH444_CTX)
  if ( p_Inp->ProfileIDC == 244 && p_Inp->symbol_mode )
  {
    snprintf(errortext, ET_SIZE, "Recompile with ENABLE_HIGH444_CTX set to one to enable the High 4:4:4 Profile with CABAC.");    
    error (errortext, 500);
  }
#endif

  // Currently to simplify things, lets copy everything (overwrites yuv_format)
  p_Inp->input_file1.format = p_Inp->source;


  if (p_Inp->idr_period && p_Inp->intra_delay && p_Inp->idr_period <= p_Inp->intra_delay)
  {
    snprintf(errortext, ET_SIZE, " IntraDelay cannot be larger than or equal to IDRPeriod.");
    error (errortext, 500);
  }

  if (p_Inp->idr_period && p_Inp->intra_delay && p_Inp->adaptive_idr_period)
  {
    snprintf(errortext, ET_SIZE, " IntraDelay can not be used with AdaptiveIDRPeriod.");
    error (errortext, 500);
  }

  // Let us set up p_Inp->jumpd from frame_skip and NumberBFrames
  p_Inp->jumpd = (p_Inp->NumberBFrames + 1) * (p_Inp->frame_skip + 1) - 1;

  updateOutFormat(p_Inp);

  if (p_Inp->no_frames == -1)
  {
    OpenFiles(&p_Inp->input_file1);
    getNumberOfFrames(p_Inp, &p_Inp->input_file1);
    CloseFiles(&p_Inp->input_file1);
  }

  storedBplus1 = (p_Inp->BRefPictures ) ? p_Inp->NumberBFrames + 1: 1;
  
  if (p_Inp->no_frames < 1)
  {      
    snprintf(errortext, ET_SIZE, "Not enough frames to encode (%d)", p_Inp->no_frames);
    error (errortext, 500);
  }

  // Direct Mode consistency check
  if(p_Inp->NumberBFrames && p_Inp->direct_spatial_mv_pred_flag != DIR_SPATIAL && p_Inp->direct_spatial_mv_pred_flag != DIR_TEMPORAL)
  {
    snprintf(errortext, ET_SIZE, "Unsupported direct mode=%d, use TEMPORAL=0 or SPATIAL=1", p_Inp->direct_spatial_mv_pred_flag);
    error (errortext, 400);
  }

  if (p_Inp->PicInterlace>0 || p_Inp->MbInterlace>0)
  {
    if (p_Inp->directInferenceFlag==0)
      printf("\nWarning: DirectInferenceFlag set to 1 due to interlace coding.");
    p_Inp->directInferenceFlag = 1;
  }

#if TRACE
  if (strlen (p_Inp->TraceFile) > 0 && (p_Enc->p_trace = fopen(p_Inp->TraceFile,"w"))==NULL)
  {
    snprintf(errortext, ET_SIZE, "Error open file %s", p_Inp->TraceFile);
    error (errortext, 500);
  }
#endif

  
  
  if ((p_Inp->slice_mode == 1)&&(p_Inp->MbInterlace != 0))
  {
    if ((p_Inp->slice_argument & 0x01)!=0)
    {
      fprintf ( stderr, "Warning: slice border within macroblock pair. ");
      if (p_Inp->slice_argument > 1)
      {
        p_Inp->slice_argument--;
      }
      else
      {
        p_Inp->slice_argument++;
      }
      fprintf ( stderr, "Using %d MBs per slice.\n", p_Inp->slice_argument);
    }
  }  

  if (p_Inp->WPMCPrecision && (p_Inp->RDPictureDecision != 1 || p_Inp->GenerateMultiplePPS != 1) )
  {
    snprintf(errortext, ET_SIZE, "WPMCPrecision requires both RDPictureDecision=1 and GenerateMultiplePPS=1.\n");
    error (errortext, 400);
  }
  if (p_Inp->WPMCPrecision && p_Inp->WPMCPrecFullRef && p_Inp->num_ref_frames < 16 )
  {
    p_Inp->num_ref_frames++;
    if ( p_Inp->P_List0_refs )
      p_Inp->P_List0_refs++;
    else
      p_Inp->P_List0_refs = p_Inp->num_ref_frames;
    if ( p_Inp->B_List0_refs )
      p_Inp->B_List0_refs++;
    else
      p_Inp->B_List0_refs = p_Inp->num_ref_frames;
    if ( p_Inp->B_List1_refs )
      p_Inp->B_List1_refs++;
    else
      p_Inp->B_List1_refs = p_Inp->num_ref_frames;
  }
  else if ( p_Inp->WPMCPrecision && p_Inp->WPMCPrecFullRef )
  {
    snprintf(errortext, ET_SIZE, "WPMCPrecFullRef requires NumberReferenceFrames < 16.\n");
    error (errortext, 400);
  }

  if (p_Inp->ReferenceReorder && p_Inp->MbInterlace )
  {
    snprintf(errortext, ET_SIZE, "ReferenceReorder not supported with MBAFF\n");
    error (errortext, 400);
  }

  if (p_Inp->SetFirstAsLongTerm && p_Inp->ReferenceReorder == 1)
  {
    printf("SetFirstAsLongTerm is set. ReferenceReorder is not supported and therefore disabled. \n");
    p_Inp->ReferenceReorder = 0;
  }

  if (p_Inp->PocMemoryManagement && p_Inp->MbInterlace )
  {
    snprintf(errortext, ET_SIZE, "PocMemoryManagement not supported with MBAFF\n");
    error (errortext, 400);
  }

  if ((!p_Inp->rdopt)&&(p_Inp->MbInterlace==2))
  {
    snprintf(errortext, ET_SIZE, "MB AFF is not compatible with non-rd-optimized coding.");
    error (errortext, 500);
  }

  // check RDoptimization mode and profile. FMD does not support Frex Profiles.
  if (p_Inp->rdopt==2 && ( p_Inp->ProfileIDC>=FREXT_HP || p_Inp->ProfileIDC==FREXT_CAVLC444 ))
  {
    snprintf(errortext, ET_SIZE, "Fast Mode Decision methods not supported in FREX Profiles");
    error (errortext, 500);
  }

  // Tian Dong: May 31, 2002
  // The number of frames in one sub-seq in enhanced layer should not exceed
  // the number of reference frame number.
  if ( p_Inp->NumFramesInELSubSeq > p_Inp->num_ref_frames || p_Inp->NumFramesInELSubSeq < 0 )
  {
    snprintf(errortext, ET_SIZE, "NumFramesInELSubSeq (%d) is out of range [0,%d).", p_Inp->NumFramesInELSubSeq, p_Inp->num_ref_frames);
    error (errortext, 500);
  }
  // Tian Dong: Enhanced GOP is not supported in bitstream mode. September, 2002
  if ( p_Inp->NumFramesInELSubSeq > 0 )
  {
    snprintf(errortext, ET_SIZE, "Enhanced GOP is not properly supported yet.");
    error (errortext, 500);
  }
  // Tian Dong (Sept 2002)
  // The AFF is not compatible with spare picture for the time being.
  if ((p_Inp->PicInterlace || p_Inp->MbInterlace) && p_Inp->SparePictureOption == TRUE)
  {
    snprintf(errortext, ET_SIZE, "AFF is not compatible with spare picture.");
    error (errortext, 500);
  }

  // Only the RTP mode is compatible with spare picture for the time being.
  if (p_Inp->of_mode != PAR_OF_RTP && p_Inp->SparePictureOption == TRUE)
  {
    snprintf(errortext, ET_SIZE, "Only RTP output mode is compatible with spare picture features.");
    error (errortext, 500);
  }

  if( (p_Inp->WeightedPrediction > 0 || p_Inp->WeightedBiprediction > 0) && (p_Inp->MbInterlace))
  {
    snprintf(errortext, ET_SIZE, "Weighted prediction coding is not supported for MB AFF currently.");
    error (errortext, 500);
  }
  if ( p_Inp->NumFramesInELSubSeq > 0 && p_Inp->WeightedPrediction > 0)
  {
    snprintf(errortext, ET_SIZE, "Enhanced GOP is not supported in weighted prediction coding mode yet.");
    error (errortext, 500);
  }

  // Rate control
  if(p_Inp->RCEnable)
  { 

    if ( p_Inp->RCUpdateMode == RC_MODE_1 && 
      !( (p_Inp->intra_period == 1 || p_Inp->idr_period == 1 || p_Inp->BRefPictures == 2 ) && !p_Inp->NumberBFrames ) )
    {
      snprintf(errortext, ET_SIZE, "Use RCUpdateMode = 1 only for all intra or all B-slice coding.");
      error (errortext, 500);
    }

    if ( p_Inp->BRefPictures == 2 && p_Inp->intra_period == 0 && p_Inp->RCUpdateMode != RC_MODE_1 )
    {
      snprintf(errortext, ET_SIZE, "Use RCUpdateMode = 1 for all B-slice coding.");
      error (errortext, 500);
    }

    if ( p_Inp->HierarchicalCoding && p_Inp->RCUpdateMode != RC_MODE_2 && p_Inp->RCUpdateMode != RC_MODE_3 )
    {
      snprintf(errortext, ET_SIZE, "Use RCUpdateMode = 2 or 3 for hierarchical B-picture coding.");
      error (errortext, 500);
    }
    if ( (p_Inp->RCUpdateMode != RC_MODE_1) && (p_Inp->intra_period == 1) )
    {
      snprintf(errortext, ET_SIZE, "Use RCUpdateMode = 1 for all intra coding.");
      error (errortext, 500);
    }
  }

  if ((p_Inp->NumberBFrames)&&(p_Inp->BRefPictures)&&(p_Inp->idr_period)&&(p_Inp->pic_order_cnt_type!=0))
  {
    error("Stored B pictures combined with IDR pictures only supported in Picture Order Count type 0\n",-1000);
  }

  if( !p_Inp->direct_spatial_mv_pred_flag && p_Inp->num_ref_frames<2 && p_Inp->NumberBFrames >0)
    error("temporal direct needs at least 2 ref frames\n",-1000);

  if (p_Inp->SearchMode == FAST_FULL_SEARCH && p_Inp->MEErrorMetric[F_PEL] > ERROR_SSE)
  {
    snprintf(errortext, ET_SIZE, "\nOnly SAD and SSE distortion computation supported with Fast Full Search.");
    error (errortext, 500);
  }

  if (p_Inp->rdopt == 0)
  {
    if (p_Inp->DisableSubpelME)
    {
      if (p_Inp->MEErrorMetric[F_PEL] != p_Inp->ModeDecisionMetric)
      {
        snprintf(errortext, ET_SIZE, "\nLast refinement level (FPel) distortion not the same as Mode decision distortion.\nPlease update MEDistortionFPel (%d) and/or  MDDistortion(%d).", p_Inp->MEErrorMetric[F_PEL], p_Inp->ModeDecisionMetric);
        error (errortext, 500);
      }
    }
    else if (p_Inp->MEErrorMetric[Q_PEL] != p_Inp->ModeDecisionMetric)
    {
      snprintf(errortext, ET_SIZE, "\nLast refinement level (QPel) distortion not the same as Mode decision distortion.\nPlease update MEDistortionQPel (%d) and/or  MDDistortion(%d).", p_Inp->MEErrorMetric[Q_PEL], p_Inp->ModeDecisionMetric);
      error (errortext, 500);
    }
  }
  // frext
  if(p_Inp->Transform8x8Mode && p_Inp->sp_periodicity /*SP-frames*/)
  {
    snprintf(errortext, ET_SIZE, "\nThe new 8x8 mode is not implemented for sp-frames.");
    error (errortext, 500);
  }

  if(p_Inp->Transform8x8Mode && ( p_Inp->ProfileIDC<FREXT_HP && p_Inp->ProfileIDC!=FREXT_CAVLC444 ))
  {
    snprintf(errortext, ET_SIZE, "\nTransform8x8Mode may be used only with ProfileIDC %d to %d.", FREXT_HP, FREXT_Hi444);
    error (errortext, 500);
  }

  if (p_Inp->DisableIntra4x4 == 1 && p_Inp->DisableIntra16x16 == 1 && p_Inp->EnableIPCM == 0 && p_Inp->Transform8x8Mode == 0)
  {
    snprintf(errortext, ET_SIZE, "\nAt least one intra prediction mode needs to be enabled.");
    error (errortext, 500);
  }

  if(p_Inp->ScalingMatrixPresentFlag && ( p_Inp->ProfileIDC<FREXT_HP && p_Inp->ProfileIDC!=FREXT_CAVLC444 ))
  {
    snprintf(errortext, ET_SIZE, "\nScalingMatrixPresentFlag may be used only with ProfileIDC %d to %d.", FREXT_HP, FREXT_Hi444);
    error (errortext, 500);
  }

  if(p_Inp->yuv_format==YUV422 && ( p_Inp->ProfileIDC < FREXT_Hi422 && p_Inp->ProfileIDC!=FREXT_CAVLC444 ))
  {
    snprintf(errortext, ET_SIZE, "\nFRExt Profile(YUV Format) Error!\nYUV422 can be used only with ProfileIDC %d or %d\n",FREXT_Hi422, FREXT_Hi444);
    error (errortext, 500);
  }
  if(p_Inp->yuv_format==YUV444 && ( p_Inp->ProfileIDC < FREXT_Hi444 && p_Inp->ProfileIDC!=FREXT_CAVLC444 ))
  {
    snprintf(errortext, ET_SIZE, "\nFRExt Profile(YUV Format) Error!\nYUV444 can be used only with ProfileIDC %d.\n",FREXT_Hi444);
    error (errortext, 500);
  }

  if (p_Inp->NumberBFrames && ((p_Inp->BiPredMotionEstimation) && (p_Inp->search_range < p_Inp->BiPredMESearchRange)))
  {
    snprintf(errortext, ET_SIZE, "\nBiPredMESearchRange must be smaller or equal SearchRange.");
    error (errortext, 500);
  }

  if (p_Inp->BiPredMotionEstimation)
  {
    p_Inp->BiPredMotionEstimation = 0;
    for (i = 0 ; i < 4; i++)
      p_Inp->BiPredMotionEstimation |= p_Inp->BiPredSearch[i];
  }
  else
  {
    for (i = 0 ; i < 4; i++)
      p_Inp->BiPredSearch[i] = 0;
  }

  // check consistency
  if ( p_Inp->ChromaMEEnable && !(p_Inp->ChromaMCBuffer) ) 
  {
    snprintf(errortext, ET_SIZE, "\nChromaMCBuffer must be set to 1 if ChromaMEEnable is set.");
    error (errortext, 500);
  }

  if ( p_Inp->ChromaMEEnable && p_Inp->yuv_format ==  YUV400) 
  {
    fprintf(stderr, "Warning: ChromaMEEnable cannot be used with monochrome color format, disabling ChromaMEEnable.\n");
    p_Inp->ChromaMEEnable = 0;
  }

  if ( (p_Inp->ChromaMCBuffer == 0) && (( p_Inp->yuv_format ==  YUV444) && (!p_Inp->separate_colour_plane_flag)) )
  {
    fprintf(stderr, "Warning: Enabling ChromaMCBuffer for 4:4:4 combined color coding.\n");
    p_Inp->ChromaMCBuffer = 1;
  }


  if (p_Inp->EnableOpenGOP)
    p_Inp->ReferenceReorder = 1;

  if (p_Inp->SearchMode != EPZS)
    p_Inp->EPZSSubPelGrid = 0;

  if (p_Inp->redundant_pic_flag)
  {
    if (p_Inp->PicInterlace || p_Inp->MbInterlace)
    {
      snprintf(errortext, ET_SIZE, "Redundant pictures cannot be used with interlaced tools.");
      error (errortext, 500);
    }
    if (p_Inp->RDPictureDecision)
    {
      snprintf(errortext, ET_SIZE, "Redundant pictures cannot be used with RDPictureDecision.");
      error (errortext, 500);
    }
    if (p_Inp->NumberBFrames)
    {
      snprintf(errortext, ET_SIZE, "Redundant pictures cannot be used with B frames.");
      error (errortext, 500);
    }
    if (p_Inp->PrimaryGOPLength < (1 << p_Inp->NumRedundantHierarchy))
    {
      snprintf(errortext, ET_SIZE, "PrimaryGOPLength must be equal or greater than 2^NumRedundantHierarchy.");
      error (errortext, 500);
    }
    if (p_Inp->num_ref_frames < p_Inp->PrimaryGOPLength)
    {
      snprintf(errortext, ET_SIZE, "NumberReferenceFrames must be greater than or equal to PrimaryGOPLength.");
      error (errortext, 500);
    }
  }

  if (p_Inp->num_ref_frames == 1 && p_Inp->NumberBFrames)
  {
    fprintf( stderr, "\nWarning: B slices used but only one reference allocated within reference buffer.\n");
    fprintf( stderr, "         Performance may be considerably compromised! \n");
    fprintf( stderr, "         2 or more references recommended for use with B slices.\n");
  }
  if ((p_Inp->HierarchicalCoding || p_Inp->BRefPictures) && p_Inp->NumberBFrames)
  {
    fprintf( stderr, "\nWarning: Hierarchical coding or Referenced B slices used.\n");
    fprintf( stderr, "         Make sure that you have allocated enough references\n");
    fprintf( stderr, "         in reference buffer to achieve best performance.\n");
  }

  if (p_Inp->FastMDEnable == 0)
  {
    p_Inp->FastIntraMD = 0;
    p_Inp->FastIntra16x16 = 0;
    p_Inp->FastIntra4x4 = 0;
    p_Inp->FastIntra8x8 = 0;
    p_Inp->FastIntraChroma = 0;
  }



  if (p_Inp->UseRDOQuant == 1)
  {
    if (p_Inp->rdopt == 0)
    {
      snprintf(errortext, ET_SIZE, "RDO Quantization not supported with low complexity RDO.");
      error (errortext, 500);
    }

    if (p_Inp->MbInterlace != 0)
    {
      printf("RDO Quantization currently not supported with MBAFF. Option disabled.\n");
      p_Inp->UseRDOQuant = 0;
      p_Inp->RDOQ_QP_Num = 1;
      p_Inp->RDOQ_CP_MV = 0;
      p_Inp->RDOQ_CP_Mode = 0;
    }
    else
    {
      p_Inp->AdaptiveRounding = 0;
      printf("AdaptiveRounding is disabled when RDO Quantization is used\n");
      if (p_Inp->RDOQ_QP_Num < 2)
      {
        p_Inp->RDOQ_CP_MV = 0;
        p_Inp->RDOQ_CP_Mode = 0;
      }
    }
  }
  else
  {
    p_Inp->RDOQ_QP_Num = 1;
    p_Inp->RDOQ_CP_MV = 0;
    p_Inp->RDOQ_CP_Mode = 0;
  }

  if(p_Inp->num_slice_groups_minus1 > 0 && (p_Inp->GenerateMultiplePPS ==1 && p_Inp->RDPictureDecision == 1))
  {
    printf("Warning: Weighted Prediction is disabled as it may not function correctly for multiple slices\n"); 
  }


  ProfileCheck(p_Inp);

}
