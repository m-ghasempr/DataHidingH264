
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
 *   configinput->values.
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

char *GetConfigFileContent (char *Filename);
static void ParseContent (char *buf, int bufsize);
static int ParameterNameToMapIndex (char *s);
static int InitEncoderParams(void);
static int TestEncoderParams(int bitdepth_qp_scale);
static int DisplayEncoderParams(void);
static void PatchInp (void);
static void ProfileCheck(void);
static void LevelCheck(void);

static int mb_width_cr[4] = {0,8, 8,16};
static int mb_height_cr[4]= {0,8,16,16};

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
int64 getVideoFileSize(void)
{
   int64 fsize;   

   lseek(p_in, 0, SEEK_END); 
   fsize = tell((int) p_in); 
   lseek(p_in, 0, SEEK_SET); 

   return fsize;
}

/*!
 ************************************************************************
 * \brief
 *    Updates the number of frames to encode based on the file size
 *
 ************************************************************************
 */
static void getNumberOfFrames (void)
{
  int64 fsize = getVideoFileSize();
  int64 isize = input->img_width * input->img_height + 2 * (input->img_width_cr * input->img_height_cr);
  if ((input->BitDepthLuma > input->BitDepthChroma) || (input->yuv_format == YUV400))
    isize <<= (input->BitDepthLuma > 8)? 1 : 0;
  else
    isize <<= (input->BitDepthChroma > 8)? 1 : 0;

  input->no_frames = (int) (((fsize - input->infile_header)/ isize) - input->start_frame - 1) / (1 + input->jumpd) + 1;
}

/*!
 ***********************************************************************
 * \brief
 *    Parse the command line parameters and read the config files.
 * \param ac
 *    number of command line parameters
 * \param av
 *    command line parameters
 ***********************************************************************
 */
void Configure (int ac, char *av[])
{
  char *content;
  int CLcount, ContentLen, NumberParams;
  char *filename=DEFAULTCONFIGFILENAME;

  memset (&configinput, 0, sizeof (InputParameters));
  //Set default parameters.
  printf ("Setting Default Parameters...\n");
  InitEncoderParams();

  // Process default config file
  CLcount = 1;

  if (ac==2)
  {
    if (0 == strncmp (av[1], "-h", 2))
    {
      JMHelpExit();
    }
  }

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
  ParseContent (content, strlen(content));
  printf ("\n");
  free (content);

  // Parse the command line

  while (CLcount < ac)
  {
    if (0 == strncmp (av[CLcount], "-h", 2))
    {
      JMHelpExit();
    }

    if (0 == strncmp (av[CLcount], "-f", 2))  // A file parameter?
    {
      content = GetConfigFileContent (av[CLcount+1]);
      if (NULL==content)
        error (errortext, 300);
      printf ("Parsing Configfile %s", av[CLcount+1]);
      ParseContent (content, strlen (content));
      printf ("\n");
      free (content);
      CLcount += 2;
    } else
    {
      if (0 == strncmp (av[CLcount], "-p", 2))  // A config change?
      {
        // Collect all data until next parameter (starting with -<x> (x is any character)),
        // put it into content, and parse content.

        CLcount++;
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
            } else
              *destin++=*source;
            source++;
          }
          *destin = '\0';
          CLcount++;
        }
        printf ("Parsing command line string '%s'", content);
        ParseContent (content, strlen(content));
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
  PatchInp();
  if (input->DisplayEncParams)
    DisplayEncoderParams();
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
  if (FileSize < 0 || FileSize > 60000)
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
 * \param buf
 *    buffer to be parsed
 * \param bufsize
 *    buffer size of buffer
 ***********************************************************************
 */
void ParseContent (char *buf, int bufsize)
{

  char *items[MAX_ITEMS_TO_PARSE];
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
        p++;
        break;
      case '#':                 // Found comment
        *p = '\0';              // Replace '#' with '\0' in case of comment immediately following integer or string
        while (*p != '\n' && p < bufend)  // Skip till EOL or EOF, whichever comes first
          p++;
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
    if (0 > (MapIdx = ParameterNameToMapIndex (items[i])))
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
        strncpy ((char *) Map[MapIdx].Place, items [i+2], FILE_NAME_SIZE);
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
  memcpy (input, &configinput, sizeof (InputParameters));
}

/*!
 ***********************************************************************
 * \brief
 *    Returns the index number from Map[] for a given parameter name.
 * \param s
 *    parameter name string
 * \return
 *    the index number if the string is a valid parameter name,         \n
 *    -1 for error
 ***********************************************************************
 */
static int ParameterNameToMapIndex (char *s)
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
static int InitEncoderParams(void)
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
static int TestEncoderParams(int bitdepth_qp_scale)
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
        if ( * (int *) (Map[i].Place) < (int) (Map[i].min_limit - bitdepth_qp_scale) || * (int *) (Map[i].Place) > (int) Map[i].max_limit )
        {
          snprintf(errortext, ET_SIZE, "Error in input parameter %s. Check configuration file. Value should be in [%d, %d] range.", Map[i].TokenName, (int) (Map[i].min_limit - bitdepth_qp_scale),(int)Map[i].max_limit );
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
static int DisplayEncoderParams(void)
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
 *    calculate Ceil(Log2(uiVal))
 ************************************************************************
 */
unsigned CeilLog2( unsigned uiVal)
{
  unsigned uiTmp = uiVal-1;
  unsigned uiRet = 0;

  while( uiTmp != 0 )
  {
    uiTmp >>= 1;
    uiRet++;
  }
  return uiRet;
}

/*!
 ************************************************************************
 * \brief
 *    read the slice group configuration file. Returns without action
 *    if type is not 0, 2 or 6
 ************************************************************************
 */
void read_slice_group_info()
{
  FILE * sgfile=NULL;
  int i;
  int ret;

  if ((input->slice_group_map_type != 0) && (input->slice_group_map_type != 2) && (input->slice_group_map_type != 6))
  {
    // nothing to do
    return;
  }

  // do we have a file name (not only NULL character)
  if (strlen (input->SliceGroupConfigFileName) <= 1)
    error ("No slice group config file name specified", 500);
    
  // open file
  sgfile = fopen(input->SliceGroupConfigFileName,"r");

  if ( NULL==sgfile )
  {
    snprintf(errortext, ET_SIZE, "Error opening slice group file %s", input->SliceGroupConfigFileName);
    error (errortext, 500);
  }

  switch (input->slice_group_map_type)
  {
  case 0:
    input->run_length_minus1=(int *)malloc(sizeof(int)*(input->num_slice_groups_minus1+1));
    if ( NULL==input->run_length_minus1 )
    {
      fclose(sgfile);
      no_mem_exit("PatchInp: input->run_length_minus1");
    }

    // each line contains one 'run_length_minus1' value
    for(i=0;i<=input->num_slice_groups_minus1;i++)
    {
      ret = fscanf(sgfile,"%d",(input->run_length_minus1+i));
      fscanf(sgfile,"%*[^\n]");
      if ( 1!=ret )
      {
        fclose(sgfile);
        snprintf(errortext, ET_SIZE, "Error while reading slice group config file (line %d)", i+1);
        error (errortext, 500);
      }
    }
    break;

  case 2:
    input->top_left=(int *)malloc(sizeof(int)*input->num_slice_groups_minus1);
    input->bottom_right=(int *)malloc(sizeof(int)*input->num_slice_groups_minus1);
    if (NULL==input->top_left)
    {
      fclose(sgfile);
      no_mem_exit("PatchInp: input->top_left");
    }
    if (NULL==input->bottom_right)
    {
      fclose(sgfile);
      no_mem_exit("PatchInp: input->bottom_right");
    }

    // every two lines contain 'top_left' and 'bottom_right' value
    for(i=0;i<input->num_slice_groups_minus1;i++)
    {
      ret = fscanf(sgfile,"%d",(input->top_left+i));
      fscanf(sgfile,"%*[^\n]");
      if ( 1!=ret )
      {
        fclose(sgfile);
        snprintf(errortext, ET_SIZE, "Error while reading slice group config file (line %d)", 2*i +1);
        error (errortext, 500);
      }
      ret = fscanf(sgfile,"%d",(input->bottom_right+i));
      fscanf(sgfile,"%*[^\n]");
      if ( 1!=ret )
      {
        fclose(sgfile);
        snprintf(errortext, ET_SIZE, "Error while reading slice group config file (line %d)", 2*i + 2);
        error (errortext, 500);
      }
    }
    break;

  case 6:
    {
      int tmp;
      int frame_mb_only;
      int mb_width, mb_height, mapunit_height;

      frame_mb_only = !(input->PicInterlace || input->MbInterlace);
      mb_width= (input->img_width+img->auto_crop_right)>>4;
      mb_height= (input->img_height+img->auto_crop_bottom)>>4;
      mapunit_height=mb_height/(2-frame_mb_only);

      input->slice_group_id=(byte * ) malloc(sizeof(byte)*mapunit_height*mb_width);
      if (NULL==input->slice_group_id)
      {
        fclose(sgfile);
        no_mem_exit("PatchInp: input->slice_group_id");
      }

      // each line contains slice_group_id for one Macroblock
      for (i=0;i<mapunit_height*mb_width;i++)
      {
        ret = fscanf(sgfile,"%d", &tmp);
        input->slice_group_id[i]= (byte) tmp;
        if ( 1!=ret )
        {
          fclose(sgfile);
          snprintf(errortext, ET_SIZE, "Error while reading slice group config file (line %d)", i + 1);
          error (errortext, 500);
        }
        if ( *(input->slice_group_id+i) > input->num_slice_groups_minus1 )
        {
          fclose(sgfile);
          snprintf(errortext, ET_SIZE, "Error while reading slice group config file: slice_group_id not allowed (line %d)", i + 1);
          error (errortext, 500);
        }
        fscanf(sgfile,"%*[^\n]");
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
static void PatchInp (void)
{
  int bitdepth_qp_scale = 6*(input->BitDepthLuma - 8);

  // These variables are added for FMO
  int i,j;
  int storedBplus1;

  TestEncoderParams(bitdepth_qp_scale);

  if (input->FrameRate == 0.0)
    input->FrameRate = INIT_FRAME_RATE;

  // Set block sizes

  // Skip/Direct16x16
  input->part_size[0][0] = 4;
  input->part_size[0][1] = 4;
  // 16x16
  input->part_size[1][0] = 4;
  input->part_size[1][1] = 4;
  // 16x8
  input->part_size[2][0] = 4;
  input->part_size[2][1] = 2;
  // 8x16
  input->part_size[3][0] = 2;
  input->part_size[3][1] = 4;
  // 8x8
  input->part_size[4][0] = 2;
  input->part_size[4][1] = 2;
  // 8x4
  input->part_size[5][0] = 2;
  input->part_size[5][1] = 1;
  // 4x8
  input->part_size[6][0] = 1;
  input->part_size[6][1] = 2;
  // 4x4
  input->part_size[7][0] = 1;
  input->part_size[7][1] = 1;

  input->blocktype_lut[0][0] = 7; // 4x4
  input->blocktype_lut[0][1] = 6; // 4x8
  input->blocktype_lut[1][0] = 5; // 8x4
  input->blocktype_lut[1][1] = 4; // 8x8
  input->blocktype_lut[1][3] = 3; // 8x16
  input->blocktype_lut[3][1] = 2; // 16x8
  input->blocktype_lut[3][3] = 1; // 16x16

  for (j = 0; j<8;j++)
  {
    for (i = 0; i<2; i++)
    {
      input->blc_size[j][i] = input->part_size[j][i] * BLOCK_SIZE;
    }
  }

  if (input->idr_period && input->intra_delay && input->idr_period <= input->intra_delay)
  {
    snprintf(errortext, ET_SIZE, " IntraDelay cannot be larger than or equal to IDRPeriod.");
    error (errortext, 500);
  }

  if (input->idr_period && input->intra_delay && input->EnableIDRGOP == 0)
  {
    snprintf(errortext, ET_SIZE, " IntraDelay can only be used with only 1 IDR or with EnableIDRGOP=1.");
    error (errortext, 500);
  }

  if (input->idr_period && input->intra_delay && input->adaptive_idr_period)
  {
    snprintf(errortext, ET_SIZE, " IntraDelay can not be used with AdaptiveIDRPeriod.");
    error (errortext, 500);
  }

  storedBplus1 = (input->BRefPictures ) ? input->successive_Bframe + 1: 1;
  
  if (input->Log2MaxFNumMinus4 == -1)
  {    
    log2_max_frame_num_minus4 = iClip3(0,12, (int) (CeilLog2(input->no_frames * storedBplus1) - 4));    
  }
  else  
    log2_max_frame_num_minus4 = input->Log2MaxFNumMinus4;
  max_frame_num = 1 << (log2_max_frame_num_minus4 + 4);

  if (log2_max_frame_num_minus4 == 0 && input->num_ref_frames == 16)
  {
    snprintf(errortext, ET_SIZE, " NumberReferenceFrames=%d and Log2MaxFNumMinus4=%d may lead to an invalid value of frame_num.", input->num_ref_frames, input-> Log2MaxFNumMinus4);
    error (errortext, 500);
  }

  // set proper log2_max_pic_order_cnt_lsb_minus4.
  if (input->Log2MaxPOCLsbMinus4 == - 1)
    log2_max_pic_order_cnt_lsb_minus4 = iClip3(0,12, (int) (CeilLog2( 2*input->no_frames * (input->jumpd + 1)) - 4));
  else
    log2_max_pic_order_cnt_lsb_minus4 = input->Log2MaxPOCLsbMinus4;
  max_pic_order_cnt_lsb = 1 << (log2_max_pic_order_cnt_lsb_minus4 + 4);

  if (((1<<(log2_max_pic_order_cnt_lsb_minus4 + 3)) < input->jumpd * 4) && input->Log2MaxPOCLsbMinus4 != -1)
    error("log2_max_pic_order_cnt_lsb_minus4 might not be sufficient for encoding. Increase value.",400);

  // B picture consistency check
  if(input->successive_Bframe > input->jumpd)
  {
    snprintf(errortext, ET_SIZE, "Number of B-frames %d can not exceed the number of frames skipped", input->successive_Bframe);
    error (errortext, 400);
  }

  // Direct Mode consistency check
  if(input->successive_Bframe && input->direct_spatial_mv_pred_flag != DIR_SPATIAL && input->direct_spatial_mv_pred_flag != DIR_TEMPORAL)
  {
    snprintf(errortext, ET_SIZE, "Unsupported direct mode=%d, use TEMPORAL=0 or SPATIAL=1", input->direct_spatial_mv_pred_flag);
    error (errortext, 400);
  }

  if (input->PicInterlace>0 || input->MbInterlace>0)
  {
    if (input->directInferenceFlag==0)
      printf("\nWarning: DirectInferenceFlag set to 1 due to interlace coding.");
    input->directInferenceFlag = 1;
  }

  // Open Files
  if ((p_in=open(input->infile, OPENFLAGS_READ))==-1)
  {
    snprintf(errortext, ET_SIZE, "Input file %s does not exist",input->infile);
    error (errortext, 500);
  }

  if (input->yuv_format != YUV400)
  {
    input->img_width_cr  = (input->img_width  * mb_width_cr [input->yuv_format]) >> 4;
    input->img_height_cr = (input->img_height * mb_height_cr[input->yuv_format]) >> 4;
  }
  else
  {
    input->img_width_cr  = 0;
    input->img_height_cr = 0;
  }

  if (input->no_frames == -1)
  {
    getNumberOfFrames();
  }

  if (input->no_frames < 1)
  {      
    snprintf(errortext, ET_SIZE, "Not enough frames to encode (%d)", input->no_frames);
    error (errortext, 500);
  }


  if (strlen (input->ReconFile) > 0 && (p_dec=open(input->ReconFile, OPENFLAGS_WRITE, OPEN_PERMISSIONS))==-1)
  {
    snprintf(errortext, ET_SIZE, "Error open file %s", input->ReconFile);
    error (errortext, 500);
  }

#if TRACE
  if (strlen (input->TraceFile) > 0 && (p_trace=fopen(input->TraceFile,"w"))==NULL)
  {
    snprintf(errortext, ET_SIZE, "Error open file %s", input->TraceFile);
    error (errortext, 500);
  }
#endif

  if (input->img_width % 16 != 0)
  {
    img->auto_crop_right = 16-(input->img_width % 16);
  }
  else
  {
    img->auto_crop_right=0;
  }
  if (input->PicInterlace || input->MbInterlace)
  {
    if (input->img_height % 2 != 0)
    {
      error ("even number of lines required for interlaced coding", 500);
    }
    if (input->img_height % 32 != 0)
    {
      img->auto_crop_bottom = 32-(input->img_height % 32);
    }
    else
    {
      img->auto_crop_bottom=0;
    }
  }
  else
  {
    if (input->img_height % 16 != 0)
    {
      img->auto_crop_bottom = 16-(input->img_height % 16);
    }
    else
    {
      img->auto_crop_bottom=0;
    }
  }
  if (img->auto_crop_bottom || img->auto_crop_right)
  {
    fprintf (stderr, "Warning: Automatic cropping activated: Coded frame Size: %dx%d\n", input->img_width+img->auto_crop_right, input->img_height+img->auto_crop_bottom);
  }

  if ((input->slice_mode==1)&&(input->MbInterlace!=0))
  {
    if ((input->slice_argument%2)!=0)
    {
      fprintf ( stderr, "Warning: slice border within macroblock pair. ");
      if (input->slice_argument > 1)
      {
        input->slice_argument--;
      }
      else
      {
        input->slice_argument++;
      }
      fprintf ( stderr, "Using %d MBs per slice.\n", input->slice_argument);
    }
  }

  // read the slice group configuration file. Only for types 0, 2 or 6
  if ( 0 != input->num_slice_groups_minus1 )
  {
    read_slice_group_info();
  }

  if (input->ReferenceReorder && (input->PicInterlace || input->MbInterlace))
  {
    snprintf(errortext, ET_SIZE, "ReferenceReorder Not supported with Interlace encoding methods\n");
    error (errortext, 400);
  }

  if (input->PocMemoryManagement && (input->PicInterlace || input->MbInterlace))
  {
    snprintf(errortext, ET_SIZE, "PocMemoryManagement not supported with Interlace encoding methods\n");
    error (errortext, 400);
  }

  // frame/field consistency check
  if (input->PicInterlace != FRAME_CODING && input->PicInterlace != ADAPTIVE_CODING && input->PicInterlace != FIELD_CODING)
  {
    snprintf (errortext, ET_SIZE, "Unsupported PicInterlace=%d, use frame based coding=0 or field based coding=1 or adaptive=2",input->PicInterlace);
    error (errortext, 400);
  }

  // frame/field consistency check
  if (input->MbInterlace != FRAME_CODING && input->MbInterlace != ADAPTIVE_CODING && input->MbInterlace != FIELD_CODING && input->MbInterlace != FRAME_MB_PAIR_CODING)
  {
    snprintf (errortext, ET_SIZE, "Unsupported MbInterlace=%d, use frame based coding=0 or field based coding=1 or adaptive=2 or frame MB pair only=3",input->MbInterlace);
    error (errortext, 400);
  }


  if ((!input->rdopt)&&(input->MbInterlace))
  {
    snprintf(errortext, ET_SIZE, "MB AFF is not compatible with non-rd-optimized coding.");
    error (errortext, 500);
  }

  /*if (input->rdopt>2)
  {
    snprintf(errortext, ET_SIZE, "RDOptimization=3 mode has been deactivated do to diverging of real and simulated decoders.");
    error (errortext, 500);
  }*/

  // check RDoptimization mode and profile. FMD does not support Frex Profiles.
  if (input->rdopt==2 && ( input->ProfileIDC>=FREXT_HP || input->ProfileIDC==FREXT_CAVLC444 ))
  {
    snprintf(errortext, ET_SIZE, "Fast Mode Decision methods not supported in FREX Profiles");
    error (errortext, 500);
  }

  if ( (input->MEErrorMetric[Q_PEL] == ERROR_SATD && input->MEErrorMetric[H_PEL] == ERROR_SAD && input->MEErrorMetric[F_PEL] == ERROR_SAD)
    && input->SearchMode > FAST_FULL_SEARCH && input->SearchMode < EPZS)
  {
    snprintf(errortext, ET_SIZE, "MEDistortionQPel=2, MEDistortionHPel=0, MEDistortionFPel=0 is not allowed when SearchMode is set to 1 or 2.");
    error (errortext, 500);
  }

  // Tian Dong: May 31, 2002
  // The number of frames in one sub-seq in enhanced layer should not exceed
  // the number of reference frame number.
  if ( input->NumFramesInELSubSeq > input->num_ref_frames || input->NumFramesInELSubSeq < 0 )
  {
    snprintf(errortext, ET_SIZE, "NumFramesInELSubSeq (%d) is out of range [0,%d).", input->NumFramesInELSubSeq, input->num_ref_frames);
    error (errortext, 500);
  }
  // Tian Dong: Enhanced GOP is not supported in bitstream mode. September, 2002
  if ( input->NumFramesInELSubSeq > 0 && input->of_mode == PAR_OF_ANNEXB )
  {
    snprintf(errortext, ET_SIZE, "Enhanced GOP is not supported in bitstream mode and RTP mode yet.");
    error (errortext, 500);
  }
  // Tian Dong (Sept 2002)
  // The AFF is not compatible with spare picture for the time being.
  if ((input->PicInterlace || input->MbInterlace) && input->SparePictureOption == TRUE)
  {
    snprintf(errortext, ET_SIZE, "AFF is not compatible with spare picture.");
    error (errortext, 500);
  }

  // Only the RTP mode is compatible with spare picture for the time being.
  if (input->of_mode != PAR_OF_RTP && input->SparePictureOption == TRUE)
  {
    snprintf(errortext, ET_SIZE, "Only RTP output mode is compatible with spare picture features.");
    error (errortext, 500);
  }

  if( (input->WeightedPrediction > 0 || input->WeightedBiprediction > 0) && (input->MbInterlace))
  {
    snprintf(errortext, ET_SIZE, "Weighted prediction coding is not supported for MB AFF currently.");
    error (errortext, 500);
  }
  if ( input->NumFramesInELSubSeq > 0 && input->WeightedPrediction > 0)
  {
    snprintf(errortext, ET_SIZE, "Enhanced GOP is not supported in weighted prediction coding mode yet.");
    error (errortext, 500);
  }

  //! the number of slice groups is forced to be 1 for slice group type 3-5
  if(input->num_slice_groups_minus1 > 0)
  {
    if( (input->slice_group_map_type >= 3) && (input->slice_group_map_type<=5) )
      input->num_slice_groups_minus1 = 1;
  }

  // Rate control
  if(input->RCEnable)
  {
    if ( ((input->img_height+img->auto_crop_bottom)*(input->img_width+img->auto_crop_right)/256) % input->basicunit != 0)
    {
      snprintf(errortext, ET_SIZE, "Frame size in macroblocks must be a multiple of BasicUnit.");
      error (errortext, 500);
    }
    // input basicunit represents the number of BUs in a frame
    // internal basicunit is the number of MBs in a BU
    input->basicunit = ((input->img_height+img->auto_crop_bottom)*(input->img_width+img->auto_crop_right)/256) / input->basicunit;

    if ( input->RCEnable && (input->successive_Bframe || input->jumpd) && input->RCUpdateMode == RC_MODE_1 )
    {
      snprintf(errortext, ET_SIZE, "Use RCUpdateMode = 1 only for all intra or all B-slice coding.");
      error (errortext, 500);
    }

    if ( input->RCEnable && input->BRefPictures == 2 && input->intra_period == 0 && input->RCUpdateMode != RC_MODE_1 )
    {
      snprintf(errortext, ET_SIZE, "Use RCUpdateMode = 1 for all B-slice coding.");
      error (errortext, 500);
    }

    if ( input->RCEnable && input->HierarchicalCoding && input->RCUpdateMode != RC_MODE_2 && input->RCUpdateMode != RC_MODE_3 )
    {
      snprintf(errortext, ET_SIZE, "Use RCUpdateMode = 2 or 3 for hierarchical B-picture coding.");
      error (errortext, 500);
    }

    if ( input->RCEnable && (input->RCUpdateMode != RC_MODE_1) && (input->intra_period == 1) )
    {
      snprintf(errortext, ET_SIZE, "Use RCUpdateMode = 1 for all intra coding.");
      error (errortext, 500);
    }
  }

  if ((input->successive_Bframe)&&(input->BRefPictures)&&(input->idr_period)&&(input->pic_order_cnt_type!=0))
  {
    error("Stored B pictures combined with IDR pictures only supported in Picture Order Count type 0\n",-1000);
  }

  if( !input->direct_spatial_mv_pred_flag && input->num_ref_frames<2 && input->successive_Bframe >0)
    error("temporal direct needs at least 2 ref frames\n",-1000);

  // frext
  if(input->Transform8x8Mode && input->sp_periodicity /*SP-frames*/)
  {
    snprintf(errortext, ET_SIZE, "\nThe new 8x8 mode is not implemented for sp-frames.");
    error (errortext, 500);
  }

  if(input->Transform8x8Mode && ( input->ProfileIDC<FREXT_HP && input->ProfileIDC!=FREXT_CAVLC444 ))
  {
    snprintf(errortext, ET_SIZE, "\nTransform8x8Mode may be used only with ProfileIDC %d to %d.", FREXT_HP, FREXT_Hi444);
    error (errortext, 500);
  }
  if(input->ScalingMatrixPresentFlag && ( input->ProfileIDC<FREXT_HP && input->ProfileIDC!=FREXT_CAVLC444 ))
  {
    snprintf(errortext, ET_SIZE, "\nScalingMatrixPresentFlag may be used only with ProfileIDC %d to %d.", FREXT_HP, FREXT_Hi444);
    error (errortext, 500);
  }

  if(input->yuv_format==YUV422 && ( input->ProfileIDC < FREXT_Hi422 && input->ProfileIDC!=FREXT_CAVLC444 ))
  {
    snprintf(errortext, ET_SIZE, "\nFRExt Profile(YUV Format) Error!\nYUV422 can be used only with ProfileIDC %d or %d\n",FREXT_Hi422, FREXT_Hi444);
    error (errortext, 500);
  }
  if(input->yuv_format==YUV444 && ( input->ProfileIDC < FREXT_Hi444 && input->ProfileIDC!=FREXT_CAVLC444 ))
  {
    snprintf(errortext, ET_SIZE, "\nFRExt Profile(YUV Format) Error!\nYUV444 can be used only with ProfileIDC %d.\n",FREXT_Hi444);
    error (errortext, 500);
  }

  if (input->successive_Bframe && ((input->BiPredMotionEstimation) && (input->search_range < input->BiPredMESearchRange)))
  {
    snprintf(errortext, ET_SIZE, "\nBiPredMESearchRange must be smaller or equal SearchRange.");
    error (errortext, 500);
  }

  // check consistency
  if ( input->ChromaMEEnable && !(input->ChromaMCBuffer) ) 
  {
    snprintf(errortext, ET_SIZE, "\nChromaMCBuffer must be set to 1 if ChromaMEEnable is set.");
    error (errortext, 500);
  }

  if ( input->ChromaMEEnable && input->yuv_format ==  YUV400) 
  {
    fprintf(stderr, "Warning: ChromaMEEnable cannot be used with YUV400 color format, disabling ChromaMEEnable.\n");
    input->ChromaMEEnable = 0;
  }

  if ( (0 == input->ChromaMCBuffer) && (( input->yuv_format ==  YUV444) && (!input->separate_colour_plane_flag)) )
  {
    fprintf(stderr, "Warning: Enabling ChromaMCBuffer for YUV444 combined color coding.\n");
    input->ChromaMCBuffer = 1;
  }

  if (input->EnableOpenGOP && input->PicInterlace)
  {
    snprintf(errortext, ET_SIZE, "Open GOP currently not supported for Field coded pictures.");
    error (errortext, 500);
  }

  if (input->EnableOpenGOP)
    input->ReferenceReorder = 1;

  input->EPZSGrid = input->EPZSSubPelGrid << 1;

  if (input->redundant_pic_flag)
  {
    if (input->PicInterlace || input->MbInterlace)
    {
      snprintf(errortext, ET_SIZE, "Redundant pictures cannot be used with interlaced tools.");
      error (errortext, 500);
    }
    if (input->RDPictureDecision)
    {
      snprintf(errortext, ET_SIZE, "Redundant pictures cannot be used with RDPictureDecision.");
      error (errortext, 500);
    }
    if (input->successive_Bframe)
    {
      snprintf(errortext, ET_SIZE, "Redundant pictures cannot be used with B frames.");
      error (errortext, 500);
    }
    if (input->PrimaryGOPLength < (1 << input->NumRedundantHierarchy))
    {
      snprintf(errortext, ET_SIZE, "PrimaryGOPLength must be equal or greater than 2^NumRedundantHierarchy.");
      error (errortext, 500);
    }
    if (input->num_ref_frames < input->PrimaryGOPLength)
    {
      snprintf(errortext, ET_SIZE, "NumberReferenceFrames must be greater than or equal to PrimaryGOPLength.");
      error (errortext, 500);
    }
  }

  if (input->num_ref_frames == 1 && input->successive_Bframe)
  {
    fprintf( stderr, "\nWarning: B slices used but only one reference allocated within reference buffer.\n");
    fprintf( stderr, "         Performance may be considerably compromised! \n");
    fprintf( stderr, "         2 or more references recommended for use with B slices.\n");
  }
  if ((input->HierarchicalCoding || input->BRefPictures) && input->successive_Bframe)
  {
    fprintf( stderr, "\nWarning: Hierarchical coding or Referenced B slices used.\n");
    fprintf( stderr, "         Make sure that you have allocated enough references\n");
    fprintf( stderr, "         in reference buffer to achieve best performance.\n");
  }

  ProfileCheck();
  LevelCheck();
}

void PatchInputNoFrames(void)
{
  // Tian Dong: May 31, 2002
  // If the frames are grouped into two layers, "FramesToBeEncoded" in the config file
  // will give the number of frames which are in the base layer. Here we let input->no_frames
  // be the total frame numbers.
  input->no_frames = 1 + (input->no_frames - 1) * (input->NumFramesInELSubSeq + 1);
}

static void ProfileCheck(void)
{
  if((input->ProfileIDC != 66 ) &&
     (input->ProfileIDC != 77 ) &&
     (input->ProfileIDC != 88 ) &&
     (input->ProfileIDC != FREXT_HP    ) &&
     (input->ProfileIDC != FREXT_Hi10P ) &&
     (input->ProfileIDC != FREXT_Hi422 ) &&
     (input->ProfileIDC != FREXT_Hi444 ) &&
     (input->ProfileIDC != FREXT_CAVLC444 ))
  {
    snprintf(errortext, ET_SIZE, "Profile must be in\n\n  66 (Baseline),\n  77 (Main),\n  88 (Extended),\n 100 (High),\n 110 (High 10 or High 10 Intra)\n"
      " 122 (High 4:2:2 or High 4:2:2 Intra),\n 244 (High 4:4:4 predictive or High 4:4:4 Intra),\n  44 (CAVLC 4:4:4 Intra)\n");
    error (errortext, 500);
  }

  if ((input->partition_mode) && (input->symbol_mode==CABAC))
  {
    snprintf(errortext, ET_SIZE, "Data partitioning and CABAC is not supported in any profile.");
    error (errortext, 500);
  }

  if (input->redundant_pic_flag)
  {
    if (input->ProfileIDC != 66)
    {
      snprintf(errortext, ET_SIZE, "Redundant pictures are only allowed in Baseline profile (ProfileIDC = 66).");
      error (errortext, 500);
    }
  }

  if (input->partition_mode)
  {
    if (input->ProfileIDC != 88)
    {
      snprintf(errortext, ET_SIZE, "Data partitioning is only allowed in Extended profile (ProfileIDC = 88).");
      error (errortext, 500);
    }
  }

  if (input->ChromaIntraDisable && input->FastCrIntraDecision)
  {
    fprintf( stderr, "\n Warning: ChromaIntraDisable and FastCrIntraDecision cannot be combined together.\n Using only Chroma Intra DC mode.\n");
    input->FastCrIntraDecision=0;
  }

  if ((input->sp_periodicity) && (input->ProfileIDC != 88 ))
  {
    snprintf(errortext, ET_SIZE, "SP pictures are only allowed in Extended profile (ProfileIDC = 88).");
    error (errortext, 500);
  }

  // baseline
  if (input->ProfileIDC == 66 )
  {
    if ((input->successive_Bframe || input->BRefPictures==2) && input->PReplaceBSlice == 0)
    {
      snprintf(errortext, ET_SIZE, "B slices are not allowed in Baseline profile (ProfileIDC = 66).");
      error (errortext, 500);
    }
    if (input->WeightedPrediction)
    {
      snprintf(errortext, ET_SIZE, "Weighted prediction is not allowed in Baseline profile (ProfileIDC = 66).");
      error (errortext, 500);
    }
    if (input->WeightedBiprediction)
    {
      snprintf(errortext, ET_SIZE, "Weighted prediction is not allowed in Baseline profile (ProfileIDC = 66).");
      error (errortext, 500);
    }
    if (input->symbol_mode == CABAC)
    {
      snprintf(errortext, ET_SIZE, "CABAC is not allowed in Baseline profile (ProfileIDC = 66).");
      error (errortext, 500);
    }
  }

  // main
  if (input->ProfileIDC == 77 )
  {
    if (input->num_slice_groups_minus1)
    {
      snprintf(errortext, ET_SIZE, "num_slice_groups_minus1>0 (FMO) is not allowed in Main profile (ProfileIDC = 77).");
      error (errortext, 500);
    }
  }

  // extended
  if (input->ProfileIDC == 88 )
  {
    if (!input->directInferenceFlag)
    {
      snprintf(errortext, ET_SIZE, "direct_8x8_inference flag must be equal to 1 in Extended profile (ProfileIDC = 88).");
      error (errortext, 500);
    }

    if (input->symbol_mode == CABAC)
    {
      snprintf(errortext, ET_SIZE, "CABAC is not allowed in Extended profile (ProfileIDC = 88).");
      error (errortext, 500);
    }
  }

  //FRExt
  if ( input->separate_colour_plane_flag )
  {
    if( input->yuv_format!=3 )
    {
      fprintf( stderr, "\nWarning: SeparateColourPlane has only effect in 4:4:4 chroma mode (YUVFormat=3),\n         disabling SeparateColourPlane.");
      input->separate_colour_plane_flag = 0;
    }

    if ( input->ChromaMEEnable )
    {
      snprintf(errortext, ET_SIZE, "\nChromaMEEnable is not allowed when SeparateColourPlane is enabled.");
      error (errortext, 500);
    }
  }

  // CAVLC 4:4:4 Intra
  if ( input->ProfileIDC == FREXT_CAVLC444 )
  {
    if ( input->symbol_mode != CAVLC )
    {
      snprintf(errortext, ET_SIZE, "\nCABAC is not allowed in CAVLC 4:4:4 Intra profile (ProfileIDC = 44).");
      error (errortext, 500);
    }
    if ( !input->IntraProfile )
    {
      fprintf (stderr, "\nWarning: ProfileIDC equal to 44 implies an Intra only profile, setting IntraProfile = 1.");
      input->IntraProfile = 1;
    }
  }

  // Intra only profiles
  if (input->IntraProfile && ( input->ProfileIDC<FREXT_HP && input->ProfileIDC!=FREXT_CAVLC444 ))
  {
    snprintf(errortext, ET_SIZE, "\nIntraProfile is allowed only with ProfileIDC %d to %d.", FREXT_HP, FREXT_Hi444);
    error (errortext, 500);
  }

  if (input->IntraProfile && !input->idr_period) 
  {
    snprintf(errortext, ET_SIZE, "\nIntraProfile requires IDRPeriod >= 1.");
    error (errortext, 500);
  }

  if (input->IntraProfile && input->intra_period != 1) 
  {
    snprintf(errortext, ET_SIZE, "\nIntraProfile requires IntraPeriod equal 1.");
    error (errortext, 500);
  }

  if (input->IntraProfile && input->num_ref_frames) 
  {
    fprintf( stderr, "\nWarning: Setting NumberReferenceFrames to 0 in IntraProfile.");
    input->num_ref_frames = 0;
  }

  if (input->IntraProfile == 0 && input->num_ref_frames == 0) 
  {
    snprintf(errortext, ET_SIZE, "\nProfiles other than IntraProfile require NumberReferenceFrames > 0.");
    error (errortext, 500);
  }
}


// Level Limit             -  -  -  -  -  -  -  -  -  1b  10  11   12   13   -  -  -  -  -  -  20   21   22    -  -  -  -  -  -  -
unsigned int  MaxFs [] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 99, 99, 396, 396, 396, 0, 0, 0, 0, 0, 0, 396, 792, 1620, 0, 0, 0, 0, 0, 0, 0,  
//                        30    31    32    -  -  -  -  -  -  -  40    41    42    -  -  -  -  -  -  -  50     51
                          1620, 3600, 5120, 0, 0, 0, 0, 0, 0, 0, 8192, 8192, 8704, 0, 0, 0, 0, 0, 0, 0, 22080, 36864 };

unsigned getMaxFs (unsigned int levelIdc)
{
  unsigned int ret;

  if ( (levelIdc < 9) || (levelIdc > 51))
    error ("getMaxFs: Unknown LevelIdc", 500);

  // in Baseline, Main and Extended: Level 1b is specified with LevelIdc==11 and constrained_set3_flag == 1

  ret = MaxFs[levelIdc];

  if ( 0 == ret )
    error ("getMaxFs: Unknown LevelIdc", 500);

  return ret;
}

static void LevelCheck(void)
{
  unsigned int PicSizeInMbs = ( (input->img_width + img->auto_crop_right) * (input->img_height + img->auto_crop_bottom) ) >> 8;

  if ( (input->LevelIDC>=30) && (input->directInferenceFlag==0))
  {
    fprintf( stderr, "\nWarning: LevelIDC 3.0 and above require direct_8x8_inference to be set to 1. Please check your settings.\n");
    input->directInferenceFlag=1;
  }
  if ( ((input->LevelIDC<21) || (input->LevelIDC>41)) && (input->PicInterlace > 0 || input->MbInterlace > 0) )
  {
    snprintf(errortext, ET_SIZE, "\nInterlace modes only supported for LevelIDC in the range of 21 and 41. Please check your settings.\n");
    error (errortext, 500);
  }

  if ( PicSizeInMbs > getMaxFs(input->LevelIDC) )
  {
    snprintf(errortext, ET_SIZE, "\nPicSizeInMbs exceeds maximum allowed size at specified LevelIdc %d\n", input->LevelIDC);
    error (errortext, 500);
  }
  
  if (input->IntraProfile && (PicSizeInMbs > 1620) && input->slice_mode != 1) 
  {
    error ("\nIntraProfile with PicSizeInMbs > 1620 requires SliceMode equal 1.", 500);
  }

  if (input->IntraProfile && (PicSizeInMbs > 1620) && ((unsigned int)input->slice_argument > (  getMaxFs(input->LevelIDC) >> 2 ) ) )
  {
    //when PicSizeInMbs is greater than 1620, the number of macroblocks in any coded slice shall not exceed MaxFS / 4
    snprintf(errortext, ET_SIZE, "\nIntraProfile requires SliceArgument smaller or equal to 1/4 MaxFs at specified LevelIdc %d.", input->LevelIDC);
    error (errortext, 500);
  }
}
