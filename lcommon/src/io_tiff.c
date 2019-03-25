/*!
 *************************************************************************************
 * \file io_tiff.c
 *
 * \brief
 *    I/O functions related to TIFF images
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Larry Luther                    <lzl@dolby.com>
 *     - Alexis Michael Tourapis         <alexismt@ieee.org>
 *     
 *************************************************************************************
 */
#include "contributors.h"

#include "global.h"
#include "report.h"
#include "io_tiff.h"

int64 bufferSize = 0;


/*!
 ************************************************************************
 * \brief
 *    Open file containing a single frame
 ************************************************************************
 */
void OpenTiffFile( VideoDataFile *input_file, InputParameters *p_Inp, int FrameNumberInFile)
{
  char infile [FILE_NAME_SIZE], in_number[16];
  int length = 0;

  if (input_file->num_digits > 0)
  {
    in_number[length]='\0';
    length = strlen(input_file->fhead);
    strncpy(infile, input_file->fhead, length);
    infile[length]='\0';
    if (input_file->zero_pad)       
      snprintf(in_number, 16, "%0*d", input_file->num_digits, FrameNumberInFile);
    else
      snprintf(in_number, 16, "%*d", input_file->num_digits, FrameNumberInFile);

    strncat(infile, in_number, sizeof(in_number));
    length += sizeof(in_number);
    infile[length]='\0';
    strncat(infile, input_file->ftail,strlen(input_file->ftail));
    length += strlen(input_file->ftail);
    infile[length]='\0';
  }
  else
  {
    strcpy(infile, input_file->fname);    
  }

  if ((input_file->f_num = open(infile, OPENFLAGS_READ)) == -1)
  {
    printf ("OpenTiffFile: cannot open file %s\n", infile);
    report_stats_on_error();
  }    
}

/*!
 ************************************************************************
 * \brief
 *    Create Frame Memory buffer
 *
 ************************************************************************
 */
int AllocateTIFFBufferMemory (unsigned char *buf, int buffersize)
{
  if (NULL == (buf = (unsigned char *) realloc (buf, buffersize)))
    return (1);
  else
    return (0);
}

/*!
 ************************************************************************
 * \brief
 *    Delete Frame Memory buffer
 *
 ************************************************************************
 */
void DeleteTiffBufferMemory (unsigned char *buf)
{
  if (buf)
  {
    free (buf);
    buf = NULL;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Reads Input File Size 
 *
 ************************************************************************
 */
static int64 ReadTIFFSize(int video_file)
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
 *    Reads Tiff File Header
 *
 ************************************************************************
 */
static void ReadTIFFHeader (int vfile, TIFFHeader *tiff)
{
  if (read(vfile, tiff, sizeof(TIFFHeader)) != sizeof(TIFFHeader))
  {
    printf ("ReadTIFFHeader: cannot read header info from input file. Exiting...\n");
    report_stats_on_error();
  }    
}

/*!
 ************************************************************************
 * \brief
 *    Reads Tiff IFD entry
 *
 ************************************************************************
 */
static void ReadTIFFIFDEntry (int vfile, TIFFHeader *tiffHeader, uint16 *ifd_count, TIFFIFDEntry **tiffIFD)
{
  int i;
  
  //printf("%s (%d) tiffHeader->IFDoff %d\n", tiffHeader->border, tiffHeader->version, tiffHeader->IFDoff);

  if (lseek (vfile, tiffHeader->IFDoff, SEEK_SET) == -1)
  {    
    error ("ReadTIFFIFDEntry: cannot lseek to first IDF entry", -1);
  }

  if (read(vfile, ifd_count, sizeof(uint16)) != sizeof(uint16))
  {
    printf ("ReadTIFFIFDEntry: cannot read number of IFD entries from input file. Exiting...\n");
    report_stats_on_error();
  }

  if ((*tiffIFD = (TIFFIFDEntry *) malloc (*ifd_count* sizeof(TIFFIFDEntry))) == NULL)
  {
    printf ("ReadTIFFIFDEntry: cannot allocate memory for IFD entries. Exiting...\n");
    report_stats_on_error();
  }
  
  for (i = 0; i < *ifd_count; i++)
  {
    if (read(vfile, &(*tiffIFD)[i], sizeof(TIFFIFDEntry)) != sizeof(TIFFIFDEntry))
    {
      printf ("ReadTIFFIFDEntry: cannot read IFD entry from input file. Exiting...\n");
      report_stats_on_error();
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Parse Tiff IFD entries
 *
 ************************************************************************
 */
static void ParseTIFFIFD (uint16 ifd_count, TIFFIFDEntry *tiffIFD, FrameFormat *source)
{
  int i;
  
  for (i = 0; i < ifd_count; i++)
  {
    switch (tiffIFD[i].tIFD_tag)
    {
    case TIFFTAG_COMPRESSION :
      if (tiffIFD[i].tIFD_count != 1)
      {
        printf ("ParseTIFFIFD: Compressed TIFF data not supported. Exiting...\n");
        report_stats_on_error();
      }
      break;
    case TIFFTAG_IMAGEWIDTH : 
      if (tiffIFD[i].tIFD_offset != (unsigned int) source->width)
      {
        printf ("ParseTIFFIFD: Tiff width (%d) different from encoder input width (%d) . Exiting...\n", tiffIFD[i].tIFD_offset, source->width);
        report_stats_on_error();
      }
      break;
    case TIFFTAG_IMAGELENGTH :
      if (tiffIFD[i].tIFD_offset != (unsigned int) source->height)
      {
        printf ("ParseTIFFIFD: Tiff height different from encoder input height. Exiting...\n");
        report_stats_on_error();
      }
      break;
    case TIFFTAG_BITSPERSAMPLE :
      if (tiffIFD[i].tIFD_count != 3) 
      {
        printf ("ParseTIFFIFD: Only 3 channel TIFF files supported. Exiting...\n");
        report_stats_on_error();
      }
      //printf("values %d %d\n", (tiffIFD[i].tIFD_offset >> 24) & 0x07, source->bit_depth[0]);
      //printf("values %d %d\n", (tiffIFD[i].tIFD_offset >> 16) & 0x07, source->bit_depth[1]);
      //printf("values %d %d\n", (tiffIFD[i].tIFD_offset >> 8) & 0x07, source->bit_depth[2]);
      //printf("values %d %d\n", (tiffIFD[i].tIFD_offset     ) & 0x07, tiffIFD[i].tIFD_offset);
      /*
      if ((tiffIFD[i].tIFD_offset >> 24) & 0x07 != source->bit_depth[0])
      {
        printf ("ParseTIFFIFD: Tiff bitdepth different from encoder input bitdepth. Exiting...\n");
        report_stats_on_error();
      }
      if ((tiffIFD[i].tIFD_offset >> 16) & 0x07 != source->bit_depth[1])
      {
        printf ("ParseTIFFIFD: Tiff bitdepth different from encoder input bitdepth. Exiting...\n");
        report_stats_on_error();
      }
      if ((tiffIFD[i].tIFD_offset >> 8) & 0x07 != source->bit_depth[2])
      {
        printf ("ParseTIFFIFD: Tiff bitdepth different from encoder input bitdepth. Exiting...\n");
        report_stats_on_error();
      }
      */
      break;
    default:
      break;
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    Reads entire tiff file from harddrive. Any processing is done
 *    in memory, reducing I/O processing
 *
 * \param p_Inp
 *    Input configuration parameters 
 * \param input_file
 *    Input file to read from
 * \param FrameNoInFile
 *    Frame number in the source file
 * \param source
 *    source file (on disk) information 
 * \param buf
 *    memory buffer
 ************************************************************************
 */
int ReadTIFFImage (InputParameters *p_Inp, VideoDataFile *input_file, int FrameNoInFile, FrameFormat *source, unsigned char *buf)
{
  int file_read = 0;
  int64 fileSize = 0; 
  int f_num;
  
  TIFFHeader   tiffHeader;
  TIFFIFDEntry *tiffIFD = NULL;
  uint16 ifd_count = 0;
  
  OpenTiffFile( input_file, p_Inp, FrameNoInFile + p_Inp->start_frame);
  f_num = input_file->f_num;

  fileSize = ReadTIFFSize(f_num);

  fileSize = source->size;
  
  ReadTIFFHeader   (f_num, &tiffHeader);
  ReadTIFFIFDEntry (f_num, &tiffHeader, &ifd_count, &tiffIFD);
  ParseTIFFIFD     (ifd_count, tiffIFD, source);
#define DEBUG_TIF
#ifdef DEBUG_TIF
  {
    int i;
    for (i = 0; i < ifd_count; i++)
      printf("Value (%d) %d %d %d %d\n", i, tiffIFD[i].tIFD_tag, tiffIFD[i].tIFD_type, tiffIFD[i].tIFD_count, tiffIFD[i].tIFD_offset);
  }
#endif

  if (fileSize != bufferSize)
  {
    //DeleteTiffBufferMemory(buf);
    AllocateTIFFBufferMemory(buf, (int) fileSize);
    bufferSize = fileSize;
  }

  // Read data
  if (read(f_num, buf, (int) fileSize) != (int) fileSize)
  {
    printf ("ReadTIFFImage: cannot read %d bytes from input file, unexpected EOF, exiting...\n", (int) fileSize);
    file_read = 0;
  }
  else
  {
    file_read = 1;
  }

  close(input_file->f_num);  

  if (tiffIFD != NULL)
    free(tiffIFD);

  return file_read;
}
