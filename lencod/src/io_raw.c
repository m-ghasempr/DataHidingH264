/*!
 *************************************************************************************
 * \file io_raw.c
 *
 * \brief
 *    I/O functions related to raw images
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Karsten Sühring                 <suehring@hhi.de>
 *     - Alexis Michael Tourapis         <alexismt@ieee.org>
 *     
 *************************************************************************************
 */
#include "contributors.h"

#include "global.h"
#include "report.h"
#include "img_io.h"

//unsigned char *buf;
#define FAST_READ 1

#if FAST_READ
static inline void ReadData (int vfile,  FrameFormat *source, unsigned char *buf)
{
  unsigned char *cur_buf = buf;
  int i, j;
  for (i = 0; i < source->height; i++)
  {
    if (read(vfile, cur_buf, source->width) != source->width)
    {
      printf ("ReadOneFrame: cannot read %d bytes from input file, unexpected EOF, exiting...\n", source->width);
      report_stats_on_error();
    }
    cur_buf += source->width;
  }
  for (j = 0; j < 2; j++)
  {
    for (i = 0; i < source->height_cr; i++)
    {
      if (read(vfile, cur_buf, source->width_cr) != source->width_cr)
      {
        printf ("ReadOneFrame: cannot read %d bytes from input file, unexpected EOF, exiting...\n", source->width_cr);
        report_stats_on_error();
      }
      cur_buf += source->width_cr;
    }
  }
}
#else
static inline void ReadData (int vfile, int framesize_in_bytes, unsigned char *buf)
{
  if (read(vfile, buf, (int) framesize_in_bytes) != (int) framesize_in_bytes)
  {
    printf ("ReadOneFrame: cannot read %d bytes from input file, unexpected EOF, exiting...\n", (int) framesize_in_bytes);
    report_stats_on_error();
  }
}
#endif


/*!
 ************************************************************************
 * \brief
 *    Reads one new frame from concatenated raw file
 *    Code does not distinguish between planar and interleaved data
 *
 * \param input_file
 *    Input file to read from
 * \param FrameNoInFile
 *    Frame number in the source file
 * \param HeaderSize
 *    Number of bytes in the source file to be skipped
 * \param source
 *    source file (on disk) information 
 * \param buf
 *    image buffer data
 ************************************************************************
 */
void ReadFrameConcatenated (VideoDataFile *input_file, int FrameNoInFile, int HeaderSize, FrameFormat *source, unsigned char *buf)
{
  int vfile = input_file->f_num;
  unsigned int symbol_size_in_bytes = img->pic_unit_size_on_disk >> 3;

  const int imgsize_y = source->width * source->height;
  const int imgsize_uv = source->width_cr * source->height_cr;

  const int bytes_y = imgsize_y * symbol_size_in_bytes;
  const int bytes_uv = imgsize_uv * symbol_size_in_bytes;

  const int64 framesize_in_bytes = bytes_y + 2*bytes_uv;
  
#if 0
  // skip Header
  if (lseek (vfile, HeaderSize, SEEK_SET) != HeaderSize)
  {
    error ("ReadOneFrame: cannot fseek to (Header size) in input file", -1);
  }

  // skip starting frames
  if (lseek (vfile, framesize_in_bytes * params->start_frame, SEEK_CUR) == -1)
  {
    snprintf(errortext, ET_SIZE, "ReadOneFrame: cannot advance file pointer in input file beyond frame %d\n", params->start_frame);
    error (errortext,-1);
  }

  // seek to current frame
  if (lseek (vfile, framesize_in_bytes * (FrameNoInFile), SEEK_CUR) == -1)
  {
    snprintf(errortext, ET_SIZE, "ReadOneFrame: cannot advance file pointer in input file beyond frame %d\n", params->start_frame + FrameNoInFile);
    error (errortext,-1);
  }
#else
  // Let us seek directly to the current frame
  if (lseek (vfile, HeaderSize + framesize_in_bytes * (FrameNoInFile + params->start_frame), SEEK_SET) == -1)
  {
    snprintf(errortext, ET_SIZE, "ReadOneFrame: cannot advance file pointer in input file beyond frame %d\n", params->start_frame + FrameNoInFile);
    error (errortext,-1);
  }
#endif
  // Here we are at the correct position for the source frame in the file.  
  // Now read it.
  if ((img->pic_unit_size_on_disk & 0x07) == 0)
  {
#if FAST_READ
    ReadData (vfile, source, buf);
#else
    ReadData (vfile, (int) framesize_in_bytes, buf);
#endif
  }
  else
  {
    printf ("ReadOneFrame (NOT IMPLEMENTED): pic unit size on disk must be divided by 8");
    exit (-1);
  }
}

/*!
 ************************************************************************
 * \brief
 *    Reads one new frame from separate data files
 *    Code does not distinguish between planar and interleaved data
 *
 * \param input_file
 *    Input file to read from
 * \param FrameNoInFile
 *    Frame number in the source file
 * \param HeaderSize
 *    Number of bytes in the source file to be skipped
 * \param source
 *    source file (on disk) information 
 * \param buf
 *    taget buffer
 ************************************************************************
 */
void ReadFrameSeparate (VideoDataFile *input_file, int FrameNoInFile, int HeaderSize, FrameFormat *source, unsigned char *buf)
{
  int vfile = input_file->f_num;

  OpenFrameFile( input_file, params, FrameNoInFile + params->start_frame);

  // skip Header
  if (lseek (vfile, HeaderSize, SEEK_SET) != HeaderSize)
  {
    error ("ReadOneFrame: cannot fseek to (Header size) in input file", -1);
  }

  // Read data
  if ((img->pic_unit_size_on_disk & 0x07) == 0)
  {
#if FAST_READ
    ReadData (vfile, source, buf);
#else
    unsigned int symbol_size_in_bytes = img->pic_unit_size_on_disk >> 3;

    const int imgsize_y = source->width * source->height;
    const int imgsize_uv = source->width_cr * source->height_cr;

    const int bytes_y = imgsize_y * symbol_size_in_bytes;
    const int bytes_uv = imgsize_uv * symbol_size_in_bytes;
    const int64 framesize_in_bytes = bytes_y + 2*bytes_uv;

    ReadData (vfile, (int) framesize_in_bytes, buf);
#endif
  }
  else
  {
    printf ("ReadOneFrame (NOT IMPLEMENTED): pic unit size on disk must be divided by 8");
    exit (-1);
  }

  if (vfile != -1)
    close(vfile);
}
