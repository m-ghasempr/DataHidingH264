
/*!
 ************************************************************************
 * \file io_raw.h
 *
 * \brief
 *    I/O functions related to raw images
 *
 * \author
 *     - Alexis Michael Tourapis         <alexismt@ieee.org>
 *
 ************************************************************************
 */

#ifndef _IO_RAW_H_
#define _IO_RAW_H_

extern void ReadFrameConcatenated  (VideoDataFile *input_file, int FrameNoInFile, int HeaderSize, FrameFormat *source, unsigned char *buf);
extern void ReadFrameSeparate      (VideoDataFile *input_file, int FrameNoInFile, int HeaderSize, FrameFormat *source, unsigned char *buf);

#endif

