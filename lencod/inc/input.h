
/*!
 ************************************************************************
 * \file input.h
 *
 * \brief
 *    Input related definitions
 *
 * \author
 *
 ************************************************************************
 */

#ifndef _INPUT_H_
#define _INPUT_H_

int testEndian(void);
void initInput(FrameFormat *source, FrameFormat *output);
void AllocateFrameMemory (ImageParameters *img, InputParameters *params, int size);
void DeleteFrameMemory (void);

void ReadOneFrame (VideoDataFile *input_file, int FrameNoInFile, int HeaderSize, FrameFormat *source, FrameFormat *output, imgpel **pImage[3]);
extern void (*buf2img) ( imgpel** imgX, unsigned char* buf, int size_x, int size_y, int osize_x, int o_size_y, int symbol_size_in_bytes, int bitshift);

#endif

