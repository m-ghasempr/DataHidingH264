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

extern int testEndian(void);
extern void initInput(ImageParameters *p_Img, FrameFormat *source, FrameFormat *output);
extern void AllocateFrameMemory (ImageParameters *p_Img, InputParameters *p_Inp, FrameFormat *source);
extern void DeleteFrameMemory (ImageParameters *p_Img);

extern int ReadOneFrame (ImageParameters *p_Img, InputParameters *p_Inp, VideoDataFile *input_file, int FrameNoInFile, int HeaderSize, FrameFormat *source, FrameFormat *output, imgpel **pImage[3]);
extern void PaddAutoCropBorders( FrameFormat output, int img_size_x, int img_size_y, int img_size_x_cr, int img_size_y_cr, imgpel **pImage[3]);

#endif

