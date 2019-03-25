/*!
 ***************************************************************************
 *
 * \file intra8x8.h
 *
 * \brief
*    prototypes of 8x8 transform functions
  *
 * \date
 *    9. October 2003
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    Yuri Vatis
 *    Alexis Michael Tourapis
 **************************************************************************/

#ifndef _INTRA8X8_H_
#define _INTRA8X8_H_

extern void set_intrapred_8x8 (Macroblock *currMB, ColorPlane pl, int img_x,int img_y, int *left_available, int *up_available, int *up_left_available);
extern void get_intrapred_8x8 (Macroblock *currMB, ColorPlane pl, int i8x8_mode, int left_available, int up_available);
extern void LowPassForIntra8x8Pred(imgpel *PredPel, int block_up_left, int block_up, int block_left);

#endif //_TRANSFORM8X8_H_
