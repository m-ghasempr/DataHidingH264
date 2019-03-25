/*!
 ***************************************************************************
 * \file
 *    q_around.h
 *
 * \author
 *    Alexis Michael Tourapis
 *
 * \brief
 *    Headerfile for Quantization Adaptive Rounding
 **************************************************************************
 */

#ifndef _Q_AROUND_H_
#define _Q_AROUND_H_

void update_adaptive_rounding_8x8(RD_8x8DATA* dataTr, int**** ARCofAdj);
void store_adaptive_rounding_4x4 (ImageParameters *img, int****ARCofAdj, int mode, int block_y, int block_x);
void update_adaptive_rounding_4x4 (ImageParameters *img, int****ARCofAdj , int mode, int block_y, int block_x);
void update_adaptive_rounding_16x16(ImageParameters *img, int****ARCofAdj , int mode);
void store_adaptive_rounding_16x16 (ImageParameters *img, int****ARCofAdj, int mode);
void reset_adaptive_rounding_direct();

void update_offset_params    (Macroblock *currMB, int mode, int luma_transform_size_8x8_flag);

#endif

