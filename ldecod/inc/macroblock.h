/*!
 ************************************************************************
 * \file macroblock.h
 *
 * \brief
 *    Arrays for macroblock encoding
 *
 * \author
 *    Inge Lille-Langoy               <inge.lille-langoy@telenor.com>
 *    Copyright (C) 1999 Telenor Satellite Services, Norway
 ************************************************************************
 */

#ifndef _MACROBLOCK_H_
#define _MACROBLOCK_H_

#include "block.h"

extern const byte QP_SCALE_CR[52];

extern void set_interpret_mb_mode(int slice_type);
extern void get_neighbors(Macroblock *currMB, 
                          PixelPos *block_a, PixelPos *block_b, PixelPos *block_c, PixelPos *block_d,
                          int mb_x, int mb_y, int blockshape_x);
extern void GetMotionVectorPredictor (Macroblock *currMB, 
                               PixelPos *block_a, PixelPos *block_b, PixelPos *block_c, 
                               short  pmv[2], char  ref_frame, char   **refPic, short  ***tmp_mv,                               
                               int mb_x, int mb_y, int blockshape_x, int blockshape_y);

extern void start_macroblock     (ImageParameters *img, Macroblock **currMB);
extern void read_one_macroblock  (ImageParameters *img, Slice *currSlice, Macroblock *currMB);
extern int  decode_one_macroblock(ImageParameters *img, Macroblock *currMB, StorablePicture *dec_picture);
extern Boolean  exit_macroblock  (ImageParameters *img, int eos_bit);
#endif

