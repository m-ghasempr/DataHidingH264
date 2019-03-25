
/*!
 ************************************************************************
 * \file block.h
 *
 * \brief
 *    definitions for block decoding functions
 *
 * \author
 *  Inge Lille-Langoy               <inge.lille-langoy@telenor.com>    \n
 *  Telenor Satellite Services                                         \n
 *  P.O.Box 6914 St.Olavs plass                                        \n
 *  N-0130 Oslo, Norway
 *
 ************************************************************************
 */

#ifndef _BLOCK_H_
#define _BLOCK_H_

#include "global.h"

#define DQ_BITS         6

extern const byte QP_SCALE_CR[52] ;
extern const int  dequant_coef[6][4][4];
extern const unsigned char subblk_offset_x[3][8][4];
extern const unsigned char subblk_offset_y[3][8][4];
extern void iMBtrans4x4(ColorPlane pl, struct img_par *img, int smb);
extern void iMBtrans8x8(ColorPlane pl, struct img_par *img);

extern void itrans_sp_cr(struct img_par *img, int uv);
extern void itrans4x4(struct img_par *img,int ioff,int joff, int yuv);
extern void itrans_sp(struct img_par *img,int ioff,int joff,int i0,int j0);
extern int  intrapred(Macroblock *currMB, ColorPlane pl, struct img_par *img,int ioff,int joff,int i4,int j4);
extern void itrans_2 (ColorPlane pl, struct img_par *img);
extern void iTransform(Macroblock *currMB, ColorPlane pl, struct img_par *img, int need_4x4_transform, int smb, int yuv);
#endif

