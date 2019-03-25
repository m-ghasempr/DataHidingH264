
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

extern const byte QP_SCALE_CR[52] ;
//! look up tables for FRExt_chroma support
extern const unsigned char subblk_offset_x[3][8][4];
extern const unsigned char subblk_offset_y[3][8][4];

extern void iMBtrans4x4(ImageParameters *img, Macroblock *currMB, ColorPlane pl, int smb);
extern void iMBtrans8x8(ImageParameters *img, Macroblock *currMB, ColorPlane pl);

extern void itrans_sp_cr(ImageParameters *img, int uv);

extern int  intrapred_luma_16x16  (ImageParameters *img, Macroblock *currMB, ColorPlane pl, int predmode);
extern void intrapred_chroma      (ImageParameters *img, Macroblock *currMB, int uv);

void (*itrans_4x4)(ImageParameters *img, ColorPlane pl, int ioff, int joff);
void (*itrans_8x8)(ImageParameters *img, Macroblock *currMB, ColorPlane pl, int ioff, int joff);

extern void Inv_Residual_trans_4x4(ImageParameters *img, ColorPlane pl, int ioff, int joff);
extern void Inv_Residual_trans_8x8(ImageParameters *img, Macroblock *currMB, ColorPlane pl, int ioff,int joff);

extern void itrans8x8   (ImageParameters *img, Macroblock *currMB, ColorPlane pl, int ioff, int joff);
extern void itrans4x4   (ImageParameters *img, ColorPlane pl, int ioff, int joff);
extern void itrans4x4_ls(ImageParameters *img, ColorPlane pl, int ioff, int joff);
extern void itrans_sp   (ImageParameters *img, ColorPlane pl, int ioff, int joff);
extern int  intrapred   (ImageParameters *img, Macroblock *currMB, ColorPlane pl, int ioff,int joff,int i4,int j4);
extern void itrans_2    (ImageParameters *img, Macroblock *currMB, ColorPlane pl);
extern void iTransform  (ImageParameters *img, Macroblock *currMB, ColorPlane pl, int need_4x4_transform, int smb);

extern int  allocate_block_mem(void);
extern void free_block_mem(void);

#endif

