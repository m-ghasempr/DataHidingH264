/*
***********************************************************************
*  COPYRIGHT  AND  WARRANTY INFORMATION
*
* Copyright 2001, International Telecommunications Union, Geneva
*
* DISCLAIMER OF WARRANTY
*
* These software programs are available to the user without any
* license fee or royalty on an "as is" basis. The ITU disclaims
* any and all warranties, whether express, implied, or
* statutory, including any implied warranties of merchantability
* or of fitness for a particular purpose.  In no event shall the
* contributor or the ITU be liable for any incidental, punitive, or
* consequential damages of any kind whatsoever arising from the
* use of these programs.
*
* This disclaimer of warranty extends to the user of these programs
* and user's customers, employees, agents, transferees, successors,
* and assigns.
*
* The ITU does not represent or warrant that the programs furnished
* hereunder are free of infringement of any third-party patents.
* Commercial implementations of ITU-T Recommendations, including
* shareware, may be subject to royalty fees to patent holders.
* Information regarding the ITU-T patent policy is available from
* the ITU Web site at http://www.itu.int.
*
* THIS IS NOT A GRANT OF PATENT RIGHTS - SEE THE ITU-T PATENT POLICY.
************************************************************************
*/

/*
 *************************************************************************************
 * \file
 *    abt.c
 *
 * \brief
 *    Description
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *     -  Mathias Wien        <wien@ient.rwth-aachen.de>
 *     -  Achim Dahlhoff      <dahlhoff@ient.rwth-aachen.de>
 *
 * \date
 *    Fri Mar 8 2002
 *
 *  copyright : (C) 2002 by Mathias Wien
 *                            Institut und Lehrstuhl für Nachrichtentechnik
 *                            RWTH Aachen University
 *                            52072 Aachen
 *                            Germany
 *************************************************************************************
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include "defines.h"
#include "global.h"
#include "abt.h"
#include "elements.h"

/*!
 ************************************************************************
 * \brief
 *    Transform 8x8 block using ABT. Is also used for transforming
 *    single sub-blocks. ATTENTION input curr_blk has to be [lines][pixels]
 ************************************************************************
 */
void transform_ABT_B8(int abt_mode,                   // block tiling mode used for ABT
                      int blk_off_x, int blk_off_y,   // block offset. WHOLE_BLK for transforming all subblocks
                      int curr_blk[B8_SIZE][B8_SIZE]  // block to be transformed.
                      )
{
  int xx, yy, kk, ox, oy, blk_x, blk_y;
  int tmp, val[B8_SIZE][B8_SIZE];
  int blk_size_x, blk_size_y;     // transform block size
  int num_blk_x, num_blk_y;       // number of transformed block in x,y direction
  int tr_idx_x, tr_idx_y;
  int bshift[2][2];
  int rshift[2][2];

  assert(abt_mode<=B4x4);

  for   (yy=0;yy<2;yy++)
    for (xx=0;xx<2;xx++)
    {
      bshift[yy][xx] = ABT_SHIFT0[abt_mode][yy][xx];
      rshift[yy][xx] = (bshift[yy][xx]>0) ? (1<<(bshift[yy][xx]-1)) : 0;
    }

  blk_size_x = ABT_TRSIZE[abt_mode][PIX];
  blk_size_y = ABT_TRSIZE[abt_mode][LIN];
  num_blk_x  = num_blk_y  = 1;
  tr_idx_x   = ABT_TRIDX[abt_mode][PIX];
  tr_idx_y   = ABT_TRIDX[abt_mode][LIN];

  if (blk_off_x == WHOLE_BLK)
  {
    assert(blk_off_y == WHOLE_BLK); // just to be shure, may be removed later
    blk_off_x = 0;
    blk_off_y = 0;
    num_blk_x  = ABT_NUMTR[abt_mode][PIX];
    num_blk_y  = ABT_NUMTR[abt_mode][LIN];
  }
  else
  {
    assert(blk_off_x+blk_size_x <= B8_SIZE);
    assert(blk_off_y+blk_size_y <= B8_SIZE);
  }

  for (blk_y=0; blk_y<num_blk_y; blk_y++)
  {
    oy = blk_off_y + blk_y*blk_size_y;
    for (blk_x=0; blk_x<num_blk_x; blk_x++)
    {
      ox = blk_off_x + blk_x*blk_size_x;

      // horizontal transform
      for(yy=0; yy<blk_size_y; yy++)
        for(xx=0; xx<blk_size_x; xx++)
        {
          val[yy][xx] = 0;
          for (kk=0; kk<blk_size_x; kk++)
            val[yy][xx] += curr_blk[oy+yy][ox+kk] * ABT_matrix[tr_idx_x][xx][kk];
        }

      // vertical transform
      for(yy=0; yy<blk_size_y; yy++)
        for(xx=0; xx<blk_size_x; xx++)
        {
          tmp = 0;
          for (kk=0; kk<blk_size_y; kk++)
            tmp += val[kk][xx] * ABT_matrix[tr_idx_y][yy][kk];

          curr_blk[oy+yy][ox+xx] = sign( ((absm(tmp) + rshift[yy&1][xx&1])>> bshift[yy&1][xx&1]),tmp);
        }
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Inverse transform 8x8 block using ABT. Is also used for transforming
 *    single sub-blocks. ATTENTION input curr_blk has to be [lines][pixels]
 ************************************************************************
 */
void inv_transform_ABT_B8(int abt_mode,                   // block tiling mode used for ABT
                          int blk_off_x, int blk_off_y,   // block offset. WHOLE_BLK for transforming all subblocks
                          int curr_blk[B8_SIZE][B8_SIZE]  // block to be inverse transformed.
                          )
{
  int xx, yy, kk, ox, oy, blk_x, blk_y;
  int tmp, val[B8_SIZE][B8_SIZE];
  int blk_size_x, blk_size_y;     // transform block size
  int num_blk_x, num_blk_y;       // number of transformed block in x,y direction
  int tr_idx_x, tr_idx_y;
  int bshift;
  int rshift;
  int th4scale[2] = {0,0}; // used to emulate the JM 4x4 inv transform
  int tv4scale[2] = {0,0}; // used to emulate the JM 4x4 inv transform

  assert(abt_mode<=B4x4);

  bshift = ABT_SHIFT1[abt_mode];
  rshift = (bshift>0) ? (1<<(bshift-1)) : 0;

  blk_size_x = ABT_TRSIZE[abt_mode][PIX];
  blk_size_y = ABT_TRSIZE[abt_mode][LIN];
  num_blk_x  = num_blk_y  = 1;
  tr_idx_x   = ABT_TRIDX[abt_mode][PIX];
  tr_idx_y   = ABT_TRIDX[abt_mode][LIN];

  if (blk_off_x == WHOLE_BLK)
  {
    assert(blk_off_y == WHOLE_BLK); // just to be shure, may be removed later
    blk_off_x = 0;
    blk_off_y = 0;
    num_blk_x  = ABT_NUMTR[abt_mode][PIX];
    num_blk_y  = ABT_NUMTR[abt_mode][LIN];
  }
  else
  {
    assert(blk_off_x+blk_size_x <= B8_SIZE);
    assert(blk_off_y+blk_size_y <= B8_SIZE);
  }

  if (blk_size_x == 4)
    th4scale[1] = 1;
  if (blk_size_y == 4)
    tv4scale[1] = 1;

  for (blk_y=0; blk_y<num_blk_y; blk_y++)
  {
    oy = blk_off_y + blk_y*blk_size_y;
    for (blk_x=0; blk_x<num_blk_x; blk_x++)
    {
      ox = blk_off_x + blk_x*blk_size_x;

      // horizontal inverse transform
      for(yy=0; yy<blk_size_y; yy++)
        for(xx=0; xx<blk_size_x; xx++)
        {
          val[yy][xx] = 0;
          for (kk=0; kk<blk_size_x; kk++)
            val[yy][xx] += (curr_blk[oy+yy][ox+kk] * ABT_matrix[tr_idx_x][kk][xx]) >> th4scale[kk&1];

          val[yy][xx] = sign( (absm(val[yy][xx]) + rshift) >> bshift,val[yy][xx]);
        }

      // vertical inverse transform
      for(yy=0; yy<blk_size_y; yy++)
        for(xx=0; xx<blk_size_x; xx++)
        {
          tmp = 0;
          for (kk=0; kk<blk_size_y; kk++)
            tmp += (val[kk][xx] * ABT_matrix[tr_idx_y][kk][yy]) >> tv4scale[kk&1];

          curr_blk[oy+yy][ox+xx] = tmp; // is reduced to 16bit in inv quantization
        }
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    Quantization of a transformed ABT block
 ************************************************************************
 */
void quant_abt_B8(int qp,                         // Quantization parameter
                  int abt_mode,                   // block tiling mode used for ABT. USED ALSO FOR SIGNALING INTRA (abt_mode+4)
                  int blk_off_x, int blk_off_y,   // block offset. WHOLE_BLK for quantizing all subblocks
                  int curr_blk[B8_SIZE][B8_SIZE]
                  )
{
  int Q[2][2];
  int xx, yy;
  int blk_size_x, blk_size_y;     // transform block size
  int val, tmp;
  int qp_const[2][2];
  int q_bits[2][2];
  int intra   = 0;
  if (abt_mode>B4x4) // abt_mode 0..3 inter, abt_mode 4..7 intra
  {
    intra     = 1;
    abt_mode -= 4;
  }

  // size of block to be quantized
  if (blk_off_x == WHOLE_BLK)
  {
    assert(blk_off_y == WHOLE_BLK); // just to be shure, may be removed later
    blk_off_x  = 0;
    blk_off_y  = 0;
    blk_size_x = B8_SIZE;
    blk_size_y = B8_SIZE;
  }
  else
  {
    blk_size_x = ABT_TRSIZE[abt_mode][PIX];
    blk_size_y = ABT_TRSIZE[abt_mode][LIN];
    assert(blk_off_x+blk_size_x <= B8_SIZE);
    assert(blk_off_y+blk_size_y <= B8_SIZE);
  }

  // get quantization values.
  // It is assumed that blk_off_x,blk_off_y are allways even.
  // Needed for proper quantization table indexing. assert?
  for   (yy=0;yy<2;yy++)
    for (xx=0;xx<2;xx++)
    {
      Q[yy][xx]      = ABT_Q[abt_mode][qp%QUANT_PERIOD][ABT_QMAP[abt_mode][yy][xx]];
      q_bits[yy][xx] = ABT_N[abt_mode][0] -ABT_SHIFT0[abt_mode][yy][xx] + qp/QUANT_PERIOD;
      if (intra)
        qp_const[yy][xx]  = (1<<q_bits[yy][xx])/3;      // intra
      else
        qp_const[yy][xx]  = (1<<q_bits[yy][xx])/6;      // inter
    }


  // quantize
  for   (yy=blk_off_y; yy<blk_off_y+blk_size_y; yy++)
    for (xx=blk_off_x; xx<blk_off_x+blk_size_x; xx++)
    {
      val              = curr_blk[yy][xx];
      tmp              = absm(val);
      curr_blk[yy][xx] = sign( ((tmp * Q[yy&1][xx&1] + qp_const[yy&1][xx&1]) >> q_bits[yy&1][xx&1]) , val);
    }
}


/*!
 ************************************************************************
 * \brief
 *    Dequantization and inverse Transform for an ABT block
 ************************************************************************
 */
void idct_dequant_abt_B8(int block8x8,
                         int qp,                         // Quantization parameter
                         int abt_mode,                   // block tiling mode used for ABT.
                         int blk_off_x, int blk_off_y,   // block offset. WHOLE_BLK for quantizing all subblocks
                         int curr_blk[B8_SIZE][B8_SIZE],
                         struct img_par *img
                         )
{
  int  xx,yy,x1,y1;
  int  blk_size_x, blk_size_y;     // transform block size
  int  boff_x, boff_y;
  int  q_bits, q_shift;
  int  qp_const;
  int  b8_y       = (block8x8 / 2) << 3;
  int  b8_x       = (block8x8 % 2) << 3;
  int  curr_val;

  assert(USEABT);   // just to be shure.
  // size of block to be quantized
  if (blk_off_x == WHOLE_BLK)
  {
    assert(blk_off_y == WHOLE_BLK); // just to be shure.
    boff_x  = 0;
    boff_y  = 0;
    blk_size_x = B8_SIZE;
    blk_size_y = B8_SIZE;
  }
  else
  {
    boff_x  = blk_off_x;
    boff_y  = blk_off_y;
    blk_size_x = ABT_TRSIZE[abt_mode][PIX];
    blk_size_y = ABT_TRSIZE[abt_mode][LIN];
    assert(blk_off_x+blk_size_x <= B8_SIZE);
    assert(blk_off_y+blk_size_y <= B8_SIZE);
  }

  // dequantization values.
  //  q_bits   = ABT_N[abt_mode][1] - ABT_SHIFT1[abt_mode] - qp/QUANT_PERIOD;
  q_bits   = ABT_N[abt_mode][1] - ABT_SHIFT1[abt_mode];
  q_shift  = qp/QUANT_PERIOD;
  qp_const = 1<<(q_bits-1);

  x1=boff_x+(blk_off_x==WHOLE_BLK? B8_SIZE : blk_size_x);
  y1=boff_y+(blk_off_y==WHOLE_BLK? B8_SIZE : blk_size_y);

  // inverse transform
  inv_transform_ABT_B8(abt_mode, blk_off_x, blk_off_y, curr_blk);

  // normalize
  for(yy=boff_y;yy<y1;yy++)
    for(xx=boff_x;xx<x1;xx++)
    {
      curr_val = ((img->mpr[b8_x+xx][b8_y+yy]<<q_bits) + curr_blk[yy][xx] + qp_const)>>q_bits;

      img->m7[xx][yy] = curr_blk[yy][xx] =
      imgY[img->pix_y+b8_y+yy][img->pix_x+b8_x+xx] = clamp(curr_val,0,255);
    }
}

/*!
 ************************************************************************
 * \brief
 *    SP-Frames: Dequantization and inverse Transform for an ABT block
 ************************************************************************
 */
void idct_dequant_abt_sp(int block8x8,
                         int abt_mode,                   // block tiling mode used for ABT.
                         int blk_off_x, int blk_off_y,   // block offset. WHOLE_BLK for quantizing all subblocks
                         int curr_blk[B8_SIZE][B8_SIZE],
                         struct img_par *img
                         )
{
  int  xx,yy,x1,y1;
  int  blk_size_x, blk_size_y;     // transform block size
  int  boff_x, boff_y;
  int  qp, dq_bits, q_shift, dqp_const;
  int  q_bits[2][2];
  int  qp_const[2][2];
  int  Q[2][2], R[2][2];
  int  qpsp, q_shift_sp;
  int  q_bits_sp[2][2];
  int  qp_const_sp[2][2];
  int  Qsp[2][2], Rsp[2][2];
  int  vq_bits, vqp_const, vq_bits_sp, vqp_const_sp, vQ, vQsp;
  int  pred_blk[B8_SIZE][B8_SIZE];
  int  b8_y          = (block8x8 / 2) << 3;
  int  b8_x          = (block8x8 % 2) << 3;
  int  curr_val, tmp;

  assert(USEABT);   // just to be shure.
  // size of block to be quantized
  if (blk_off_x == WHOLE_BLK)
  {
    assert(blk_off_y == WHOLE_BLK); // just to be shure.
    boff_x  = 0;
    boff_y  = 0;
    blk_size_x = B8_SIZE;
    blk_size_y = B8_SIZE;
  }
  else
  {
    boff_x  = blk_off_x;
    boff_y  = blk_off_y;
    blk_size_x = ABT_TRSIZE[abt_mode][PIX];
    blk_size_y = ABT_TRSIZE[abt_mode][LIN];
    assert(blk_off_x+blk_size_x <= B8_SIZE);
    assert(blk_off_y+blk_size_y <= B8_SIZE);
  }

  // quantization values.
  qp        = img->qp+QP_OFS-MIN_QP;
  q_shift   = qp/QUANT_PERIOD;
  dq_bits   = ABT_N[abt_mode][1] - ABT_SHIFT1[abt_mode];
  dqp_const = 1<<(dq_bits-1);

  qpsp        = img->qpsp+QP_OFS-MIN_QP;
  q_shift_sp  = qpsp/QUANT_PERIOD;

  get_quant_consts(abt_mode,qp,2,Q,q_bits,qp_const);
  get_dequant_consts(abt_mode,qp,R);
  get_quant_consts(abt_mode,qpsp,2,Qsp,q_bits_sp,qp_const_sp);
  get_dequant_consts(abt_mode,qpsp,Rsp);

  if (img->sp_switch)
  {
    qp      = qpsp;
    q_shift = qp/QUANT_PERIOD;
    get_quant_consts(abt_mode,qp,2,Q,q_bits,qp_const);
    get_dequant_consts(abt_mode,qp,R);
  }

  x1=boff_x+(blk_off_x==WHOLE_BLK? B8_SIZE : blk_size_x);
  y1=boff_y+(blk_off_y==WHOLE_BLK? B8_SIZE : blk_size_y);


  // get prediction block
  for(yy=boff_y;yy<y1;yy++)
    for(xx=boff_x;xx<x1;xx++)
      pred_blk[yy][xx] = img->mpr[b8_x+xx][b8_y+yy];

  // transform prediction block
  transform_ABT_B8(abt_mode, blk_off_x, blk_off_y, pred_blk);

  for(yy=boff_y;yy<y1;yy++)
    for(xx=boff_x;xx<x1;xx++)
    {
      // quantization values:
      vqp_const    = qp_const   [yy&1][xx&1];
      vqp_const_sp = qp_const_sp[yy&1][xx&1];
      vq_bits      = q_bits     [yy&1][xx&1];
      vq_bits_sp   = q_bits_sp  [yy&1][xx&1];
      vQ           = Q  [yy&1][xx&1];
      vQsp         = Qsp[yy&1][xx&1];

      // recovering coefficient since they are already dequantized earlier
      tmp      = (curr_blk[yy][xx]>>q_shift) / R[yy&1][xx&1];
      curr_val = sign((( (abs(tmp)<<vq_bits) + vQ/2) / vQ),tmp) + pred_blk[yy][xx];
      curr_blk[yy][xx] = sign(( abs(curr_val)*vQsp + vqp_const_sp)>>vq_bits_sp, curr_val) * Rsp[yy&1][xx&1] << q_shift_sp;
    }

  // inverse transform
  inv_transform_ABT_B8(abt_mode, blk_off_x, blk_off_y, curr_blk);


  // clamp and copy to imgY
  for(yy=boff_y;yy<y1;yy++)
    for(xx=boff_x;xx<x1;xx++)
    {
      curr_val = (curr_blk[yy][xx] + dqp_const) >> dq_bits;
      img->m7[xx][yy] =
        curr_blk[yy][xx] =
        imgY[img->pix_y+b8_y+yy][img->pix_x+b8_x+xx] = clamp(curr_val,0,255);
    }
}







/*!
 ************************************************************************
 * \brief
 *    Read ABT coefficients of one 8x8 block
 ************************************************************************
*/
void readLumaCoeffABT_B8(int block8x8, struct inp_par *inp, struct img_par *img)
{
  int mb_nr          = img->current_mb_nr;
  Macroblock *currMB = &img->mb_data[mb_nr];
  const int cbp      = currMB->cbp;
  Slice *currSlice   = img->currentSlice;
  DataPartition *dP;
  int *partMap       = assignSE2partition[inp->partition_mode];
  SyntaxElement currSE;
  unsigned int is_last_levelrun;
  int
    intra,
    abt_mode,
    inumblk,                    /* number of blocks per CBP*/
    inumcoeff,                  /* number of coeffs per block */
    iblk,                       /* current block */
    icoef,                      /* current coefficient */
    ipos,
    run, level,
    len,        // obsolete?
    ii,jj,
    sbx, sby;
  int boff_x, boff_y;
  int any_coeff;
  int vlc_numcoef;
  int  cbp_blk_mask;

  int symbol2D,mode2D,Golomb_se_type;
  const char (*table2D)[8];
  char (*table2D_dec)[2];

  int qp, q_shift, R[2][2]; // dequantization parameters
//  static const int blk_table[4][4] = { {0,-1,-1,-1}, {0,2,-1,-1}, {0,1,-1,-1}, {0,1,2,3}};
  static const int blkmode2ctx [4] = {LUMA_8x8, LUMA_8x4, LUMA_4x8, LUMA_4x4};
  int  frm_fld; // indicate progressive/interlaced scan use.

#ifdef _ALT_SCAN_
  frm_fld = (img->structure>FRAME); // img->structure==FRAME: progressive, else interlaced
#else
  frm_fld = 0;
#endif

  // this has to be done for each subblock seperately
  intra     = (currMB->b8mode[block8x8]==IBLOCK)||IS_INTRA(currMB);
  abt_mode  = currMB->abt_mode[block8x8];
  inumblk   = ABT_NUMTR [abt_mode][0] * ABT_NUMTR [abt_mode][1];
  inumcoeff = ABT_TRSIZE[abt_mode][0] * ABT_TRSIZE[abt_mode][1] + 1; //  all positions + EOB

  assert(USEABT);   // just to be shure.
  // ========= dequantization values ===========
  qp       = currMB->qp+QP_OFS-MIN_QP; // using old style qp in ABT routines.
  q_shift  = qp/QUANT_PERIOD;
  get_dequant_consts(abt_mode,qp,R);
  // ============================================

  //make decoder table for 2D code
  if(ABT_2D_VLC_dec[0][0][1]<0)                                                          // Don't need to set this every time. rewrite later.
  {
    memset(ABT_2D_VLC_dec,-1,sizeof(ABT_2D_VLC_dec));
    for(iblk=0;iblk<4;iblk++)
    {
      table2D=ABT_2D_VLC[iblk];
      for(run=0;run<16;run++)
        for(level=0;level<8;level++)
        {
          ipos=table2D[run][level];
          assert(ipos<64);
          if(ipos>=0)
          {
            ABT_2D_VLC_dec[iblk][ipos][0]=level;
            ABT_2D_VLC_dec[iblk][ipos][1]=run;
            if(level)
            {
              ABT_2D_VLC_dec[iblk][ipos+1][0]=-level;
              ABT_2D_VLC_dec[iblk][ipos+1][1]=run;
            }
          }
        }
    }
    assert(ABT_2D_VLC_dec[0][0][1]>=0);        //otherwise, tables are bad.
  }

  if( !intra ) // bug fix mwi 020505
  {
    mode2D=2;
    ii=0;        //first table is the inter table.
  }else{
    mode2D=2;
    if(abt_mode>B8x8)mode2D--;
    if(abt_mode>=B4x4)mode2D--;
    ii=(qp-6)>>3;        //find table to use for intra ( 6..13 -> 0 , 14..21 -> 1 , 22..29 -> 2 )
    if(ii<0)ii=0;
    if(ii>2)ii=2;
    ii+=1;        //first table is the inter table. Intra tables start off 1.
  }
  table2D    =ABT_2D_VLC[ii];
  table2D_dec=ABT_2D_VLC_dec[ii];

  //clear cbp_blk bits of thie 8x8 block (and not all 4!)
  cbp_blk_mask = cbp_blk_masks[0][0] ;
  if(block8x8&1)cbp_blk_mask<<=2;
  if(block8x8&2)cbp_blk_mask<<=8;
  currMB->cbp_blk&=~cbp_blk_mask;

  vlc_numcoef=-1;

  if (inp->symbol_mode == UVLC)
  {
    Golomb_se_type=SE_LUM_AC_INTER;
    if( intra )
    {
      vlc_numcoef=0;        //this means 'use numcoeffs symbol'.
      Golomb_se_type=SE_LUM_AC_INTRA;
    }
  }

  // === decoding ===
  if ( cbp & (1<<block8x8) )
  {
    for (iblk=0; iblk<inumblk; iblk++)
    {
      // === set offset in current macroblock ===
      boff_x = ( (block8x8%2)<<3 ) + (iblk%ABT_NUMTR[abt_mode][0])*ABT_TRSIZE[abt_mode][0];
      boff_y = ( (block8x8/2)<<3 ) + (iblk/ABT_NUMTR[abt_mode][0])*ABT_TRSIZE[abt_mode][1];
      img->subblock_x = boff_x>>2;
      img->subblock_y = boff_y>>2;

      ipos  = -1;
      any_coeff=0;

      if( vlc_numcoef>=0 && inp->symbol_mode==UVLC )
      {
        //get numcoeffs symbol
        currSE.type=Golomb_se_type;
        currSE.golomb_grad=0;
        if( intra )
          currSE.golomb_grad=2;
        currSE.golomb_maxlevels=31;
        dP = &(currSlice->partArr[partMap[currSE.type]]);
        dP->readSyntaxElement(&currSE,img,inp,dP);
        vlc_numcoef=currSE.value1;
      }

      is_last_levelrun=0;

      for (icoef=0; icoef<inumcoeff; icoef++)        //will break if level is zero or if EOB was sent
      {
        if (inp->symbol_mode == UVLC)
        {
          //decode from VLC
          if( vlc_numcoef>=0 && icoef>=vlc_numcoef )
            break;

          //read 2D symbol
          currSE.type=Golomb_se_type;
          currSE.golomb_grad=mode2D;
          currSE.golomb_maxlevels=6-mode2D;        //make code with little less than 64 symbols.
          dP = &(currSlice->partArr[partMap[currSE.type]]);
          dP->readSyntaxElement(&currSE,img,inp,dP);
          symbol2D=currSE.value1;

          if(vlc_numcoef>=0)
            symbol2D++;

          if(symbol2D<CODE2D_ESCAPE_SYMBOL)
          {
            level=table2D_dec[symbol2D][0];
            run=table2D_dec[symbol2D][1];
            if(!level)break;

          }else{        //if out-of-table symbol, do normal read.

            //decode level
            currSE.type=Golomb_se_type;
            currSE.golomb_grad=3;
            currSE.golomb_maxlevels=28;
            dP = &(currSlice->partArr[partMap[currSE.type]]);
            dP->readSyntaxElement(&currSE,img,inp,dP);

            level=(currSE.value1>>1)+1;
            if( currSE.value1 & 1 )                //sign bit
              level=-level;

            //decode run
            currSE.type=Golomb_se_type;
            currSE.golomb_grad=2;
            currSE.golomb_maxlevels=29;
            dP = &(currSlice->partArr[partMap[currSE.type]]);
            dP->readSyntaxElement(&currSE,img,inp,dP);
            run=currSE.value1;

          }
        }else{
          //decode from CABAC
          currSE.reading = readRunLevelFromBuffer_CABAC;
          currSE.type         = (intra ? (icoef==0 ? SE_LUM_DC_INTRA : SE_LUM_AC_INTRA) : (icoef==0 ? SE_LUM_DC_INTER : SE_LUM_AC_INTER));
          currSE.context      = blkmode2ctx[abt_mode];
          img->is_intra_block = intra;

#if TRACE
          sprintf(currSE.tracestring, " Luma sng ");
#endif
          dP = &(currSlice->partArr[partMap[currSE.type]]);
          dP->readSyntaxElement(&currSE,img,inp,dP);
          level = currSE.value1;
          run =  currSE.value2;
          len = currSE.len;
        }

        if(!level)
          break;
        any_coeff=1;
        ipos += (run+1);

        assert(ipos<inumcoeff);

        ii = ABT_SCAN[frm_fld][abt_mode][ipos][0];
        jj = ABT_SCAN[frm_fld][abt_mode][ipos][1];

        img->m7[boff_x + ii][boff_y + jj] = level*R[jj&1][ii&1]<<q_shift;

        if(is_last_levelrun)
          break;
      }//end loop icoef

      // CAVLC store number of coefficients
      sbx = img->subblock_x;
      sby = img->subblock_y;
      switch (abt_mode)
      {
      case 0: // 8x8
        img->nz_coeff [img->mb_x ][img->mb_y ][sbx  ][sby  ] = icoef/4;
        img->nz_coeff [img->mb_x ][img->mb_y ][sbx+1][sby  ] = icoef/4;
        img->nz_coeff [img->mb_x ][img->mb_y ][sbx  ][sby+1] = icoef/4;
        img->nz_coeff [img->mb_x ][img->mb_y ][sbx+1][sby+1] = icoef/4;
        break;
      case 1: // 8x4
        img->nz_coeff [img->mb_x ][img->mb_y ][sbx  ][sby] = icoef/2;
        img->nz_coeff [img->mb_x ][img->mb_y ][sbx+1][sby] = icoef/2;
        break;
      case 2: // 4x8
        img->nz_coeff [img->mb_x ][img->mb_y ][sbx][sby  ] = icoef/2;
        img->nz_coeff [img->mb_x ][img->mb_y ][sbx][sby+1] = icoef/2;
        break;
      case 3: // 4x4
        img->nz_coeff [img->mb_x ][img->mb_y ][sbx][sby] = icoef;
        break;
      }


      if(any_coeff)
      {
        cbp_blk_mask = cbp_blk_masks[ abt_mode&3 ][ ((boff_y&4)>>1)|((boff_x&4)>>2) ] ;
        if(block8x8&1)cbp_blk_mask<<=2;
        if(block8x8&2)cbp_blk_mask<<=8;

        currMB->cbp_blk |= cbp_blk_mask;
      }
    }//end loop iblk
  }//end if ( cbp & (1<<icbp) )
}




/*!
 ************************************************************************
 * \brief
 *    Save block mode for the current 8x8 block. Writes to colB8mode.
 ************************************************************************
*/
void setDirectModeABT(int block8x8, struct img_par *img)
{
  int dirmode;
  int blk_x          = (img->pix_x>>3) + block8x8%2;
  int blk_y          = (img->pix_y>>3) + block8x8/2;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int pstruct        = img->structure; // distinguish frame/top/bottom array. 020603 mwi

  assert(pstruct<3); // needs revision for other values. 020603 mwi

  if      (currMB->mb_type==I16MB)
  {
    dirmode = I16MB;
  }
  else if (currMB->mb_type==0)
  {
    dirmode = B8x8+4; // 8x8 MB mode
  }
  else if( currMB->mb_type==I4MB || (currMB->mb_type==P8x8&&currMB->b8mode[block8x8]==IBLOCK) )
    dirmode=4+(currMB->abt_mode[block8x8]&3);
  else
    dirmode = currMB->b8mode[block8x8];

  colB8mode[pstruct][blk_y][blk_x] = dirmode;
};





/*!
 ************************************************************************
 * \brief
 *    Get Direct Mode block tiling for the current B-Frame 8x8 block.
 *    Writes to currMB->abt_mode,abt_pred_mode and returns abt_mode.
 ************************************************************************
*/
int getDirectModeABT(int block8x8, struct img_par *img)
{
  int dirmode;
  int blk_x          = (img->pix_x>>3) + block8x8%2;
  int blk_y          = (img->pix_y>>3) + block8x8/2;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  static const int  abt_mode_table[12] = {-2, 0, 0, 0, 0, 1, 2, 3, -2, 3, 0, 3};  // DO NOT CHANGE ORDER !!!   // b8mode->abt_mode
  int pstruct        = img->structure; // distinguish frame/top/bottom array. 020603 mwi

  assert(pstruct<3); // needs revision for other values. 020603 mwi

  dirmode = abt_mode_table[colB8mode[pstruct][blk_y][blk_x]];

  assert (dirmode != -2); // just to be sure. P8x8 or '0' are not valid for an ABT 8x8 block mode.

  return (currMB->abt_mode[block8x8] = currMB->abt_pred_mode[block8x8] = dirmode);
};




/*!
 ************************************************************************
 * \brief
 *    Copy region of img->m7 corresponding to block8x8 to curr_blk[][].
 *    Attention: img->m7 is [x][y] and curr_blk is [y][x].
 ************************************************************************
*/
void get_curr_blk( int block8x8,struct img_par *img, int curr_blk[B8_SIZE][B8_SIZE])
{
  int  xx, yy;
  int  mb_y       = (block8x8 / 2) << 3;
  int  mb_x       = (block8x8 % 2) << 3;

  assert(USEABT);   // just to be shure.
  for   (yy=0; yy<B8_SIZE; yy++)
    for (xx=0; xx<B8_SIZE; xx++)
      curr_blk[yy][xx] = img->m7[mb_x+xx][mb_y+yy];
}



/*!
 ************************************************************************
 * \brief
 *    Make Intra prediction for all 9 modes for larger blocks than 4*4,
 *    that is for 4*8, 8*4 and 8*8 blocks.
 *    bs_x and bs_y may be only 4 and 8.
 *    img_x and img_y are pixels offsets in the picture.
 ************************************************************************
*/
int intrapred_ABT(struct img_par *img,int img_x,int img_y,int bs_x,int bs_y)
{
  unsigned char edgepixels[40];
#define EP (edgepixels+20)
  int x,y,sum,number,last_pix,new_pix;
  unsigned int x_off,y_off;
  int i,predmode;
  int b4_x,b4_y,b4_ox,b4_oy;
  int b4,b8,mb,ob4,ob8,omb;
  int block_available_up,block_available_up_right;
  int block_available_left,block_available_left_down;

    //This intrapred uses stronger filtering for the larger blocks than the normal intraprediction
    assert( (bs_x==4||bs_x==8) && (bs_y==4||bs_y==8) ); //check blocksize
    assert( !(img_x&(bs_x-1)) && !(img_y&(bs_y-1)) );   //check position must be multiple of blocksize

    x_off=img_x&0x0FU;
    y_off=img_y&0x0FU;
    if( bs_x==4 && bs_y==4 )
    {
      if (intrapred(img,x_off,y_off,img_x>>2,img_y>>2)==SEARCH_SYNC)  /* make JM 4x4 prediction */
      {
        return SEARCH_SYNC;                   /* bit error */
      }
      else
      {
        return DECODING_OK;
      }
    }
#ifdef USE_6_INTRA_MODES
    assert(0);        //not supported (could be added)
#endif

    predmode = img->ipredmode[img_x/BLOCK_SIZE+1][img_y/BLOCK_SIZE+1];
    assert(predmode>=0);

    b4_x=img_x>>2;
    b4_y=img_y>>2;
    mb=(b4_x>>2)+(b4_y>>2)*(img->width>>2);  b8=((b4_x&2)>>1)|(b4_y&2);  b4=(b4_x&1)|((b4_y&1)<<1);
/*
    block_available_up       =( img->ipredmode[img_x/BLOCK_SIZE+1][img_y/BLOCK_SIZE] >= 0 );
    block_available_up_right =( img->ipredmode[img_x/BLOCK_SIZE+1+(bs_x>>2)][img_y/BLOCK_SIZE] >= 0 );
    block_available_left     =( img->ipredmode[img_x/BLOCK_SIZE][img_y/BLOCK_SIZE+1] >= 0 );
    block_available_left_down=( img->ipredmode[img_x/BLOCK_SIZE][img_y/BLOCK_SIZE+1+(bs_y>>2)] >= 0 );
*/
    //check block up
    b4_ox=b4_x  ;b4_oy=b4_y-1;        //other block
    omb=(b4_ox>>2)+(b4_oy>>2)*(img->width>>2);  ob8=((b4_ox&2)>>1)|(b4_oy&2);  ob4=(b4_ox&1)|((b4_oy&1)<<1);
    block_available_up=( b4_oy>=0 && img->ipredmode[1+b4_ox][1+b4_oy]>=0 && ((omb<<4)|(ob8<<2)|ob4)<((mb<<4)|(b8<<2)|b4) );        //available and earlier in sequence?

    //check block up right
    b4_ox=b4_x+(bs_x>>2);b4_oy=b4_y-1;        //other block
    omb=(b4_ox>>2)+(b4_oy>>2)*(img->width>>2);  ob8=((b4_ox&2)>>1)|(b4_oy&2);  ob4=(b4_ox&1)|((b4_oy&1)<<1);
    block_available_up_right=( b4_oy>=0 && img->ipredmode[1+b4_ox][1+b4_oy]>=0 && ((omb<<4)|(ob8<<2)|ob4)<((mb<<4)|(b8<<2)|b4) );        //available and earlier in sequence?

    //check block left
    b4_ox=b4_x-1;b4_oy=b4_y  ;        //other block
    omb=(b4_ox>>2)+(b4_oy>>2)*(img->width>>2);  ob8=((b4_ox&2)>>1)|(b4_oy&2);  ob4=(b4_ox&1)|((b4_oy&1)<<1);
    block_available_left=( b4_ox>=0 && img->ipredmode[1+b4_ox][1+b4_oy]>=0 && ((omb<<4)|(ob8<<2)|ob4)<((mb<<4)|(b8<<2)|b4) );        //available and earlier in sequence?

    //check block left down
    b4_ox=b4_x-1;b4_oy=b4_y+(bs_y>>2);        //other block
    omb=(b4_ox>>2)+(b4_oy>>2)*(img->width>>2);  ob8=((b4_ox&2)>>1)|(b4_oy&2);  ob4=(b4_ox&1)|((b4_oy&1)<<1);
    block_available_left_down=( b4_ox>=0 && img->ipredmode[1+b4_ox][1+b4_oy]>=0 && ((omb<<4)|(ob8<<2)|ob4)<((mb<<4)|(b8<<2)|b4) );        //available and earlier in sequence?

    sum  =  ( ((img_x+(bs_x>>3))>>2) & 3 )  |  (      img_y        & 12 ) ;
    if( (1<<sum) & 0x0515U )block_available_up_right=0;
    sum  =  ( (    img_x        >>2) & 3 )  |  ( (img_y+(bs_y>>3)) & 12 ) ;
    if( (1<<sum) & 0xA8A0U )block_available_left_down=0;

    //get prediciton pixels

    if(block_available_up)
    {
      for(x=0;x<bs_x;x++)
        EP[x+1]=imgY[img_y-1][img_x+x];
      if(block_available_up_right)
      {
        for(x=0;x<bs_x;x++)
          EP[1+x+bs_x]=imgY[img_y-1][img_x+bs_x+x];
        for(;x<bs_y;x++)
          EP[1+x+bs_x]=imgY[img_y-1][img_x+bs_x+bs_x-1];
      }else
        for(x=0;x<bs_y;x++)
          EP[1+x+bs_x]=EP[bs_x];
      for(;x<bs_y+2;x++)
        EP[1+x+bs_x]=EP[bs_x+x];
      EP[0]=imgY[img_y-1][img_x];
    }
    if(block_available_left)
    {
      for(y=0;y<bs_y;y++)
        EP[-1-y]=imgY[img_y+y][img_x-1];
      if(block_available_left_down)
      {
        for(y=0;y<bs_y;y++)
          EP[-1-y-bs_y]=imgY[img_y+bs_y+y][img_x-1];
        for(;y<bs_x;y++)
          EP[-1-y-bs_y]=imgY[img_y+bs_y+bs_y-1][img_x-1];
      }else
        for(y=0;y<bs_x;y++)
          EP[-1-y-bs_y]=EP[-bs_y];
      for(;y<bs_x+2;y++)
        EP[-1-y-bs_y]=EP[-y-bs_y];
      EP[0]=imgY[img_y][img_x-1];
    }
    if(block_available_up&&block_available_left)
      EP[0]=imgY[img_y-1][img_x-1];

    //lowpass (Those emlements that are not needed will not disturb)
    last_pix=EP[-(bs_x+bs_y)];
    for(i=-(bs_x+bs_y);i<=(bs_x+bs_y);i++)
    {
      new_pix=( last_pix + (EP[i]<<1) + EP[i+1] + 2 )>>2;
      last_pix=EP[i];
      EP[i]=(unsigned char)new_pix;
    }

    switch(predmode)
    {
    case DC_PRED:// 0 DC
      number=0UL;sum=0UL;
      if(block_available_up)
      {
        for(x=0UL;x<bs_x;x++)sum+=EP[1+x];
        number+=bs_x;
      }
      if(block_available_left)
      {
        for(y=0UL;y<bs_y;y++)sum+=EP[-1-y];
        number+=bs_y;
      }
      if(number>0)
        sum=(sum+(number>>1))/number;        //division could be replaced by a table lookup
      else
        sum=128UL;
      for(y=0UL;y<bs_y;y++)
        for(x=0UL;x<bs_x;x++)
          img->mpr[x+x_off][y+y_off]=(int)sum;
      break;

    case VERT_PRED:// 1 vertical
      if(!block_available_up)printf("!!invalid intrapred mode %d!!\n",predmode);
      for(y=0UL;y<bs_y;y++)
        for(x=0UL;x<bs_x;x++)
          img->mpr[x+x_off][y+y_off]=EP[1+x];
      break;

    case HOR_PRED:// 2 horizontal
      if(!block_available_left)printf("!!invalid intrapred mode %d!!\n",predmode);
      for(y=0UL;y<bs_y;y++)
        for(x=0UL;x<bs_x;x++)
          img->mpr[x+x_off][y+y_off]=EP[-1-y];
      break;

    case DIAG_PRED_SE:// 3 down-right
      if(!block_available_left||!block_available_up)printf("!!invalid intrapred mode %d!!\n",predmode);
      for(y=0UL;y<bs_y;y++)
        for(x=0UL;x<bs_x;x++)
          img->mpr[x+x_off][y+y_off]=EP[(int)x-(int)y];
      break;

    case DIAG_PRED_NE:// 4 up-right bidirectional
      if(!block_available_left||!block_available_up)printf("!!invalid intrapred mode %d!!\n",predmode);
      for(y=0UL;y<bs_y;y++)
        for(x=0UL;x<bs_x;x++)
          img->mpr[x+x_off][y+y_off]=(EP[2+x+y]+EP[-2-(int)(x+y)])>>1;
      break;

    case DIAG_PRED_SSE:// 5 down-right-down
      if(!block_available_left||!block_available_up)printf("!!invalid intrapred mode %d!!\n",predmode);
      for(y=0UL;y<bs_y;y+=2)//even lines
        for(x=0UL;x<bs_x;x++)
          if( (i=x-(int)(y>>1)) >= 0 )
            img->mpr[x+x_off][y+y_off]=(EP[i]+EP[1+i])>>1;
          else
            img->mpr[x+x_off][y+y_off]=EP[1+2*(int)x-(int)y];
      for(y=1UL;y<bs_y;y+=2)//odd lines
        for(x=0UL;x<bs_x;x++)
          if( (i=x-(int)(y>>1)) >= 0 )
            img->mpr[x+x_off][y+y_off]=EP[i];
          else
            img->mpr[x+x_off][y+y_off]=EP[1+2*(int)x-(int)y];
      break;

    case DIAG_PRED_NNE:// 6 down-left-down
      if(!block_available_up)printf("!!invalid intrapred mode %d!!\n",predmode);
      for(y=0UL;y<bs_y;y+=2)//even lines
        for(x=0UL;x<bs_x;x++)
          img->mpr[x+x_off][y+y_off]=(EP[1+x+(y>>1)]+EP[2+x+(y>>1)])>>1;
      for(y=1UL;y<bs_y;y+=2)//odd lines
        for(x=0UL;x<bs_x;x++)
          img->mpr[x+x_off][y+y_off]=EP[2+x+(y>>1)];
      break;

    case DIAG_PRED_ENE:// 7 right-up-right
      if(!block_available_left)printf("!!invalid intrapred mode %d!!\n",predmode);
      for(y=0UL;y<bs_y;y++)//even columns
        for(x=0UL;x<bs_x;x+=2)
          img->mpr[x+x_off][y+y_off]=(EP[-1-(int)(y+(x>>1))]+EP[-2-(int)(y+(x>>1))])>>1;
      for(y=0UL;y<bs_y;y++)//odd columns
        for(x=1UL;x<bs_x;x+=2)
          img->mpr[x+x_off][y+y_off]=EP[-2-(int)(y+(x>>1))];
      break;

    case DIAG_PRED_ESE:// 8 right-down-right
      if(!block_available_left||!block_available_up)printf("!!invalid intrapred mode %d!!\n",predmode);
      for(y=0UL;y<bs_y;y++)//even columns
        for(x=0UL;x<bs_x;x+=2)
          if( (i=-(int)y+(x>>1)) <= 0 )
            img->mpr[x+x_off][y+y_off]=(EP[i]+EP[i-1])>>1;
          else
            img->mpr[x+x_off][y+y_off]=EP[-1-2*(int)y+(int)x];
      for(y=0UL;y<bs_y;y++)//odd columns
        for(x=1UL;x<bs_x;x+=2)
          if( (i=-(int)y+(x>>1)) <= 0 )
            img->mpr[x+x_off][y+y_off]=EP[i];
          else
            img->mpr[x+x_off][y+y_off]=EP[-1-2*(int)y+(int)x];
      break;
    default:
      printf("Error: illegal prediction mode input: %d\n",predmode);
      return SEARCH_SYNC;
      break;
    }//end switch predmode

    return DECODING_OK;
#undef EP
}

/*!
 ************************************************************************
 * \brief
 *    Get constants needed for quantization
 ************************************************************************
*/
void get_quant_consts(int abt_mode,int qp,int intra,int Q[2][2],int Qrshift[2][2],int qp_const[2][2])
{
 int x,y;
 int divfac;

  divfac=( intra ? ((intra==1)?3:2) : 6 );     // inter=0:6, intra=1:3, sp=2:2
  for(y=0;y<2;y++)
    for(x=0;x<2;x++)
    {
      Q[y][x]      = ABT_Q[abt_mode][qp%QUANT_PERIOD][ABT_QMAP[abt_mode][y][x]] ;
      Qrshift[y][x] = ABT_N[abt_mode][0] - ABT_SHIFT0[abt_mode][y][x] + qp/QUANT_PERIOD ;
      qp_const[y][x] = (1<<Qrshift[y][x])/divfac;
    }
}

/*!
 ************************************************************************
 * \brief
 *    Get constants needed for de-quantization
 ************************************************************************
*/
void get_dequant_consts(int abt_mode,int qp,int R[2][2])
{
 int x,y;
  for(y=0;y<2;y++)
    for (x=0;x<2;x++)
      R[y][x]=ABT_R[abt_mode][qp%QUANT_PERIOD][ABT_QMAP[abt_mode][y][x]];
}



