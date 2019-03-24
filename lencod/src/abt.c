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
#include "rdopt_coding_state.h"

extern const byte mapTab[9];

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

  assert(input->abt); // just to be shure.
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

  assert(input->abt); // just to be shure.
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
 *    inverse DCT and transform of ABT block
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
  assert(input->abt); // just to be shure.
  assert(qp>=0); // remove after redesign of quantization
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
  // It is assumed that blk_off_x,blk_off_y are always even.
  // Needed for proper quantization table indexing. assert?
  get_quant_consts(abt_mode,qp,intra,Q,q_bits,qp_const);

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
 *    Quantization, scan, and reconstruction of a transformed ABT block
 *    Return coeff_cost of the block.
 ************************************************************************
 */
int scanquant_ABT_B8(int qp, int abt_mode, int block8x8, // block8x8: number of current 8x8 block
                     int blk_off_x,int blk_off_y,
                     int curr_blk[B8_SIZE][B8_SIZE],
                     int scrFlag,
                     int *cbp,
                     int *cbp_blk)
{
  int  run;
  int  xx, yy,x1,y1;
  int  intra;
  int  iblk_x, iblk_y, icoef, ipos;
  int  blk_size_x, blk_size_y;
  int  boff_x, boff_y, sblk_x, sblk_y;
  int  b8_y       = (block8x8 / 2) << 3;
  int  b8_x       = (block8x8 % 2) << 3;
  int  num_blk_x,  num_blk_y;
  int  num_coeff;
  int  blk_pos;                       // index into the JM 4x4-block raster for each 8x8 block
  int  coeff_cost    = 0;
  int  cbp_mask      = 1 << block8x8;
  int  cbp_blk_mask;
  int  q_bits, q_shift;
  int  qp_const;
  int  R[2][2];
  int* ACLevel;
  int* ACRun;
  int  frm_fld; // indicate progressive/interlaced scan use.
  int  curr_val;

#ifdef _ALT_SCAN_
  frm_fld = (img->pstruct>0); // img->pstruct==0: progressive, else interlaced
#else
  frm_fld = 0;
#endif

  assert(input->abt);
  assert(abt_mode>=0);
  // =========================================================
  // quantize
  // =========================================================
  quant_abt_B8(qp, abt_mode, blk_off_x, blk_off_y, curr_blk);

  // =========================================================
  // scan
  // =========================================================
  intra=0;
  if(abt_mode>3)intra=1;
  abt_mode %= 4; // abt_mode+=4 is used to signal intra coding to quant_abt_B8().

  // dequantization values.
  //q_bits   = ABT_N[abt_mode][1] - ABT_SHIFT1[abt_mode] - qp/QUANT_PERIOD;
  q_bits   = ABT_N[abt_mode][1] - ABT_SHIFT1[abt_mode];
  q_shift  = qp/QUANT_PERIOD;
  qp_const = 1<<(q_bits-1);

  get_dequant_consts(abt_mode,qp,R);

  // general block information
  blk_size_x = ABT_TRSIZE[abt_mode][PIX];
  blk_size_y = ABT_TRSIZE[abt_mode][LIN];
  num_coeff  = blk_size_x * blk_size_y;
  num_blk_x  = num_blk_y  = 1;

  if (blk_off_x == WHOLE_BLK)
  {
    assert(blk_off_y == WHOLE_BLK); // just to be shure, may be removed later
    boff_x = 0;
    boff_y = 0;
    blk_pos   = 0;
    num_blk_x = ABT_NUMTR[abt_mode][PIX];
    num_blk_y = ABT_NUMTR[abt_mode][LIN];
  }
  else
  {
    assert(blk_off_x+blk_size_x <= B8_SIZE);
    assert(blk_off_y+blk_size_y <= B8_SIZE);
    boff_x   = blk_off_x;
    boff_y   = blk_off_y;
    blk_pos  = ((blk_off_y&4)>>1) | ((blk_off_x&4)>>2) ;
  }

  for (iblk_y=0; iblk_y<num_blk_y; iblk_y++)
  {
    sblk_y = iblk_y * blk_size_y;
    for (iblk_x=0; iblk_x<num_blk_x; iblk_x++)
    {
      sblk_x = iblk_x * blk_size_x;

      ACLevel  = img->cofAC[block8x8][blk_pos+(iblk_y<<1)+iblk_x][0];
      ACRun    = img->cofAC[block8x8][blk_pos+(iblk_y<<1)+iblk_x][1];
      for (xx=0; xx<65; xx++)
        ACRun[xx] = ACLevel[xx] = 0; // just to be shure.

      run  = -1;
      ipos = 0;

      for (icoef=0; icoef<num_coeff; icoef++)
      {
        run++;
        xx = ABT_SCAN[frm_fld][abt_mode][icoef][0];
        yy = ABT_SCAN[frm_fld][abt_mode][icoef][1];
        curr_val = curr_blk[boff_y + sblk_y + yy][boff_x + sblk_x + xx];

        if (curr_val != 0)
        {
          ACLevel[ipos] = curr_val;
          ACRun[ipos]   = run;

          if (scrFlag && absm(ACLevel[ipos])==1)
            coeff_cost += (scrFlag==1)? 1 : ABT_COEFF_COST[abt_mode][run];
          else
            coeff_cost += MAX_VALUE; // block has to be saved

          curr_blk[boff_y + sblk_y + yy][boff_x + sblk_x + xx] = (curr_val*R[yy&1][xx&1])<<q_shift; // dequantization

          run = -1;
          ipos++;
        }
      }

      if (ipos>0) // there are coefficients
      {
        cbp_blk_mask = cbp_blk_masks[ abt_mode&3 ][ (iblk_y<<1)|iblk_x ] ;
        if(block8x8&1)cbp_blk_mask<<=2;
        if(block8x8&2)cbp_blk_mask<<=8;

        (*cbp)     |= cbp_mask;
        (*cbp_blk) |= cbp_blk_mask;
      }

    }
  }

  // =========================================================
  // reconstruct current 8x8 block and add to prediction block
  // =========================================================


  // add prediction and normalize (normalize here to avoid problem of rounding of negative values)
  x1=boff_x+num_blk_x*ABT_TRSIZE[abt_mode][0];
  y1=boff_y+num_blk_y*ABT_TRSIZE[abt_mode][1];

  // inverse transform
  inv_transform_ABT_B8(abt_mode, blk_off_x, blk_off_y, curr_blk);   // here: abt_mode 0..3  .  blk_off_x,y can be WHOLE_BLK

  //normal reconstruction
  for(yy=boff_y;yy<y1;yy++)
    for(xx=boff_x;xx<x1;xx++)
    {
      curr_val = ((img->mpr[b8_x+xx][b8_y+yy]<<q_bits) + curr_blk[yy][xx] + qp_const) >> q_bits;
      img->m7[xx][yy] = curr_blk[yy][xx] = clamp(curr_val,0,255);

      imgY[img->pix_y+b8_y+yy][img->pix_x+b8_x+xx] = curr_blk[yy][xx];
    }
  return coeff_cost;
}

int trans_scanquant_ABT_sp(int abt_mode, int block8x8, // block8x8: number of current 8x8 block
                           int blk_off_x,int blk_off_y,
                           int curr_blk[B8_SIZE][B8_SIZE],
                           int scrFlag,
                           int *cbp,
                           int *cbp_blk)
{
  int  run;
  int  xx, yy, x1, y1;
  int  ixx, iyy;
  int  iblk_x, iblk_y, icoef, ipos;
  int  blk_size_x, blk_size_y;
  int  boff_x, boff_y, sblk_x, sblk_y;
  int  b8_y       = (block8x8 / 2) << 3;
  int  b8_x       = (block8x8 % 2) << 3;
  int  num_blk_x,  num_blk_y;
  int  num_coeff;
  int  blk_pos;                       // index into the JM 4x4-block raster for each 8x8 block
  int  coeff_cost    = 0;
  int  cbp_mask      = 1 << block8x8;
  int  cbp_blk_mask;
  int  qp, qpsp;
  int  q_bits[2][2];
  int  qp_const[2][2];
  int  q_shift, dq_bits, dqp_const;
  int  Q[2][2];
  int  R[2][2];
  int  q_bits_sp[2][2];
  int  qp_const_sp[2][2];
  int  q_shift_sp;
  int  Qsp[2][2];
  int  Rsp[2][2];
  int  vq_bits, vqp_const,vq_bits_sp, vqp_const_sp, vQ, vQsp;
  int  pred_blk[B8_SIZE][B8_SIZE];
  int* ACLevel;
  int* ACRun;
  int  frm_fld; // indicate progressive/interlaced scan use.
  int  curr_val;

  int    level;
  int    c_err, c_err1, c_err2, level1, level2;
  double D_dis1, D_dis2;
  int    len, info;
  double lambda_mode   = 0.85 * pow (2, img->qp/3.0) * 4;

#ifdef _ALT_SCAN_
  frm_fld = (img->pstruct>0); // img->pstruct==0: progressive, else interlaced
#else
  frm_fld = 0;
#endif

  assert(input->abt);
  assert(abt_mode>=0);
  assert(abt_mode<4);

  qp        = img->qp+QP_OFS-MIN_QP;
  qpsp      = img->qpsp+QP_OFS-MIN_QP;

  dq_bits   = ABT_N[abt_mode][1] - ABT_SHIFT1[abt_mode];
  dqp_const = 1<<(dq_bits-1);

  q_shift    = qp/QUANT_PERIOD;
  q_shift_sp = qpsp/QUANT_PERIOD;

  get_quant_consts(abt_mode,qp,2,Q,q_bits,qp_const);
  get_dequant_consts(abt_mode,qp,R);
  get_quant_consts(abt_mode,qpsp,2,Qsp,q_bits_sp,qp_const_sp);
  get_dequant_consts(abt_mode,qpsp,Rsp);


  // general block information
  blk_size_x = ABT_TRSIZE[abt_mode][PIX];
  blk_size_y = ABT_TRSIZE[abt_mode][LIN];
  num_coeff  = blk_size_x * blk_size_y;
  num_blk_x  = num_blk_y  = 1;

  if (blk_off_x == WHOLE_BLK)
  {
    assert(blk_off_y == WHOLE_BLK); // just to be shure, may be removed later
    boff_x = 0;
    boff_y = 0;
    blk_pos   = 0;
    num_blk_x = ABT_NUMTR[abt_mode][PIX];
    num_blk_y = ABT_NUMTR[abt_mode][LIN];
  }
  else
  {
    assert(blk_off_x+blk_size_x <= B8_SIZE);
    assert(blk_off_y+blk_size_y <= B8_SIZE);
    boff_x   = blk_off_x;
    boff_y   = blk_off_y;
    blk_pos  = ((blk_off_y&4)>>1) | ((blk_off_x&4)>>2) ;
  }
  x1=boff_x+num_blk_x*ABT_TRSIZE[abt_mode][0];
  y1=boff_y+num_blk_y*ABT_TRSIZE[abt_mode][1];

  // get prediction block
  for(yy=boff_y;yy<y1;yy++)
    for(xx=boff_x;xx<x1;xx++)
    {
      curr_blk[yy][xx] += img->mpr[b8_x+xx][b8_y+yy]; // MC+pred.error =>the original block
      pred_blk[yy][xx]  = img->mpr[b8_x+xx][b8_y+yy]; // the MC block
    }

  transform_ABT_B8(abt_mode,WHOLE_BLK,WHOLE_BLK,curr_blk);
  transform_ABT_B8(abt_mode,WHOLE_BLK,WHOLE_BLK,pred_blk);

  // do scan
  for (iblk_y=0; iblk_y<num_blk_y; iblk_y++)
  {
    sblk_y = iblk_y * blk_size_y;
    for (iblk_x=0; iblk_x<num_blk_x; iblk_x++)
    {
      sblk_x = iblk_x * blk_size_x;

      ACLevel  = img->cofAC[block8x8][blk_pos+(iblk_y<<1)+iblk_x][0];
      ACRun    = img->cofAC[block8x8][blk_pos+(iblk_y<<1)+iblk_x][1];
      for (xx=0; xx<65; xx++)
        ACRun[xx] = ACLevel[xx] = 0; // just to be shure.

      run  = -1;
      ipos = 0;

      for (icoef=0; icoef<num_coeff; icoef++)
      {
        run++;
        curr_val = 0;
        xx = ABT_SCAN[frm_fld][abt_mode][icoef][0];
        yy = ABT_SCAN[frm_fld][abt_mode][icoef][1];

        ixx = boff_x + sblk_x + xx;
        iyy = boff_y + sblk_y + yy;

        // quantization
        //curr_val = curr_blk[iyy][ixx];

        // decide prediction

        // img->m7[i][j]               -> curr_blk[iyy][ixx]
        // quant_coef[qp_rem][i][j]    -> Q[yy&1][xx&1]
        // quant_coef[qp_rem_sp][i][j] -> Qsp[yy&1][xx&1]

        // quantization values:
        vqp_const    = qp_const   [yy&1][xx&1];
        vqp_const_sp = qp_const_sp[yy&1][xx&1];
        vq_bits      = q_bits     [yy&1][xx&1];
        vq_bits_sp   = q_bits_sp  [yy&1][xx&1];
        vQ           = Q  [yy&1][xx&1];
        vQsp         = Qsp[yy&1][xx&1];

        // case 1
        level1 = (absm (pred_blk[iyy][ixx]) * vQsp + vqp_const_sp) >> vq_bits_sp;
        level1 = (level1 << vq_bits_sp) / vQsp;
        c_err1 = curr_blk[iyy][ixx]-sign(level1, pred_blk[iyy][ixx]);
        level1 = (abs (c_err1) * vQ + vqp_const) >> vq_bits;

        // case 2
        c_err2 = curr_blk[iyy][ixx]-pred_blk[iyy][ixx];
        level2 = (abs (c_err2) * vQ + vqp_const) >> vq_bits;

        // select prediction
        if ((level1 != level2) && (level1 != 0) && (level2 != 0))
        {
          D_dis1 = curr_blk[iyy][ixx] - sign(((level1 << vq_bits) + vQ / 2) / vQ, c_err1) - pred_blk[iyy][ixx];
          len = level1+run; info=0;// levrun_linfo_inter(level1, run, &len, &info); // workaround. need routine to get golomb-word length
          D_dis1 = D_dis1*D_dis1 + lambda_mode * len;

          D_dis2 = curr_blk[iyy][ixx] - sign(((level2 << vq_bits) + vQ / 2) / vQ, c_err2) - pred_blk[iyy][ixx];
          len = level2+run; info=0;// levrun_linfo_inter(level2, run, &len, &info); // workaround. need routine to get golomb-word length
          D_dis2 = D_dis2 * D_dis2 + lambda_mode * len;

          if (D_dis1 == D_dis2)
            level = (abs(level1) < abs(level2)) ? level1 : level2;
          else
          {
            if (D_dis1 < D_dis2)
              level = level1;
            else
              level = level2;
          }
          c_err = (level == level1) ? c_err1 : c_err2;
        }
        else if (level1 == level2)
        {
          level = level1;
          c_err = c_err1;
        }
        else
        {
          level = (level1 == 0) ? level1 : level2;
          c_err = (level1 == 0) ? c_err1 : c_err2;
        }

        if (level != 0)
        {
          ACLevel[ipos] = sign(level,c_err);
          ACRun[ipos]   = run;

          if (level==1)
            coeff_cost += (scrFlag==1)? 1 : ABT_COEFF_COST[abt_mode][run];
          else
            coeff_cost += MAX_VALUE; // block has to be saved

          run = -1;
          ipos++;
          curr_val = level;
        }

        curr_val = sign(((curr_val << vq_bits) + vQ / 2) / vQ, c_err) + pred_blk[iyy][ixx];
        curr_blk[iyy][ixx] = sign((abs(curr_val) * vQsp + vqp_const_sp)>> vq_bits_sp, curr_val) * Rsp[yy&1][xx&1] << q_shift_sp;
      }

      if (ipos>0) // there are coefficients
      {
        cbp_blk_mask = cbp_blk_masks[ abt_mode&3 ][ (iblk_y<<1)|iblk_x ] ;
        if(block8x8&1)cbp_blk_mask<<=2;
        if(block8x8&2)cbp_blk_mask<<=8;

        (*cbp)     |= cbp_mask;
        (*cbp_blk) |= cbp_blk_mask;
      }

    }
  }

  // inverse transform
  inv_transform_ABT_B8(abt_mode, blk_off_x, blk_off_y, curr_blk);

  for(yy=boff_y;yy<y1;yy++)
    for(xx=boff_x;xx<x1;xx++)
    {
      level = (curr_blk[yy][xx] + dqp_const) >> dq_bits;
      img->m7[xx][yy] = curr_blk[yy][xx] = clamp(level,0,255);

      imgY[img->pix_y+b8_y+yy][img->pix_x+b8_x+xx] = curr_blk[yy][xx];
    }

  return coeff_cost;



}


/*!
 ************************************************************************
 * \brief
 *    calculate SAD or ABT-SATD for a prediction error block of size
 *    iSizeX x iSizeY.
 ************************************************************************
 */
int find_sad_abt(int iMode, int iSizeX, int iSizeY, int iOffX, int iOffY, int m7[MB_BLOCK_SIZE][MB_BLOCK_SIZE])
{
  int
    i,j,
    bmode,
    ishift = 0,
    sad    = 0;

    assert( (iSizeX+iOffX) <= MB_BLOCK_SIZE );
    assert( (iSizeY+iOffY) <= MB_BLOCK_SIZE );

  assert(input->abt);
  // m7[y,j,line][x,i,pixel]
  switch (iMode)
    {
    case 0 : // ---------- SAD ----------
      for (j=iOffY; j < iOffY+iSizeY; j++)
        for (i=iOffX; i < iOffX+iSizeX; i++)
          sad += absm(m7[j][i]);
      break;

    case 1 : // --------- SATD ----------
      bmode = iSizeX+iSizeY;
      if (bmode<24) // 4x4-8x8
        {
          sad = sad_hadamard(iSizeY,iSizeX,iOffY,iOffX,m7); // Attention: sad_hadamard() is X/Y flipped
          switch (bmode)
            {
            case 8 :                // 4x4       : sad/1.75
              ishift = 4;
              sad   *= 9;
              break;
            case 12 :               // 4x8 8x4   : sad/sqrt(8)
              ishift = 9;
              sad   *= 181;
              break;
            case 16 :               // 8x8       : sad/4
              ishift = 2;
              break;
            default :
              assert(0==1);
            }
          sad = (sad + (1<<(ishift-1)))>>ishift;
        }
      else // 8x16-16x16
        {
          switch (bmode)
            {
            case 24 :               // 16x8 8x16
              sad  = sad_hadamard(8,8,iOffY,iOffX,m7);
              sad += sad_hadamard(8,8,iOffY+((iSizeY==16)?8:0) , iOffX+((iSizeX==16)?8:0) ,m7);
              ishift = 2;
              break;
            case 32 :               // 16x16
              sad  = sad_hadamard(8,8,0,0,m7);
              sad += sad_hadamard(8,8,8,0,m7);
              sad += sad_hadamard(8,8,0,8,m7);
              sad += sad_hadamard(8,8,8,8,m7);
              ishift = 2;
              break;
            default :
              assert(0==1);
            }
          sad = (sad + (1<<(ishift-1)))>>ishift;
        }
      break;

    default :
      assert(0==1);                // more switches may be added here later
    }

  return sad;
}



/*!
 ************************************************************************
 * \brief
 *    calculates the SAD of the Hadamard transformed block of
 *    size iSizeX*iSizeY. Block may have an offset of (iOffX,iOffY).
 *    If offset!=0 then iSizeX/Y has to be <=8.
 ************************************************************************
 */
int sad_hadamard(int iSizeX, int iSizeY, int iOffX, int iOffY, int m7[MB_BLOCK_SIZE][MB_BLOCK_SIZE])
{
  int
    i,j,
    ii,

    m1[MB_BLOCK_SIZE][MB_BLOCK_SIZE],
    m2[MB_BLOCK_SIZE][MB_BLOCK_SIZE],
    m3[MB_BLOCK_SIZE][MB_BLOCK_SIZE],

    iy[MB_BLOCK_SIZE] =
    {
        iOffY,iOffY,iOffY,iOffY,
        iOffY,iOffY,iOffY,iOffY,
        iOffY,iOffY,iOffY,iOffY,
        iOffY,iOffY,iOffY,iOffY
    },

    sad = 0;

  // ==============================================
  // in this routine, cols are j,y and rows are i,x
  // ==============================================

  assert( ((iOffX==0)||(iSizeX<=8)) && ((iOffY==0)||(iSizeY<=8)) );

  for (j=1; j<iSizeY; j++)
    iy[j] += j;

  // vertical transform
  if (iSizeY == 4)
    for (i=0; i < iSizeX; i++)
    {
      ii = i+iOffX;

      m1[i][0] = m7[ii][iy[0]] + m7[ii][iy[3]];
      m1[i][1] = m7[ii][iy[1]] + m7[ii][iy[2]];
      m1[i][2] = m7[ii][iy[1]] - m7[ii][iy[2]];
      m1[i][3] = m7[ii][iy[0]] - m7[ii][iy[3]];

      m3[i][0] = m1[i][0] + m1[i][1];
      m3[i][1] = m1[i][0] - m1[i][1];
      m3[i][2] = m1[i][2] + m1[i][3];
      m3[i][3] = m1[i][3] - m1[i][2];
    }
  else
    for (i=0; i < iSizeX; i++)
    {
      ii = i+iOffX;

      m1[i][0] = m7[ii][iy[0]] + m7[ii][iy[4]];
      m1[i][1] = m7[ii][iy[1]] + m7[ii][iy[5]];
      m1[i][2] = m7[ii][iy[2]] + m7[ii][iy[6]];
      m1[i][3] = m7[ii][iy[3]] + m7[ii][iy[7]];
      m1[i][4] = m7[ii][iy[0]] - m7[ii][iy[4]];
      m1[i][5] = m7[ii][iy[1]] - m7[ii][iy[5]];
      m1[i][6] = m7[ii][iy[2]] - m7[ii][iy[6]];
      m1[i][7] = m7[ii][iy[3]] - m7[ii][iy[7]];

      m2[i][0] = m1[i][0] + m1[i][2];
      m2[i][1] = m1[i][1] + m1[i][3];
      m2[i][2] = m1[i][0] - m1[i][2];
      m2[i][3] = m1[i][1] - m1[i][3];
      m2[i][4] = m1[i][4] + m1[i][6];
      m2[i][5] = m1[i][5] + m1[i][7];
      m2[i][6] = m1[i][4] - m1[i][6];
      m2[i][7] = m1[i][5] - m1[i][7];

      m3[i][0] = m2[i][0] + m2[i][1];
      m3[i][1] = m2[i][0] - m2[i][1];
      m3[i][2] = m2[i][2] + m2[i][3];
      m3[i][3] = m2[i][2] - m2[i][3];
      m3[i][4] = m2[i][4] + m2[i][5];
      m3[i][5] = m2[i][4] - m2[i][5];
      m3[i][6] = m2[i][6] + m2[i][7];
      m3[i][7] = m2[i][6] - m2[i][7];
    }


  // horizontal transform
  if (iSizeX == 4)
    for (j=0; j < iSizeY; j++)
    {
      m1[0][j]=m3[0][j]+m3[3][j];
      m1[1][j]=m3[1][j]+m3[2][j];
      m1[2][j]=m3[1][j]-m3[2][j];
      m1[3][j]=m3[0][j]-m3[3][j];

      m2[0][j]=m1[0][j]+m1[1][j];
      m2[1][j]=m1[0][j]-m1[1][j];
      m2[2][j]=m1[2][j]+m1[3][j];
      m2[3][j]=m1[3][j]-m1[2][j];

      for (i=0; i < iSizeX; i++)
        sad += absm(m2[i][j]);
    }
  else
    for (j=0; j < iSizeY; j++)
    {
      m2[0][j] = m3[0][j] + m3[4][j];
      m2[1][j] = m3[1][j] + m3[5][j];
      m2[2][j] = m3[2][j] + m3[6][j];
      m2[3][j] = m3[3][j] + m3[7][j];
      m2[4][j] = m3[0][j] - m3[4][j];
      m2[5][j] = m3[1][j] - m3[5][j];
      m2[6][j] = m3[2][j] - m3[6][j];
      m2[7][j] = m3[3][j] - m3[7][j];

      m1[0][j] = m2[0][j] + m2[2][j];
      m1[1][j] = m2[1][j] + m2[3][j];
      m1[2][j] = m2[0][j] - m2[2][j];
      m1[3][j] = m2[1][j] - m2[3][j];
      m1[4][j] = m2[4][j] + m2[6][j];
      m1[5][j] = m2[5][j] + m2[7][j];
      m1[6][j] = m2[4][j] - m2[6][j];
      m1[7][j] = m2[5][j] - m2[7][j];

      m2[0][j] = m1[0][j] + m1[1][j];
      m2[1][j] = m1[0][j] - m1[1][j];
      m2[2][j] = m1[2][j] + m1[3][j];
      m2[3][j] = m1[2][j] - m1[3][j];
      m2[4][j] = m1[4][j] + m1[5][j];
      m2[5][j] = m1[4][j] - m1[5][j];
      m2[6][j] = m1[6][j] + m1[7][j];
      m2[7][j] = m1[6][j] - m1[7][j];

      for (i=0; i < iSizeX; i++)
        sad += absm(m2[i][j]);
    }

  return(sad);
}




/*!
 ************************************************************************
 * \brief
 *    Write ABT coefficients of one 8x8 block
 ************************************************************************
*/
int writeLumaCoeffABT_B8(int b8,int intra,int blk_off_x,int blk_off_y)
{
  int no_bits        = 0;
  int mb_nr          = img->current_mb_nr;
  Macroblock *currMB = &img->mb_data[mb_nr];
  const int cbp      = currMB->cbp;
  int *bitCount      = currMB->bitcounter;
  Slice *currSlice   = img->currentSlice;
  int *partMap       = assignSE2partition[input->partition_mode];
  DataPartition *dataPart;
  SyntaxElement *currSE = &img->MB_SyntaxElements[currMB->currSEnr];

  int
    abt_mode,
    inumblk,                    /* number of blocks per CBP*/
    inumcoeff,                  /* number of coeffs per block */
    iblk,                       /* current block */
    icoef,                      /* current coefficient */
    ipos,                       /* current position in cof_abt */
    run, level,
    blk_max_x,blk_max_y,i,j,
    sbx, sby;
  int use_vlc_numcoef,abs_level;
  int b4;
  int* ACLevel;
  int* ACRun;

  int symbol2D,mode2D,Golomb_se_type;
  const char (*table2D)[8];

  static const int blk_table[4][4] = { {0,-1,-1,-1}, {0,2,-1,-1}, {0,1,-1,-1}, {0,1,2,3}};
  static const int blkmode2ctx [4] = {LUMA_8x8, LUMA_8x4, LUMA_4x8, LUMA_4x4};

  assert(input->abt);

  abt_mode  = currMB->abt_mode[b8];
  inumblk   = ABT_NUMTR [abt_mode][0] * ABT_NUMTR [abt_mode][1];
  inumcoeff = ABT_TRSIZE[abt_mode][0] * ABT_TRSIZE[abt_mode][1] + 1; // all positions + EOB

  assert(abt_mode>=0);

  if(blk_off_x==WHOLE_BLK)
  {
    blk_off_x=blk_off_y=0;
    blk_max_x=blk_max_y=8;
  }else{
    blk_max_x=blk_off_x+ABT_TRSIZE[abt_mode][0];
    blk_max_y=blk_off_y+ABT_TRSIZE[abt_mode][1];
  }

  use_vlc_numcoef=0;

  if (input->symbol_mode == UVLC)
  {
    Golomb_se_type=SE_LUM_AC_INTER;
    if( intra ) // bug fix mwi 020726
    {
      use_vlc_numcoef=1;
      Golomb_se_type=SE_LUM_AC_INTRA;
    }
  }

  if( !intra ) // bug fix: the ELSE part is intra. mwi 020505
  {
    mode2D=2;
    i=0;        //first table is the inter table.
  }else{
    mode2D=2;
    if(abt_mode>B8x8)mode2D--;
    if(abt_mode>=B4x4)mode2D--;
    i=(currMB->qp+QP_OFS-MIN_QP-6)>>3;        //find table to use for intra ( 6..13 -> 0 , 14..21 -> 1 , 22..29 -> 2 )
    if(i<0)i=0;
    if(i>2)i=2;
    i+=1;        //first table is the inter table. Intra tables start off 1.
  }
  table2D=ABT_2D_VLC[i];

  if ( cbp & (1<<b8) )
  {
    //code all symbols
    for (iblk=0; iblk<inumblk; iblk++)
    {
      i=iblk*(ABT_TRSIZE[abt_mode][0]>>2);        // =b4
      j=(i&~1)<<1; //off_x
      i=(i&1)<<2;  //off_y
      if( (i<blk_off_x) || (i>=blk_max_x) || (j<blk_off_y) || (j>=blk_max_y) )
        continue;        //used in ABT Intra RDopt. code one one subblock.

      b4 = blk_table[abt_mode][iblk];



      level = 1; // get inside loop
      ipos  = 0;
      ACLevel = img->cofAC[b8][b4][0];
      ACRun   = img->cofAC[b8][b4][1];

      img->subblock_x = ((b8&0x1)==0)?(((b4&0x1)==0)?0:1):(((b4&0x1)==0)?2:3); // horiz. position for coeff_count context
      img->subblock_y = (b8<2)?((b4<2)?0:1):((b4<2)?2:3); // vert.  position for coeff_count context


      for(icoef=0;icoef<inumcoeff;icoef++) //count coeffs
        if(!ACLevel[icoef])
          break;

      if( use_vlc_numcoef && input->symbol_mode==UVLC )
      {
        //code numcoeffs symbol
        currSE->type=SE_LUM_AC_INTER;
        currSE->golomb_grad=0;
        currSE->golomb_maxlevels=31;
        if( intra )        //Intra blocks contain more. use higher grade here.
        {
          currSE->type=SE_LUM_AC_INTRA; // Bug Fix: INTER->INTRA. mwi 020726
          currSE->golomb_grad=2;
        }
        currSE->value1=icoef;
        currSE->value2=0;
        // choose the appropriate data partition
        if (img->type != B_IMG && img->type != BS_IMG) dataPart = &(currSlice->partArr[partMap[currSE->type]]);
        else                                           dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
        dataPart->writeSyntaxElement(currSE, dataPart);
        bitCount[BITS_COEFF_Y_MB]+=currSE->len;
        no_bits+=currSE->len;
      }

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
      

      for (icoef=0; icoef<inumcoeff && level!=0; icoef++)
      {
        level = ACLevel[icoef];
        run   = ACRun[icoef];

        if (input->symbol_mode == UVLC)
        {
          //VLC coding
          symbol2D=CODE2D_ESCAPE_SYMBOL;        //symbol for out-of-table
          if( level>-8 && level<8 && run<16 )
          {
            symbol2D=table2D[run][level-((level+level)&(level>>31))];
            if(symbol2D>=0)
            {
              if(level<0)symbol2D++;
              if(use_vlc_numcoef)
              {
                if(!level)break;
                assert(symbol2D);        //when using numcoeff, there cannot be the EOB symbol.
                symbol2D--;
              }
            }else
              symbol2D=CODE2D_ESCAPE_SYMBOL;        //symbol for out-of-table
          }
          if(symbol2D==CODE2D_ESCAPE_SYMBOL)
            img+=0;
          //code symbol2D
          currSE->type=Golomb_se_type;
          currSE->value1=symbol2D;
          currSE->value2=0;
          currSE->golomb_grad=mode2D;
          currSE->golomb_maxlevels=6-mode2D;        //make code with little less than 64 symbols.
          // choose the appropriate data partition
          if (img->type != B_IMG && img->type != BS_IMG) dataPart = &(currSlice->partArr[partMap[currSE->type]]);
          else                                           dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
          //code
          dataPart->writeSyntaxElement(currSE, dataPart);
          bitCount[BITS_COEFF_Y_MB]+=currSE->len;//count
          no_bits+=currSE->len;
          currSE++;// proceed to next SE
          currMB->currSEnr++;


          if(!level)break;       //EOB was just coded.


          if(symbol2D>=CODE2D_ESCAPE_SYMBOL)      //was out-of-table ?
          {
            //code level
            currSE->type=Golomb_se_type;
            currSE->golomb_grad=3;
            currSE->golomb_maxlevels=28;
            abs_level=level-((level+level)&(level>>31));        //abs_level=abs(level)
            currSE->value1=(abs_level-1)<<1;
            if(level<0)currSE->value1|=1;           //sign bit
            // choose the appropriate data partition
            if (img->type != B_IMG && img->type != BS_IMG) dataPart = &(currSlice->partArr[partMap[currSE->type]]);
            else                                           dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
            //code
            dataPart->writeSyntaxElement(currSE, dataPart);
            bitCount[BITS_COEFF_Y_MB]+=currSE->len;//count
            no_bits+=currSE->len;
            currSE++;// proceed to next SE
            currMB->currSEnr++;


            //code run with signbit
            currSE->type=Golomb_se_type;
            currSE->value1=run;
            currSE->value2=0;
            currSE->golomb_grad=2;
            currSE->golomb_maxlevels=29;
            // choose the appropriate data partition
            if (img->type != B_IMG && img->type != BS_IMG) dataPart = &(currSlice->partArr[partMap[currSE->type]]);
            else                                           dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);            //code
            dataPart->writeSyntaxElement(currSE, dataPart);
            bitCount[BITS_COEFF_Y_MB]+=currSE->len;//count
            no_bits+=currSE->len;
            currSE++;// proceed to next SE
            currMB->currSEnr++;

          }
        }else{
          //CABAC coding
          currSE->value1 = level;
          currSE->value2 = run;

          currSE->writing = writeRunLevel2Buffer_CABAC;

          currSE->type        = (intra ? (icoef==0 && run==0 ? SE_LUM_DC_INTRA : SE_LUM_AC_INTRA) : (icoef==0 && run==0 ? SE_LUM_DC_INTER : SE_LUM_AC_INTER));
          currSE->context     = blkmode2ctx[abt_mode];
          img->is_intra_block = (intra ? 1 : 0);

          // choose the appropriate data partition
          if (img->type != B_IMG && img->type != BS_IMG) dataPart = &(currSlice->partArr[partMap[currSE->type]]);
          else                                           dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);

          dataPart->writeSyntaxElement(currSE, dataPart);
          bitCount[BITS_COEFF_Y_MB]+=currSE->len;
          no_bits+=currSE->len;
#if TRACE
          snprintf(currSE->tracestring, TRACESTRING_SIZE, "Luma sng(%2d) level =%3d run =%2d", icoef, level,run);
#endif
          // proceed to next SE
          currSE++;
          currMB->currSEnr++;
        }//end if UVLC or CABAC
      }//end loop coeffs
    }
  }

  return no_bits;
}






/*!
 ************************************************************************
 * \brief
 *    Save block mode for the current 8x8 block. Writes to colB8mode.
 ************************************************************************
*/
void setDirectModeABT(int block8x8)
{
  int dirmode;
  int blk_x          = (img->pix_x>>3) + block8x8%2;
  int blk_y          = (img->pix_y>>3) + block8x8/2;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int pstruct        = img->pstruct; // distinguish frame/top/bottom array. 020603 mwi

  assert(pstruct<4); // needs revision for other pstruct values. 020603 mwi

  if      (currMB->mb_type==I16MB)
  {
    dirmode = I16MB;
  }
  else if (currMB->mb_type==0)
  {
    dirmode = B8x8+4; // 8x8 MB mode
  }
  else if( currMB->mb_type==I4MB || (currMB->mb_type==P8x8&&currMB->b8mode[block8x8]==IBLOCK) )
    dirmode=4+(currMB->abt_mode[block8x8]&3);   // gives B4x4 if abt_mode[block8x8] = -1
  else
    dirmode = currMB->b8mode[block8x8];

  colB8mode[pstruct][blk_y][blk_x] = dirmode;
};





/*!
 ************************************************************************
 * \brief
 *    Get Direct MRüdiger Krabbeode block tiling for the current B-Frame 8x8 block.
 *    Writes to currMB->abt_mode,abt_pred_mode and returns abt_mode.
 ************************************************************************
*/
int getDirectModeABT(int block8x8)
{
  int dirmode;
  int blk_x          = (img->pix_x>>3) + block8x8%2;
  int blk_y          = (img->pix_y>>3) + block8x8/2;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  static const int  abt_mode_table[MAXMODE] = {-2, 0, 0, 0, 0, 1, 2, 3, -2, 3, 0, 3, -1};  // DO NOT CHANGE ORDER !!!   // b8mode->abt_mode
  int pstruct        = img->pstruct; // distinguish frame/top/bottom array. 020603 mwi

  assert(pstruct<4); // needs revision for other values. 020603 mwi

  dirmode = abt_mode_table[colB8mode[pstruct][blk_y][blk_x]];

  assert (dirmode != -2); // just to be shure. P8x8 or '0' are not valid for an ABT 8x8 block mode.

  return (currMB->abt_mode[block8x8] = currMB->abt_pred_mode[block8x8] = dirmode);
};

/*!
 ************************************************************************
 * \brief
 *    Encode SP Copy Block with ABT
 ************************************************************************
*/
void copyblock_SP_ABT(int abt_mode,int b8,int blk_off_x, int blk_off_y)
{
  int xx,yy,x1,y1;
  int boff_x,boff_y,num_blk_x,num_blk_y;
  int qp_const,q_bits,q_shift;
  int b8_x,b8_y;
  int pred[8][8];
  int Q1_qpsp[2][2],Qrs_qpsp[2][2],Qc_qpsp[2][2];
  int Q2_qpsp[2][2];
  int val;
  Macroblock *currMB;

  currMB=img->mb_data+img->current_mb_nr;

  // dequantization values.
  q_bits   = ABT_N[abt_mode][1] - ABT_SHIFT1[abt_mode];
  q_shift  = (img->qpsp+QP_OFS-MIN_QP)/QUANT_PERIOD;
  qp_const = 1<<(q_bits-1);

  get_quant_consts(abt_mode,img->qpsp+QP_OFS-MIN_QP,2,Q1_qpsp,Qrs_qpsp,Qc_qpsp);        //Get quant factors for the predicted block
  get_dequant_consts(abt_mode,img->qpsp+QP_OFS-MIN_QP,Q2_qpsp);

  num_blk_x  = num_blk_y  = 1;
  if (blk_off_x == WHOLE_BLK)
  {
    assert(blk_off_y == WHOLE_BLK); // just to be shure, may be removed later
    boff_x    = 0;
    boff_y    = 0;
    num_blk_x = ABT_NUMTR[abt_mode][PIX];
    num_blk_y = ABT_NUMTR[abt_mode][LIN];
  }else{
    boff_x    = blk_off_x;
    boff_y    = blk_off_y;
  }


  b8_x=(b8&1)<<3;
  b8_y=(b8&2)<<2;
  x1=boff_x+num_blk_x*ABT_TRSIZE[abt_mode][PIX];
  y1=boff_y+num_blk_y*ABT_TRSIZE[abt_mode][LIN];


  for(yy=boff_y;yy<y1;yy++)                //copy prediciton
    for(xx=boff_x;xx<x1;xx++)
      pred[yy][xx]=img->mpr[b8_x+xx][b8_y+yy];

  transform_ABT_B8(abt_mode,blk_off_x,blk_off_y,pred);        //transform and quant

  // quant and dequant
  for(yy=boff_y;yy<y1;yy++)
    for(xx=boff_x;xx<x1;xx++)
    {
      val          = pred[yy][xx];
      pred[yy][xx] = sign(
                      (((absm(val) * Q1_qpsp[yy&1][xx&1] + Qc_qpsp[yy&1][xx&1]) >> Qrs_qpsp[yy&1][xx&1]) * Q2_qpsp[yy&1][xx&1])<<q_shift,
                      val);
    }

  inv_transform_ABT_B8(abt_mode,blk_off_x,blk_off_y,pred);

  for(yy=boff_y;yy<y1;yy++)
    for(xx=boff_x;xx<x1;xx++)
    {
      val=(pred[yy][xx]+qp_const)>>q_bits;
      imgY[img->pix_y+b8_y+yy][img->pix_x+b8_x+xx]=clamp(val,0,255);
    }
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
      Qrshift[y][x] = ABT_N[abt_mode][0] - ABT_SHIFT0[abt_mode][y][x] + qp/QUANT_PERIOD;
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
      R[y][x] = ABT_R[abt_mode][qp%QUANT_PERIOD][ABT_QMAP[abt_mode][y][x]];
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
void intrapred_luma_ABT(int img_x,int img_y,int bs_x,int bs_y)
{
  unsigned char edgepixels[40];
#define EP (edgepixels+20)
  int x,y,sum,number,last_pix,new_pix;
  int b4_x,b4_y,b4_ox,b4_oy;
  int b4,b8,mb,ob4,ob8,omb;
  int i;
  int block_available_up,block_available_up_right;
  int block_available_left,block_available_left_down;

    //This intrapred uses stronger filtering for the larger blocks than the normal intraprediction
    assert( (bs_x==4||bs_x==8) && (bs_y==4||bs_y==8) ); //check blocksize
    assert( !(img_x&(bs_x-1)) && !(img_y&(bs_y-1)) );   //check position must be multiple of blocksize
    if( bs_x==4 && bs_y==4 )
    {
      intrapred_luma(img_x,img_y);
      return;
    }
#ifdef USE_6_INTRA_MODES
    assert(0);        //not supported (could be added)
#endif

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

    for(i=0;i<9;i++)
      img->mprr[i][0][0]=-1;

    // 0 DC
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
        img->mprr[DC_PRED][y][x]=(int)sum;

    // 1 vertical
    if(block_available_up)
      for(y=0UL;y<bs_y;y++)
        for(x=0UL;x<bs_x;x++)
          img->mprr[VERT_PRED][y][x]=EP[1+x];

    // 2 horizontal
    if(block_available_left)
      for(y=0UL;y<bs_y;y++)
        for(x=0UL;x<bs_x;x++)
          img->mprr[HOR_PRED][y][x]=EP[-1-y];

    // 3 down-right
    if(block_available_left&&block_available_up)
      for(y=0UL;y<bs_y;y++)
        for(x=0UL;x<bs_x;x++)
          img->mprr[DIAG_DOWN_RIGHT_PRED][y][x]=EP[x-y];

    // 4 up-right bidirectional
    if(block_available_left&&block_available_up)
      for(y=0UL;y<bs_y;y++)
        for(x=0UL;x<bs_x;x++)
          img->mprr[DIAG_DOWN_LEFT_PRED][y][x]=(EP[2+x+y]+EP[-2-(x+y)])>>1;

    // 5 down-right-down
    if(block_available_left&&block_available_up)
    {
      for(y=0UL;y<bs_y;y+=2)//even lines
        for(x=0UL;x<bs_x;x++)
          if( (i=x-(int)(y>>1)) >= 0 )
            img->mprr[VERT_RIGHT_PRED][y][x]=(EP[i]+EP[1+i])>>1;
          else
            img->mprr[VERT_RIGHT_PRED][y][x]=EP[1+2*(int)x-(int)y];
      for(y=1UL;y<bs_y;y+=2)//odd lines
        for(x=0UL;x<bs_x;x++)
          if( (i=x-(int)(y>>1)) >= 0 )
            img->mprr[VERT_RIGHT_PRED][y][x]=EP[i];
          else
            img->mprr[VERT_RIGHT_PRED][y][x]=EP[1+2*(int)x-(int)y];
    }

    // 6 down-left-down
    if(block_available_left&&block_available_up)
    {
      for(y=0UL;y<bs_y;y+=2)//even lines
        for(x=0UL;x<bs_x;x++)
          img->mprr[VERT_LEFT_PRED][y][x]=(EP[1+x+(y>>1)]+EP[2+x+(y>>1)])>>1;
      for(y=1UL;y<bs_y;y+=2)//odd lines
        for(x=0UL;x<bs_x;x++)
          img->mprr[VERT_LEFT_PRED][y][x]=EP[2+x+(y>>1)];
    }

    // 7 right-up-right
    if(block_available_left&&block_available_up)
    {
      for(y=0UL;y<bs_y;y++)//even columns
        for(x=0UL;x<bs_x;x+=2)
          img->mprr[HOR_UP_PRED][y][x]=(EP[-1-(int)(y+(x>>1))]+EP[-2-(int)(y+(x>>1))])>>1;
      for(y=0UL;y<bs_y;y++)//odd columns
        for(x=1UL;x<bs_x;x+=2)
          img->mprr[HOR_UP_PRED][y][x]=EP[-2-(int)(y+(x>>1))];
    }

    // 8 right-down-right
    if(block_available_left&&block_available_up)
    {
      for(y=0UL;y<bs_y;y++)//even columns
        for(x=0UL;x<bs_x;x+=2)
          if( (i=-(int)y+(x>>1)) <= 0 )
            img->mprr[HOR_DOWN_PRED][y][x]=(EP[i]+EP[i-1])>>1;
          else
            img->mprr[HOR_DOWN_PRED][y][x]=EP[-1-2*(int)y+(int)x];
      for(y=0UL;y<bs_y;y++)//odd columns
        for(x=1UL;x<bs_x;x+=2)
          if( (i=-(int)y+(x>>1)) <= 0 )
            img->mprr[HOR_DOWN_PRED][y][x]=EP[i];
          else
            img->mprr[HOR_DOWN_PRED][y][x]=EP[-1-2*(int)y+(int)x];
    }

#undef EP
}



//This might also be placed in rdopt.c behind Mode_Decision_for_4x4IntraBlocks().

/*!
 ************************************************************************
 * \brief
 *    Mode Decision for ABT intra blocks
 *    This might also be placed in rdopt.c behind Mode_Decision_for_4x4IntraBlocks().
 ************************************************************************
*/
int Mode_Decision_for_ABT_IntraBlocks(int b8,int b4,double lambda,int *min_cost,int bs_x,int bs_y)
{
  int ipmode,best_ipmode,i,j,x,y,cost,loc_cbp,loc_cbp_blk;
  int tmp_x,tmp_y;
  int ABTmode;
  int c_nz, nonzero, rec4x4[8][8], diff[16][16];
  double rdcost;
  int block_x;
  int block_y;
  int pic_pix_x;
  int pic_pix_y;
  int pic_block_x;
  int pic_block_y;
  int tmp_block_88[8][8];
  Macroblock *currMB;
  double min_rdcost ;
  int     upMode;
  int     leftMode;
  int     mostProbableMode;

  currMB=img->mb_data+img->current_mb_nr;

  block_x    =8*(b8&1)+4*(b4&1);
  block_y    =8*(b8>>1)+4*(b4>>1);
  pic_pix_x  =img->pix_x+block_x;
  pic_pix_y  =img->pix_y+block_y;
  pic_block_x=pic_pix_x>>2;
  pic_block_y=pic_pix_y>>2;
  min_rdcost =1e30;
  ABTmode=(bs_x>4?( bs_y>4 ? B8x8 : B8x4 ):( bs_y>4 ? B4x8 : B4x4 ));
  currMB->abt_mode[b8]=ABTmode;

  *min_cost = (1<<20);

  upMode           = img->ipredmode[pic_block_x+1][pic_block_y  ];
  leftMode         = img->ipredmode[pic_block_x  ][pic_block_y+1];
  mostProbableMode = (upMode < 0 || leftMode < 0) ? DC_PRED : mapTab[upMode] < mapTab[leftMode] ? upMode : leftMode;


  //===== INTRA PREDICTION FOR 4x4 BLOCK =====
  intrapred_luma_ABT(pic_pix_x,pic_pix_y,bs_x,bs_y);


  //===== LOOP OVER ALL INTRA PREDICTION MODES =====
  for (ipmode=0;ipmode<NO_INTRA_PMODE;ipmode++)
  {
    if(img->mprr[ipmode][0][0]>=0)        //changed this. the intrapred function marks the invalid modes. At least one is always valid (the DC).
    {
      if(!input->rdopt)
      {
        for (j=0;j<bs_y;j++)
          for (i=0;i<bs_x;i++)
            diff[j][i] = imgY_org[pic_pix_y+j][pic_pix_x+i] - img->mprr[ipmode][j][i]; // bug fix: index of diff was flipped. mwi 020701
        cost  = (ipmode == mostProbableMode) ? 0 : (int)floor(4 * lambda );

        cost+=find_sad_abt(input->hadamard,bs_x,bs_y,0,0,diff);

        if (cost < *min_cost)
        {
          best_ipmode = ipmode;
          *min_cost   = cost;
        }
      }else{
        // get prediction and prediction error
        for (j=0;j<bs_y;j++)
          for (i=0;i<bs_x;i++)
          {
            img->mpr[block_x+i][block_y+j]=img->mprr[ipmode][j][i];
            img->m7[i][j] = imgY_org[pic_pix_y+j][pic_pix_x+i] - img->mpr[block_x+i][block_y+j] ;
          }

        //===== store the coding state =====
        store_coding_state (cs_cm);

        // get and check rate-distortion cost
        if ((rdcost = RDCost_for_ABTIntraBlocks(&c_nz,b8,b4,ipmode,lambda,min_rdcost,bs_x,bs_y, mostProbableMode)) < min_rdcost)
        {
          //--- set coefficients ---
          for(j=0;j<2;j++)
            for(i=0;i<65;i++)
              cofAC4x4[j][i]=img->cofAC[b8][b4][j][i];

          //--- set reconstruction ---
          for(y=0; y<bs_y; y++)
            for(x=0; x<bs_y; x++)
              rec4x4[y][x]=imgY[pic_pix_y+y][pic_pix_x+x];

          //--- flag if dct-coefficients must be coded ---
          nonzero = c_nz;

          //--- set best mode update minimum cost ---
          min_rdcost  = rdcost;
          *min_cost=(int)min_rdcost;
          best_ipmode = ipmode;
        }
        reset_coding_state (cs_cm);
      }
    }
  }

//  assert( *min_cost < (1<<20) ); //make shure we found one.

  //===== set intra mode prediction =====
  for(j=0;j<(bs_y>>2);j++)
    for(i=0;i<(bs_x>>2);i++)
    {
      img->ipredmode[pic_block_x+i+1][pic_block_y+j+1] = best_ipmode;
      currMB->intra_pred_modes[(b8<<2)+b4+(j<<1)+i]=0;
    }

  currMB->intra_pred_modes[(b8<<2)+b4] = mostProbableMode == best_ipmode ? -1 : mapTab[best_ipmode] < mapTab[mostProbableMode] ? mapTab[best_ipmode] : mapTab[best_ipmode]-1;

  // get prediction and prediction error
  tmp_x=(b4&1)<<2;tmp_y=(b4&2)<<1;
  for (j=0;j<bs_y; j++)
    for (i=0;i<bs_x; i++)
    {
      img->mpr[block_x+i][block_y+j] = img->mprr[best_ipmode][j][i];
      tmp_block_88[tmp_y+j][tmp_x+i] = imgY_org[pic_pix_y+j][pic_pix_x+i] - img->mprr[best_ipmode][j][i];
    }
  //nonzero = dct_luma (block_x, block_y, &dummy, 1);
  transform_ABT_B8(ABTmode,(b4&1)<<2,(b4&2)<<1,tmp_block_88);
  loc_cbp=loc_cbp_blk=0;
  scanquant_ABT_B8(img->qp+QP_OFS-MIN_QP,ABTmode|4,b8,(b4&1)<<2,(b4&2)<<1,tmp_block_88,0,&loc_cbp,&loc_cbp_blk); // '|4' indicate intra for quantization
  //scanquant returns a SCR value!.....
  nonzero=(loc_cbp!=0);
  currMB->cbp|=loc_cbp;
  currMB->cbp_blk|=loc_cbp_blk;
  if(nonzero)nonzero=1;

  return nonzero;
}

/*!
 ************************************************************************
 * \brief
 *    Get RD cost for ABT intra block
 ************************************************************************
*/
double RDCost_for_ABTIntraBlocks(int *nonzero,
                                 int b8,int b4,
                                 int ipmode,
                                 double lambda,
                                 double  min_rdcost,
                                 int bs_x,int bs_y,
                                 int mostProbableMode)
{
 int x,y,rate,ABTmode,tmp_cbp,tmp_cbp_blk;
 int tmp_x,tmp_y;
 int distortion;
 int block_x;
 int block_y;
 int pic_pix_x;
 int pic_pix_y;
 int pic_block_x;
 int pic_block_y;
 int even_block;
 Slice   *currSlice;
 Macroblock  *currMB;
 SyntaxElement *currSE;
 const int   *partMap;
 DataPartition *dataPart;
 int tmp_block_88_inv[8][8];

  block_x    =8*(b8%2)+4*(b4%2);
  block_y    =8*(b8/2)+4*(b4/2);
  pic_pix_x  =img->pix_x+block_x;
  pic_pix_y  =img->pix_y+block_y;
  pic_block_x=pic_pix_x/4;
  pic_block_y=pic_pix_y/4;
  even_block =b4%2;
  currSlice=img->currentSlice;
  currMB   =&img->mb_data[img->current_mb_nr];
  currSE   =&img->MB_SyntaxElements[currMB->currSEnr];
  partMap  =assignSE2partition[input->partition_mode];


  //===== perform DCT, Q, IQ, IDCT, Reconstruction =====
  ABTmode=(bs_x>4?( bs_y>4 ? B8x8 : B8x4 ):( bs_y>4 ? B4x8 : B4x4 ));
  tmp_x=(b4&1)<<2;tmp_y=(b4&2)<<1;
  for(y=0;y<bs_y;y++)
    for(x=0;x<bs_x;x++)
      tmp_block_88_inv[tmp_y+y][tmp_x+x]=img->m7[x][y];        //the subblock to be processed is in the top left corner of img->m7[][].
  transform_ABT_B8(ABTmode,(b4&1)<<2,(b4&2)<<1,tmp_block_88_inv);
  tmp_cbp=tmp_cbp_blk=0;
  scanquant_ABT_B8(img->qp+QP_OFS-MIN_QP,ABTmode|4,b8,(b4&1)<<2,(b4&2)<<1,tmp_block_88_inv,0,&tmp_cbp,&tmp_cbp_blk); // '|4' indicate intra for quantization
  //scanquant returns a SCR value!.....
  *nonzero=(tmp_cbp!=0);

  //===== get distortion (SSD) of 4x4 block =====
  distortion =0;
  for (y=pic_pix_y; y<pic_pix_y+bs_y; y++)
    for (x=pic_pix_x; x<pic_pix_x+bs_x; x++)
      distortion += img->quad [imgY_org[y][x] - imgY[y][x]];

  //===== RATE for INTRA PREDICTION MODE  (SYMBOL MODE MUST BE SET TO UVLC) =====
  currSE->value1 = (mostProbableMode == ipmode) ? -1 : ipmode < mostProbableMode ? ipmode : ipmode-1;

  //--- set position and type ---
  currSE->context = 4*b8 + b4;
  currSE->type    = SE_INTRAPREDMODE;

  //--- set function pointer ----
  if (input->symbol_mode != UVLC)    currSE->writing = writeIntraPredMode2Buffer_CABAC;

  //--- choose data partition ---
  if (img->type != B_IMG && img->type != BS_IMG)   dataPart = &(currSlice->partArr[partMap[currSE->type]]);
  else                    dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);

  //--- encode and update rate ---
  if (input->symbol_mode == UVLC)    writeSyntaxElement_Intra4x4PredictionMode(currSE, dataPart);
  else                               dataPart->writeSyntaxElement (currSE, dataPart);
  rate = currSE->len;
  currSE++;
  currMB->currSEnr++;

  //===== RATE for LUMINANCE COEFFICIENTS =====
  tmp_x=(b4&1)<<2;tmp_y=(b4&2)<<1;
    x=currMB->cbp;
  currMB->cbp=tmp_cbp;//writeLumaCoeffABT_B8 needs a correct Macroblock->cbp .
  rate+=writeLumaCoeffABT_B8(b8,1,tmp_x,tmp_y);
  currMB->cbp=x;

  //calc RD and return it.
  return (double)distortion+lambda*(double)rate;
}


