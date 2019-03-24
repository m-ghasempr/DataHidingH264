/*
***********************************************************************
* COPYRIGHT AND WARRANTY INFORMATION
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

/*!
 *************************************************************************************
 * \file b_frame.c
 *
 * \brief
 *    B picture decoding
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Byeong-Moon Jeon                <jeonbm@lge.com>
 *    - Yoon-Seong Soh                  <yunsung@lge.com>
 *    - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
 *************************************************************************************
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "global.h"
#include "mbuffer.h"
#include "b_frame.h"
#include "elements.h"

#define POS 0

/*!
 ************************************************************************
 * \brief
 *    Write previous decoded P frame to output file
 ************************************************************************
 */
void write_prev_Pframe(struct img_par *img, FILE *p_out)
{
  int i,j;

  for(i=0;i<img->height;i++)
    for(j=0;j<img->width;j++)
      fputc(imgY_prev[i][j],p_out);

  for(i=0;i<img->height_cr;i++)
    for(j=0;j<img->width_cr;j++)
      fputc(imgUV_prev[0][i][j],p_out);

  for(i=0;i<img->height_cr;i++)
    for(j=0;j<img->width_cr;j++)
      fputc(imgUV_prev[1][i][j],p_out);
}



/*!
 ************************************************************************
 * \brief
 *    Copy decoded P frame to temporary image array
 ************************************************************************
 */
void copy_Pframe(struct img_par *img, int postfilter)
{
  int i,j;

  /*
   * the mmin, mmax macros are taken out, because it makes no sense due to limited range of data type
   */

  if(postfilter)
  {
    for(i=0;i<img->height;i++)
      for(j=0;j<img->width;j++)
      {
        imgY_prev[i][j] = imgY_pf[i][j];
      }
    for(i=0;i<img->height_cr;i++)
      for(j=0;j<img->width_cr;j++)
      {
        imgUV_prev[0][i][j] = imgUV_pf[0][i][j];
      }
    for(i=0;i<img->height_cr;i++)
      for(j=0;j<img->width_cr;j++)
      {
        imgUV_prev[1][i][j] = imgUV_pf[1][i][j];
      }
  }
  else
  {
    for(i=0;i<img->height;i++)
      for(j=0;j<img->width;j++)
      {
        imgY_prev[i][j] = imgY[i][j];
      }
    for(i=0;i<img->height_cr;i++)
      for(j=0;j<img->width_cr;j++)
      {
        imgUV_prev[0][i][j] = imgUV[0][i][j];
      }
    for(i=0;i<img->height_cr;i++)
      for(j=0;j<img->width_cr;j++)
      {
        imgUV_prev[1][i][j] = imgUV[1][i][j];
      }
  }
}




/*!
 ************************************************************************
 * \brief
 *    init macroblock B frames
 ************************************************************************
 */
void init_macroblock_Bframe(struct img_par *img)
{
  int i,j;
  int fw_predframe_no=0;
    Macroblock *currMB = &img->mb_data[img->current_mb_nr];

    currMB->ref_frame = 0;
    currMB->predframe_no = 0;


  // reset vectors and pred. modes
  for (i=0;i<BLOCK_SIZE;i++)
  {
    for(j=0;j<BLOCK_SIZE;j++)
    {
      img->fw_mv[img->block_x+i+4][img->block_y+j][0]=img->fw_mv[img->block_x+i+4][img->block_y+j][1]=0;
      img->bw_mv[img->block_x+i+4][img->block_y+j][0]=img->bw_mv[img->block_x+i+4][img->block_y+j][1]=0;
      img->dfMV[img->block_x+i+4][img->block_y+j][0]=img->dfMV[img->block_x+i+4][img->block_y+j][1]=0;
      img->dbMV[img->block_x+i+4][img->block_y+j][0]=img->dbMV[img->block_x+i+4][img->block_y+j][1]=0;
      img->ipredmode[img->block_x+i+1][img->block_y+j+1] = 0;
    }
  }

  // Set the reference frame information for motion vector prediction
  if (img->imod==INTRA_MB_OLD || img->imod==INTRA_MB_NEW)
  {
    for (j = 0;j < BLOCK_SIZE;j++)
    {
      for (i = 0;i < BLOCK_SIZE;i++)
      {
        img->fw_refFrArr[img->block_y+j][img->block_x+i] =
        img->bw_refFrArr[img->block_y+j][img->block_x+i] = -1;
      }
    }
  }
  else if(img->imod == B_Forward)
  {
    for (j = 0;j < BLOCK_SIZE;j++)
    {
      for (i = 0;i < BLOCK_SIZE;i++)
      {
        img->fw_refFrArr[img->block_y+j][img->block_x+i] = fw_predframe_no; // previous P
        img->bw_refFrArr[img->block_y+j][img->block_x+i] = -1;
      }
    }
  }
  else if(img->imod == B_Backward)
  {
    for (j = 0;j < BLOCK_SIZE;j++)
    {
      for (i = 0;i < BLOCK_SIZE;i++)
      {
        img->fw_refFrArr[img->block_y+j][img->block_x+i] = -1;
        img->bw_refFrArr[img->block_y+j][img->block_x+i] = 0; // next P
      }
    }
  }
  else if(img->imod == B_Bidirect)
  {
    for (j = 0;j < BLOCK_SIZE;j++)
    {
      for (i = 0;i < BLOCK_SIZE;i++)
      {
        img->fw_refFrArr[img->block_y+j][img->block_x+i] = fw_predframe_no; // previous P
        img->bw_refFrArr[img->block_y+j][img->block_x+i] = 0; // next P
      }
    }
  }
  else if(img->imod == B_Direct)
  {
    for (j = 0;j < BLOCK_SIZE;j++)
    {
      for (i = 0;i < BLOCK_SIZE;i++)
      {
        img->fw_refFrArr[img->block_y+j][img->block_x+i] = -1; // previous P
        img->bw_refFrArr[img->block_y+j][img->block_x+i] = -1; // next P
      }
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Get for a given MB of a B picture the reference frame
 *    parameter and the motion vectors from the NAL
 *
 ************************************************************************
 */
void readMotionInfoFromNAL_Bframe(struct img_par *img,struct inp_par *inp)
{
  int i,j,ii,jj,ie,j4,i4,k,pred_vec=0;
  int vec;
  int ref_frame;
  int fw_predframe_no=0, bw_predframe_no;     // frame number which current MB is predicted from
  int fw_blocktype=0, bw_blocktype=0, blocktype=0;
  int l,m;

  // keep track of neighbouring macroblocks available for prediction
  int mb_nr = img->current_mb_nr;
  int mb_width = img->width/16;
  int mb_available_up = (img->mb_y == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width].slice_nr);
  int mb_available_left = (img->mb_x == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-1].slice_nr);
  int mb_available_upleft  = (img->mb_x == 0 || img->mb_y == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width-1].slice_nr);
  int mb_available_upright = (img->mb_x >= mb_width-1 || img->mb_y == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width+1].slice_nr);

  // keep track of neighbouring blocks available for motion vector prediction
  int block_available_up, block_available_left, block_available_upright, block_available_upleft;
  int mv_a, mv_b, mv_c, mv_d;
  int mvPredType, rFrameL, rFrameU, rFrameUR;

  SyntaxElement THEcurrSE;
  SyntaxElement *currSE = &THEcurrSE;

  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  Slice *currSlice = img->currentSlice;
  DataPartition *dp;
  int *partMap = assignSE2partition[currSlice->dp_mode];

  bw_predframe_no=0;


  // /////////////////////////////////////////////////////////////
  // find fw_predfame_no
  // /////////////////////////////////////////////////////////////
  if(img->imod==B_Forward || img->imod==B_Bidirect)
  {
    if(img->type==B_IMG_MULT)
    {

#if TRACE
      snprintf(currSE->tracestring, TRACESTRING_SIZE, "B_fw_Reference frame no  ");
#endif
      currSE->type = SE_REFFRAME;
      dp = &(currSlice->partArr[partMap[SE_BFRAME]]);

      if (inp->symbol_mode == UVLC)
        currSE->mapping = linfo;
      else
        currSE->reading = readRefFrameFromBuffer_CABAC;

      dp->readSyntaxElement(  currSE, img, inp, dp);
      fw_predframe_no = currSE->value1;


      if (fw_predframe_no > img->buf_cycle)
      {
        set_ec_flag(SE_REFFRAME);
        fw_predframe_no = 1;
      }

      ref_frame=fw_predframe_no;

 /*     if (ref_frame > img->number)
      {
        ref_frame = 0;
      }*/

      currMB->predframe_no = fw_predframe_no;
      currMB->ref_frame = img->ref_frame = ref_frame;

      // Update the reference frame information for motion vector prediction
      for (j = 0;j < BLOCK_SIZE;j++)
      {
        for (i = 0;i < BLOCK_SIZE;i++)
        {
          img->fw_refFrArr[img->block_y+j][img->block_x+i] = fw_predframe_no;
        }
      }
    } 
    else
      fw_predframe_no = currMB->predframe_no;
  }

  // /////////////////////////////////////////////////////////////
  // find blk_size
  // /////////////////////////////////////////////////////////////
  if(img->imod==B_Bidirect)
  {

#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "B_bidiret_fw_blk");
#endif

    currSE->type = SE_BFRAME;

    dp = &(currSlice->partArr[partMap[SE_BFRAME]]);

    if (inp->symbol_mode == UVLC)
      currSE->mapping = linfo;
    else
      currSE->reading = readBiDirBlkSize2Buffer_CABAC;

    dp->readSyntaxElement(  currSE, img, inp, dp);

    fw_blocktype = currSE->value1+1;


#if TRACE
    snprintf(currSE->tracestring, TRACESTRING_SIZE, "B_bidiret_bw_blk");
#endif

    currSE->type = SE_BFRAME;

    dp = &(currSlice->partArr[partMap[SE_BFRAME]]);

    if (inp->symbol_mode == UVLC)
      currSE->mapping = linfo;
    else
      currSE->reading = readBiDirBlkSize2Buffer_CABAC;


    dp->readSyntaxElement(  currSE, img, inp, dp);

    bw_blocktype = currSE->value1+1;

  }

  // /////////////////////////////////////////////////////////////
  // find MVDFW
  // /////////////////////////////////////////////////////////////
  if(img->imod==B_Forward || img->imod==B_Bidirect)
  {
    // forward : note realtion between blocktype and img->mb_mode
    // bidirect : after reading fw_blk_size, fw_pmv is obtained
    if(img->mb_mode==1)
      blocktype=1;
    else if(img->mb_mode>3)
      blocktype=img->mb_mode/2;
    else if(img->mb_mode==3)
      blocktype=fw_blocktype;

    img->fw_blc_size_h = BLOCK_STEP[blocktype][0]*4;
    img->fw_blc_size_v = BLOCK_STEP[blocktype][1]*4;

    ie=4-BLOCK_STEP[blocktype][0];
    for(j=0;j<4;j=j+BLOCK_STEP[blocktype][1])   // Y position in 4-pel resolution
    {
      block_available_up = mb_available_up || (j > 0);
      j4=img->block_y+j;
      for(i=0;i<4;i=i+BLOCK_STEP[blocktype][0]) // X position in 4-pel resolution
      {
        i4=img->block_x+i;
        block_available_left = mb_available_left || (i > 0);

        if (j > 0)
          block_available_upright = i != ie ? 1 : 0;
        else if (i != ie)
          block_available_upright = block_available_up;
        else
          block_available_upright = mb_available_upright;

        if (i > 0)
          block_available_upleft = j > 0 ? 1 : block_available_up;
        else if (j > 0)
          block_available_upleft = block_available_left;
        else
          block_available_upleft = mb_available_upleft;

        mvPredType = MVPRED_MEDIAN;

        rFrameL    = block_available_left    ? img->fw_refFrArr[j4][i4-1]   : -1;
        rFrameU    = block_available_up      ? img->fw_refFrArr[j4-1][i4]   : -1;
        rFrameUR   = block_available_upright ? img->fw_refFrArr[j4-1][i4+BLOCK_STEP[blocktype][0]] :
                     block_available_upleft  ? img->fw_refFrArr[j4-1][i4-1] : -1;

        // Prediction if only one of the neighbors uses the selected reference frame
        if(rFrameL == fw_predframe_no && rFrameU != fw_predframe_no && rFrameUR != fw_predframe_no)
          mvPredType = MVPRED_L;
        else if(rFrameL != fw_predframe_no && rFrameU == fw_predframe_no && rFrameUR != fw_predframe_no)
          mvPredType = MVPRED_U;
        else if(rFrameL != fw_predframe_no && rFrameU != fw_predframe_no && rFrameUR == fw_predframe_no)
          mvPredType = MVPRED_UR;

        // Directional predictions
        else if(blocktype == 3)
        {
          if(i == 0)
          {
            if(rFrameL == fw_predframe_no)
              mvPredType = MVPRED_L;
          }
          else
          {
            if(rFrameUR == fw_predframe_no)
              mvPredType = MVPRED_UR;
          }
        }
        else if(blocktype == 2)
        {
          if(j == 0)
          {
            if(rFrameU == fw_predframe_no)
              mvPredType = MVPRED_U;
          }
          else
          {
            if(rFrameL == fw_predframe_no)
              mvPredType = MVPRED_L;
          }
        }
        else if(blocktype == 5 && i == 2)
          mvPredType = MVPRED_L;
        else if(blocktype == 6 && j == 2)
          mvPredType = MVPRED_U;

        for (k=0;k<2;k++)                       // vector prediction. find first koordinate set for pred
        {
          mv_a = block_available_left ? img->fw_mv[i4-1+BLOCK_SIZE][j4][k] : 0;
          mv_b = block_available_up   ? img->fw_mv[i4+BLOCK_SIZE][j4-1][k]   : 0;
          mv_d = block_available_upleft  ? img->fw_mv[i4-1+BLOCK_SIZE][j4-1][k] : 0;
          mv_c = block_available_upright ? img->fw_mv[i4+BLOCK_STEP[blocktype][0]+BLOCK_SIZE][j4-1][k] : mv_d;

          switch (mvPredType)
          {
          case MVPRED_MEDIAN:
            if(!(block_available_upleft || block_available_up || block_available_upright))
              pred_vec = mv_a;
            else
              pred_vec =mv_a+mv_b+mv_c-min(mv_a,min(mv_b,mv_c))-max(mv_a,max(mv_b,mv_c));
            break;
          case MVPRED_L:
            pred_vec = mv_a;
            break;
          case MVPRED_U:
            pred_vec = mv_b;
            break;
          case MVPRED_UR:
            pred_vec = mv_c;
            break;
          default:
            break;
          }

#if TRACE
          snprintf(currSE->tracestring, TRACESTRING_SIZE, " MVD(%d) ",k);
#endif

          currSE->type = SE_MVD;
          dp = &(currSlice->partArr[partMap[SE_BFRAME]]);

          if (inp->symbol_mode == UVLC)
            currSE->mapping = linfo_mvd; // linfo_2
          else
          {
            img->subblock_x = i; // position used for context determination
            img->subblock_y = j; // position used for context determination
                currSE->value2 = 2*k; // identifies the component and the direction; only used for context determination
                        currSE->reading = readBiMVD2Buffer_CABAC;
          }

          dp->readSyntaxElement(  currSE, img, inp, dp);

          vec = currSE->value1;

          for (l=0; l < BLOCK_STEP[blocktype][1]; l++)
            for (m=0; m < BLOCK_STEP[blocktype][0]; m++)
              currMB->mvd[0][j+l][i+m][k] =  vec;

          vec=vec+pred_vec;

          for(ii=0;ii<BLOCK_STEP[blocktype][0];ii++)
          {
            for(jj=0;jj<BLOCK_STEP[blocktype][1];jj++)
            {
              img->fw_mv[i4+ii+BLOCK_SIZE][j4+jj][k]=vec;
            }
          }

        }
      }
    }
  }

  // /////////////////////////////////////////////////////////////
  // find MVDBW
  // /////////////////////////////////////////////////////////////
  if(img->imod==B_Backward || img->imod==B_Bidirect)
  {
    // backward : note realtion between blocktype and img->mb_mode
    // bidirect : after reading bw_blk_size, bw_pmv is obtained

    if(img->mb_mode==2)
            blocktype=1;
    else if(img->mb_mode>4)
            blocktype=(img->mb_mode-1)/2;
    else if(img->mb_mode==3)
            blocktype=bw_blocktype;

        img->bw_blc_size_h = BLOCK_STEP[blocktype][0]*4;
        img->bw_blc_size_v = BLOCK_STEP[blocktype][1]*4;

    ie=4-BLOCK_STEP[blocktype][0];
    for(j=0;j<4;j=j+BLOCK_STEP[blocktype][1])   // Y position in 4-pel resolution
    {
      block_available_up = mb_available_up || (j > 0);
      j4=img->block_y+j;
      for(i=0;i<4;i=i+BLOCK_STEP[blocktype][0]) // X position in 4-pel resolution
      {
        i4=img->block_x+i;
        block_available_left = mb_available_left || (i > 0);

        if (j > 0)
          block_available_upright = i != ie ? 1 : 0;
        else if (i != ie)
          block_available_upright = block_available_up;
        else
          block_available_upright = mb_available_upright;

        if (i > 0)
          block_available_upleft = j > 0 ? 1 : block_available_up;
        else if (j > 0)
          block_available_upleft = block_available_left;
        else
          block_available_upleft = mb_available_upleft;

        mvPredType = MVPRED_MEDIAN;

        rFrameL    = block_available_left    ? img->bw_refFrArr[j4][i4-1]   : -1;
        rFrameU    = block_available_up      ? img->bw_refFrArr[j4-1][i4]   : -1;
        rFrameUR   = block_available_upright ? img->bw_refFrArr[j4-1][i4+BLOCK_STEP[blocktype][0]] :
                     block_available_upleft  ? img->bw_refFrArr[j4-1][i4-1] : -1;

        // Prediction if only one of the neighbors uses the selected reference frame
        if(rFrameL == bw_predframe_no && rFrameU != bw_predframe_no && rFrameUR != bw_predframe_no)
          mvPredType = MVPRED_L;
        else if(rFrameL != bw_predframe_no && rFrameU == bw_predframe_no && rFrameUR != bw_predframe_no)
          mvPredType = MVPRED_U;
        else if(rFrameL != bw_predframe_no && rFrameU != bw_predframe_no && rFrameUR == bw_predframe_no)
          mvPredType = MVPRED_UR;

        // Directional predictions
        else if(blocktype == 3)
        {
          if(i == 0)
          {
            if(rFrameL == bw_predframe_no)
              mvPredType = MVPRED_L;
          }
          else
          {
            if(rFrameUR == bw_predframe_no)
              mvPredType = MVPRED_UR;
          }
        }
        else if(blocktype == 2)
        {
          if(j == 0)
          {
            if(rFrameU == bw_predframe_no)
              mvPredType = MVPRED_U;
          }
          else
          {
            if(rFrameL == bw_predframe_no)
              mvPredType = MVPRED_L;
          }
        }
        else if(blocktype == 5 && i == 2)
          mvPredType = MVPRED_L;
        else if(blocktype == 6 && j == 2)
          mvPredType = MVPRED_U;

        for (k=0;k<2;k++)                       // vector prediction. find first koordinate set for pred
        {
          mv_a = block_available_left ? img->bw_mv[i4-1+BLOCK_SIZE][j4][k] : 0;
          mv_b = block_available_up   ? img->bw_mv[i4+BLOCK_SIZE][j4-1][k]   : 0;
          mv_d = block_available_upleft  ? img->bw_mv[i4-1+BLOCK_SIZE][j4-1][k] : 0;
          mv_c = block_available_upright ? img->bw_mv[i4+BLOCK_STEP[blocktype][0]+BLOCK_SIZE][j4-1][k] : mv_d;

          switch (mvPredType)
          {
          case MVPRED_MEDIAN:
            if(!(block_available_upleft || block_available_up || block_available_upright))
              pred_vec = mv_a;
            else
              pred_vec =mv_a+mv_b+mv_c-min(mv_a,min(mv_b,mv_c))-max(mv_a,max(mv_b,mv_c));
            break;
          case MVPRED_L:
            pred_vec = mv_a;
            break;
          case MVPRED_U:
            pred_vec = mv_b;
            break;
          case MVPRED_UR:
            pred_vec = mv_c;
            break;
          default:
            break;
          }

#if TRACE
          snprintf(currSE->tracestring, TRACESTRING_SIZE, " MVD(%d) ",k);
#endif

          currSE->type = SE_MVD;
          dp = &(currSlice->partArr[partMap[SE_BFRAME]]);

          if (inp->symbol_mode == UVLC)
            currSE->mapping = linfo_mvd; // linfo_2
          else
          {
            img->subblock_x = i; // position used for context determination
            img->subblock_y = j; // position used for context determination
            currSE->value2 = 2*k+1; // identifies the component and the direction; only used for context determination
                        currSE->reading = readBiMVD2Buffer_CABAC;
          }

          dp->readSyntaxElement(  currSE, img, inp, dp);

          vec = currSE->value1;

          for (l=0; l < BLOCK_STEP[blocktype][1]; l++)
            for (m=0; m < BLOCK_STEP[blocktype][0]; m++)
              currMB->mvd[1][j+l][i+m][k] =  vec;

          vec=vec+pred_vec;

          for(ii=0;ii<BLOCK_STEP[blocktype][0];ii++)
          {
            for(jj=0;jj<BLOCK_STEP[blocktype][1];jj++)
            {
              img->bw_mv[i4+ii+BLOCK_SIZE][j4+jj][k]=vec;
            }
          }

        }
      }
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    decode one macroblock for a B frame
 ************************************************************************
 */
int decode_one_macroblock_Bframe(struct img_par *img)
{
  int tmp_block[BLOCK_SIZE][BLOCK_SIZE];
  int tmp_blockbw[BLOCK_SIZE][BLOCK_SIZE];
  int js[2][2];
  int i,j,ii,jj,i1,j1,j4,i4;
  int js0,js1,js2,js3,jf,ifx;
  int uv;
  int vec1_x,vec1_y,vec2_x,vec2_y, vec1_xx, vec1_yy;
  int fw_pred, bw_pred;

  int ioff,joff;
  int ii0,jj0,ii1,jj1,if1,jf1,if0,jf0;

  int hv;
  byte refP_tr, TRb, TRp;
  // keep track of neighbouring macroblocks available for prediction

  int mb_nr = img->current_mb_nr;
  int mb_width = img->width/16;
  int mb_available_up = (img->mb_y == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width].slice_nr);
  int mb_available_left = (img->mb_x == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-1].slice_nr);
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  int ref_frame_bw = 0;
  int ref_frame_fw = currMB->ref_frame+1;
  int ref_frame = ref_frame_fw;
  int mv_mul, f1, f2, f3, f4;



#ifdef _ADAPT_LAST_GROUP_
  extern int *last_P_no;
#endif

#if POS
  int color;
  if( mb_nr == 122)
    color = 0;
  else
    color = 0xff;
#endif

  // set variables depending on mv_res
  if(img->mv_res)
  {
    mv_mul=8;
    f1=16;
    f2=15;
  }
  else
  {
    mv_mul=4;
    f1=8;
    f2=7;
  }

  f3=f1*f1;
  f4=f3/2;

  // ////////////////////////////////////////////////////////////
  // start_scan decoding
  // ////////////////////////////////////////////////////////////
  /**********************************
   *    luma                        *
   *********************************/

  if (img->imod==INTRA_MB_NEW)
    intrapred_luma_2( img, currMB->intra_pred_modes[0]);

  for(j=0;j<MB_BLOCK_SIZE/BLOCK_SIZE;j++)
  {
    joff=j*4;
    j4=img->block_y+j;
    for(i=0;i<MB_BLOCK_SIZE/BLOCK_SIZE;i++)
    {
      ioff=i*4;
      i4=img->block_x+i;
      if(img->imod==INTRA_MB_OLD)
      {
        if (intrapred(img,ioff,joff,i4,j4)==SEARCH_SYNC)  // make 4x4 prediction block mpr from given prediction img->mb_mode
        return SEARCH_SYNC;                   // bit error
      }

      // //////////////////////////////////
      // Forward MC using img->fw_MV
      // //////////////////////////////////
      else if(img->imod==B_Forward)
      {
        vec1_x = i4*4*mv_mul + img->fw_mv[i4+BLOCK_SIZE][j4][0];
        vec1_y = j4*4*mv_mul + img->fw_mv[i4+BLOCK_SIZE][j4][1];

        get_block(ref_frame_fw,vec1_x,vec1_y,img,tmp_block);

        for(ii=0;ii<BLOCK_SIZE;ii++)
          for(jj=0;jj<BLOCK_SIZE;jj++)
            img->mpr[ii+ioff][jj+joff] = tmp_block[ii][jj];
      }

      // //////////////////////////////////
      // Backward MC using img->bw_MV
      // //////////////////////////////////
      else if(img->imod==B_Backward)
      {
        vec1_x = i4*4*mv_mul + img->bw_mv[i4+BLOCK_SIZE][j4][0];
        vec1_y = j4*4*mv_mul + img->bw_mv[i4+BLOCK_SIZE][j4][1];

        get_block(ref_frame_bw,vec1_x,vec1_y,img,tmp_blockbw);

        for(ii=0;ii<BLOCK_SIZE;ii++)
          for(jj=0;jj<BLOCK_SIZE;jj++)
            img->mpr[ii+ioff][jj+joff] = tmp_blockbw[ii][jj];
      }

      // //////////////////////////////////
      // Bidirect MC using img->fw_MV, bw_MV
      // //////////////////////////////////
      else if(img->imod==B_Bidirect)
      {
        vec2_x=i4*4*mv_mul;
        vec1_x=vec2_x+img->fw_mv[i4+BLOCK_SIZE][j4][0];
        vec1_xx=vec2_x+img->bw_mv[i4+BLOCK_SIZE][j4][0];
        vec2_y=j4*4*mv_mul;
        vec1_y=vec2_y+img->fw_mv[i4+BLOCK_SIZE][j4][1];
        vec1_yy=vec2_y+img->bw_mv[i4+BLOCK_SIZE][j4][1];

        get_block(ref_frame_fw,vec1_x,vec1_y,img,tmp_block);
        get_block(ref_frame_bw,vec1_xx,vec1_yy,img,tmp_blockbw);

        for(ii=0;ii<BLOCK_SIZE;ii++)
          for(jj=0;jj<BLOCK_SIZE;jj++)
            img->mpr[ii+ioff][jj+joff] = (tmp_block[ii][jj]+tmp_blockbw[ii][jj]+1)/2;
      }

      // //////////////////////////////////
      // Direct MC using img->mv
      // //////////////////////////////////
      else if(img->imod==B_Direct)
      {
        // next P is intra mode
        if(refFrArr[j4][i4]==-1)
        {
          for(hv=0; hv<2; hv++)
          {
            img->dfMV[i4+BLOCK_SIZE][j4][hv]=img->dbMV[i4+BLOCK_SIZE][j4][hv]=0;
          }
          ref_frame = 1;
        }
        // next P is skip or inter mode
        else
        {
#ifdef _ADAPT_LAST_GROUP_
          refP_tr = last_P_no[refFrArr[j4][i4]];
#else
          refP_tr = nextP_tr-((refFrArr[j4][i4]+1)*P_interval);
#endif
          TRb = img->tr-refP_tr;
          TRp = nextP_tr-refP_tr;

          img->dfMV[i4+BLOCK_SIZE][j4][0]=TRb*img->mv[i4+BLOCK_SIZE][j4][0]/TRp;
          img->dfMV[i4+BLOCK_SIZE][j4][1]=TRb*img->mv[i4+BLOCK_SIZE][j4][1]/TRp;
          img->dbMV[i4+BLOCK_SIZE][j4][0]=(TRb-TRp)*img->mv[i4+BLOCK_SIZE][j4][0]/TRp;
          img->dbMV[i4+BLOCK_SIZE][j4][1]=(TRb-TRp)*img->mv[i4+BLOCK_SIZE][j4][1]/TRp;
          ref_frame=1 + refFrArr[j4][i4];
        }

        vec2_x=i4*4*mv_mul;
        vec1_x=vec2_x+img->dfMV[i4+BLOCK_SIZE][j4][0];
        vec1_xx=vec2_x+img->dbMV[i4+BLOCK_SIZE][j4][0];
        vec2_y=j4*4*mv_mul;
        vec1_y=vec2_y+img->dfMV[i4+BLOCK_SIZE][j4][1];
        vec1_yy=vec2_y+img->dbMV[i4+BLOCK_SIZE][j4][1];

        get_block(ref_frame,vec1_x,vec1_y,img,tmp_block);
        get_block(ref_frame_bw,vec1_xx,vec1_yy,img,tmp_blockbw);

        for(ii=0;ii<BLOCK_SIZE;ii++)
          for(jj=0;jj<BLOCK_SIZE;jj++)
            img->mpr[ii+ioff][jj+joff] = (tmp_block[ii][jj]+tmp_blockbw[ii][jj]+1)/2;
      }

      itrans(img,ioff,joff,i,j);      // use DCT transform and make 4x4 block m7 from prediction block mpr

      for(ii=0;ii<BLOCK_SIZE;ii++)
      {
        for(jj=0;jj<BLOCK_SIZE;jj++)
        {
          imgY[j4*BLOCK_SIZE+jj][i4*BLOCK_SIZE+ii]=img->m7[ii][jj]; // contruct picture from 4x4 blocks
        }
      }
    }
  }




#if POS
  imgY[img->block_y*BLOCK_SIZE][img->block_x*BLOCK_SIZE]= color;
#endif

  /**********************************
   *    chroma                      *
   *********************************/
  for(uv=0;uv<2;uv++)
  {
    if (img->imod==INTRA_MB_OLD || img->imod==INTRA_MB_NEW)// intra mode
    {
      js0=0;
      js1=0;
      js2=0;
      js3=0;
      for(i=0;i<4;i++)
      {
        if(mb_available_up)
        {
          js0=js0+imgUV[uv][img->pix_c_y-1][img->pix_c_x+i];
          js1=js1+imgUV[uv][img->pix_c_y-1][img->pix_c_x+i+4];
        }
        if(mb_available_left)
        {
          js2=js2+imgUV[uv][img->pix_c_y+i][img->pix_c_x-1];
          js3=js3+imgUV[uv][img->pix_c_y+i+4][img->pix_c_x-1];
        }
      }
      if(mb_available_up && mb_available_left)
      {
        js[0][0]=(js0+js2+4)/8;
        js[1][0]=(js1+2)/4;
        js[0][1]=(js3+2)/4;
        js[1][1]=(js1+js3+4)/8;
      }
      if(mb_available_up && !mb_available_left)
      {
        js[0][0]=(js0+2)/4;
        js[1][0]=(js1+2)/4;
        js[0][1]=(js0+2)/4;
        js[1][1]=(js1+2)/4;
      }
      if(mb_available_left && !mb_available_up)
      {
        js[0][0]=(js2+2)/4;
        js[1][0]=(js2+2)/4;
        js[0][1]=(js3+2)/4;
        js[1][1]=(js3+2)/4;
      }
      if(!mb_available_up && !mb_available_left)
      {
        js[0][0]=128;
        js[1][0]=128;
        js[0][1]=128;
        js[1][1]=128;
      }
    }

    for (j=4;j<6;j++)
    {
      joff=(j-4)*4;
      j4=img->pix_c_y+joff;
      for(i=0;i<2;i++)
      {
        ioff=i*4;
        i4=img->pix_c_x+ioff;
        // make pred
        if(img->imod==INTRA_MB_OLD|| img->imod==INTRA_MB_NEW)// intra
        {
          for(ii=0;ii<4;ii++)
            for(jj=0;jj<4;jj++)
            {
              img->mpr[ii+ioff][jj+joff]=js[i][j-4];
            }
        }
        // //////////////////////////////////
        // Forward Reconstruction
        // //////////////////////////////////
        else if(img->imod==B_Forward)
        {
          for(jj=0;jj<4;jj++)
          {
            jf=(j4+jj)/2;
            for(ii=0;ii<4;ii++)
            {
              if1=(i4+ii)/2;
              i1=(img->pix_c_x+ii+ioff)*f1+img->fw_mv[if1+4][jf][0];
              j1=(img->pix_c_y+jj+joff)*f1+img->fw_mv[if1+4][jf][1];

              ii0=max (0, min (i1/f1, img->width_cr-1));
              jj0=max (0, min (j1/f1, img->height_cr-1));
              ii1=max (0, min ((i1+f2)/f1, img->width_cr-1));
              jj1=max (0, min ((j1+f2)/f1, img->height_cr-1));

              if1=(i1 & f2);
              jf1=(j1 & f2);
              if0=f1-if1;
              jf0=f1-jf1;
              img->mpr[ii+ioff][jj+joff]=(if0*jf0*mcef[ref_frame_fw][uv][jj0][ii0]+
                                          if1*jf0*mcef[ref_frame_fw][uv][jj0][ii1]+
                                          if0*jf1*mcef[ref_frame_fw][uv][jj1][ii0]+
                                          if1*jf1*mcef[ref_frame_fw][uv][jj1][ii1]+f4)/f3;
            }
          }
        }
        // //////////////////////////////////
        // Backward Reconstruction
        // //////////////////////////////////
        else if(img->imod==B_Backward)
        {
          for(jj=0;jj<4;jj++)
          {
            jf=(j4+jj)/2;
            for(ii=0;ii<4;ii++)
            {
              if1=(i4+ii)/2;
              i1=(img->pix_c_x+ii+ioff)*f1+img->bw_mv[if1+4][jf][0];
              j1=(img->pix_c_y+jj+joff)*f1+img->bw_mv[if1+4][jf][1];

              ii0=max (0, min (i1/f1, img->width_cr-1));
              jj0=max (0, min (j1/f1, img->height_cr-1));
              ii1=max (0, min ((i1+f2)/f1, img->width_cr-1));
              jj1=max (0, min ((j1+f2)/f1, img->height_cr-1));

              if1=(i1 & f2);
              jf1=(j1 & f2);
              if0=f1-if1;
              jf0=f1-jf1;

              img->mpr[ii+ioff][jj+joff]=(if0*jf0*mcef[ref_frame_bw][uv][jj0][ii0]+
                                          if1*jf0*mcef[ref_frame_bw][uv][jj0][ii1]+
                                          if0*jf1*mcef[ref_frame_bw][uv][jj1][ii0]+
                                          if1*jf1*mcef[ref_frame_bw][uv][jj1][ii1]+f4)/f3;
            }
          }
        }
        // //////////////////////////////////
        // Bidirect Reconstruction
        // //////////////////////////////////
        else if(img->imod==B_Bidirect)
        {
          for(jj=0;jj<4;jj++)
          {
            jf=(j4+jj)/2;
            for(ii=0;ii<4;ii++)
            {
              ifx=(i4+ii)/2;
              i1=(img->pix_c_x+ii+ioff)*f1+img->fw_mv[ifx+4][jf][0];
              j1=(img->pix_c_y+jj+joff)*f1+img->fw_mv[ifx+4][jf][1];

              ii0=max (0, min (i1/f1, img->width_cr-1));
              jj0=max (0, min (j1/f1, img->height_cr-1));
              ii1=max (0, min ((i1+f2)/f1, img->width_cr-1));
              jj1=max (0, min ((j1+f2)/f1, img->height_cr-1));
              if1=(i1 & f2);
              jf1=(j1 & f2);
              if0=f1-if1;
              jf0=f1-jf1;

              fw_pred=(if0*jf0*mcef[ref_frame_fw][uv][jj0][ii0]+
                       if1*jf0*mcef[ref_frame_fw][uv][jj0][ii1]+
                       if0*jf1*mcef[ref_frame_fw][uv][jj1][ii0]+
                       if1*jf1*mcef[ref_frame_fw][uv][jj1][ii1]+f4)/f3;

              i1=(img->pix_c_x+ii+ioff)*f1+img->bw_mv[ifx+4][jf][0];
              j1=(img->pix_c_y+jj+joff)*f1+img->bw_mv[ifx+4][jf][1];

              ii0=max (0, min (i1/f1, img->width_cr-1));
              jj0=max (0, min (j1/f1, img->height_cr-1));
              ii1=max (0, min ((i1+f2)/f1, img->width_cr-1));
              jj1=max (0, min ((j1+f2)/f1, img->height_cr-1));
              if1=(i1 & f2);
              jf1=(j1 & f2);
              if0=f1-if1;
              jf0=f1-jf1;
              bw_pred=(if0*jf0*mcef[ref_frame_bw][uv][jj0][ii0]+
                       if1*jf0*mcef[ref_frame_bw][uv][jj0][ii1]+
                       if0*jf1*mcef[ref_frame_bw][uv][jj1][ii0]+
                       if1*jf1*mcef[ref_frame_bw][uv][jj1][ii1]+f4)/f3;

              img->mpr[ii+ioff][jj+joff]=(int)((fw_pred+bw_pred)/2.+.5);
            }
          }
        }
        // //////////////////////////////////
        // Direct Reconstruction
        // //////////////////////////////////
        else if(img->imod==B_Direct)
        {
          for(jj=0;jj<4;jj++)
          {
            jf=(j4+jj)/2;
            for(ii=0;ii<4;ii++)
            {
              ifx=(i4+ii)/2;

              i1=(img->pix_c_x+ii+ioff)*f1+img->dfMV[ifx+4][jf][0];
              j1=(img->pix_c_y+jj+joff)*f1+img->dfMV[ifx+4][jf][1];

              ii0=max (0, min (i1/f1, img->width_cr-1));
              jj0=max (0, min (j1/f1, img->height_cr-1));
              ii1=max (0, min ((i1+f2)/f1, img->width_cr-1));
              jj1=max (0, min ((j1+f2)/f1, img->height_cr-1));

              if1=(i1 & f2);
              jf1=(j1 & f2);
              if0=f1-if1;
              jf0=f1-jf1;

              if(refFrArr[jf][ifx]==-1)
                ref_frame=1;
              else
                ref_frame=1 + refFrArr[jf][ifx];

              fw_pred=(if0*jf0*mcef[ref_frame][uv][jj0][ii0]+
                       if1*jf0*mcef[ref_frame][uv][jj0][ii1]+
                       if0*jf1*mcef[ref_frame][uv][jj1][ii0]+
                       if1*jf1*mcef[ref_frame][uv][jj1][ii1]+f4)/f3;


              i1=(img->pix_c_x+ii+ioff)*f1+img->dbMV[ifx+4][jf][0];
              j1=(img->pix_c_y+jj+joff)*f1+img->dbMV[ifx+4][jf][1];

              ii0=max (0, min (i1/f1, img->width_cr-1));
              jj0=max (0, min (j1/f1, img->height_cr-1));
              ii1=max (0, min ((i1+f2)/f1, img->width_cr-1));
              jj1=max (0, min ((j1+f2)/f1, img->height_cr-1));

              if1=(i1 & f2);
              jf1=(j1 & f2);
              if0=f1-if1;
              jf0=f1-jf1;

              bw_pred=(if0*jf0*mcef[ref_frame_bw][uv][jj0][ii0]+
                       if1*jf0*mcef[ref_frame_bw][uv][jj0][ii1]+
                       if0*jf1*mcef[ref_frame_bw][uv][jj1][ii0]+
                       if1*jf1*mcef[ref_frame_bw][uv][jj1][ii1]+f4)/f3;
              // LG : direct residual coding
              img->mpr[ii+ioff][jj+joff]=(int)((fw_pred+bw_pred)/2.+.5);
            }
          }
        }

        itrans(img,ioff,joff,2*uv+i,j);
        for(ii=0;ii<4;ii++)
          for(jj=0;jj<4;jj++)
          {
            imgUV[uv][j4+jj][i4+ii]=img->m7[ii][jj];
          }
      } // for(i=0;i<2;i++)
    } // for (j=4;j<6;j++)
  } // for(uv=0;uv<2;uv++)


#if POS
  imgUV[0][img->pix_c_y][img->pix_c_x]= color;
  imgUV[1][img->pix_c_y][img->pix_c_x]= color;
#endif

  return DECODING_OK;
}

