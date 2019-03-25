#define ABIPRED 1
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
 ***************************************************************************
 * \file rdopt.c
 *
 * \brief
 *    Rate-Distortion optimized mode decision
 *
 * \author
 *    Heiko Schwarz <hschwarz@hhi.de>
 *
 * \date
 *    12. April 2001
 **************************************************************************
 */

#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <assert.h>
#include "rdopt_coding_state.h"
#include "elements.h"
#include "refbuf.h"
#include "intrarefresh.h"
#include "vlc.h"
#include "mbuffer.h"
#include "image.h"
#include "mb_access.h"

extern       int  QP2QUANT  [40];

//==== MODULE PARAMETERS ====
int   best_mode;
int   rec_mbY[16][16], rec_mbU[8][8], rec_mbV[8][8], rec_mbY8x8[16][16];    // reconstruction values
int   mpr8x8[16][16];
int   ****cofAC=NULL, ****cofAC8x8=NULL;        // [8x8block][4x4block][level/run][scan_pos]
int   ***cofDC=NULL;                       // [yuv][level/run][scan_pos]
int   **cofAC4x4=NULL, ****cofAC4x4intern=NULL; // [level/run][scan_pos]
int   cbp, cbp8x8, cnt_nonz_8x8;
int   cbp_blk, cbp_blk8x8;
//int   frefframe[4], brefframe[4], b8mode[4], b8pdir[4];
int   frefframe[4][4], brefframe[4][4], b8mode[4], b8pdir[4];
int   best8x8mode [4];                // [block]
int   best8x8pdir [MAXMODE][4];       // [mode][block]
int   best8x8ref  [MAXMODE][4];       // [mode][block]
int   b8_ipredmode[16], b8_intra_pred_modes[16];
CSptr cs_mb=NULL, cs_b8=NULL, cs_cm=NULL, cs_imb=NULL, cs_ib8=NULL, cs_ib4=NULL, cs_pc=NULL;
int   best_c_imode;
int   best_i16offset;

int   best8x8bwref     [MAXMODE][4];       // [mode][block]
int   best8x8abp_type  [MAXMODE][4];       // [mode][block]
int   abp_typeframe[4][4];

/*!
 ************************************************************************
 * \brief
 *    delete structure for RD-optimized mode decision
 ************************************************************************
 */
void clear_rdopt ()
{
  free_mem_DCcoeff (cofDC);
  free_mem_ACcoeff (cofAC);
  free_mem_ACcoeff (cofAC8x8);
  free_mem_ACcoeff (cofAC4x4intern);

  // structure for saving the coding state
  delete_coding_state (cs_mb);
  delete_coding_state (cs_b8);
  delete_coding_state (cs_cm);
  delete_coding_state (cs_imb);
  delete_coding_state (cs_ib8);
  delete_coding_state (cs_ib4);
  delete_coding_state (cs_pc);
}


/*!
 ************************************************************************
 * \brief
 *    create structure for RD-optimized mode decision
 ************************************************************************
 */
void init_rdopt ()
{
  get_mem_DCcoeff (&cofDC);
  get_mem_ACcoeff (&cofAC);
  get_mem_ACcoeff (&cofAC8x8);
  get_mem_ACcoeff (&cofAC4x4intern);
  cofAC4x4 = cofAC4x4intern[0][0];

  // structure for saving the coding state
  cs_mb  = create_coding_state ();
  cs_b8  = create_coding_state ();
  cs_cm  = create_coding_state ();
  cs_imb = create_coding_state ();
  cs_ib8 = create_coding_state ();
  cs_ib4 = create_coding_state ();
  cs_pc  = create_coding_state ();
}



/*! 
 *************************************************************************************
 * \brief
 *    Updates the pixel map that shows, which reference frames are reliable for
 *    each MB-area of the picture.
 *
 * \note
 *    The new values of the pixel_map are taken from the temporary buffer refresh_map
 *
 *************************************************************************************
 */
void UpdatePixelMap()
{
  int mx,my,y,x,i,j;
  if (img->type==I_SLICE)
  {
    for (y=0; y<img->height; y++)
    for (x=0; x<img->width; x++)
    {
      pixel_map[y][x]=1;
    }
  }
  else
  {
    for (my=0; my<img->height/8; my++)
    for (mx=0; mx<img->width/8;  mx++)
    {
      j = my*8 + 8;
      i = mx*8 + 8;
      if (refresh_map[my][mx])
      {
        for (y=my*8; y<j; y++)
        for (x=mx*8; x<i; x++)  pixel_map[y][x] = 1;
      }
      else
      {
        for (y=my*8; y<j; y++)
        for (x=mx*8; x<i; x++)  pixel_map[y][x] = min(pixel_map[y][x]+1, input->num_reference_frames+1);
      }
    }
  }
}

/*! 
 *************************************************************************************
 * \brief
 *    Checks if a given reference frame is reliable for the current 
 *    macroblock, given the motion vectors that the motion search has 
 *    returned.
 *
 * \param ref_frame
 *    The number of the reference frame that we want to check
 *
 * \return
 *    If the return value is 1, the reference frame is reliable. If it 
 *    is 0, then it is not reliable.
 *
 * \note
 *    A specific area in each reference frame is assumed to be unreliable
 *    if the same area has been intra-refreshed in a subsequent frame.
 *    The information about intra-refreshed areas is kept in the pixel_map.
 *
 *************************************************************************************
 */
int CheckReliabilityOfRef (int block, int ref, int mode)
{
  int y,x, block_y, block_x, dy, dx, y_pos, x_pos, yy, xx, pres_x, pres_y;
  int maxold_x  = img->width-1;
  int maxold_y  = img->height-1;
  int ref_frame = ref+1;

  int by0 = (mode>=4?2*(block/2):mode==2?2*block:0);
  int by1 = by0 + (mode>=4||mode==2?2:4);
  int bx0 = (mode>=4?2*(block%2):mode==3?2*block:0);
  int bx1 = bx0 + (mode>=4||mode==3?2:4);

  for (block_y=by0; block_y<by1; block_y++)
    for (block_x=bx0; block_x<bx1; block_x++)
    {
      y_pos  = img->all_mv[block_x][block_y][ref][mode][1];
      y_pos += (img->block_y+block_y) * BLOCK_SIZE * 4;
      x_pos  = img->all_mv[block_x][block_y][ref][mode][0];
      x_pos += (img->block_x+block_x) * BLOCK_SIZE * 4;
      
      /* Here we specify which pixels of the reference frame influence
         the reference values and check their reliability. This is
         based on the function Get_Reference_Pixel */
      
      dy = y_pos & 3;
      dx = x_pos & 3;

      y_pos = (y_pos-dy)/4;
      x_pos = (x_pos-dx)/4;

      if (dy==0 && dx==0) //full-pel
      {
        for (y=0 ; y < BLOCK_SIZE ; y++)
          for (x=0 ; x < BLOCK_SIZE ; x++)
            if (pixel_map[max(0,min(maxold_y,y_pos+y))][max(0,min(maxold_x,x_pos+x))] < ref_frame)
              return 0;
      }
      else  /* other positions */
      {
        if (dy == 0) 
        {
          for (y=0 ; y < BLOCK_SIZE ; y++)
            for (x=0 ; x < BLOCK_SIZE ; x++)
            {
              pres_y = max(0,min(maxold_y,y_pos+y));
              for(xx=-2;xx<4;xx++) {
                pres_x = max(0,min(maxold_x,x_pos+x+xx));
                if (pixel_map[pres_y][pres_x] < ref_frame)
                  return 0;
              }
            }
        }

        else if (dx == 0) 
        {
          for (y=0 ; y < BLOCK_SIZE ; y++)
            for (x=0 ; x < BLOCK_SIZE ; x++)
            {
              pres_x = max(0,min(maxold_x,x_pos+x));
              for(yy=-2;yy<4;yy++) {
                pres_y = max(0,min(maxold_y,y_pos+yy+y));
                if (pixel_map[pres_y][pres_x] < ref_frame)
                  return 0;
              }
            }
        }
        else if (dx == 2) 
        {
          for (y=0 ; y < BLOCK_SIZE ; y++)
            for (x=0 ; x < BLOCK_SIZE ; x++)
            {
              for(yy=-2;yy<4;yy++) {
                pres_y = max(0,min(maxold_y,y_pos+yy+y));
                for(xx=-2;xx<4;xx++) {
                  pres_x = max(0,min(maxold_x,x_pos+xx+x));
                  if (pixel_map[pres_y][pres_x] < ref_frame)
                    return 0;
                }
              }
            }
        }
        else if (dy == 2) 
        {
          for (y=0 ; y < BLOCK_SIZE ; y++)
            for (x=0 ; x < BLOCK_SIZE ; x++)
            {
              for(xx=-2;xx<4;xx++) {
                pres_x = max(0,min(maxold_x,x_pos+xx+x));
                for(yy=-2;yy<4;yy++) {
                  pres_y = max(0,min(maxold_y,y_pos+yy+y));
                  if (pixel_map[pres_y][pres_x] < ref_frame)
                    return 0;
                }
              }
            }
        }
        else 
        {
          for (y=0 ; y < BLOCK_SIZE ; y++)
            for (x=0 ; x < BLOCK_SIZE ; x++)
            {
              pres_y = dy == 1 ? y_pos+y : y_pos+y+1;
              pres_y = max(0,min(maxold_y,pres_y));

              for(xx=-2;xx<4;xx++) {
                pres_x = max(0,min(maxold_x,x_pos+xx+x));
                if (pixel_map[pres_y][pres_x] < ref_frame)
                  return 0;
              }

              pres_x = dx == 1 ? x_pos+x : x_pos+x+1;
              pres_x = max(0,min(maxold_x,pres_x));

              for(yy=-2;yy<4;yy++) {
                pres_y = max(0,min(maxold_y,y_pos+yy+y));
                if (pixel_map[pres_y][pres_x] < ref_frame)
                  return 0;
              }
            }
        }

      }
    }

  return 1;
}



/*! 
 *************************************************************************************
 * \brief
 *    R-D Cost for an 4x4 Intra block
 *************************************************************************************
 */
double RDCost_for_4x4IntraBlocks (int*    nonzero,
                           int     b8,
                           int     b4,
                           int    ipmode,
                           double  lambda,
                           double  min_rdcost,
                           int mostProbableMode)
{
  double  rdcost;
  int     dummy, x, y, rate;
  int     distortion  = 0;
  int     block_x     = 8*(b8%2)+4*(b4%2);
  int     block_y     = 8*(b8/2)+4*(b4/2);
  int     pic_pix_x   = img->pix_x+block_x;
  int     pic_pix_y   = img->pix_y+block_y;
  int     frame_pic_pix_y = pic_pix_y;
  int     pic_block_y = pic_pix_y/4;
  byte    **imgY_orig  = imgY_org;
  byte    **imgY       = enc_picture->imgY;

  Slice          *currSlice    =  img->currentSlice;
  Macroblock     *currMB       = &img->mb_data[img->current_mb_nr];
  SyntaxElement  *currSE       = &img->MB_SyntaxElements[currMB->currSEnr];
  const int      *partMap      = assignSE2partition[input->partition_mode];
  DataPartition  *dataPart;

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    pic_pix_y   = img->field_pix_y + block_y;
    pic_block_y = pic_pix_y/4;
    imgY_orig   = img->top_field ? imgY_org_top : imgY_org_bot;
  }

  //===== perform DCT, Q, IQ, IDCT, Reconstruction =====
  dummy = 0;
  *nonzero = dct_luma (block_x, block_y, &dummy, 1);

  //===== get distortion (SSD) of 4x4 block =====
  for (y=0; y<4; y++)
  for (x=pic_pix_x; x<pic_pix_x+4; x++)  
    distortion += img->quad [imgY_orig[pic_pix_y+y][x] - imgY[frame_pic_pix_y+y][x]];

  //===== RATE for INTRA PREDICTION MODE  (SYMBOL MODE MUST BE SET TO UVLC) =====
  currSE->value1 = (mostProbableMode == ipmode) ? -1 : ipmode < mostProbableMode ? ipmode : ipmode-1;

  //--- set position and type ---
  currSE->context = 4*b8 + b4;
  currSE->type    = SE_INTRAPREDMODE;

  //--- set function pointer ----
  if (input->symbol_mode != UVLC)    currSE->writing = writeIntraPredMode_CABAC;

  //--- choose data partition ---
  if (img->type!=B_SLICE && img->type!=BS_IMG)   dataPart = &(currSlice->partArr[partMap[SE_INTRAPREDMODE]]);
  else                                         dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);

  //--- encode and update rate ---
  if (input->symbol_mode == UVLC)    writeSyntaxElement_Intra4x4PredictionMode(currSE, dataPart);
  else                               dataPart->writeSyntaxElement (currSE, dataPart);
  rate = currSE->len;
  currSE++;
  currMB->currSEnr++;

  //===== RATE for LUMINANCE COEFFICIENTS =====
  if (input->symbol_mode == UVLC)
  {
    rate  += writeCoeff4x4_CAVLC (LUMA, b8, b4, 0);
  }
  else
  {
    rate  += writeLumaCoeff4x4_CABAC (b8, b4, 1);
  }
  rdcost = (double)distortion + lambda*(double)rate;

  return rdcost;
}

/*! 
 *************************************************************************************
 * \brief
 *    Mode Decision for an 4x4 Intra block
 *************************************************************************************
 */
int Mode_Decision_for_4x4IntraBlocks (int  b8,  int  b4,  double  lambda,  int*  min_cost)
{
  int     ipmode, best_ipmode = 0, i, j, k, x, y, cost, dummy;
  int     c_nz, nonzero = 0, rec4x4[4][4], diff[16];
  double  rdcost;
  int     block_x     = 8*(b8%2)+4*(b4%2);
  int     block_y     = 8*(b8/2)+4*(b4/2);
  int     pic_pix_x   = img->pix_x+block_x;
  int     pic_pix_y   = img->pix_y+block_y;
  int     frame_pic_pix_y = pic_pix_y;
  int     pic_block_x = pic_pix_x/4;
  int     pic_block_y = pic_pix_y/4;
  double  min_rdcost  = 1e30;
  byte    **imgY_orig  = imgY_org;

  int left_available, up_available, all_available;

  int     upMode;
  int     leftMode;
  int     mostProbableMode;

  PixelPos left_block;
  PixelPos top_block;

  getLuma4x4Neighbour(img->current_mb_nr, block_x/4, block_y/4, -1,  0, &left_block);
  getLuma4x4Neighbour(img->current_mb_nr, block_x/4, block_y/4,  0, -1, &top_block);


  if(img->MbaffFrameFlag)
  {
    imgY_orig   = (img->top_field ? imgY_org_top:imgY_org_bot);
  }

  upMode            = top_block.available ? img->ipredmode[top_block.pos_x ][top_block.pos_y ] : -1;
  leftMode          = left_block.available ? img->ipredmode[left_block.pos_x][left_block.pos_y] : -1;
  
  mostProbableMode  = (upMode < 0 || leftMode < 0) ? DC_PRED : upMode < leftMode ? upMode : leftMode;

  *min_cost = (1<<20);


  //===== INTRA PREDICTION FOR 4x4 BLOCK =====
  intrapred_luma (pic_pix_x, pic_pix_y, &left_available, &up_available, &all_available);

  //===== LOOP OVER ALL 4x4 INTRA PREDICTION MODES =====
  for (ipmode=0; ipmode<NO_INTRA_PMODE; ipmode++)
  {
    if( (ipmode==DC_PRED) ||
        ((ipmode==VERT_PRED||ipmode==VERT_LEFT_PRED||ipmode==DIAG_DOWN_LEFT_PRED) && up_available ) ||
        ((ipmode==HOR_PRED||ipmode==HOR_UP_PRED) && left_available ) ||
        (all_available) )
    {
      if (!input->rdopt)
      {
        for (k=j=0; j<4; j++)
          for (i=0; i<4; i++, k++)
          {
            diff[k] = imgY_orig[pic_pix_y+j][pic_pix_x+i] - img->mprr[ipmode][j][i];
          }
        cost  = (ipmode == mostProbableMode) ? 0 : (int)floor(4 * lambda );
        cost += SATD (diff, input->hadamard);
        if (cost < *min_cost)
        {
          best_ipmode = ipmode;
          *min_cost   = cost;
        }
      }
      else
      {
        // get prediction and prediction error
        for (j=0; j<4; j++)
        for (i=0; i<4; i++)
        {
          img->mpr[block_x+i][block_y+j]  = img->mprr[ipmode][j][i];
          img->m7[i][j]                   = imgY_orig[pic_pix_y+j][pic_pix_x+i] - img->mprr[ipmode][j][i];
        }

        //===== store the coding state =====
        store_coding_state (cs_cm);
        // get and check rate-distortion cost
        if ((rdcost = RDCost_for_4x4IntraBlocks (&c_nz, b8, b4, ipmode, lambda, min_rdcost, mostProbableMode)) < min_rdcost)
        {
          //--- set coefficients ---
          for (j=0; j<2; j++)
          for (i=0; i<18;i++)  cofAC4x4[j][i]=img->cofAC[b8][b4][j][i];

          //--- set reconstruction ---
          for (y=0; y<4; y++)
          for (x=0; x<4; x++)  rec4x4[y][x] = enc_picture->imgY[frame_pic_pix_y+y][pic_pix_x+x];

          //--- flag if dct-coefficients must be coded ---
          nonzero = c_nz;

          //--- set best mode update minimum cost ---
          min_rdcost  = rdcost;
          best_ipmode = ipmode;
        }
        reset_coding_state (cs_cm);
      }
    }
  }

  //===== set intra mode prediction =====
  img->ipredmode[pic_block_x][pic_block_y] = best_ipmode;
  img->mb_data[img->current_mb_nr].intra_pred_modes[4*b8+b4] = mostProbableMode == best_ipmode ? -1 : best_ipmode < mostProbableMode ? best_ipmode : best_ipmode-1;

  if (!input->rdopt)
  {
    // get prediction and prediction error
    for (j=0; j<4; j++)
      for (i=0; i<4; i++)
      {
        img->mpr[block_x+i][block_y+j]  = img->mprr[best_ipmode][j][i];
        img->m7[i][j]                   = imgY_orig[pic_pix_y+j][pic_pix_x+i] - img->mprr[best_ipmode][j][i];
      }
    nonzero = dct_luma (block_x, block_y, &dummy, 1);
/*    if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
    {
      for (y=0; y<4; y++)
        for (x=0; x<4; x++)
        {
          if(img->top_field)
            imgY_top[pic_pix_y +y][pic_pix_x+x] = imgY[frame_pic_pix_y+y][pic_pix_x+x];
          else
            imgY_bot[pic_pix_y +y][pic_pix_x+x] = imgY[frame_pic_pix_y+y][pic_pix_x+x];
        }
    }*/
  }
  else
  {
    //===== restore coefficients =====
    for (j=0; j<2; j++)
    for (i=0; i<18;i++)  img->cofAC[b8][b4][j][i]=cofAC4x4[j][i];
  
    //===== restore reconstruction and prediction (needed if single coeffs are removed) =====
    for (y=0; y<4; y++)
    for (x=0; x<4; x++)
    {
      enc_picture->imgY[frame_pic_pix_y+y][pic_pix_x+x] = rec4x4[y][x];
      img->mpr[block_x+x][block_y+y] = img->mprr[best_ipmode][y][x];
      /*
      if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
      {
        if(img->top_field)
          imgY_top[pic_pix_y +y][pic_pix_x+x] = rec4x4[y][x];
        else
          imgY_bot[pic_pix_y +y][pic_pix_x+x] = rec4x4[y][x];
      }
      */
    }
  }

  return nonzero;
}


/*! 
 *************************************************************************************
 * \brief
 *    Mode Decision for an 8x8 Intra block
 *************************************************************************************
 */
int Mode_Decision_for_8x8IntraBlocks(int b8,double lambda,int *cost)
{
  int  nonzero=0, b4;
  int  cost4x4;
  
  *cost = (int)floor(6.0 * lambda + 0.4999);

  for (b4=0; b4<4; b4++)
  {
    if (Mode_Decision_for_4x4IntraBlocks (b8, b4, lambda, &cost4x4))
    {
      nonzero        = 1;
    }
    *cost += cost4x4;
  }

  return nonzero;
}

/*! 
 *************************************************************************************
 * \brief
 *    4x4 Intra mode decision for an macroblock
 *************************************************************************************
 */
int Mode_Decision_for_Intra4x4Macroblock (double lambda,  int* cost)

{
  int  cbp=0, b8, cost8x8;

  for (*cost=0, b8=0; b8<4; b8++)
  {
    if (Mode_Decision_for_8x8IntraBlocks (b8, lambda, &cost8x8))
    {
      cbp |= (1<<b8);
    }
    *cost += cost8x8;
  }

  return cbp;
}


/*! 
 *************************************************************************************
 * \brief
 *    R-D Cost for an 8x8 Partition
 *************************************************************************************
 */
double RDCost_for_8x8blocks (int*    cnt_nonz,   // --> number of nonzero coefficients
                      int*    cbp_blk,    // --> cbp blk
                      double  lambda,     // <-- lagrange multiplier
                      int     block,      // <-- 8x8 block number
                      int     mode,       // <-- partitioning mode
                      int     pdir,       // <-- prediction direction
                      int     ref,        // <-- reference frame
                      int     bwd_ref)   // <-- abp type
{
  int  i, j, k;
  int  rate=0, distortion=0;
  int  dummy, mrate;
  int  fw_mode, bw_mode;
  int  cbp     = 0;
  int  pax     = 8*(block%2);
  int  pay     = 8*(block/2);
  int  i0      = pax/4;
  int  j0      = pay/4;
  int  bframe  = (img->type==B_SLICE || img->type==BS_IMG);
  int  direct  = (bframe && mode==0);
  int  b8value = B8Mode2Value (mode, pdir);

  Macroblock    *currMB    = &img->mb_data[img->current_mb_nr];
  SyntaxElement *currSE    = &img->MB_SyntaxElements[currMB->currSEnr];
  Slice         *currSlice = img->currentSlice;
  DataPartition *dataPart;
  const int     *partMap   = assignSE2partition[input->partition_mode];

  int  **frefarr = refFrArr;  // For MB level field/frame 
  int  block_y = img->block_y;
  int pix_y = img->pix_y;
  byte **imgY_original  =  imgY_org;

  EncodingEnvironmentPtr eep_dp;

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    block_y = img->field_block_y;
    pix_y   = img->field_pix_y;
    if(img->top_field)
    {
      imgY_original = imgY_org_top;
      frefarr   = refFrArr_top;

    }
    else
    {
      imgY_original = imgY_org_bot;
      frefarr   = refFrArr_bot;
    }
  }

  //=====
  //=====  GET COEFFICIENTS, RECONSTRUCTIONS, CBP
  //=====
  if (direct)
  {
    if (input->direct_type)
      *cnt_nonz = LumaResidualCoding8x8 (&cbp, cbp_blk, block, 0, 0, max(0,fwdir_refFrArr[block_y+j0][img->block_x+i0]), 0);
    else        
      *cnt_nonz = LumaResidualCoding8x8 (&cbp, cbp_blk, block, 0, 0, max(0,frefarr[block_y+j0][img->block_x+i0]), 0);      
  }
  else
  {
    fw_mode   = (pdir==0||pdir==2 ? mode : 0);
    bw_mode   = (pdir==1||pdir==2 ? mode : 0);
    *cnt_nonz = LumaResidualCoding8x8 (&cbp, cbp_blk, block, fw_mode, bw_mode, ref, bwd_ref);
  }

  //===== get residue =====
  if (input->rdopt==2 && img->type!=B_SLICE)
  {
    // We need the reconstructed prediction residue for the simulated decoders.
    compute_residue_b8block (block, -1);
  }

  //=====
  //=====   GET DISTORTION
  //=====
  if (input->rdopt==2 && img->type!=B_SLICE)
  {
    for (k=0; k<input->NoOfDecoders ;k++)
    {
      decode_one_b8block (k, P8x8, block, mode, ref);
      for (j=img->pix_y+pay; j<img->pix_y+pay+8; j++)
      for (i=img->pix_x+pax; i<img->pix_x+pax+8; i++)
      {
        distortion += img->quad[imgY_org[j][i] - decs->decY[k][j][i]];
      }
    }
    distortion /= input->NoOfDecoders;
  }
  else
  {
    for (j=pay; j<pay+8; j++)
    for (i=img->pix_x+pax; i<img->pix_x+pax+8; i++)
    {
      distortion += img->quad [imgY_original[pix_y+j][i] - enc_picture->imgY[img->pix_y+j][i]];
    }
  }

  //=====
  //=====   GET RATE
  //=====
  //----- block 8x8 mode -----
  if (input->symbol_mode == UVLC)
  {
    ue_linfo (b8value, dummy, &mrate, &dummy);
    rate += mrate;
  }
  else
  {
    currSE->value1  = b8value;
    currSE->writing = writeB8_typeInfo_CABAC;
    currSE->type    = SE_MBTYPE;
    if (img->type==B_SLICE && img->type!=BS_IMG)   dataPart = &(currSlice->partArr[partMap[currSE->type]]);
    else                                         dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);
    dataPart->writeSyntaxElement (currSE, dataPart);
    rate += currSE->len;
    currSE++;
    currMB->currSEnr++;
  }

  //----- motion information -----
  if (!direct)
  {
    if ((input->num_reference_frames>1 || input->add_ref_frame>0) && (pdir==0 || pdir==2))
      rate  += writeReferenceFrame (mode, i0, j0, 1, ref);
    if(input->StoredBPictures > 0)
    {
      if (pdir==1 || pdir==2)
      {
        rate  += writeReferenceFrame (mode, i0, j0, 0, bwd_ref);
      }
    }

    if (pdir==0 || pdir==2)
    {
      if(img->type==BS_IMG && (pdir==1 || pdir==2) )
        rate  += writeMotionVector8x8 (i0, j0, i0+2, j0+2, ref, 1/*DMV*/, 1, mode);
      else
        rate  += writeMotionVector8x8 (i0, j0, i0+2, j0+2, ref, 0, 1, mode);
    }
    if (pdir==1 || pdir==2)
    {
      if(img->type==BS_IMG && (pdir==0 || pdir==2) )
        rate  += writeMotionVector8x8 (i0, j0, i0+2, j0+2, bwd_ref, 1/*DMV*/, 0, mode);
      else
        rate  += writeMotionVector8x8 (i0, j0, i0+2, j0+2, bwd_ref, 0, 0, mode);
    }
  }

  //----- coded block pattern (for CABAC only) -----
  if (input->symbol_mode == CABAC)
  {
    dataPart = &(currSlice->partArr[partMap[SE_CBP_INTER]]);
    eep_dp   = &(dataPart->ee_cabac);
    mrate    = arienco_bits_written (eep_dp);
    writeCBP_BIT_CABAC (block, ((*cnt_nonz>0)?1:0), cbp8x8, currMB, 1, eep_dp);
    mrate    = arienco_bits_written (eep_dp) - mrate;
    rate    += mrate;
  }

  //----- luminance coefficients -----
  if (*cnt_nonz)
  {
    rate += writeLumaCoeff8x8 (block, 0);
  }

  return (double)distortion + lambda * (double)rate;
}


/*! 
 *************************************************************************************
 * \brief
 *    Gets mode offset for intra16x16 mode
 *************************************************************************************
 */
int I16Offset (int cbp, int i16mode)
{
  return (cbp&15?13:1) + i16mode + ((cbp&0x30)>>2);
}


/*! 
 *************************************************************************************
 * \brief
 *    Sets modes and reference frames for an macroblock
 *************************************************************************************
 */
void SetModesAndRefframeForBlocks (int mode)
{
  int i,j,k,l;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int  bframe  = (img->type==B_SLICE || img->type==BS_IMG);
  int  **abp_type_arr = abp_type_FrArr;
  int    **fwrefarr = fw_refFrArr;
  int    **bwrefarr = bw_refFrArr;
  int        **frefarr = refFrArr;
  int        block_y = img->block_y;

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    block_y = img->field_block_y;
    if(img->top_field)
    {
      frefarr  = refFrArr_top;
      fwrefarr = fw_refFrArr_top;
      bwrefarr = bw_refFrArr_top;
      abp_type_arr = abp_type_FrArr_top;
    }
    else
    {
      frefarr = refFrArr_bot;
      fwrefarr = fw_refFrArr_bot;
      bwrefarr = bw_refFrArr_bot;
      abp_type_arr = abp_type_FrArr_bot;
    }
  }

  //--- macroblock type ---
  currMB->mb_type = mode;

  //--- block 8x8 mode and prediction direction ---
  switch (mode)
  {
    case 0:
      for(i=0;i<4;i++)
      {
        currMB->b8mode[i] = 0;
        currMB->b8pdir[i] = (bframe?2:0);
      }
      break;
    case 1:
    case 2:
    case 3:
      for(i=0;i<4;i++)
      {
        currMB->b8mode[i] = mode;
        currMB->b8pdir[i] = best8x8pdir[mode][i];
      }
      break;
    case P8x8:
      for(i=0;i<4;i++)
      {
        currMB->b8mode[i]   = best8x8mode[i];
        currMB->b8pdir[i]   = best8x8pdir[mode][i];
      }
      break;
    case I4MB:
      for(i=0;i<4;i++)
      {
        currMB->b8mode[i] = IBLOCK;
        currMB->b8pdir[i] = -1;
     }
      break;
    case I16MB:
      for(i=0;i<4;i++)
      {
        currMB->b8mode[i] = 0;
        currMB->b8pdir[i] = -1;
      }
      break;
    default:
      printf ("Unsupported mode in SetModesAndRefframeForBlocks!\n");
      exit (1);
  }

#define IS_FW ((best8x8pdir[mode][k]==0 || best8x8pdir[mode][k]==2) && (mode!=P8x8 || best8x8mode[k]!=0 || !bframe))
#define IS_BW ((best8x8pdir[mode][k]==1 || best8x8pdir[mode][k]==2) && (mode!=P8x8 || best8x8mode[k]!=0))
  //--- reference frame arrays ---
  if (mode==0 || mode==I4MB || mode==I16MB)
  {
    if (bframe)
    { 
      for (j=0;j<4;j++)
      for (i=0;i<4;i++)
      {
        if(!mode)
        {
          if (img->type==BS_IMG)
          {
            fwrefarr[block_y+j][img->block_x+i] = fwdir_refFrArr[block_y+j][img->block_x+i];
            bwrefarr[block_y+j][img->block_x+i] = bwdir_refFrArr[block_y+j][img->block_x+i];            
            frefarr [block_y+j][img->block_x+i] = fwdir_refFrArr[block_y+j][img->block_x+i];
            abp_type_arr[block_y+j][img->block_x+i] = best8x8abp_type[mode][2*(j/2)+(i/2)];
          }
          else 
          {
            if (input->direct_type)
            {
              fwrefarr[block_y+j][img->block_x+i] = fwdir_refFrArr[block_y+j][img->block_x+i];
              bwrefarr[block_y+j][img->block_x+i] = bwdir_refFrArr[block_y+j][img->block_x+i];            
            }
            else
            {
              fwrefarr[block_y+j][img->block_x+i] = frefarr[block_y+j][img->block_x+i];
              bwrefarr[block_y+j][img->block_x+i] = min(frefarr[block_y+j][img->block_x+i],0);
            }
          }
        }
        else
        {
          fwrefarr[block_y+j][img->block_x+i] = -1;
          bwrefarr[block_y+j][img->block_x+i] = -1;          
          if (img->type==BS_IMG)
          {
            frefarr [block_y+j][img->block_x+i] = -1;
            abp_type_arr[block_y+j][img->block_x+i] = best8x8abp_type[mode][2*(j/2)+(i/2)];
          }
        }
      }
    }
    else
    {
      for (j=0;j<4;j++)
      for (i=0;i<4;i++)
      {
        frefarr [block_y+j][img->block_x+i] = (mode==0?0:-1);
      }
    }
  }
  else
  {
    if (bframe)
    {
      for (j=0;j<4;j++)
      for (i=0;i<4;i++)
      {
        k = 2*(j/2)+(i/2);
        l = 2*(j%2)+(i%2);

        if (img->type==BS_IMG)
        {
          if(mode == P8x8 && best8x8mode[k]==0)
          {
            fwrefarr[block_y+j][img->block_x+i] = fwdir_refFrArr[block_y+j][img->block_x+i];
            bwrefarr[block_y+j][img->block_x+i] = bwdir_refFrArr[block_y+j][img->block_x+i];            
          }
          else
          {
            fwrefarr[block_y+j][img->block_x+i] = (IS_FW ? best8x8ref[mode][k]   : -1);
            bwrefarr[block_y+j][img->block_x+i] = (IS_BW ? best8x8bwref[mode][k] : -1);
          }
          frefarr [block_y+j][img->block_x+i] = fwrefarr[block_y+j][img->block_x+i];
          abp_type_arr[block_y+j][img->block_x+i] = (IS_BW ? best8x8abp_type[mode][k] : 0);
        }
        else
        {
          if(mode == P8x8 && best8x8mode[k]==0)
          {           
            if (input->direct_type)
            {
              fwrefarr[block_y+j][img->block_x+i] = fwdir_refFrArr[block_y+j][img->block_x+i];
              bwrefarr[block_y+j][img->block_x+i] = bwdir_refFrArr[block_y+j][img->block_x+i];            
            }
            else
            {
              fwrefarr[block_y+j][img->block_x+i] = frefarr[block_y+j][img->block_x+i];
              bwrefarr[block_y+j][img->block_x+i] = min(frefarr[block_y+j][img->block_x+i],0);
            }
          }
          else
          {
            fwrefarr[block_y+j][img->block_x+i] = (IS_FW ? best8x8ref[mode][k] : -1);
            bwrefarr[block_y+j][img->block_x+i] = (IS_BW ?                   0 : -1);
          }
        }
      }
    }
    else
    {
      for (j=0;j<4;j++)
      for (i=0;i<4;i++)
      {
        k = 2*(j/2)+(i/2);
        l = 2*(j%2)+(i%2);
        frefarr [block_y+j][img->block_x+i] = (IS_FW ? best8x8ref[mode][k] : -1);
      }
    }
  }
#undef IS_FW
#undef IS_BW
}


/*! 
 *************************************************************************************
 * \brief
 *    Intra 16x16 mode decision
 *************************************************************************************
 */
void
Intra16x16_Mode_Decision (Macroblock* currMB, int* i16mode)
{
  intrapred_luma_16x16 ();   /* make intra pred for all 4 new modes */
  find_sad_16x16 (i16mode);   /* get best new intra mode */
  currMB->cbp = dct_luma_16x16 (*i16mode);
}



/*! 
 *************************************************************************************
 * \brief
 *    Sets Coefficients and reconstruction for an 8x8 block
 *************************************************************************************
 */
void SetCoeffAndReconstruction8x8 (Macroblock* currMB)
{
  int block, k, j, i;
  int block_y = img->block_y;
  int **ipredmodes = img->ipredmode;
 
  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    block_y   = img->field_block_y;
  }

  //--- restore coefficients ---
  for (block=0; block<6; block++)
  for (k=0; k<4; k++)
  for (j=0; j<2; j++)
  for (i=0; i<65; i++)
    img->cofAC[block][k][j][i] = cofAC8x8[block][k][j][i];

  if (cnt_nonz_8x8<=5 && img->type!=SP_SLICE)
  {
    currMB->cbp     = 0;
    currMB->cbp_blk = 0;
    for (j=0; j<16; j++)
    for (i=0; i<16; i++)  enc_picture->imgY[img->pix_y+j][img->pix_x+i] = mpr8x8[j][i];
  }
  else
  {
    currMB->cbp     = cbp8x8;
    currMB->cbp_blk = cbp_blk8x8;
    for (j=0; j<16; j++)
    for (i=0; i<16; i++)  enc_picture->imgY[img->pix_y+j][img->pix_x+i] = rec_mbY8x8[j][i];
  }

  //===== restore intra prediction modes for 8x8+ macroblock mode =====
  for (k=0, j=block_y; j<block_y+4; j++)
  for (     i=img->block_x; i<img->block_x+4; i++, k++)
  {
    ipredmodes    [i][j] = b8_ipredmode       [k];
    currMB->intra_pred_modes[k] = b8_intra_pred_modes[k];
  }
}


/*! 
 *************************************************************************************
 * \brief
 *    Sets motion vectors for an macroblock
 *************************************************************************************
 */
void SetMotionVectorsMB (Macroblock* currMB, int bframe)
{
  int i, j, k, l, mode8, pdir8, ref, by, bx, bxr, dref;
  int ***mv = tmp_mv;
  int *****all_mv = img->all_mv;
  int *****all_bmv = img->all_bmv;
  int ***fwMV    = tmp_fwMV;
  int ***bwMV    = tmp_bwMV;
  int *****imgmv   = img->mv;
  int *****p_fwMV  = img->p_fwMV;
  int *****p_bwMV  = img->p_bwMV;
  int **refar    = refFrArr;
  int **frefar     = fw_refFrArr;
  int block_y    = img->block_y;
  int **brefar     = bw_refFrArr;
  int  bw_ref;

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    if(img->top_field)
    {
      mv=tmp_mv_top;
      all_mv = img->all_mv_top;
      all_bmv = img->all_bmv_top;
      fwMV = tmp_fwMV_top;
      bwMV = tmp_bwMV_top;
      imgmv = img->mv_top;
      p_fwMV = img->p_fwMV_top;
      p_bwMV = img->p_bwMV_top;
      block_y = img->block_y/2;
      refar    = refFrArr_top;
      frefar     = fw_refFrArr_top;
      brefar     = bw_refFrArr_top;
    }
    else
    {
      mv=tmp_mv_bot;
      all_mv = img->all_mv_bot;
      all_bmv = img->all_bmv_bot;
      fwMV = tmp_fwMV_bot;
      bwMV = tmp_bwMV_bot;
      imgmv = img->mv_bot;
      p_fwMV = img->p_fwMV_bot;
      p_bwMV = img->p_bwMV_bot;
      block_y = (img->block_y-4)/2;
      refar    = refFrArr_bot;
      frefar     = fw_refFrArr_bot;
      brefar     = bw_refFrArr_bot;
    }
  }
  for (j=0; j<4; j++)
  for (i=0; i<4; i++)
  {
    mode8 = currMB->b8mode[k=2*(j/2)+(i/2)];
    pdir8 = currMB->b8pdir[k];
    l     = 2*(j%2)+(i%2);
    by    = block_y+j;
    bxr   = img->block_x+i;
    bx    = img->block_x+i+4;
    ref   = (bframe?frefar:refar)[by][bxr];
    bw_ref = (bframe?brefar:refar)[by][bxr];

    if (!bframe)
    {
      if (mode8!=IBLOCK && ref != -1)
      {
        mv   [0][by][bx] = all_mv [i][j][ ref][mode8][0];
        mv   [1][by][bx] = all_mv [i][j][ ref][mode8][1];
      }
      else
      {
        mv   [0][by][bx] = 0;
        mv   [1][by][bx] = 0;
      }

    if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
    {
      rdopt->tmp_mv[0][j][i]  = mv[0][by][bx];
      rdopt->tmp_mv[1][j][i]  = mv[1][by][bx];
    }

    }
    else
    {
      if (pdir8==-1) // intra
      {
        fwMV [0][by][bx] = 0;
        fwMV [1][by][bx] = 0;
        bwMV [0][by][bx] = 0;
        bwMV [1][by][bx] = 0;
        dfMV     [0][by][bx] = 0;
        dfMV     [1][by][bx] = 0;
        dbMV     [0][by][bx] = 0;
        dbMV     [1][by][bx] = 0;
        if (img->type==BS_IMG)
        {
          mv   [0][by][bx] = 0;
          mv   [1][by][bx] = 0;
        }
      }
      else if (pdir8==0) // forward
      {
        fwMV [0][by][bx] = all_mv [i][j][ ref][mode8][0];
        fwMV [1][by][bx] = all_mv [i][j][ ref][mode8][1];
        bwMV [0][by][bx] = 0;
        bwMV [1][by][bx] = 0;
        dfMV     [0][by][bx] = 0;
        dfMV     [1][by][bx] = 0;
        dbMV     [0][by][bx] = 0;
        dbMV     [1][by][bx] = 0;
        if (img->type==BS_IMG)
        {
          mv   [0][by][bx] = all_mv [i][j][ ref][mode8][0];
          mv   [1][by][bx] = all_mv [i][j][ ref][mode8][1];
        }
      }
      else if (pdir8==1) // backward
      {
        fwMV [0][by][bx] = 0;
        fwMV [1][by][bx] = 0;
        if (img->type==BS_IMG)
        {
          bwMV [0][by][bx] = all_bmv[i][j][bw_ref][mode8][0];
          bwMV [1][by][bx] = all_bmv[i][j][bw_ref][mode8][1];
          mv   [0][by][bx] = 0;
          mv   [1][by][bx] = 0;
        }
        else
        {
          bwMV [0][by][bx] = all_bmv[i][j][   0][mode8][0];
          bwMV [1][by][bx] = all_bmv[i][j][   0][mode8][1];
        }
        dfMV     [0][by][bx] = 0;
        dfMV     [1][by][bx] = 0;
        dbMV     [0][by][bx] = 0;
        dbMV     [1][by][bx] = 0;
      }
      else if (mode8!=0) // bidirect
      {
        fwMV [0][by][bx] = all_mv [i][j][ ref][mode8][0];
        fwMV [1][by][bx] = all_mv [i][j][ ref][mode8][1];
        if (img->type==BS_IMG)
        {
          bwMV [0][by][bx] = all_bmv[i][j][bw_ref][mode8][0];
          bwMV [1][by][bx] = all_bmv[i][j][bw_ref][mode8][1];
          mv   [0][by][bx] = all_mv [i][j][ref][mode8][0];
          mv   [1][by][bx] = all_mv [i][j][ref][mode8][1];
        }
        else
        {
          bwMV [0][by][bx] = all_bmv[i][j][   0][mode8][0];
          bwMV [1][by][bx] = all_bmv[i][j][   0][mode8][1];
        }
        dfMV     [0][by][bx] = 0;
        dfMV     [1][by][bx] = 0;
        dbMV     [0][by][bx] = 0;
        dbMV     [1][by][bx] = 0;
      }
      else // direct
      {
        if(img->type==BS_IMG)
        {
          fwMV [0][by][bx] = dfMV     [0][by][bx] = all_mv [i][j][max(0,fwdir_refFrArr[by][bxr])][    0][0];
          fwMV [1][by][bx] = dfMV     [1][by][bx] = all_mv [i][j][max(0,fwdir_refFrArr[by][bxr])][    0][1];
          if(bwdir_refFrArr[by][bxr]==-1)
          {
            bwMV [0][by][bx] = dbMV     [0][by][bx] = all_bmv[i][j][                             1][    0][0];
            bwMV [1][by][bx] = dbMV     [1][by][bx] = all_bmv[i][j][                             1][    0][1];
          }
          else
          {
            bwMV [0][by][bx] = dbMV     [0][by][bx] = all_bmv[i][j][bwdir_refFrArr[by][bxr]][    0][0];
            bwMV [1][by][bx] = dbMV     [1][by][bx] = all_bmv[i][j][bwdir_refFrArr[by][bxr]][    0][1];
          }
          mv   [0][by][bx] = fwMV [0][by][bx];
          mv   [1][by][bx] = fwMV [1][by][bx];
        }
        else
        {
          if (input->direct_type)
            dref = max(0,fwdir_refFrArr[by][bxr]);        
          else
            dref = max(0,refar[by][bxr]);        

          fwMV [0][by][bx] = dfMV     [0][by][bx] = all_mv [i][j][dref][    0][0];
          fwMV [1][by][bx] = dfMV     [1][by][bx] = all_mv [i][j][dref][    0][1];
          bwMV [0][by][bx] = dbMV     [0][by][bx] = all_bmv[i][j][   0][    0][0];
          bwMV [1][by][bx] = dbMV     [1][by][bx] = all_bmv[i][j][   0][    0][1];
        }
      }

      if(img->type==BS_IMG)
      {
        if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
        {
          rdopt->tmp_mv[0][j][i]  = mv[0][by][bx];
          rdopt->tmp_mv[1][j][i]  = mv[1][by][bx];
        }
      }
      if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
      {
        rdopt->tmp_fwMV[0][j][i]  = fwMV[0][by][bx];
        rdopt->tmp_fwMV[1][j][i]  = fwMV[1][by][bx];
        rdopt->tmp_bwMV[0][j][i]  = bwMV[0][by][bx];
        rdopt->tmp_bwMV[1][j][i]  = bwMV[1][by][bx];
        rdopt->dfMV[0][j][i] = dfMV     [0][by][bx];
        rdopt->dfMV[1][j][i] = dfMV     [1][by][bx];                  
        rdopt->dbMV[0][j][i] = dbMV     [0][by][bx];
        rdopt->dbMV[1][j][i] = dbMV     [1][by][bx];                  
      }
    }
  }

// copy all the motion vectors into rdopt structure
// Can simplify this by copying the MV's of the best mode (TBD)

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
  {
    for(i=0;i<4;i++)
      for(j=0;j<4;j++)
        for(k=0;k<img->buf_cycle;k++)
          for(l=0;l<9;l++)
          {
            rdopt->all_mv[i][j][k][l][0] = all_mv[i][j][k][l][0];
            rdopt->all_bmv[i][j][k][l][0] = all_bmv[i][j][k][l][0];
            rdopt->p_fwMV[i][j][k][l][0] = p_fwMV[i][j][k][l][0];
            rdopt->p_bwMV[i][j][k][l][0] = p_bwMV[i][j][k][l][0];
            rdopt->mv[i][j][k][l][0]    = imgmv[i][j][k][l][0];
          
            rdopt->all_mv[i][j][k][l][1] = all_mv[i][j][k][l][1];
            rdopt->all_bmv[i][j][k][l][1] =all_bmv[i][j][k][l][1];
            rdopt->p_fwMV[i][j][k][l][1] = p_fwMV[i][j][k][l][1];
            rdopt->p_bwMV[i][j][k][l][1] = p_bwMV[i][j][k][l][1];
            rdopt->mv[i][j][k][l][1]    = imgmv[i][j][k][l][1];
          }
  }

}



/*! 
 *************************************************************************************
 * \brief
 *    R-D Cost for a macroblock
 *************************************************************************************
 */
int RDCost_for_macroblocks (double   lambda,      // <-- lagrange multiplier
                        int      mode,        // <-- modus (0-COPY/DIRECT, 1-16x16, 2-16x8, 3-8x16, 4-8x8(+), 5-Intra4x4, 6-Intra16x16)
                        double*  min_rdcost)  // <-> minimum rate-distortion cost
{
  int         i, j, k; //, k, ****ip4;
  int         i16mode, rate=0, distortion=0;
  double      rdcost;
  Macroblock  *currMB   = &img->mb_data[img->current_mb_nr];
  int         bframe    = (img->type==B_SLICE || img->type==BS_IMG);
  int         tmp_cc;
  int         use_of_cc =  (img->type!=I_SLICE &&  input->symbol_mode!=CABAC);
  int         cc_rate, dummy;
  byte        **imgY_orig = imgY_org;
  byte        ***imgUV_orig = imgUV_org;
  int         pix_y = img->pix_y;
  int         pix_c_y = img->pix_c_y;

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    pix_y = img->field_pix_y;
    pix_c_y = img->field_pix_c_y;
    imgY_orig = img->top_field ? imgY_org_top:imgY_org_bot;
    imgUV_orig = img->top_field ? imgUV_org_top:imgUV_org_bot;
  }


  //=====
  //=====  SET REFERENCE FRAMES AND BLOCK MODES
  //=====
  SetModesAndRefframeForBlocks (mode);

  //=====
  //=====  GET COEFFICIENTS, RECONSTRUCTIONS, CBP
  //=====
  if (mode<P8x8)
  {
    LumaResidualCoding ();
  }
  else if (mode==P8x8)
  {
    SetCoeffAndReconstruction8x8 (currMB);
  }
  else if (mode==I4MB)
  {
    currMB->cbp = Mode_Decision_for_Intra4x4Macroblock (lambda, &dummy);
  }
  else if (mode==I16MB)
  {
    Intra16x16_Mode_Decision  (currMB, &i16mode);
  }

  if (input->rdopt==2 && img->type!=B_SLICE)
  {
    // We need the reconstructed prediction residue for the simulated decoders.
    compute_residue_mb (mode==I16MB?i16mode:-1);
  }

  img->i16offset = 0;
  dummy = 0;
  ChromaResidualCoding (&dummy);
  if (mode==I16MB)     img->i16offset = I16Offset  (currMB->cbp, i16mode);


  //=====
  //=====   GET DISTORTION
  //=====
  // LUMA
  if (input->rdopt==2 && img->type!=B_SLICE)
  {
    for (k=0; k<input->NoOfDecoders ;k++)
    {
      decode_one_mb (k, currMB);
      for (j=0; j<MB_BLOCK_SIZE; j++)
      for (i=img->pix_x; i<img->pix_x+MB_BLOCK_SIZE; i++)
      {
        distortion += img->quad [imgY_orig[pix_y+j][i] - decs->decY[k][img->pix_y+j][i]];
      }
    }
    distortion /= input->NoOfDecoders;
  }
  else
  {
    for (j=0; j<16; j++)
    for (i=img->pix_x; i<img->pix_x+16; i++)
    {
      distortion += img->quad [imgY_orig[j+pix_y][i] - enc_picture->imgY[img->pix_y+j][i]];
    }
  }

  // CHROMA
  for (j=0; j<8; j++)
  for (i=img->pix_c_x; i<img->pix_c_x+8; i++)
  {
    distortion += img->quad [imgUV_orig[0][j+pix_c_y][i] - enc_picture->imgUV[0][img->pix_c_y+j][i]];
    distortion += img->quad [imgUV_orig[1][j+pix_c_y][i] - enc_picture->imgUV[1][img->pix_c_y+j][i]];
  }


  //=====   S T O R E   C O D I N G   S T A T E   =====
  //---------------------------------------------------
  store_coding_state (cs_cm);
  

  //=====
  //=====   GET RATE
  //=====
  //----- macroblock header -----
  if (use_of_cc)
  {
    if (currMB->mb_type!=0 || (bframe && currMB->cbp!=0))
    {
      // cod counter and macroblock mode are written ==> do not consider code counter
      tmp_cc = img->cod_counter;
      rate   = writeMBHeader (1); 
      ue_linfo (tmp_cc, dummy, &cc_rate, &dummy);
      rate  -= cc_rate;
      img->cod_counter = tmp_cc;
    }
    else
    {
      // cod counter is just increased  ==> get additional rate
      ue_linfo (img->cod_counter+1, dummy, &rate,    &dummy);
      ue_linfo (img->cod_counter,   dummy, &cc_rate, &dummy);
      rate -= cc_rate;
    }
  }
  else
  {
    rate = writeMBHeader (1); 
  }
  if (mode)
  {
    //----- motion information -----
    rate  += writeMotionInfo2NAL  ();
  }
  if (mode || (bframe && (currMB->cbp!=0 || input->symbol_mode==CABAC)))
  {
    rate  += writeCBPandLumaCoeff ();

    rate  += writeChromaCoeff     ();
  }


  //=====   R E S T O R E   C O D I N G   S T A T E   =====
  //-------------------------------------------------------
  reset_coding_state (cs_cm);


  rdcost = (double)distortion + lambda * (double)rate;
 
  if (rdcost >= *min_rdcost)
  {
    return 0;
  }


  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode && mode == 0)
  {
    if(!bframe && TopFieldIsSkipped)
      return 0;
    else if(bframe && TopFieldIsSkipped && currMB->cbp == 0)
      return 0;
  }


  //=====   U P D A T E   M I N I M U M   C O S T   =====
  //-----------------------------------------------------
  *min_rdcost = rdcost;
  return 1;
}





/*! 
 *************************************************************************************
 * \brief
 *    Store macroblock parameters
 *************************************************************************************
 */
void store_macroblock_parameters (int mode)
{
  int  i, j, k, ****i4p, ***i3p;
  Macroblock *currMB  = &img->mb_data[img->current_mb_nr];
  int        bframe   = (img->type==B_SLICE || img->type==BS_IMG);
  int        **frefar = ((img->type==B_SLICE || img->type==BS_IMG)? fw_refFrArr : refFrArr);
  int        **abp_type_arr = abp_type_FrArr;
  int        **brefar = bw_refFrArr;
  int    block_y  = img->block_y;
  int    pix_y    = img->pix_y;
  int    pix_c_y  = img->pix_c_y;

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    block_y = img->field_block_y;
    pix_y   = img->field_pix_y;
    pix_c_y = img->field_pix_c_y;
    if(img->top_field)
    {
      frefar = ((img->type==B_SLICE || img->type==BS_IMG) ? fw_refFrArr_top : refFrArr_top);
      abp_type_arr = abp_type_FrArr_top;
      brefar = bw_refFrArr_top;
    }
    else
    {
      frefar = ((img->type==B_SLICE || img->type==BS_IMG) ? fw_refFrArr_bot : refFrArr_bot);
      abp_type_arr = abp_type_FrArr_bot;
      brefar = bw_refFrArr_bot;
    }

  }

  //--- store best mode ---
  best_mode = mode;
  best_c_imode = currMB->c_ipred_mode;
  best_i16offset = img->i16offset;
  for (i=0; i<4; i++)
  {
    b8mode[i]   = currMB->b8mode[i];
    b8pdir[i]   = currMB->b8pdir[i];
  }

  //--- reconstructed blocks ----
  for (j=0; j<16; j++)
  for (i=0; i<16; i++)
  {
    rec_mbY[j][i] = enc_picture->imgY[img->pix_y+j][img->pix_x+i];
  }
  for (j=0; j<8; j++)
  for (i=0; i<8; i++)
  {
    rec_mbU[j][i] = enc_picture->imgUV[0][img->pix_c_y+j][img->pix_c_x+i];
    rec_mbV[j][i] = enc_picture->imgUV[1][img->pix_c_y+j][img->pix_c_x+i];
  }


  //--- store results of decoders ---
  if (input->rdopt==2 && img->type!=B_SLICE)
  {
    for (k=0; k<input->NoOfDecoders; k++)
    {
      for (j=img->pix_y; j<img->pix_y+16; j++)
      for (i=img->pix_x; i<img->pix_x+16; i++)
      {
        // Keep the decoded values of each MB for updating the ref frames
        decs->decY_best[k][j][i] = decs->decY[k][j][i];
      }
    }
  }

  //--- coeff, cbp, kac ---
  if (mode || bframe)
  {
    i4p=cofAC; cofAC=img->cofAC; img->cofAC=i4p;
    i3p=cofDC; cofDC=img->cofDC; img->cofDC=i3p;
    cbp     = currMB->cbp;
    cbp_blk = currMB->cbp_blk;
  }
  else
  {
    cbp = cbp_blk = 0;
  }

  for (j=0; j<4; j++)
  for (i=0; i<4; i++)
  {  
    frefframe[j][i] = frefar[block_y+j  ][img->block_x +i ];
    if (bframe)
    {
      brefframe[j][i] = brefar[block_y +j ][img->block_x+i  ];
      abp_typeframe[j][i] = abp_type_arr[block_y + j][img->block_x + i];     
    }
  }

/*
  //--- references ---
  frefframe[0] = frefar[block_y  ][img->block_x  ];
  frefframe[1] = frefar[block_y  ][img->block_x+2];
  frefframe[2] = frefar[block_y+2][img->block_x  ];
  frefframe[3] = frefar[block_y+2][img->block_x+2];
  if (bframe)
  {
    brefframe[0] = brefar[block_y  ][img->block_x  ]; 
    brefframe[1] = brefar[block_y  ][img->block_x+2];
    brefframe[2] = brefar[block_y+2][img->block_x  ];
    brefframe[3] = brefar[block_y+2][img->block_x+2];
  }
  */
}


/*! 
 *************************************************************************************
 * \brief
 *    Set stored macroblock parameters
 *************************************************************************************
 */
void set_stored_macroblock_parameters ()
{
  int  i, j, k, ****i4p, ***i3p,l;
  Macroblock  *currMB  = &img->mb_data[img->current_mb_nr];
  int         mode     = best_mode;
  int         bframe   = (img->type==B_SLICE || img->type==BS_IMG);
  int         **frefar = ((img->type==B_SLICE || img->type==BS_IMG)? fw_refFrArr : refFrArr);
  int         **abp_type_arr = abp_type_FrArr;
  int         **brefar = bw_refFrArr;
  int     **ipredmodes = img->ipredmode;
  int     block_y = img->block_y; 
  
  byte    **imgY=enc_picture->imgY;
  byte    ***imgUV=enc_picture->imgUV;

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    block_y = img->field_block_y;
    if(img->top_field)
    {
      frefar = (img->type==B_SLICE || img->type==BS_IMG) ? fw_refFrArr_top : refFrArr_top;
      abp_type_arr = abp_type_FrArr_top;
      brefar = bw_refFrArr_top;
    }
    else
    {
      frefar = (img->type==B_SLICE || img->type==BS_IMG) ? fw_refFrArr_bot : refFrArr_bot;
      abp_type_arr = abp_type_FrArr_bot;
      brefar = bw_refFrArr_bot;

    }
  }

  //===== reconstruction values =====
  for (j=0; j<16; j++)
  for (i=0; i<16; i++)
  {
    imgY[img->pix_y+j][img->pix_x+i] = rec_mbY[j][i];
    if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
      rdopt->rec_mbY[j][i]       = rec_mbY[j][i]; 
  }

  for (j=0; j<8; j++)
  for (i=0; i<8; i++)
  {
    imgUV[0][img->pix_c_y+j][img->pix_c_x+i] = rec_mbU[j][i];
    imgUV[1][img->pix_c_y+j][img->pix_c_x+i] = rec_mbV[j][i];
    if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
    {
      rdopt->rec_mbU[j][i]           = rec_mbU[j][i];   
      rdopt->rec_mbV[j][i]           = rec_mbV[j][i];   
    }
  }

  //===== coefficients and cbp =====
  i4p=cofAC; cofAC=img->cofAC; img->cofAC=i4p;
  i3p=cofDC; cofDC=img->cofDC; img->cofDC=i3p;
  currMB->cbp      = cbp;
  currMB->cbp_blk = cbp_blk;
  //==== macroblock type ====
  currMB->mb_type = mode;
  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
  {
    rdopt->mode = mode;
    rdopt->i16offset = img->i16offset;  // For MBINTLC  -Rajeev
    rdopt->cbp = cbp;
    rdopt->cbp_blk = cbp_blk;
    rdopt->mb_type  = mode;


    for(i=0;i<6;i++)
      for(j=0;j<4;j++)
        for(k=0;k<2;k++)
          for(l=0;l<18;l++)
            rdopt->cofAC[i][j][k][l] = img->cofAC[i][j][k][l];

    for(i=0;i<3;i++)
        for(k=0;k<2;k++)
          for(l=0;l<18;l++)
            rdopt->cofDC[i][k][l] = img->cofDC[i][k][l];
  }


  for (i=0; i<4; i++)
  {
    currMB->b8mode[i]   = b8mode[i];
    currMB->b8pdir[i]   = b8pdir[i];
    if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
    {
      rdopt->b8mode[i]  = b8mode[i];                  
      rdopt->b8pdir[i]  = b8pdir[i];                  
    }
  }
  if (input->rdopt==2 && img->type!=B_SLICE)
  {
    //! save the MB Mode of every macroblock
    decs->dec_mb_mode[img->mb_x][img->mb_y] = mode;
  }

  //==== reference frames =====
  for (j=0; j<4; j++)
  for (i=0; i<4; i++)
  {
    frefar[block_y+j][img->block_x+i] = frefframe[j][i];
//    frefar[block_y+j][img->block_x+i] = frefframe[2*(j/2) + (i/2)];
    if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
      rdopt->frefar[j][i]            = frefar[block_y+j][img->block_x+i];
  }
  if (bframe)
  {
    for (j=0; j<4; j++)
    for (i=0; i<4; i++)
    {
//      brefar[block_y+j][img->block_x+i] = brefframe[2*(j/2) + (i/2)];
      brefar[block_y+j][img->block_x+i] = brefframe[j][i];
      if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
        rdopt->brefar[j][i]            = brefar[block_y+j][img->block_x+i];
      if (img->type==BS_IMG)
      {
        refFrArr [block_y+j][img->block_x+i] = frefar[block_y+j][img->block_x+i];
        abp_type_arr[block_y+j][img->block_x+i] = abp_typeframe[j][i];
      }

    }
  }

  //==== intra prediction modes ====
  currMB->c_ipred_mode = best_c_imode;
  img->i16offset = best_i16offset;
  if (mode==P8x8)
  {
    for (k=0, j=block_y; j<block_y+4; j++)
    for (     i=img->block_x; i<img->block_x+4; i++, k++)
    {
      ipredmodes    [i][j] = b8_ipredmode       [k];
      currMB->intra_pred_modes[k] = b8_intra_pred_modes[k];
    }
  }
  else if (mode!=I4MB)
  {
    for (k=0, j=block_y; j<block_y+4; j++)
      for (     i=img->block_x; i<img->block_x+4; i++, k++)
      {
        ipredmodes    [i][j] = DC_PRED;
        currMB->intra_pred_modes[k] = DC_PRED;
      }
  }

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
  {
    for (k=0, j=block_y; j<block_y+4; j++)
      for (     i=img->block_x; i<img->block_x+4; i++, k++)
      {
        rdopt->ipredmode[i][j] = ipredmodes[i][j];
        rdopt->intra_pred_modes[k] = currMB->intra_pred_modes[k];
      }
    rdopt->c_ipred_mode = currMB->c_ipred_mode;
    rdopt->i16offset = img->i16offset;  // DH
  }

  //==== motion vectors =====
  SetMotionVectorsMB (currMB, bframe);
}



/*! 
 *************************************************************************************
 * \brief
 *    Set reference frames and motion vectors
 *************************************************************************************
 */
void SetRefAndMotionVectors (int block, int mode, int ref, int pdir, int abp_type)
{
  int i, j;
  int     bframe  = (img->type==B_SLICE || img->type==BS_IMG);
  int**   abp_type_arr = abp_type_FrArr;
  int     pmode   = (mode==1||mode==2||mode==3?mode:4);
  int     j0      = ((block/2)<<1);
  int     i0      = ((block%2)<<1);
  int     j1      = j0 + (input->blc_size[pmode][1]>>2);
  int     i1      = i0 + (input->blc_size[pmode][0]>>2);
  int**   frefArr = (bframe ? fw_refFrArr : refFrArr);
  int**   brefArr = bw_refFrArr;
  int***  fmvArr  = (bframe ? tmp_fwMV    : tmp_mv);
  int***  bmvArr  = tmp_bwMV;
  int***** all_mv_arr = img->all_mv;
  int***** all_bmv_arr = img->all_bmv;
  int   block_y = img->block_y;

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    block_y = img->field_block_y;
    if(img->top_field)
    {
      frefArr = (bframe ? fw_refFrArr_top : refFrArr_top);
      brefArr = bw_refFrArr_top;
      fmvArr  = (bframe ? tmp_fwMV_top  : tmp_mv_top);
      bmvArr  = tmp_bwMV_top;
      all_mv_arr = img->all_mv_top;
      all_bmv_arr= img->all_bmv_top;
      abp_type_arr = abp_type_FrArr_top;
    }
    else
    {
      frefArr = (bframe ? fw_refFrArr_bot : refFrArr_bot);
      brefArr = bw_refFrArr_bot;
      fmvArr  = (bframe ? tmp_fwMV_bot  : tmp_mv_bot);
      bmvArr  = tmp_bwMV_bot;
      all_mv_arr = img->all_mv_bot;
      all_bmv_arr= img->all_bmv_bot;
      abp_type_arr = abp_type_FrArr_bot;
    }
  }   
  if ((pdir==0 || pdir==2) && (mode!=IBLOCK && mode!=0))
  {
    for (j=j0; j<j1; j++)
    for (i=i0; i<i1; i++)
    {
      fmvArr[0][block_y+j][img->block_x+i+4] = all_mv_arr[i][j][ref][mode][0];
      fmvArr[1][block_y+j][img->block_x+i+4] = all_mv_arr[i][j][ref][mode][1];
      frefArr  [block_y+j][img->block_x+i  ] = ref;
    }
  }
  else
  {
    if(!mode&&bframe)
    {
      for (j=j0; j<j0+2; j++)
        for (i=i0; i<i0+2; i++)
        {
          if (input->direct_type)
            ref = fwdir_refFrArr[block_y+j][img->block_x+i];        
          else
          {
            ref = refFrArr[block_y+j][img->block_x+i];             
            if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
              ref = img->top_field ? refFrArr_top[block_y+j][img->block_x+i]:refFrArr_bot[block_y+j][img->block_x+i];
          }
          abp_type_arr[block_y+j][img->block_x+i] = abp_type;
          if(ref==-1)
          {
            fmvArr[0][block_y+j][img->block_x+i+4] = 0;
            fmvArr[1][block_y+j][img->block_x+i+4] = 0;
            frefArr  [block_y+j][img->block_x+i  ] = -1;
          }
          else
          {
            fmvArr[0][block_y+j][img->block_x+i+4] = all_mv_arr[i][j][ref][mode][0];
            fmvArr[1][block_y+j][img->block_x+i+4] = all_mv_arr[i][j][ref][mode][1];
            frefArr  [block_y+j][img->block_x+i  ] = ref;
          }
        }
    }
    else    
      for (j=j0; j<j0+2; j++)
        for (i=i0; i<i0+2; i++)
        {
          fmvArr[0][block_y+j][img->block_x+i+4] = 0;
          fmvArr[1][block_y+j][img->block_x+i+4] = 0;
          frefArr  [block_y+j][img->block_x+i  ] = -1;
          abp_type_arr[block_y+j][img->block_x+i] = abp_type;
        }
  }
  if ((pdir==1 || pdir==2) && (mode!=IBLOCK && mode!=0))
  {
    for (j=j0; j<j0+2; j++)
    for (i=i0; i<i0+2; i++)
    {
      if (abp_type==1 && img->type==BS_IMG)
      {
        bmvArr[0][block_y+j][img->block_x+i+4] = all_bmv_arr[i][j][0][mode][0];
        bmvArr[1][block_y+j][img->block_x+i+4] = all_bmv_arr[i][j][0][mode][1];
        brefArr  [block_y+j][img->block_x+i  ] = 0;
        abp_type_arr[block_y+j][img->block_x+i] = 1;
      }
      else if (abp_type==2 && img->type==BS_IMG)
      {
        bmvArr[0][block_y+j][img->block_x+i+4] = all_bmv_arr[i][j][1][mode][0];
        bmvArr[1][block_y+j][img->block_x+i+4] = all_bmv_arr[i][j][1][mode][1];
        brefArr  [block_y+j][img->block_x+i  ] = 1;
        abp_type_arr[block_y+j][img->block_x+i] = 2;
      }
      else
      {
        bmvArr[0][block_y+j][img->block_x+i+4] = all_bmv_arr[i][j][0][mode][0];
        bmvArr[1][block_y+j][img->block_x+i+4] = all_bmv_arr[i][j][0][mode][1];
        brefArr  [block_y+j][img->block_x+i  ] = 0;
        abp_type_arr[block_y+j][img->block_x+i] = 0;
      }
    }
  }
  else if (bframe)
  {
    if(!mode)
      for (j=j0; j<j0+2; j++)
      for (i=i0; i<i0+2; i++)
      {
        if(img->type==BS_IMG)
        {
          if (bwdir_refFrArr[block_y+j][img->block_x+i]==-1)
          {
            bmvArr[0][block_y+j][img->block_x+i+4] = all_bmv_arr[i][j][1][mode][0];
            bmvArr[1][block_y+j][img->block_x+i+4] = all_bmv_arr[i][j][1][mode][1];
          }
          else
          {
            bmvArr[0][block_y+j][img->block_x+i+4] = all_bmv_arr[i][j][bwdir_refFrArr[block_y+j][img->block_x+i]][mode][0];
            bmvArr[1][block_y+j][img->block_x+i+4] = all_bmv_arr[i][j][bwdir_refFrArr[block_y+j][img->block_x+i]][mode][1];
          }
          brefArr  [block_y+j][img->block_x+i  ] = bwdir_refFrArr[block_y+j][img->block_x+i];
        }
        else
        {
          bmvArr[0][block_y+j][img->block_x+i+4] = all_bmv_arr[i][j][0][mode][0];
          bmvArr[1][block_y+j][img->block_x+i+4] = all_bmv_arr[i][j][0][mode][1];
          if (input->direct_type)
            brefArr  [block_y+j][img->block_x+i  ] = bwdir_refFrArr[block_y+j][img->block_x+i];
          else
            brefArr  [block_y+j][img->block_x+i  ] = min(ref,0); 
        }
        abp_type_arr[block_y+j][img->block_x+i] = abp_type;
      }   
    else
      for (j=j0; j<j0+2; j++)
      for (i=i0; i<i0+2; i++)
      {
        bmvArr[0][block_y+j][img->block_x+i+4] = 0;
        bmvArr[1][block_y+j][img->block_x+i+4] = 0;
        brefArr  [block_y+j][img->block_x+i  ] = -1;
        abp_type_arr[block_y+j][img->block_x+i] = abp_type;
      }
  }
}




/*! 
 *************************************************************************************
 * \brief
 *    Mode Decision for a macroblock
 *************************************************************************************
 */
void encode_one_macroblock ()
{
  static const int  b8_mode_table[6]  = {0, 4, 5, 6, 7};         // DO NOT CHANGE ORDER !!!
  static const int  mb_mode_table[7]  = {0, 1, 2, 3, P8x8, I16MB, I4MB}; // DO NOT CHANGE ORDER !!!

  int         valid[MAXMODE];
  int         rerun, block, index, mode, i0, i1, j0, j1, pdir, ref, i, j, k, ctr16x16, dummy;
  double      qp, lambda_mode, lambda_motion, min_rdcost, rdcost = 0, max_rdcost=1e30;
  int         lambda_motion_factor;
  int         fw_mcost, bw_mcost, bid_mcost, mcost, max_mcost=(1<<30);
  int         curr_cbp_blk, cnt_nonz = 0, best_cnt_nonz = 0, best_fw_ref = 0, best_pdir;
  int         cost, min_cost, cost8x8, cost_direct=0, have_direct=0, i16mode;
  int         intra1;

  int         intra       = (((img->type==P_SLICE||img->type==SP_SLICE) && img->mb_y==img->mb_y_upd && img->mb_y_upd!=img->mb_y_intra) || img->type==I_SLICE);
  int         spframe     = (img->type==SP_SLICE);
  int         siframe     = (img->type==SI_SLICE);
  int         bframe      = (img->type==B_SLICE);
  int         write_ref   = (input->num_reference_frames>1 || input->add_ref_frame>0);
  int         runs        = (input->RestrictRef==1 && input->rdopt==2 && (img->type==P_SLICE || img->type==SP_SLICE || img->type==BS_IMG) ? 2 : 1);
  
  int         max_ref     = listXsize[0];

  int         checkref    = (input->rdopt && input->RestrictRef && (img->type==P_SLICE || img->type==SP_SLICE));
  Macroblock* currMB      = &img->mb_data[img->current_mb_nr];

  int     **ipredmodes = img->ipredmode;
  int     block_y   = img->block_y;
  int     best_bw_ref;
  int     best_abp_type;
  int     ***tmpmvs = tmp_mv;
  int     *****allmvs = img->all_mv;
  int     **refar     = refFrArr;

#ifdef ABIPRED
  int tmp_fw_ref, tmp_bw_ref;
#endif
  
  intra |= RandomIntra (img->current_mb_nr);    // Forced Pseudo-Random Intra

  if(img->MbaffFrameFlag)
  {
    if (img->fld_flag)
      max_ref = listXsize[2];

    block_y   = img->field_block_y;
    if(img->top_field)
    {
      tmpmvs = tmp_mv_top;
      allmvs = img->all_mv_top;
      refar  = refFrArr_top;
    }
    else
    {
      tmpmvs = tmp_mv_bot;
      allmvs = img->all_mv_bot;
      refar  = refFrArr_bot;
    }
  }

  //===== SET VALID MODES =====
  valid[I4MB]   = 1;
  valid[I16MB]  = 0;
//  valid[I16MB]  = 1;

  // HS: I'm not sure when the Intra Mode on 8x8 basis should be unvalid
  //     Is it o.k. for data partitioning? (where the syntax elements have to written to?)
  valid[0]      = (!intra );
  valid[1]      = (!intra && input->InterSearch16x16);
  valid[2]      = (!intra && input->InterSearch16x8);
  valid[3]      = (!intra && input->InterSearch8x16);
  valid[4]      = (!intra && input->InterSearch8x8);
  valid[5]      = (!intra && input->InterSearch8x4);
  valid[6]      = (!intra && input->InterSearch4x8);
  valid[7]      = (!intra && input->InterSearch4x4);
  valid[P8x8]   = (valid[4] || valid[5] || valid[6] || valid[7]);
  valid[12]     = (siframe);

  //===== SET LAGRANGE PARAMETERS =====
  if (input->rdopt)
  {
    qp            = (double)img->qp - SHIFT_QP;
//    lambda_mode   = 0.85 * pow (2, qp/3.0) * (bframe? 4.0:spframe?max(1.4,min(3.0,(qp / 12.0))):1.0);  // ????? WHY FACTOR 4 FOR SP-FRAME ?????
    if (input->successive_Bframe>0)
      lambda_mode   = 0.68 * pow (2, qp/3.0) * (img->type==B_SLICE? max(2.00,min(4.00,(qp / 6.0))):spframe?max(1.4,min(3.0,(qp / 12.0))):1.0);  
    else
      lambda_mode   = 0.85 * pow (2, qp/3.0) * (img->type==B_SLICE? 4.0:spframe?max(1.4,min(3.0,(qp / 12.0))):1.0);  
    lambda_motion = sqrt (lambda_mode);
  }
  else
  {
    lambda_mode = lambda_motion = QP2QUANT[max(0,img->qp-SHIFT_QP)];
  }
  lambda_motion_factor = LAMBDA_FACTOR (lambda_motion);


  for (rerun=0; rerun<runs; rerun++)
  {
    if (runs==2)
    {
      if (rerun==0)   input->rdopt=1;
      else            input->rdopt=2;
    }

    // reset chroma intra predictor to default
    currMB->c_ipred_mode = DC_PRED_8;

    if (!intra)
    {
      //===== set direct motion vectors =====
      if (bframe)
      {
        Get_Direct_Motion_Vectors ();
      }

      if (valid[P8x8])
      {
        cost8x8 = 0;

        //===== store coding state of macroblock =====
        store_coding_state (cs_mb);

        //=====  LOOP OVER 8x8 SUB-PARTITIONS  (Motion Estimation & Mode Decision) =====
        for (cbp8x8=cbp_blk8x8=cnt_nonz_8x8=0, block=0; block<4; block++)
        {
          //--- set coordinates ---
          j0 = ((block/2)<<3);    j1 = (j0>>2);
          i0 = ((block%2)<<3);    i1 = (i0>>2);

          //=====  LOOP OVER POSSIBLE CODING MODES FOR 8x8 SUB-PARTITION  =====
          for (min_cost=(1<<20), min_rdcost=1e30, index=(bframe?0:1); index<5; index++)
          {
            if (valid[mode=b8_mode_table[index]])
            {
              curr_cbp_blk = 0;

              if (mode==0)
              {
               //--- Direct Mode ---
                if (!input->rdopt)
                {
                  cost_direct += (cost = Get_Direct_Cost8x8 ( block, lambda_mode ));
                  have_direct ++;
                }
                best_fw_ref = -1;
                best_pdir   =  2;
#if 0
                if (input->WeightedBiprediction==0 && input->StoredBPictures && img->type==BS_IMG) 
                {
                  if ((input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode) || input->InterlaceCodingOption == FIELD_CODING) 
                  {
                    best_fw_ref = 3;
                    best_bw_ref = 1;
                  }
                  else
                  {
                    best_fw_ref = 1;
                    best_bw_ref = 0;
                  }
                }
                else if (input->WeightedBiprediction==2 && input->StoredBPictures&& img->type==BS_IMG) 
                {
                  if ((input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode) || input->InterlaceCodingOption == FIELD_CODING) 
                  {
                    best_fw_ref = 1;
                    best_bw_ref = 3;
                  }
                  else
                  {
                    best_fw_ref = 0;
                    best_bw_ref = 1;
                  }
                }
                else if (img->type==B_SLICE)
                {
                  best_bw_ref = 0;
                }
                else
                {
                  best_bw_ref = 0;
                }
#endif
              } // if (mode==0)
              else
              {
                //--- motion estimation for all reference frames ---
                PartitionMotionSearch (mode, block, lambda_motion);

                //--- get cost and reference frame for forward prediction ---
                for (fw_mcost=max_mcost, ref=0; ref<max_ref; ref++) // Tian Dong. PLUS1. -adjust_ref. June 06, 2002
                {
                  if (!checkref || ref==0 || CheckReliabilityOfRef (block, ref, mode))
                  {
                    mcost  = (input->rdopt ? write_ref ? REF_COST_FWD (lambda_motion_factor, ref) : 0 : (int)(2*lambda_motion*min(ref,1)));

                    mcost += motion_cost[mode][ref+1][block];
                    if (mcost < fw_mcost)
                    {
                      fw_mcost    = mcost;
                      best_fw_ref = ref;
                    }
                  }
                }

                if (bframe)
                {
                  if (img->type==BS_IMG)
                  {
                    //--- get cost for bidirectional prediction ---
#if BIPRED_SIMPLE
                    tmp_fw_ref = best_fw_ref;
#endif
                    bid_mcost =  ABIDPartitionCost (mode, block, &tmp_fw_ref, &tmp_bw_ref, lambda_motion_factor, &best_abp_type);
 //                   bid_mcost =  BBIDPartitionCost (mode, block, &tmp_fw_ref, &tmp_bw_ref, lambda_motion_factor, &best_abp_type,  useABT, lambda_motion);

                    best_bw_ref = 0;
                    bw_mcost = 1<<20;
                  }
                  else
                  {
                    bid_mcost  = (input->rdopt ? write_ref ? REF_COST (lambda_motion_factor, best_fw_ref) : 0 : (int)(2*lambda_motion*min(best_fw_ref,1)));
                    bw_mcost   = motion_cost[mode][0][block];
                    bid_mcost += BIDPartitionCost (mode, block, best_fw_ref, lambda_motion_factor );
                  }

                  //--- get prediction direction ----
                  if      (fw_mcost<=bw_mcost && fw_mcost<=bid_mcost)
                  {
                    best_pdir = 0;
                    cost = fw_mcost;
                    best_bw_ref = 0;
                  }
                  else if (bw_mcost<=fw_mcost && bw_mcost<=bid_mcost)
                  {
                    best_pdir = 1;
                    cost = bw_mcost;
                    best_bw_ref = 0;
                  }
                  else
                  {
                    best_pdir = 2;
                    cost = bid_mcost;
                    if (img->type==BS_IMG)
                    {
                      if ((input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode) || input->InterlaceCodingOption == FIELD_CODING) 
                      {
                          best_fw_ref = tmp_fw_ref;
                          best_bw_ref = tmp_bw_ref;
                      }
                      else
                      {
                          best_fw_ref = tmp_fw_ref;
                          best_bw_ref = tmp_bw_ref;
                      }
 
                    }
                    else
                    {
                      best_bw_ref = 0;
                    }
                  }
                } // if (bframe)
                else
                {
                  best_pdir = 0;
                  cost      = fw_mcost;
                }
              } // if (mode!=0)

              //--- store coding state before coding with current mode ---
              store_coding_state (cs_cm);

              if (input->rdopt)
              {
                //--- get and check rate-distortion cost ---
                if (img->type == BS_IMG)
                { 
                  rdcost = RDCost_for_8x8blocks (&cnt_nonz, &curr_cbp_blk, lambda_mode,
                                                 block, mode, best_pdir, best_fw_ref, best_bw_ref);
                }
                else
                {
                  rdcost = RDCost_for_8x8blocks (&cnt_nonz, &curr_cbp_blk, lambda_mode,
                                                 block, mode, best_pdir, best_fw_ref, 0);
                }
              }
              else
              {
                cost += (REF_COST (lambda_motion_factor, B8Mode2Value (mode, best_pdir)) - 1);
              }

              //--- set variables if best mode has changed ---
              if (( input->rdopt && rdcost < min_rdcost) ||
                  (!input->rdopt && cost   < min_cost  )   )
              {
                min_cost                   = cost;
                min_rdcost                 = rdcost;
                best8x8mode        [block] = mode;
                best8x8pdir  [P8x8][block] = best_pdir;
                best8x8ref   [P8x8][block] = best_fw_ref;
                best8x8bwref   [P8x8][block] = best_bw_ref;
                best8x8abp_type[P8x8][block] = best_abp_type;

                //--- store number of nonzero coefficients ---
                best_cnt_nonz  = cnt_nonz;

                if (input->rdopt)
                {
                  //--- store block cbp ---
                  cbp_blk8x8    &= (~(0x33 << (((block>>1)<<3)+((block%2)<<1)))); // delete bits for block
                  cbp_blk8x8    |= curr_cbp_blk;

                  //--- store coefficients ---
                  for (k=0; k< 4; k++)
                  for (j=0; j< 2; j++)
                  for (i=0; i<65; i++)  cofAC8x8[block][k][j][i] = img->cofAC[block][k][j][i]; // 18->65 for ABT

                  //--- store reconstruction and prediction ---
                  for (j=j0; j<j0+8; j++)
                  for (i=i0; i<i0+8; i++)
                  {
                    rec_mbY8x8[j][i] = enc_picture->imgY[img->pix_y+j][img->pix_x+i];
                    mpr8x8    [j][i] = img->mpr[i][j];
                  }
                }

                //--- store coding state ---
                store_coding_state (cs_b8);
              } // if (rdcost <= min_rdcost)

              //--- re-set coding state as it was before coding with current mode was performed ---
              reset_coding_state (cs_cm);
            } // if (valid[mode=b8_mode_table[index]])
          } // for (min_rdcost=1e30, index=(bframe?0:1); index<6; index++)

          cost8x8 += min_cost;

          if (!input->rdopt)
          {
            mode = best8x8mode[block];
            pdir = best8x8pdir[P8x8][block];

            curr_cbp_blk  = 0;
            best_cnt_nonz = LumaResidualCoding8x8 (&dummy, &curr_cbp_blk, block,
                                                   (pdir==0||pdir==2?mode:0),
                                                   (pdir==1||pdir==2?mode:0),
                                                   (mode!=0?best8x8ref[P8x8][block]:max(0,refar[block_y+j1][img->block_x+i1])),
                                                   (mode!=0?best8x8bwref[P8x8][block]:0));
            cbp_blk8x8   &= (~(0x33 << (((block>>1)<<3)+((block%2)<<1)))); // delete bits for block
            cbp_blk8x8   |= curr_cbp_blk;

            //--- store coefficients ---
            for (k=0; k< 4; k++)
            for (j=0; j< 2; j++)
            for (i=0; i<65; i++)  cofAC8x8[block][k][j][i] = img->cofAC[block][k][j][i]; // 18->65 for ABT

            //--- store reconstruction and prediction ---
            for (j=j0; j<j0+8; j++)
            for (i=i0; i<i0+8; i++)
            {
              rec_mbY8x8[j][i] = enc_picture->imgY[img->pix_y+j][img->pix_x+i];
              mpr8x8    [j][i] = img->mpr[i][j];
            }
          }

          //----- set cbp and count of nonzero coefficients ---
          if (best_cnt_nonz)
          {
            cbp8x8        |= (1<<block);
            cnt_nonz_8x8  += best_cnt_nonz;
          }
          
          mode=best8x8mode[block];
          //===== reset intra prediction modes (needed for prediction, must be stored after 8x8 mode dec.) =====
          j0 = block_y+2*(block/2);
          i0 = img->block_x+2*(block%2);
          for (j=j0; j<j0+2; j++)
            for (i=i0; i<i0+2; i++) 
              ipredmodes[i][j]         = DC_PRED;
          i0 = 4*block;
          for (i=i0; i<i0+4; i++)    currMB->intra_pred_modes[i]  = DC_PRED;

          if (block<3)
          {
            //===== set motion vectors and reference frames (prediction) =====
            SetRefAndMotionVectors (block, mode, best8x8ref[P8x8][block], best8x8pdir[P8x8][block], best8x8abp_type[P8x8][block]);

            //===== re-set reconstructed block =====
            j0   = 8*(block/2);
            i0   = 8*(block%2);
            for (j=j0; j<j0+8; j++)
            for (i=i0; i<i0+8; i++)  enc_picture->imgY[img->pix_y+j][img->pix_x+i] = rec_mbY8x8[j][i];
            /*
            if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
            {
              if(img->top_field)
              {
                for (j=j0; j<j0+8; j++)
                for (i=i0; i<i0+8; i++)  imgY_top[img->field_pix_y+j][img->pix_x+i] = rec_mbY8x8[j][i];
              }
              else
              {
                for (j=j0; j<j0+8; j++)
                for (i=i0; i<i0+8; i++)  imgY_bot[img->field_pix_y+j][img->pix_x+i] = rec_mbY8x8[j][i];
              }
            }*/
          } // if (block<3)

          //===== set the coding state after current block =====
          reset_coding_state (cs_b8);
        } // for (cbp8x8=cbp_blk8x8=cnt_nonz_8x8=0, block=0; block<4; block++)

        //===== store intra prediction modes for 8x8+ macroblock mode =====
        for (k=0, j=block_y; j<block_y+4; j++)
        for (     i=img->block_x; i<img->block_x+4; i++, k++)
        {
          b8_ipredmode       [k] = ipredmodes    [i][j];
          b8_intra_pred_modes[k] = currMB->intra_pred_modes[k];
        }

        //--- re-set coding state (as it was before 8x8 block coding) ---
        reset_coding_state (cs_mb);
      }
      else // if (valid[P8x8])
      {
        cost8x8 = (1<<20);
      }



      //===== MOTION ESTIMATION FOR 16x16, 16x8, 8x16 BLOCKS =====
      for (min_cost=cost8x8, best_mode=P8x8, mode=3; mode>0; mode--)
      {
        if (valid[mode])
        {
          for (cost=0, block=0; block<(mode==1?1:2); block++)
          {
            PartitionMotionSearch (mode, block, lambda_motion);

            //--- set 4x4 block indizes (for getting MV) ---
            j = (block==1 && mode==2 ? 2 : 0);
            i = (block==1 && mode==3 ? 2 : 0);

            //--- get cost and reference frame for forward prediction ---
            for (fw_mcost=max_mcost, ref=0; ref<max_ref; ref++)
            {
              if (!checkref || ref==0 || CheckReliabilityOfRef (block, ref, mode))
              {
                mcost  = (input->rdopt ? write_ref ? REF_COST_FWD (lambda_motion_factor, ref) : 0 : (int)(2*lambda_motion*min(ref,1)));

                mcost += motion_cost[mode][ref+1][block];
                if (mcost < fw_mcost)
                {
                  fw_mcost    = mcost;
                  best_fw_ref = ref;
                }
              }
            }

            if (bframe)
            {
              //--- get cost for bidirectional prediction ---
              if (img->type==BS_IMG)
              {
#if BIPRED_SIMPLE
                tmp_fw_ref = best_fw_ref;
#endif
                bid_mcost =  BBIDPartitionCost (mode, block, &tmp_fw_ref, &tmp_bw_ref, lambda_motion_factor, &best_abp_type, (int)lambda_motion);

                best_bw_ref = 0;
                bw_mcost = 1<<20;
              }
              else
              {
                bw_mcost   = motion_cost[mode][0][block];
                bid_mcost  = (input->rdopt ? write_ref ? REF_COST (lambda_motion_factor, best_fw_ref) : 0 : (int)(2*lambda_motion*min(best_fw_ref,1)));
                bid_mcost += BIDPartitionCost (mode, block, best_fw_ref, lambda_motion_factor);
                best_abp_type = 1;
              }

              //--- get prediction direction ----
              if (fw_mcost<=bw_mcost && fw_mcost<=bid_mcost)
              {
                best_pdir = 0;
                best_bw_ref = 0;
                cost += fw_mcost;
              }
              else if (bw_mcost<=fw_mcost && bw_mcost<=bid_mcost)
              {
                best_pdir = 1;
                cost += bw_mcost;
                best_bw_ref = 0;
              }
              else
              {
                best_pdir = 2;
                cost += bid_mcost;
                if (img->type==BS_IMG)
                {
                  if ((input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode) || input->InterlaceCodingOption == FIELD_CODING) 
                  {
                      best_fw_ref = tmp_fw_ref;
                      best_bw_ref = tmp_bw_ref;
                  }
                  else
                  {
                      best_fw_ref = tmp_fw_ref;
                      best_bw_ref = tmp_bw_ref;
                  }
                }
                else
                {
                  best_abp_type = 1;
                  best_bw_ref = 0;
                }
              }
            }
            else // if (bframe)
            {
              best_pdir  = 0;
              cost      += fw_mcost;
            }


            //----- set reference frame and direction parameters -----
            if (mode==3)
            {
              best8x8ref [3][  block] = best8x8ref [3][  block+2] = best_fw_ref;
              best8x8pdir[3][  block] = best8x8pdir[3][  block+2] = best_pdir;
              best8x8bwref   [3][block] = best8x8bwref   [3][block+2] = best_bw_ref;
              best8x8abp_type[3][block] = best8x8abp_type[3][block+2] = best_abp_type;
            }
            else if (mode==2)
            {
              best8x8ref [2][2*block] = best8x8ref [2][2*block+1] = best_fw_ref;
              best8x8pdir[2][2*block] = best8x8pdir[2][2*block+1] = best_pdir;
              best8x8bwref   [2][2*block] = best8x8bwref   [2][2*block+1] = best_bw_ref;
              best8x8abp_type[2][2*block] = best8x8abp_type[2][2*block+1] = best_abp_type;
            }
            else
            {
              best8x8ref [1][0] = best8x8ref [1][1] = best8x8ref [1][2] = best8x8ref [1][3] = best_fw_ref;
              best8x8pdir[1][0] = best8x8pdir[1][1] = best8x8pdir[1][2] = best8x8pdir[1][3] = best_pdir;
              best8x8bwref   [1][0] = best8x8bwref   [1][1] = best8x8bwref   [1][2] = best8x8bwref   [1][3] = best_bw_ref;
              best8x8abp_type[1][0] = best8x8abp_type[1][1] = best8x8abp_type[1][2] = best8x8abp_type[1][3] = best_abp_type;
            }

            //--- set reference frames and motion vectors ---
            if (mode>1 && block==0)   SetRefAndMotionVectors (block, mode, best_fw_ref, best_pdir, best_abp_type);

          } // for (block=0; block<(mode==1?1:2); block++)

          if (cost < min_cost)
          {
            best_mode = mode;
            min_cost  = cost;
          }
        } // if (valid[mode])
      } // for (mode=3; mode>0; mode--)

      // Find a motion vector for the Skip mode
      if((img->type == P_SLICE)||(img->type == SP_SLICE))
        FindSkipModeMotionVector ();
    }
    else // if (img->type!=I_SLICE)
    {
      min_cost = (1<<20);
    }

    if (input->rdopt)
    {
      int mb_available_up;
      int mb_available_left;
      int mb_available_up_left;
      
      min_rdcost = max_rdcost;
      
      // precompute all new chroma intra prediction modes
      IntraChromaPrediction8x8(&mb_available_up, &mb_available_left, &mb_available_up_left);

      for (currMB->c_ipred_mode=DC_PRED_8; currMB->c_ipred_mode<=PLANE_8; currMB->c_ipred_mode++)
      {
        
        // bypass if c_ipred_mode is not allowed
        if ((currMB->c_ipred_mode==VERT_PRED_8 && !mb_available_up) ||
            (currMB->c_ipred_mode==HOR_PRED_8 && !mb_available_left) ||
            (currMB->c_ipred_mode==PLANE_8 && (!mb_available_left || !mb_available_up || !mb_available_up_left)))
           continue;
           

        //===== GET BEST MACROBLOCK MODE =====
        for (ctr16x16=0, index=0; index<7; index++)
        {
          mode = mb_mode_table[index];

          //--- for INTER16x16 check all prediction directions ---
          if (mode==1 && img->type==B_SLICE)
          {
            best8x8pdir[1][0] = best8x8pdir[1][1] = best8x8pdir[1][2] = best8x8pdir[1][3] = ctr16x16;
            if (ctr16x16 < 2) index--;
            ctr16x16++;
          }

          img->NoResidueDirect = 0;

          if (valid[mode])
          {
            // bypass if c_ipred_mode not used
            SetModesAndRefframeForBlocks (mode);
            if (currMB->c_ipred_mode == DC_PRED_8 ||
              (IS_INTRA(currMB) ))
            {
              if (RDCost_for_macroblocks (lambda_mode, mode, &min_rdcost))
              {
                store_macroblock_parameters (mode);
              }
            }
          }
          if (bframe && mode == 0 && currMB->cbp && (currMB->cbp&15) != 15)
          {
            img->NoResidueDirect = 1;
            if (RDCost_for_macroblocks (lambda_mode, mode, &min_rdcost))
            {
              store_macroblock_parameters (mode);
            }
          }
        }
      }
    }
    else
    {
      if (valid[0] && bframe) // check DIRECT MODE
      {
        cost  = (have_direct?cost_direct:Get_Direct_CostMB (lambda_mode));
        cost -= (int)floor(16*lambda_motion+0.4999);
        if (cost <= min_cost)
        {
          min_cost  = cost;
          best_mode = 0;
        }
      }
      if (valid[I4MB]) // check INTRA4x4
      {
        currMB->cbp = Mode_Decision_for_Intra4x4Macroblock (lambda_mode, &cost);
        if (cost <= min_cost)
        {
          min_cost  = cost;
          best_mode = I4MB;
        }
      }
      if (valid[I16MB]) // check INTRA16x16
      {
        intrapred_luma_16x16 ();
        cost = find_sad_16x16 (&i16mode);
        if (cost < min_cost)
        {
          best_mode   = I16MB;
          currMB->cbp = dct_luma_16x16 (i16mode);
        }
      }
    }

    if (rerun==0)
    {
      intra1 = (currMB->mb_type==I16MB || currMB->mb_type==I4MB ? 1 : 0);
    }
  } // for (rerun=0; rerun<runs; rerun++)
  
  if (input->rdopt)
  {
    set_stored_macroblock_parameters ();
  }
  else
  {
    //===== set parameters for chosen mode =====
    SetModesAndRefframeForBlocks (best_mode);
    if (best_mode==P8x8)
    {
      SetCoeffAndReconstruction8x8 (currMB);
    }
    else
    {
      if (best_mode!=I4MB)
      {
        for (k=0, j=block_y; j<block_y+4; j++)
        for (     i=img->block_x; i<img->block_x+4; i++, k++)
        {
          ipredmodes    [i][j] = DC_PRED;
          currMB->intra_pred_modes[k] = DC_PRED;
        }
        if (best_mode!=I16MB)
        {
          LumaResidualCoding ();
        }
      }
    }
    // precompute all chroma intra prediction modes
    IntraChromaPrediction8x8(NULL, NULL, NULL);
    img->i16offset = 0;
    dummy = 0;
    ChromaResidualCoding (&dummy);
    if (best_mode==I16MB)
    {
      img->i16offset = I16Offset  (currMB->cbp, i16mode);
    }
    SetMotionVectorsMB (currMB, bframe);

    //===== check for SKIP mode =====
    if ((img->type==P_SLICE || img->type==SP_SLICE) && best_mode==1 && currMB->cbp==0 &&
        refFrArr [block_y][img->block_x  ]==0 &&
        tmpmvs[0][block_y][img->block_x+4]==allmvs[0][0][0][0][0] &&
        tmpmvs[1][block_y][img->block_x+4]==allmvs[0][0][0][0][1]               )
    {
      if(!mb_adaptive)
        currMB->mb_type=currMB->b8mode[0]=currMB->b8mode[1]=currMB->b8mode[2]=currMB->b8mode[3]=0;
      else if(mb_adaptive && tmpmvs[0][block_y][img->block_x+4] == 0 && tmpmvs[1][block_y][img->block_x+4]== 0)
      {
        currMB->mb_type=currMB->b8mode[0]=currMB->b8mode[1]=currMB->b8mode[2]=currMB->b8mode[3]=0;
        SetModesAndRefframeForBlocks(0);
      }
    }

    if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
      set_mbaff_parameters();
  }

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
  {
    if(input->rdopt)
      rdopt->min_rdcost = min_rdcost;
    else
      rdopt->min_rdcost = min_cost;
    
    
    if(img->field_mode && img->top_field && (img->type != B_SLICE && img->type != BS_IMG) && best_mode == 0)
      TopFieldIsSkipped=1;  //set top field MB skipped to 1 for skipped top field MBs
    else if(img->field_mode && img->top_field && (img->type == B_SLICE || img->type == BS_IMG) && best_mode == 0)
      TopFieldIsSkipped = currMB->cbp ? 0:1;  // direct mode may or may not be skipped
    
    if(!img->field_mode && (img->type != B_SLICE && img->type != BS_IMG) && best_mode == 0 && !img->mb_y%2)
      TopFrameIsSkipped=1;  
    else if(!img->field_mode && (img->type == B_SLICE || img->type == BS_IMG) && best_mode == 0 && !img->mb_y%2)
      TopFrameIsSkipped = currMB->cbp ? 0:1;  // direct mode may or may not be skipped
    
    if(!input->rdopt && img->field_mode && !img->top_field && TopFieldIsSkipped)
      if(img->type == B_SLICE  && (best_mode == 0 && currMB->cbp == 0 && input->symbol_mode!=CABAC) )
        rdopt->min_rdcost = 1e30;  // don't allow coding of an MB pair as skip in field mode
  }
/*
  if (input->rdopt==2 && !bframe)
  {
    for (j=0 ;j<input->NoOfDecoders; j++)  DeblockMb(img, decs->decY_best[j], NULL);
  }
*/

  //===== init and update number of intra macroblocks =====
  if (img->current_mb_nr==0)                      intras=0;
  if ((img->type==P_SLICE || img->type==SP_SLICE || img->type==BS_IMG) && IS_INTRA(currMB))   intras++;

  //===== Decide if this MB will restrict the reference frames =====
  if (input->RestrictRef==1)
  {
    if (input->rdopt==1)
    {
      refresh_map[2*img->mb_y  ][2*img->mb_x  ] = (intra ? 1 : 0);
      refresh_map[2*img->mb_y  ][2*img->mb_x+1] = (intra ? 1 : 0);
      refresh_map[2*img->mb_y+1][2*img->mb_x  ] = (intra ? 1 : 0);
      refresh_map[2*img->mb_y+1][2*img->mb_x+1] = (intra ? 1 : 0);
    }
    else if (input->rdopt==2)
    {
      refresh_map[2*img->mb_y  ][2*img->mb_x  ] = (intra1==0 && (currMB->mb_type==I16MB || currMB->mb_type==I4MB) ? 1 : 0);
      refresh_map[2*img->mb_y  ][2*img->mb_x+1] = (intra1==0 && (currMB->mb_type==I16MB || currMB->mb_type==I4MB) ? 1 : 0);
      refresh_map[2*img->mb_y+1][2*img->mb_x  ] = (intra1==0 && (currMB->mb_type==I16MB || currMB->mb_type==I4MB) ? 1 : 0);
      refresh_map[2*img->mb_y+1][2*img->mb_x+1] = (intra1==0 && (currMB->mb_type==I16MB || currMB->mb_type==I4MB) ? 1 : 0);
    }
  }
  else if (input->RestrictRef==2)
  {
    refresh_map[2*img->mb_y  ][2*img->mb_x  ] = (currMB->mb_type==I16MB || currMB->mb_type==I4MB ? 1 : 0);
    refresh_map[2*img->mb_y  ][2*img->mb_x+1] = (currMB->mb_type==I16MB || currMB->mb_type==I4MB ? 1 : 0);
    refresh_map[2*img->mb_y+1][2*img->mb_x  ] = (currMB->mb_type==I16MB || currMB->mb_type==I4MB ? 1 : 0);
    refresh_map[2*img->mb_y+1][2*img->mb_x+1] = (currMB->mb_type==I16MB || currMB->mb_type==I4MB ? 1 : 0);
  }
}


void set_mbaff_parameters()
{
  int  i, j, k, l;
  Macroblock  *currMB  = &img->mb_data[img->current_mb_nr];
  int         mode     = best_mode;
  int         bframe   = (img->type==B_SLICE);
  int         **frefar = (img->type==B_SLICE ? fw_refFrArr : refFrArr);
  int         **brefar = bw_refFrArr;
  int     **ipredmodes = img->ipredmode;
  int     block_y = img->block_y; 

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    block_y = img->field_block_y;
    if(img->top_field)
    {
      frefar = (img->type == B_SLICE) ? fw_refFrArr_top : refFrArr_top;
      brefar = bw_refFrArr_top;
    }
    else
    {
      frefar = (img->type == B_SLICE) ? fw_refFrArr_bot : refFrArr_bot;
      brefar = bw_refFrArr_bot;

    }
  }

  //===== reconstruction values =====
  for (j=0; j<16; j++)
  for (i=0; i<16; i++)
    rdopt->rec_mbY[j][i]       = enc_picture->imgY[img->pix_y+j][img->pix_x+i]; 

  for (j=0; j<8; j++)
  for (i=0; i<8; i++)
  {
    rdopt->rec_mbU[j][i]           = enc_picture->imgUV[0][img->pix_c_y+j][img->pix_c_x+i];   
    rdopt->rec_mbV[j][i]           = enc_picture->imgUV[1][img->pix_c_y+j][img->pix_c_x+i];    
  }

  //===== coefficients and cbp =====
  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
  {
    rdopt->mode = mode;
    rdopt->i16offset = img->i16offset;  // For MBINTLC  -Rajeev
    rdopt->cbp = currMB->cbp;
    rdopt->cbp_blk = currMB->cbp_blk;
    rdopt->mb_type  =currMB->mb_type;
    if(rdopt->mb_type == 0 && mode != 0)
    {
      mode=0;
      rdopt->mode=0;
    }


    for(i=0;i<6;i++)
      for(j=0;j<4;j++)
        for(k=0;k<2;k++)
          for(l=0;l<18;l++)
            rdopt->cofAC[i][j][k][l] = img->cofAC[i][j][k][l];

    for(i=0;i<3;i++)
        for(k=0;k<2;k++)
          for(l=0;l<18;l++)
            rdopt->cofDC[i][k][l] = img->cofDC[i][k][l];
  }


  for (i=0; i<4; i++)
  {
  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
  {
    rdopt->b8mode[i]  = currMB->b8mode[i];                  
    rdopt->b8pdir[i]  = currMB->b8pdir[i];                  
  }
  }

  //==== reference frames =====
  for (j=0; j<4; j++)
  for (i=0; i<4; i++)
  {
    rdopt->frefar[j][i]            = frefar[block_y+j][img->block_x+i];
  }
  if (bframe)
  {
    for (j=0; j<4; j++)
    for (i=0; i<4; i++)
    {
      rdopt->brefar[j][i]            = brefar[block_y+j][img->block_x+i];

    }
  }


  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
  {
    for (k=0, j=block_y; j<block_y+4; j++)
      for (     i=img->block_x; i<img->block_x+4; i++, k++)
      {
        rdopt->ipredmode[i][j] = ipredmodes[i][j];
        rdopt->intra_pred_modes[k] = currMB->intra_pred_modes[k];
      }
  }
}
