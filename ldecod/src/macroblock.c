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
 ***********************************************************************
 * \file macroblock.c
 *
 * \brief
 *     Decode a Macroblock
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Inge Lille-Langøy               <inge.lille-langoy@telenor.com>
 *    - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *    - Jani Lainema                     <jani.lainema@nokia.com>
 *    - Sebastian Purreiter             <sebastian.purreiter@mch.siemens.de>
 *    - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
 *    - Detlev Marpe                      <marpe@hhi.de>
 *    - Gabi Blaettermann              <blaetter@hhi.de>
 *    - Ye-Kui Wang                      <wangy@cs.tut.fi>
 ***********************************************************************
*/

#include "contributors.h"

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "global.h"
#include "mbuffer.h"
#include "elements.h"
#include "macroblock.h"


/*!
 ************************************************************************
 * \brief
 *    Checks the availability of neighboring macroblocks of
 *    the current macroblock for prediction and context determination;
 *    marks the unavailable MBs for intra prediction in the
 *    ipredmode-array by -1. Only neighboring MBs in the causal
 *    past of the current MB are checked.
 ************************************************************************
 */
void CheckAvailabilityOfNeighbors(struct img_par *img)
{
  int i,j;
  const int mb_width = img->width/MB_BLOCK_SIZE;
  const int mb_nr = img->current_mb_nr;
  Macroblock *currMB = &img->mb_data[mb_nr];

  // mark all neighbors as unavailable
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      img->mb_data[mb_nr].mb_available[i][j]=NULL;
  img->mb_data[mb_nr].mb_available[1][1]=currMB; // current MB

  // Check MB to the left
  if(img->pix_x >= MB_BLOCK_SIZE)
  {
    int remove_prediction = currMB->slice_nr != img->mb_data[mb_nr-1].slice_nr;
    if(img->UseConstrainedIntraPred)
      remove_prediction = (remove_prediction || img->intra_mb[mb_nr-1] ==0);
    if(remove_prediction)
    {
      img->ipredmode[img->block_x][img->block_y+1] = -1;
      img->ipredmode[img->block_x][img->block_y+2] = -1;
      img->ipredmode[img->block_x][img->block_y+3] = -1;
      img->ipredmode[img->block_x][img->block_y+4] = -1;
    } else
      currMB->mb_available[1][0]=&(img->mb_data[mb_nr-1]);
  }


  // Check MB above
  if(img->pix_y >= MB_BLOCK_SIZE)
  {
    int remove_prediction = currMB->slice_nr != img->mb_data[mb_nr-mb_width].slice_nr;
    if(img->UseConstrainedIntraPred)
      remove_prediction = (remove_prediction || img->intra_mb[mb_nr-mb_width] ==0);
    if(remove_prediction)
    {
      img->ipredmode[img->block_x+1][img->block_y] = -1;
      img->ipredmode[img->block_x+2][img->block_y] = -1;
      img->ipredmode[img->block_x+3][img->block_y] = -1;
      img->ipredmode[img->block_x+4][img->block_y] = -1;
    } else
      currMB->mb_available[0][1]=&(img->mb_data[mb_nr-mb_width]);
  }

  // Check MB left above
  if(img->pix_x >= MB_BLOCK_SIZE && img->pix_y  >= MB_BLOCK_SIZE )
  {
    if(currMB->slice_nr == img->mb_data[mb_nr-mb_width-1].slice_nr)
      img->mb_data[mb_nr].mb_available[0][0]=&(img->mb_data[mb_nr-mb_width-1]);
  }

  // Check MB right above
  if(img->pix_y >= MB_BLOCK_SIZE && img->pix_x < (img->width-MB_BLOCK_SIZE ))
  {
    if(currMB->slice_nr == img->mb_data[mb_nr-mb_width+1].slice_nr)
      currMB->mb_available[0][2]=&(img->mb_data[mb_nr-mb_width+1]);
  }
}

/*!
 ************************************************************************
 * \brief
 *    initializes the current macroblock
 ************************************************************************
 */
void start_macroblock(struct img_par *img,struct inp_par *inp)
{
  int i,j,k,l;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  // WYK: Oct. 8, 2001, start ...
  // The following is moved and modified from exit_macroblock(), 
  // to make the decoding process correct when some macroblocks are lost
  /* Update coordinates of the current macroblock */
  img->mb_x = (img->current_mb_nr)%(img->width/MB_BLOCK_SIZE);
  img->mb_y = (img->current_mb_nr)/(img->width/MB_BLOCK_SIZE);
  
  /* Define vertical positions */
  img->block_y = img->mb_y * BLOCK_SIZE;      /* luma block position */
  img->pix_y   = img->mb_y * MB_BLOCK_SIZE;   /* luma macroblock position */
  img->pix_c_y = img->mb_y * MB_BLOCK_SIZE/2; /* chroma macroblock position */
  
  /* Define horizontal positions */
  img->block_x = img->mb_x * BLOCK_SIZE;      /* luma block position */
  img->pix_x   = img->mb_x * MB_BLOCK_SIZE;   /* luma pixel position */
  img->pix_c_x = img->mb_x * MB_BLOCK_SIZE/2; /* chroma pixel position */
  //WYK: Oct. 8, 2001, ... end

  // Save the slice number of this macroblock. When the macroblock below
  // is coded it will use this to decide if prediction for above is possible
  currMB->slice_nr = img->current_slice_nr;
  
  // If MB is next to a slice boundary, mark neighboring blocks unavailable for prediction
  CheckAvailabilityOfNeighbors(img);

  // Reset syntax element entries in MB struct
  currMB->mb_type = 0;
  currMB->ref_frame = 0;
  currMB->predframe_no = 0;
  currMB->delta_quant = 0;

  currMB->cbp     = 0;
  currMB->cbp_blk = 0;

  for (l=0; l < 2; l++)
    for (j=0; j < BLOCK_MULTIPLE; j++)
      for (i=0; i < BLOCK_MULTIPLE; i++)
        for (k=0; k < 2; k++)
          currMB->mvd[l][j][i][k] = 0;

  for (i=0; i < (BLOCK_MULTIPLE*BLOCK_MULTIPLE); i++)
    currMB->intra_pred_modes[i] = 0;
}

/*!
 ************************************************************************
 * \brief
 *    set coordinates of the next macroblock
 *    check end_of_slice condition (have to implement)
 ************************************************************************
 */
int exit_macroblock(struct img_par *img,struct inp_par *inp)
{
  const int number_mb_per_row = img->width / MB_BLOCK_SIZE ;
  Slice *currSlice = img->currentSlice;

  // Update coordinates of the next macroblock
  img->mb_x++;
  if (img->mb_x == number_mb_per_row) // next row of MBs
  {
    img->mb_x = 0; // start processing of next row
    img->mb_y++;
  }
  img->current_mb_nr++;

  // Define vertical positions
  img->block_y = img->mb_y * BLOCK_SIZE;      // luma block position
  img->pix_y   = img->mb_y * MB_BLOCK_SIZE;   // luma macroblock position
  img->pix_c_y = img->mb_y * MB_BLOCK_SIZE/2; // chroma macroblock position

  // Define horizontal positions
  img->block_x = img->mb_x * BLOCK_SIZE;      // luma block position
  img->pix_x   = img->mb_x * MB_BLOCK_SIZE;   // luma pixel position
  img->pix_c_x = img->mb_x * MB_BLOCK_SIZE/2; // chroma pixel position

  if (img->current_mb_nr == img->max_mb_nr)
  {
    if (currSlice->next_header != EOS)
      currSlice->next_header = SOP;
    return TRUE;
  }
  // ask for last mb in the slice  UVLC
  else
  {
    if(nal_startcode_follows(img, inp) == FALSE)
      return FALSE;
    if(img->type == INTRA_IMG || img->type == SP_IMG_1|| img->type == SP_IMG_MULT || inp->symbol_mode == CABAC)
      return TRUE;
    if(img->cod_counter<=0)
      return TRUE;
    return FALSE;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Interpret the mb mode for P-Frames
 ************************************************************************
 */
void interpret_mb_mode_P(struct img_par *img)
{
  const int ICBPTAB[6] = {0,16,32,15,31,47};
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];


  currMB->intraOrInter           = INTER_MB ;

  if (img->mb_mode == INTRA_MB)   // 4x4 intra
  {
    currMB->intraOrInter         = INTRA_MB_4x4 ;
    currMB->mb_imode = img->imod = INTRA_MB_OLD ;
  }
  if (img->mb_mode > INTRA_MB)    // 16x16 intra
  {
    currMB->intraOrInter         = INTRA_MB_16x16 ;
    img->imod = currMB->mb_imode = INTRA_MB_NEW;     // mod0=img->mb_mode-1;kmod=mod0 & 3;cbp = ICBPTAB[mod0/4];
    currMB->intra_pred_modes[0] = (img->mb_mode - INTRA_MB-1) & 3;
    currMB->cbp = ICBPTAB[(img->mb_mode - INTRA_MB-1)>>2];
  }
  if (img->mb_mode < INTRA_MB)    // inter prediction mode (block shape)
    img->imod = currMB->mb_imode = INTRA_MB_INTER;   // intra in inter frame
}

/*!
 ************************************************************************
 * \brief
 *    Interpret the mb mode for I-Frames
 ************************************************************************
 */
void interpret_mb_mode_I(struct img_par *img)
{
  const int ICBPTAB[6] = {0,16,32,15,31,47};
  Macroblock *currMB   = &img->mb_data[img->current_mb_nr];

  if (img->mb_mode == 0)
  {
    currMB->intraOrInter         = INTRA_MB_4x4 ;
    img->imod = currMB->mb_imode = INTRA_MB_OLD;     // 4x4 intra
  }
  else
  {
    currMB->intraOrInter         = INTRA_MB_16x16 ;
    img->imod = currMB->mb_imode = INTRA_MB_NEW;     // 16x16 intra */ /* mod0=img->mb_mode-1;kmod=mod0 & 3;cbp = ICBPTAB[mod0/4];
    currMB->intra_pred_modes[0] = (img->mb_mode - 1) & 3;
    currMB->cbp = ICBPTAB[(img->mb_mode - 1)>>2];
  }
}

/*!
 ************************************************************************
 * \brief
 *    Interpret the mb mode for B-Frames
 ************************************************************************
 */
void interpret_mb_mode_B(struct img_par *img)
{
  const int ICBPTAB[6] = {0,16,32,15,31,47};
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  currMB->intraOrInter         = INTER_MB ;
  if (img->mb_mode == INTRA_MB_B) // 4x4 intra
  {
    currMB->intraOrInter         = INTRA_MB_4x4 ;
    img->imod = currMB->mb_imode = INTRA_MB_OLD;
  }
  if (img->mb_mode > INTRA_MB_B)  // 16x16 intra
  {
    currMB->intraOrInter         = INTRA_MB_16x16 ;
    img->imod = currMB->mb_imode = INTRA_MB_NEW;
    currMB->intra_pred_modes[0] = (img->mb_mode - INTRA_MB_B-1) & 3;
    currMB->cbp = ICBPTAB[(img->mb_mode - INTRA_MB_B-1)>>2];
  }

  if (img->mb_mode < INTRA_MB_B)
  {
    if(img->mb_mode == 0)
      img->imod = currMB->mb_imode = B_Direct;
    else if(img->mb_mode == 3)
      img->imod = currMB->mb_imode = B_Bidirect;
    else if(img->mb_mode==1 || (img->mb_mode>3 && img->mb_mode%2==0))
      img->imod = currMB->mb_imode = B_Forward;
    else if(img->mb_mode==2 || (img->mb_mode>4 && img->mb_mode%2==1))
      img->imod = currMB->mb_imode = B_Backward;
  }
}

/*!
 ************************************************************************
 * \brief
 *    init macroblock I and P frames
 ************************************************************************
 */
void init_macroblock(struct img_par *img)
{
  int i,j;
  int predframe_no;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  img->mv[img->block_x+4][img->block_y][2]=img->number;

  for (i=0;i<BLOCK_SIZE;i++)
  {                           // reset vectors and pred. modes
    for(j=0;j<BLOCK_SIZE;j++)
    {
      img->mv[img->block_x+i+4][img->block_y+j][0] = 0;
      img->mv[img->block_x+i+4][img->block_y+j][1] = 0;
      img->ipredmode[img->block_x+i+1][img->block_y+j+1] = 0;
    }
  }

  currMB->ref_frame = 0;
  currMB->predframe_no = predframe_no = 0;// g.b.1;

  // Set the reference frame information for motion vector prediction
  if (img->imod == INTRA_MB_OLD || img->imod == INTRA_MB_NEW)
    for (j = 0;j < BLOCK_SIZE;j++)
      for (i = 0;i < BLOCK_SIZE;i++)
        refFrArr[img->block_y+j][img->block_x+i] = -1;
  else
    for (j = 0;j < BLOCK_SIZE;j++)
      for (i = 0;i < BLOCK_SIZE;i++)
        refFrArr[img->block_y+j][img->block_x+i] = predframe_no;
}

/*!
 ************************************************************************
 * \brief
 *    Get the syntax elements from the NAL
 ************************************************************************
 */
int read_one_macroblock(struct img_par *img,struct inp_par *inp)
{
  int i, i1, j1;

  SyntaxElement currSE;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  Slice *currSlice = img->currentSlice;
  DataPartition *dP;
  int *partMap = assignSE2partition[currSlice->dp_mode];

  int dbl_ipred_word;

  currMB->qp           = img->qp ;

  //  read MB mode *****************************************************************
  currSE.type = SE_MBTYPE;
  
  if(img->type == B_IMG_1 || img->type == B_IMG_MULT)
    dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
  else
    dP = &(currSlice->partArr[partMap[currSE.type]]);

  if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag) 
    currSE.mapping = linfo;
  else
    currSE.reading = readMB_typeInfoFromBuffer_CABAC;

  if(inp->symbol_mode == CABAC || (img->type != INTER_IMG_1 && img->type != INTER_IMG_MULT && img->type != B_IMG_1 && img->type != B_IMG_MULT))
  {
    //  read MB mode
#if TRACE
    strncpy(currSE.tracestring, "MB Type", TRACESTRING_SIZE);
#endif
    
    dP->readSyntaxElement(&currSE,img,inp,dP);
    img->mb_mode = currMB->mb_type = currSE.value1;
  } 
  else
  {
    if(img->cod_counter == -1)
    {
#if TRACE
      strncpy(currSE.tracestring, "MB runlength", TRACESTRING_SIZE);
#endif
      dP->readSyntaxElement(&currSE,img,inp,dP);
      img->cod_counter = currSE.value1;
    }
    if (img->cod_counter==0)
    {
#if TRACE
      strncpy(currSE.tracestring, "MB Type", TRACESTRING_SIZE);
#endif
      dP->readSyntaxElement(&currSE,img,inp,dP);
      if(img->type == INTER_IMG_1 || img->type == INTER_IMG_MULT)
        currSE.value1++;
      img->mb_mode = currMB->mb_type = currSE.value1;
      img->cod_counter--;
    } 
    else
    {
      img->cod_counter--;
      img->mb_mode = 0;
    }
  }

  if(img->UseConstrainedIntraPred)
  {
    if (img->type==INTER_IMG_1 || (img->type==INTER_IMG_MULT))        // inter frame
      if (img->mb_mode < 8)
        img->intra_mb[img->current_mb_nr] = 0;
  }
  
  
  //! TO for Error Concelament
  //! If we have an INTRA Macroblock and we lost the partition
  //! which contains the intra coefficients Copy MB would be better 
  //! than just a grey block.
  //! Seems to be a bit at the wrong place to do this right here, but for this case 
  //! up to now there is no other way.
  dP = &(currSlice->partArr[partMap[SE_CBP_INTRA]]);
  if(img->mb_mode == INTRA_MB && dP->bitstream->ei_flag && img->number)
  {
    img->mb_mode = COPY_MB;
    img->imod = INTRA_MB_INTER;
  }
  if(img->type == B_IMG_1 || img->type == B_IMG_MULT)
    dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
  else
    dP = &(currSlice->partArr[partMap[currSE.type]]);
  //! End TO
  
  if ((img->type==INTER_IMG_1) || (img->type==INTER_IMG_MULT))    // inter frame
    interpret_mb_mode_P(img);
  else if (img->type==INTRA_IMG)                                  // intra frame
    interpret_mb_mode_I(img);
  else if ((img->type==B_IMG_1) || (img->type==B_IMG_MULT))       // B frame
    interpret_mb_mode_B(img);
  else if ((img->type==SP_IMG_1) || (img->type==SP_IMG_MULT))    // SP frame
    interpret_mb_mode_P(img);

  if ((img->type==B_IMG_1) || (img->type==B_IMG_MULT))
    init_macroblock_Bframe(img);
  else
    init_macroblock(img);

  if (inp->symbol_mode != CABAC && img->imod==B_Direct && img->mb_mode==COPY_MB && img->cod_counter >= 0)
  {
    int i, j, iii, jjj;
    currMB->cbp = 0;
    for (i=0;i<BLOCK_SIZE;i++)
    { // reset luma coeffs
      for (j=0;j<BLOCK_SIZE;j++)
        for(iii=0;iii<BLOCK_SIZE;iii++)
          for(jjj=0;jjj<BLOCK_SIZE;jjj++)
            img->cof[i][j][iii][jjj]=0;
    }
    for (j=4;j<6;j++)
    { // reset chroma coeffs
      for (i=0;i<4;i++)
        for (iii=0;iii<4;iii++)
          for (jjj=0;jjj<4;jjj++)
            img->cof[i][j][iii][jjj]=0;
    }
    return DECODE_MB_BFRAME;
  }

  if (img->imod==INTRA_MB_INTER && img->mb_mode==COPY_MB) //keep last macroblock
  {
    return DECODE_COPY_MB;
  }

  // intra prediction modes for a macroblock 4x4 **********************************************
  if (img->imod==INTRA_MB_OLD)
  {
    currSE.type = SE_INTRAPREDMODE;
    if(img->type == B_IMG_1 || img->type == B_IMG_MULT)
      dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
    else
      dP = &(currSlice->partArr[partMap[currSE.type]]);

    if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag)
        currSE.mapping = linfo;
    else
        currSE.reading = readIntraPredModeFromBuffer_CABAC;

    for(i=0;i<MB_BLOCK_SIZE/2;i++)
    {
#if TRACE
      snprintf(currSE.tracestring, TRACESTRING_SIZE, "Intra mode ");
#endif

      dP->readSyntaxElement(&currSE,img,inp,dP);

      i1=img->block_x + 2*(i&0x01);
      j1=img->block_y + i/2;

      if (inp->symbol_mode == UVLC)
      {
        dbl_ipred_word = currSE.value1;
        // find intra prediction mode for two blocks
        img->ipredmode[i1+1][j1+1] = PRED_IPRED[img->ipredmode[i1+1][j1]+1][img->ipredmode[i1][j1+1]+1][IPRED_ORDER[dbl_ipred_word][0]];
        img->ipredmode[i1+2][j1+1] = PRED_IPRED[img->ipredmode[i1+2][j1]+1][img->ipredmode[i1+1][j1+1]+1][IPRED_ORDER[dbl_ipred_word][1]];
      }
      else
      {
        img->ipredmode[i1+1][j1+1] = PRED_IPRED[img->ipredmode[i1+1][j1]+1][img->ipredmode[i1][j1+1]+1][currSE.value1];
        img->ipredmode[i1+2][j1+1] = PRED_IPRED[img->ipredmode[i1+2][j1]+1][img->ipredmode[i1+1][j1+1]+1][currSE.value2];
      }
    }
  }

  // read inter frame vector data ********************************************************
  if ((img->type==B_IMG_1) || (img->type==B_IMG_MULT))
    readMotionInfoFromNAL_Bframe(img,inp);
  else if(img->imod==INTRA_MB_INTER)
    readMotionInfoFromNAL_Pframe(img,inp);

  // read CBP and Coeffs  ***************************************************************
  readCBPandCoeffsFromNAL(img,inp);

  return (((img->type==B_IMG_1) || (img->type==B_IMG_MULT)) ? DECODE_MB_BFRAME : DECODE_MB);
}

/*!
 ************************************************************************
 * \brief
 *    Get for a given MB of a P picture the reference frame
 *    parameter and the motion vectors from the NAL
 ************************************************************************
 */
void readMotionInfoFromNAL_Pframe(struct img_par *img,struct inp_par *inp)
{
  int i,j,k,l,m;
  int step_h,step_v;
  int curr_mvd;
  int mb_nr = img->current_mb_nr;
  Macroblock *currMB = &img->mb_data[mb_nr];
  SyntaxElement currSE;
  Slice *currSlice = img->currentSlice;
  DataPartition *dP;
  int *partMap = assignSE2partition[currSlice->dp_mode];
  int ref_frame = currMB->ref_frame;
  int predframe_no = currMB->predframe_no;


  // keep track of neighbouring macroblocks available for prediction
  int mb_width = img->width/16;
  int mb_available_up = (img->mb_y == 0) ? 0 : (currMB->slice_nr == img->mb_data[mb_nr-mb_width].slice_nr);
  int mb_available_left = (img->mb_x == 0) ? 0 : (currMB->slice_nr == img->mb_data[mb_nr-1].slice_nr);
  int mb_available_upleft  = (img->mb_x == 0 || img->mb_y == 0) ? 0 : (currMB->slice_nr == img->mb_data[mb_nr-mb_width-1].slice_nr);
  int mb_available_upright = (img->mb_x >= mb_width-1 || img->mb_y == 0) ? 0 : (currMB->slice_nr == img->mb_data[mb_nr-mb_width+1].slice_nr);

  // keep track of neighbouring blocks available for motion vector prediction
  int block_available_up, block_available_left, block_available_upright, block_available_upleft;
  int mv_a, mv_b, mv_c, mv_d;
  int mvPredType, rFrameL, rFrameU, rFrameUR;
  int ie, j4, i4, ii,jj;
  int pred_vec=0, vec;


  //  If multiple ref. frames, read reference frame for the MB *********************************
  if(img->type==INTER_IMG_MULT || img->type == SP_IMG_MULT)
  {
#if TRACE
    strncpy(currSE.tracestring,  "Reference frame no ", TRACESTRING_SIZE);
#endif
    currSE.type = SE_REFFRAME;
    dP = &(currSlice->partArr[partMap[currSE.type]]);

    if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag)
      currSE.mapping = linfo;
    else
      currSE.reading = readRefFrameFromBuffer_CABAC;

    dP->readSyntaxElement(&currSE,img,inp,dP);
    predframe_no = currMB->predframe_no = currSE.value1;
    ref_frame = currMB->ref_frame = predframe_no;

    /*!
    * \note
    * if frame lost occurs within img->buf_cycle frames and buffer of previous
    * decoded pictures is still empty, set ref_frame to last decoded
    * picture to avoid prediction from unexistent frames.
    */
//    if (ref_frame > img->number) ref_frame = 0;

    // Update the reference frame information for motion vector prediction
    for (j = 0;j < BLOCK_SIZE;j++)
      for (i = 0;i < BLOCK_SIZE;i++)
        refFrArr[img->block_y+j][img->block_x+i] = predframe_no;
  }

  // read motion vectors **********************************************************************

  currSE.type = SE_MVD;
  dP = &(currSlice->partArr[partMap[currSE.type]]);

  if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag)
    currSE.mapping = linfo_mvd;
  else
    currSE.reading = readMVDFromBuffer_CABAC;

  step_h=BLOCK_STEP[img->mb_mode][0];     // horizontal stepsize
  step_v=BLOCK_STEP[img->mb_mode][1];     // vertical stepsize

  ie=4-BLOCK_STEP[img->mb_mode][0];
  for (j=0; j < BLOCK_SIZE; j += step_v)
  {
    block_available_up = mb_available_up || (j > 0);
    j4=img->block_y+j;

    for (i=0;i < BLOCK_SIZE; i += step_h)
    {
      i4=img->block_x+i;
      // first make mv-prediction

      // D B C
      // A X

      // 1 A, B, D are set to 0 if unavailable
      // 2 If C is not available it is replaced by D

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

      rFrameL    = block_available_left    ? refFrArr[j4][i4-1]   : -1;
      rFrameU    = block_available_up      ? refFrArr[j4-1][i4]   : -1;
      rFrameUR   = block_available_upright ? refFrArr[j4-1][i4+BLOCK_STEP[img->mb_mode][0]] :
      block_available_upleft  ? refFrArr[j4-1][i4-1] : -1;

      // Prediction if only one of the neighbors uses the selected reference frame

      if(rFrameL == predframe_no && rFrameU != predframe_no && rFrameUR != predframe_no)
        mvPredType = MVPRED_L;
      else if(rFrameL != predframe_no && rFrameU == predframe_no && rFrameUR != predframe_no)
        mvPredType = MVPRED_U;
      else if(rFrameL != predframe_no && rFrameU != predframe_no && rFrameUR == predframe_no)
        mvPredType = MVPRED_UR;

      // Directional predictions

      else if(img->mb_mode == 3)
      {
        if(i == 0)
        {
          if(rFrameL == predframe_no)
            mvPredType = MVPRED_L;
        }
        else
        {
          if(rFrameUR == predframe_no)
            mvPredType = MVPRED_UR;
        }
      }
      else if(img->mb_mode == 2)
      {
        if(j == 0)
        {
          if(rFrameU == predframe_no)
            mvPredType = MVPRED_U;
        }
        else
        {
          if(rFrameL == predframe_no)
            mvPredType = MVPRED_L;
        }
      }
      else if(img->mb_mode == 5 && i == 2)
        mvPredType = MVPRED_L;
      else if(img->mb_mode == 6 && j == 2)
        mvPredType = MVPRED_U;

      for (k=0; k < 2; k++)
      {

        mv_a = block_available_left ? img->mv[i4-1+BLOCK_SIZE][j4][k] : 0;
        mv_b = block_available_up      ? img->mv[i4+BLOCK_SIZE][j4-1][k]   : 0;
        mv_d = block_available_upleft  ? img->mv[i4-1+BLOCK_SIZE][j4-1][k] : 0;
        mv_c = block_available_upright ? img->mv[i4+BLOCK_STEP[img->mb_mode][0]+BLOCK_SIZE][j4-1][k] : mv_d;

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
        snprintf(currSE.tracestring, TRACESTRING_SIZE, " MVD");
#endif
        img->subblock_x = i; // position used for context determination
        img->subblock_y = j; // position used for context determination
        currSE.value2 = k; // identifies the component; only used for context determination

        dP->readSyntaxElement(&currSE,img,inp,dP);
        curr_mvd = currSE.value1;

        vec=curr_mvd+pred_vec;           // find motion vector
        for(ii=0;ii<BLOCK_STEP[img->mb_mode][0];ii++)
          for(jj=0;jj<BLOCK_STEP[img->mb_mode][1];jj++)
            img->mv[i4+ii+BLOCK_SIZE][j4+jj][k]=vec;

          // store (oversampled) mvd
          for (l=0; l < step_v; l++)
            for (m=0; m < step_h; m++)
              currMB->mvd[0][j+l][i+m][k] =  curr_mvd;
      }
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Get coded block pattern and coefficients (run/level)
 *    from the NAL
 ************************************************************************
 */
void readCBPandCoeffsFromNAL(struct img_par *img,struct inp_par *inp)
{
  int i,j,k;
  int level, run;
  int mb_nr = img->current_mb_nr;
  int ii,jj;
  int i1,j1, m2,jg2;
  Macroblock *currMB = &img->mb_data[mb_nr];
  int cbp;
  SyntaxElement currSE;
  Slice *currSlice = img->currentSlice;
  DataPartition *dP;
  int *partMap = assignSE2partition[currSlice->dp_mode];
  int iii,jjj;
  int coef_ctr, len, i0,j0;

  int ll;
  int scan_loop_ctr;

  int block_x,block_y;
  int scan_mode, start_scan;

  // read CBP if not new intra mode
  if (img->imod != INTRA_MB_NEW)
  {
    if (img->imod == INTRA_MB_OLD)
      currSE.type = SE_CBP_INTRA;
    else
      currSE.type = SE_CBP_INTER;

    if(img->type == B_IMG_1 || img->type == B_IMG_MULT)
      dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
    else
      dP = &(currSlice->partArr[partMap[currSE.type]]);

    if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag)
    {
      if (img->imod == INTRA_MB_OLD)
        currSE.mapping = linfo_cbp_intra;
      else
        currSE.mapping = linfo_cbp_inter;
    }
    else
      currSE.reading = readCBPFromBuffer_CABAC;

#if TRACE
    snprintf(currSE.tracestring, TRACESTRING_SIZE, " CBP ");
#endif
   

    dP->readSyntaxElement(&currSE,img,inp,dP);
    currMB->cbp = cbp = currSE.value1;
    // Delta quant only if nonzero coeffs
    if (cbp !=0)
    {
      if (currMB->intraOrInter == INTER_MB)
        currSE.type = SE_DELTA_QUANT_INTER;
      else
        currSE.type = SE_DELTA_QUANT_INTRA;

      if(img->type == B_IMG_1 || img->type == B_IMG_MULT)
        dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
      else
        dP = &(currSlice->partArr[partMap[currSE.type]]);
  
      if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag)
        currSE.mapping = linfo_dquant;
      else
        currSE.reading= readDquantFromBuffer_CABAC;

#if TRACE
      snprintf(currSE.tracestring, TRACESTRING_SIZE, " Delta quant ");
#endif

      dP->readSyntaxElement(&currSE,img,inp,dP);
      currMB->delta_quant = currSE.value1;
      img->qp= (img->qp+currMB->delta_quant+32)%32;
    }
  }
  else
    cbp = currMB->cbp;

  for (i=0;i<BLOCK_SIZE;i++)
    for (j=0;j<BLOCK_SIZE;j++)
      for(iii=0;iii<BLOCK_SIZE;iii++)
        for(jjj=0;jjj<BLOCK_SIZE;jjj++)
          img->cof[i][j][iii][jjj]=0;// reset luma coeffs

  if(img->imod==INTRA_MB_NEW) // read DC coeffs for new intra modes
  {

    if (currMB->intraOrInter == INTER_MB)
      currSE.type = SE_DELTA_QUANT_INTER;
    else
      currSE.type = SE_DELTA_QUANT_INTRA;

    if(img->type == B_IMG_1 || img->type == B_IMG_MULT)
      dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
    else
      dP = &(currSlice->partArr[partMap[currSE.type]]);

    if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag)
      currSE.mapping = linfo_dquant;
    else
      currSE.reading= readDquantFromBuffer_CABAC;

#if TRACE
      snprintf(currSE.tracestring, TRACESTRING_SIZE, " Delta quant ");
#endif

    dP->readSyntaxElement(&currSE,img,inp,dP);
    currMB->delta_quant = currSE.value1;
    img->qp= (img->qp+currMB->delta_quant+32)%32;

    for (i=0;i<BLOCK_SIZE;i++)
      for (j=0;j<BLOCK_SIZE;j++)
        img->ipredmode[img->block_x+i+1][img->block_y+j+1]=0;

            // //////

    currSE.type = SE_LUM_DC_INTRA;
    if(img->type == B_IMG_1 || img->type == B_IMG_MULT)
      dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
    else
      dP = &(currSlice->partArr[partMap[currSE.type]]);

    if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag)
      currSE.mapping = linfo_levrun_inter;
    else
    {
      currSE.reading = readRunLevelFromBuffer_CABAC;
      currSE.context = 3; // for choosing context model
    }

    // ////////////
    coef_ctr=-1;
    level = 1;                            // just to get inside the loop
    for(k=0;(k<17) && (level!=0);k++)
    {
#if TRACE
      snprintf(currSE.tracestring, TRACESTRING_SIZE, "DC luma 16x16 ");
#endif
      dP->readSyntaxElement(&currSE,img,inp,dP);
      level = currSE.value1;
      run = currSE.value2;
      len = currSE.len;

      if (level != 0)                     // leave if len=1
      {
        coef_ctr=coef_ctr+run+1;

        i0=SNGL_SCAN[coef_ctr][0];
        j0=SNGL_SCAN[coef_ctr][1];

        img->cof[i0][j0][0][0]=level;// add new intra DC coeff
      }
    }
    itrans_2(img);// transform new intra DC
  }


  if (img->imod == INTRA_MB_OLD && img->qp < 24)
    scan_mode=DOUBLE_SCAN;
  else
    scan_mode=SINGLE_SCAN;


  // luma coefficients
  for (block_y=0; block_y < 4; block_y += 2) // all modes
  {
    for (block_x=0; block_x < 4; block_x += 2)
    {
      for (j=block_y; j < block_y+2; j++)
      {
        jj=j/2;
        for (i=block_x; i < block_x+2; i++)
        {
          ii=i/2;
          if (img->imod == INTRA_MB_NEW)
            start_scan = 1; // skip DC coeff
          else
            start_scan = 0; // take all coeffs

          if((cbp & (int)pow(2,(ii+2*jj))) != 0)  // are there any coeff in current block at all
          {
            if (scan_mode==SINGLE_SCAN)
            {
              coef_ctr=start_scan-1;
              level = 1;
              for(k=start_scan;(k<17) && (level!=0);k++)
              {
                  /*
                  * make distinction between INTRA and INTER coded
                  * luminance coefficients
                */
                if (k == 0)
                {
                  if (img->imod == INTRA_MB_OLD || img->imod == INTRA_MB_NEW)
                  {
                    currSE.context = 2; // for choosing context model
                    currSE.type  = SE_LUM_DC_INTRA;
                  }
                  else
                  {
                    currSE.context = 1; // for choosing context model
                    currSE.type  = SE_LUM_DC_INTER;
                  }
                }
                else
                {
                  if (img->imod == INTRA_MB_OLD /*|| img->imod == INTRA_MB_NEW*/)
                  {
                    currSE.context = 2; // for choosing context model
                    currSE.type  = SE_LUM_AC_INTRA;
                  }
                  else if ( img->imod == INTRA_MB_NEW )
                  {
                    currSE.context = 4; // for choosing context model
                    currSE.type  = SE_LUM_AC_INTRA;
                  }
                  else
                  {
                    currSE.context = 1; // for choosing context model
                    currSE.type  = SE_LUM_AC_INTER;
                  }
                }

#if TRACE
                snprintf(currSE.tracestring, TRACESTRING_SIZE, " Luma sng ");
#endif
                
                if(img->type == B_IMG_1 || img->type == B_IMG_MULT)
                  dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
                else
                  dP = &(currSlice->partArr[partMap[currSE.type]]);

                if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag)
                  currSE.mapping = linfo_levrun_inter;
                else
                  currSE.reading = readRunLevelFromBuffer_CABAC;

                dP->readSyntaxElement(&currSE,img,inp,dP);
                level = currSE.value1;
                run =  currSE.value2;
                len = currSE.len;

                if (level != 0)   // leave if len=1
                {
                  coef_ctr += run+1;
                  currMB->cbp_blk |= 1 << ((j<<2) + i) ;
                  i0=SNGL_SCAN[coef_ctr][0];
                  j0=SNGL_SCAN[coef_ctr][1];


                  img->cof[i][j][i0][j0]=level*JQ1[img->qp];
                }
              }
            }
            else    // double scan (old intra with QP<24
            {
              for(scan_loop_ctr=0;scan_loop_ctr<2;scan_loop_ctr++)
              {
                coef_ctr=start_scan-1;
                level=1;                          // just to get inside the loop
                for(k=0; k<9 && level!=0;k++)
                {
                  if (k == 0)
                    currSE.type  = SE_LUM_DC_INTRA; // element is of type DC
                  else
                    currSE.type  = SE_LUM_AC_INTRA;   // element is of type AC
#if TRACE
                  snprintf(currSE.tracestring, TRACESTRING_SIZE, "Luma dbl(%2d,%2d)  ",scan_loop_ctr,k);
#endif
                  if(img->type == B_IMG_1 || img->type == B_IMG_MULT)
                    dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
                  else
                    dP = &(currSlice->partArr[partMap[currSE.type]]);

                  if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag)
                    currSE.mapping = linfo_levrun_intra;
                  else
                  {
                    currSE.context = 0; // for choosing context model
                    currSE.reading = readRunLevelFromBuffer_CABAC;
                  }

                  dP->readSyntaxElement(&currSE,img,inp,dP);
                  level = currSE.value1;
                  run = currSE.value2;
                  len = currSE.len;

                  if (level != 0)   // leave if len=1
                  {
                    coef_ctr=coef_ctr+run+1;
                    currMB->cbp_blk |= 1 << ((j<<2) + i) ;

                    i0=DBL_SCAN[coef_ctr][0][scan_loop_ctr];
                    j0=DBL_SCAN[coef_ctr][1][scan_loop_ctr];

                    img->cof[i][j][i0][j0]=level*JQ1[img->qp];
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  for (j=4;j<6;j++) // reset all chroma coeffs before read
    for (i=0;i<4;i++)
      for (iii=0;iii<4;iii++)
        for (jjj=0;jjj<4;jjj++)
          img->cof[i][j][iii][jjj]=0;

  m2=img->mb_x*2;
  jg2=img->mb_y*2;

  // chroma 2x2 DC coeff
  if(cbp>15)
  {
    for (ll=0;ll<3;ll+=2)
    {
      for (i=0;i<4;i++)
        img->cofu[i]=0;

      coef_ctr=-1;
      level=1;
      for(k=0;(k<5)&&(level!=0);k++)
      {
        if ( img->imod == INTRA_MB_OLD || img->imod == INTRA_MB_NEW)
        {
          currSE.context = 6; // for choosing context model
          currSE.type  = SE_CHR_DC_INTRA;
        }
        else
        {
          currSE.context = 5; // for choosing context model
          currSE.type  = SE_CHR_DC_INTER;
        }

#if TRACE
        snprintf(currSE.tracestring, TRACESTRING_SIZE, " 2x2 DC Chroma ");
#endif
        if(img->type == B_IMG_1 || img->type == B_IMG_MULT)
          dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
        else
          dP = &(currSlice->partArr[partMap[currSE.type]]);

        if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag)
          currSE.mapping = linfo_levrun_c2x2;
        else
          currSE.reading = readRunLevelFromBuffer_CABAC;

        dP->readSyntaxElement(&currSE,img,inp,dP);
        level = currSE.value1;
        run = currSE.value2;
        len = currSE.len;
        if (level != 0)
        {
          currMB->cbp_blk |= 0xf0000 << (ll<<1) ;
          coef_ctr=coef_ctr+run+1;
          // Bug: img->cofu has only 4 entries, hence coef_ctr MUST be <4 (which is
          // caught by the assert().  If it is bigger than 4, it starts patching the
          // img->predmode pointer, which leads to bugs later on.
          //
          // This assert() should be left in the code, because it captures a very likely
          // bug early when testing in error prone environments (or when testing NAL
          // functionality).
          assert (coef_ctr < 4);
          img->cofu[coef_ctr]=level*JQ1[QP_SCALE_CR[img->qp]];
        }
      }

      if ((img->type==SP_IMG_1 || img->type==SP_IMG_MULT) && (currMB->mb_imode == INTRA_MB_INTER))
      {
        img->cof[0+ll][4][0][0]=img->cofu[0]/JQ1[QP_SCALE_CR[img->qp]];
        img->cof[1+ll][4][0][0]=img->cofu[1]/JQ1[QP_SCALE_CR[img->qp]];
        img->cof[0+ll][5][0][0]=img->cofu[2]/JQ1[QP_SCALE_CR[img->qp]];
        img->cof[1+ll][5][0][0]=img->cofu[3]/JQ1[QP_SCALE_CR[img->qp]];
      }
      else
      {
        img->cof[0+ll][4][0][0]=(img->cofu[0]+img->cofu[1]+img->cofu[2]+img->cofu[3])/2;
        img->cof[1+ll][4][0][0]=(img->cofu[0]-img->cofu[1]+img->cofu[2]-img->cofu[3])/2;
        img->cof[0+ll][5][0][0]=(img->cofu[0]+img->cofu[1]-img->cofu[2]-img->cofu[3])/2;
        img->cof[1+ll][5][0][0]=(img->cofu[0]-img->cofu[1]-img->cofu[2]+img->cofu[3])/2;
      }
    }
  }

  // chroma AC coeff, all zero fram start_scan
  if (cbp>31)
  {
    block_y=4;
    for (block_x=0; block_x < 4; block_x += 2)
    {
      for (j=block_y; j < block_y+2; j++)
      {
        jj=j/2;
        j1=j-4;
        for (i=block_x; i < block_x+2; i++)
        {
          ii=i/2;
          i1=i%2;
          {
            coef_ctr=0;
            level=1;
            for(k=0;(k<16)&&(level!=0);k++)
            {
              if ( img->imod == INTRA_MB_OLD || img->imod == INTRA_MB_NEW)
              {
                currSE.context = 8; // for choosing context model
                currSE.type  = SE_CHR_AC_INTRA;
              }
              else
              {
                currSE.context = 7; // for choosing context model
                currSE.type  = SE_CHR_AC_INTER;
              }
#if TRACE
              snprintf(currSE.tracestring, TRACESTRING_SIZE, " AC Chroma ");
#endif
              if(img->type == B_IMG_1 || img->type == B_IMG_MULT)
                dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
              else
                dP = &(currSlice->partArr[partMap[currSE.type]]);
              
              if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag)
                currSE.mapping = linfo_levrun_inter;
              else
                currSE.reading = readRunLevelFromBuffer_CABAC;

              dP->readSyntaxElement(&currSE,img,inp,dP);
              level = currSE.value1;
              run = currSE.value2;
              len = currSE.len;

              if (level != 0)
              {
                currMB->cbp_blk |= 1 << (16 + (j1<<1) + i1 + (block_x<<1) ) ;
                coef_ctr=coef_ctr+run+1;
                i0=SNGL_SCAN[coef_ctr][0];
                j0=SNGL_SCAN[coef_ctr][1];
                img->cof[i][j][i0][j0]=level*JQ1[QP_SCALE_CR[img->qp]];
              }
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
 *    copy current MB from last MB
 ************************************************************************
 */
void decode_one_CopyMB(struct img_par *img,struct inp_par *inp)
{
  int tmp_block[BLOCK_SIZE][BLOCK_SIZE];
  int i, j, ii, jj, uv;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int ref_frame = currMB->ref_frame;
  int mv_mul;

  if(img->mv_res)
    mv_mul=8;
  else
    mv_mul=4;

  // get luma pixel *************************************************
  for(j=0;j<MB_BLOCK_SIZE;j+=BLOCK_SIZE)
  {
    for(i=0;i<MB_BLOCK_SIZE;i+=BLOCK_SIZE)
    {
      get_block(ref_frame,(img->pix_x+i)*mv_mul,(img->pix_y+j)*mv_mul,img,tmp_block);

      for(ii=0;ii<BLOCK_SIZE;ii++)
        for(jj=0;jj<BLOCK_SIZE;jj++)
          imgY[img->pix_y+j+jj][img->pix_x+i+ii]=tmp_block[ii][jj];
    }
  }
  if (img->type==SP_IMG_1 || img->type==SP_IMG_MULT)
  {
    for(j=0;j<MB_BLOCK_SIZE;j++)
    {
      jj=img->pix_y+j;
      for(i=0;i<MB_BLOCK_SIZE;i++)
      {
        ii=img->pix_x+i;
        img->mpr[i][j]=imgY[jj][ii];
      }
    }
    for (i=0; i < MB_BLOCK_SIZE; i+=BLOCK_SIZE)
      for (j=0; j < MB_BLOCK_SIZE; j+=BLOCK_SIZE)
        copyblock_sp(img,i,j);
  }
  // get chroma pixel **********************************************
  for(uv=0;uv<2;uv++)
  {
    for(j=0;j<MB_BLOCK_SIZE/2;j++)
    {
      jj=img->pix_c_y+j;
      for(i=0;i<MB_BLOCK_SIZE/2;i++)
      {
        ii=img->pix_c_x+i;
        imgUV[uv][jj][ii]=mcef[ref_frame][uv][jj][ii];
      }
    }
    if(img->type==SP_IMG_1 || img->type==SP_IMG_MULT)
    {
      for(j=0;j<MB_BLOCK_SIZE/2;j++)
      {
        jj=img->pix_c_y+j;
        for(i=0;i<MB_BLOCK_SIZE/2;i++)
        {
          ii=img->pix_c_x+i;
          img->mpr[i][j]=imgUV[uv][jj][ii];
        }
      }
      for (j=4;j<6;j++)
        for(i=0;i<4;i++)
          for(ii=0;ii<4;ii++)
            for(jj=0;jj<4;jj++)
              img->cof[i][j][ii][jj]=0;

      itrans_sp_chroma(img,2*uv);

      for (j=4;j<6;j++)
      {
        for(i=0;i<2;i++)
        {
          itrans(img,i*4,(j-4)*4,2*uv+i,j);

          for(ii=0;ii<4;ii++)
            for(jj=0;jj<4;jj++)
            {
              imgUV[uv][img->pix_c_y+(j-4)*4+jj][img->pix_c_x+i*4+ii]=img->m7[ii][jj];
            }
        }
      }
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    decode one macroblock
 ************************************************************************
 */
int decode_one_macroblock(struct img_par *img,struct inp_par *inp)
{
  int tmp_block[BLOCK_SIZE][BLOCK_SIZE];
  int js[2][2];
  int i=0,j=0,ii=0,jj=0,i1=0,j1=0,j4=0,i4=0;
  int js0=0,js1=0,js2=0,js3=0,jf=0;
  int uv;
  int vec1_x=0,vec1_y=0,vec2_x=0,vec2_y=0;
  int ioff,joff;

  int ii0,jj0,ii1,jj1,if1,jf1,if0,jf0;
  int mv_mul,f1,f2,f3,f4;

  Macroblock *currMB   = &img->mb_data[img->current_mb_nr];
  int ref_frame        = currMB->ref_frame;

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

  // luma decoding **************************************************

  // get prediction for INTRA_MB_16x16
  if (currMB->mb_imode == INTRA_MB_NEW)
    intrapred_luma_2(img,currMB->intra_pred_modes[0]);

  for(j=0;j<MB_BLOCK_SIZE/BLOCK_SIZE;j++)
  {
    joff=j*4;
    j4=img->block_y+j;
    for(i=0;i<MB_BLOCK_SIZE/BLOCK_SIZE;i++)
    {
      ioff=i*4;
      i4=img->block_x+i;
      // get prediction for INTRA_MB_4x4
      if(currMB->mb_imode == INTRA_MB_OLD)
      {
        if (intrapred(img,ioff,joff,i4,j4)==SEARCH_SYNC)  // make 4x4 prediction block mpr from given prediction img->mb_mode
          return SEARCH_SYNC;                   // bit error
      }
      // get motion prediction for INTER_MB
      else if(currMB->mb_imode == INTRA_MB_INTER)
      {
        vec1_x = i4*4*mv_mul + img->mv[i4+BLOCK_SIZE][j4][0];
        vec1_y = j4*4*mv_mul + img->mv[i4+4][j4][1];

        get_block(ref_frame,vec1_x,vec1_y,img,tmp_block);

        for(ii=0;ii<BLOCK_SIZE;ii++)
          for(jj=0;jj<BLOCK_SIZE;jj++)
            img->mpr[ii+ioff][jj+joff] = tmp_block[ii][jj];
      }

      if ((img->type==SP_IMG_1 || img->type==SP_IMG_MULT) && (currMB->mb_imode == INTRA_MB_INTER))
        itrans_sp(img,ioff,joff,i,j);
      else itrans(img,ioff,joff,i,j);      // use DCT transform and make 4x4 block m7 from prediction block mpr

      for(ii=0;ii<BLOCK_SIZE;ii++)
      {
        for(jj=0;jj<BLOCK_SIZE;jj++)
        {
          imgY[j4*BLOCK_SIZE+jj][i4*BLOCK_SIZE+ii]=img->m7[ii][jj]; // contruct picture from 4x4 blocks
        }
      }
    }
  }
  // chroma decoding *******************************************************

  for(uv=0;uv<2;uv++)
  {
    if (img->imod==INTRA_MB_OLD || img->imod==INTRA_MB_NEW)// intra mode
    {
      int mb_available_up = (currMB->mb_available[0][1] != NULL);
      int mb_available_left = (currMB->mb_available[1][0] != NULL);
      if(img->UseConstrainedIntraPred)
      {
        const int mb_nr = img->current_mb_nr;
        const int mb_width = img->width/MB_BLOCK_SIZE;
        if (mb_available_up && (img->intra_mb[mb_nr-mb_width] ==0))
          mb_available_up = 0;
        if (mb_available_left && (img->intra_mb[mb_nr-1] ==0))
          mb_available_left = 0;
      }

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
        else
        {
          for(jj=0;jj<4;jj++)
          {
            jf=(j4+jj)/2;
            for(ii=0;ii<4;ii++)
            {
              if1=(i4+ii)/2;
              i1=(img->pix_c_x+ii+ioff)*f1+img->mv[if1+4][jf][0];
              j1=(img->pix_c_y+jj+joff)*f1+img->mv[if1+4][jf][1];

              ii0=max (0, min (i1/f1, img->width_cr-1));
              jj0=max (0, min (j1/f1, img->height_cr-1));
              ii1=max (0, min ((i1+f2)/f1, img->width_cr-1));
              jj1=max (0, min ((j1+f2)/f1, img->height_cr-1));

              if1=(i1 & f2);
              jf1=(j1 & f2);
              if0=f1-if1;
              jf0=f1-jf1;
              img->mpr[ii+ioff][jj+joff]=(if0*jf0*mcef[ref_frame][uv][jj0][ii0]+
                if1*jf0*mcef[ref_frame][uv][jj0][ii1]+
                if0*jf1*mcef[ref_frame][uv][jj1][ii0]+
                if1*jf1*mcef[ref_frame][uv][jj1][ii1]+f4)/f3;

            }
          }
        }
        if ((img->type!=SP_IMG_1 && img->type!=SP_IMG_MULT) || (currMB->mb_imode != INTRA_MB_INTER))
        {
          itrans(img,ioff,joff,2*uv+i,j);
          for(ii=0;ii<4;ii++)
            for(jj=0;jj<4;jj++)
            {
              imgUV[uv][j4+jj][i4+ii]=img->m7[ii][jj];
            }
        }
      }
    }
    if((img->type==SP_IMG_1 || img->type==SP_IMG_MULT) && (currMB->mb_imode == INTRA_MB_INTER))
    {
      itrans_sp_chroma(img,2*uv);
      for (j=4;j<6;j++)
      {
        joff=(j-4)*4;
        j4=img->pix_c_y+joff;
        for(i=0;i<2;i++)
        {
          ioff=i*4;
          i4=img->pix_c_x+ioff;
          itrans(img,ioff,joff,2*uv+i,j);

          for(ii=0;ii<4;ii++)
            for(jj=0;jj<4;jj++)
            {
              imgUV[uv][j4+jj][i4+ii]=img->m7[ii][jj];
            }
        }
      }
    }
  }

  return 0;
}
