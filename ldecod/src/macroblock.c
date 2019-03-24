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
#include "errorconcealment.h"
#include "macroblock.h"
#include "decodeiff.h"

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
    // upper blocks
    if (remove_prediction || (img->UseConstrainedIntraPred && img->intra_block[mb_nr-1][1]==0))
    {
      img->ipredmode[img->block_x][img->block_y+1] = -1;
      img->ipredmode[img->block_x][img->block_y+2] = -1;
    }
    // lower blocks
    if (remove_prediction || (img->UseConstrainedIntraPred && img->intra_block[mb_nr-1][3]==0))
    {
      img->ipredmode[img->block_x][img->block_y+3] = -1;
      img->ipredmode[img->block_x][img->block_y+4] = -1;
    }
    if (!remove_prediction)
    {
      currMB->mb_available[1][0]=&(img->mb_data[mb_nr-1]);
    }
  }


  // Check MB above
  if(img->pix_y >= MB_BLOCK_SIZE)
  {
    int remove_prediction = currMB->slice_nr != img->mb_data[mb_nr-mb_width].slice_nr;
    // upper blocks
    if (remove_prediction || (img->UseConstrainedIntraPred && img->intra_block[mb_nr-mb_width][2]==0))
    {
      img->ipredmode[img->block_x+1][img->block_y] = -1;
      img->ipredmode[img->block_x+2][img->block_y] = -1;
    }
    // lower blocks
    if (remove_prediction || (img->UseConstrainedIntraPred && img->intra_block[mb_nr-mb_width][3]==0))
    {
      img->ipredmode[img->block_x+3][img->block_y] = -1;
      img->ipredmode[img->block_x+4][img->block_y] = -1;
    }
    if (!remove_prediction)
    {
      currMB->mb_available[0][1]=&(img->mb_data[mb_nr-mb_width]);
    }
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
  currMB->qp          = img->qp ;
  currMB->mb_type     = 0;
  currMB->delta_quant = 0;
  currMB->cbp         = 0;
  currMB->cbp_blk     = 0;

  for (l=0; l < 2; l++)
    for (j=0; j < BLOCK_MULTIPLE; j++)
      for (i=0; i < BLOCK_MULTIPLE; i++)
        for (k=0; k < 2; k++)
          currMB->mvd[l][j][i][k] = 0;

  for (i=0; i < (BLOCK_MULTIPLE*BLOCK_MULTIPLE); i++)
    currMB->intra_pred_modes[i] = 0;

  for (j=0; j < BLOCK_MULTIPLE; j++)
    for (i=0; i < BLOCK_MULTIPLE; i++)
      currMB->coeffs_count[j][i] = 0;
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
    if(nal_startcode_follows(img, inp) == FALSE) return FALSE;
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
  int i;
  const int ICBPTAB[6] = {0,16,32,15,31,47};
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];
  int         mbmode = currMB->mb_type;

#define ZERO_P8x8     (mbmode==5)
#define MODE_IS_P8x8  (mbmode==4 || mbmode==5)
#define MODE_IS_I4x4  (mbmode==6)
#define I16OFFSET     (mbmode-7)

  if(mbmode <4)
  {
    currMB->mb_type = mbmode;
    for (i=0;i<4;i++) {currMB->b8mode[i]=mbmode; currMB->b8pdir[i]=0; }
  }
  else if(MODE_IS_P8x8)
  {
    currMB->mb_type = P8x8;
    img->allrefzero = ZERO_P8x8;
    // b8mode and pdir are read and set later
  }
  else if(MODE_IS_I4x4)
  {
    currMB->mb_type = I4MB;
    for (i=0;i<4;i++) {currMB->b8mode[i]=IBLOCK; currMB->b8pdir[i]=-1; }
  }
  else
  {
    currMB->mb_type = I16MB;
    for (i=0;i<4;i++) {currMB->b8mode[i]=0; currMB->b8pdir[i]=-1; }
    currMB->cbp= ICBPTAB[(I16OFFSET)>>2];
    currMB->i16mode = (I16OFFSET) & 0x03;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Interpret the mb mode for I-Frames
 ************************************************************************
 */
void interpret_mb_mode_I(struct img_par *img)
{
  int i;
  const int ICBPTAB[6] = {0,16,32,15,31,47};
  Macroblock *currMB   = &img->mb_data[img->current_mb_nr];
  int         mbmode   = currMB->mb_type;

  if (mbmode==0)
  {
    currMB->mb_type = I4MB;
    for (i=0;i<4;i++) {currMB->b8mode[i]=IBLOCK; currMB->b8pdir[i]=-1; }
  }
  else
  {
    currMB->mb_type = I16MB;
    for (i=0;i<4;i++) {currMB->b8mode[i]=0; currMB->b8pdir[i]=-1; }
    currMB->cbp= ICBPTAB[(mbmode-1)>>2];
    currMB->i16mode = (mbmode-1) & 0x03;
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
  static const int offset2pdir16x16[12]   = {0, 0, 1, 2, 0,0,0,0,0,0,0,0};
  static const int offset2pdir16x8[22][2] = {{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{1,1},{0,0},{0,1},{0,0},{1,0},
                                             {0,0},{0,2},{0,0},{1,2},{0,0},{2,0},{0,0},{2,1},{0,0},{2,2},{0,0}};
  static const int offset2pdir8x16[22][2] = {{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{1,1},{0,0},{0,1},{0,0},
                                             {1,0},{0,0},{0,2},{0,0},{1,2},{0,0},{2,0},{0,0},{2,1},{0,0},{2,2}};

  const int ICBPTAB[6] = {0,16,32,15,31,47};
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  int i, mbmode;
  int mbtype  = currMB->mb_type;
  int *b8mode = currMB->b8mode;
  int *b8pdir = currMB->b8pdir;

  //--- set mbtype, b8type, and b8pdir ---
  if (mbtype==0)       // direct
  {
    mbmode=0;       for(i=0;i<4;i++) {b8mode[i]=0;          b8pdir[i]=2;}
  }
  else if (mbtype==23) // intra4x4
  {
    mbmode=I4MB;    for(i=0;i<4;i++) {b8mode[i]=IBLOCK;     b8pdir[i]=-1;}
  }
  else if (mbtype>23) // intra16x16
  {
    mbmode=I16MB;   for(i=0;i<4;i++) {b8mode[i]=0;          b8pdir[i]=-1;}
    currMB->cbp     = ICBPTAB[(mbtype-24)>>2];
    currMB->i16mode = (mbtype-24) & 0x03;
  }
  else if (mbtype==22) // 8x8(+split)
  {
    mbmode=P8x8;       // b8mode and pdir is transmitted in additional codewords
  }
  else if (mbtype<4)   // 16x16
  {
    mbmode=1;       for(i=0;i<4;i++) {b8mode[i]=1;          b8pdir[i]=offset2pdir16x16[mbtype];}
  }
  else if (mbtype%2==0) // 16x8
  {
    mbmode=2;       for(i=0;i<4;i++) {b8mode[i]=2;          b8pdir[i]=offset2pdir16x8 [mbtype][i/2];}
  }
  else
  {
    mbmode=3;       for(i=0;i<4;i++) {b8mode[i]=3;          b8pdir[i]=offset2pdir8x16 [mbtype][i%2];}
  }
  currMB->mb_type = mbmode;
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

  predframe_no = 0;

  // Set the reference frame information for motion vector prediction
  if (IS_INTRA (currMB))
  {
    for (j=0; j<4; j++)
    for (i=0; i<4; i++)
    {
      refFrArr[img->block_y+j][img->block_x+i] = -1;
    }
  }
  else if (!IS_P8x8 (currMB))
  {
    for (j=0; j<4; j++)
    for (i=0; i<4; i++)
    {
      refFrArr[img->block_y+j][img->block_x+i] = 0;
    }
  }
  else
  {
    for (j=0; j<4; j++)
    for (i=0; i<4; i++)
    {
      refFrArr[img->block_y+j][img->block_x+i] = (currMB->b8mode[2*(j/2)+(i/2)]==IBLOCK ? -1 : 0);
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    Sets mode for 8x8 block
 ************************************************************************
 */
void
SetB8Mode (struct img_par* img, Macroblock* currMB, int value, int i)
{
  static const int p_v2b8 [ 5] = {4, 5, 6, 7, IBLOCK};
  static const int p_v2pd [ 5] = {0, 0, 0, 0, -1};
  static const int b_v2b8 [14] = {0, 4, 4, 4, 5, 6, 5, 6, 5, 6, 7, 7, 7, IBLOCK};
  static const int b_v2pd [14] = {2, 0, 1, 2, 0, 0, 1, 1, 2, 2, 0, 1, 2, -1};

  if (img->type==B_IMG_1 || img->type==B_IMG_MULT)
  {
    currMB->b8mode[i] = b_v2b8[value];
    currMB->b8pdir[i] = b_v2pd[value];
  }
  else
  {
    currMB->b8mode[i] = p_v2b8[value];
    currMB->b8pdir[i] = p_v2pd[value];
  }
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

  currMB->qp = img->qp ;


  //  read MB mode *****************************************************************
  currSE.type = SE_MBTYPE;

  if(img->type == B_IMG_1 || img->type == B_IMG_MULT) dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
  else                                                dP = &(currSlice->partArr[partMap[currSE.type]]);

  if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag)   currSE.mapping = linfo;
  else                                                      currSE.reading = readMB_typeInfoFromBuffer_CABAC;
  
  if(inp->symbol_mode == CABAC || (img->type != INTER_IMG_1 && img->type != INTER_IMG_MULT && img->type != B_IMG_1 && img->type != B_IMG_MULT))
  {
    //  read MB mode
#if TRACE
    strncpy(currSE.tracestring, "MB Type", TRACESTRING_SIZE);
#endif
    dP->readSyntaxElement(&currSE,img,inp,dP);
    currMB->mb_type = currSE.value1;
    if(!dP->bitstream->ei_flag)
      currMB->ei_flag = 0;
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
      currMB->mb_type = currSE.value1;
      if(!dP->bitstream->ei_flag)
        currMB->ei_flag = 0;
      img->cod_counter--;
    } 
    else
    {
      img->cod_counter--;
      currMB->mb_type = 0;
      currMB->ei_flag = 0;
    }
  }


  if ((img->type==INTER_IMG_1) || (img->type==INTER_IMG_MULT))    // inter frame
    interpret_mb_mode_P(img);
  else if (img->type==INTRA_IMG)                                  // intra frame
    interpret_mb_mode_I(img);
  else if ((img->type==B_IMG_1) || (img->type==B_IMG_MULT))       // B frame
    interpret_mb_mode_B(img);
  else if ((img->type==SP_IMG_1) || (img->type==SP_IMG_MULT))     // SP frame
    interpret_mb_mode_P(img);


  //====== READ 8x8 SUB-PARTITION MODES (modes of 8x8 blocks) ======
  if (IS_P8x8 (currMB))
  {
    currSE.type    = SE_MBTYPE;
    if (img->type==B_IMG_1 || img->type==B_IMG_MULT)      dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
    else                                                  dP = &(currSlice->partArr[partMap[SE_MBTYPE]]);
    if (inp->symbol_mode==UVLC || dP->bitstream->ei_flag) currSE.mapping = linfo;
    else                                                  currSE.reading = readB8_typeInfoFromBuffer_CABAC;

    for (i=0; i<4; i++)
    {
      dP->readSyntaxElement (&currSE, img, inp, dP);
      SetB8Mode (img, currMB, currSE.value1, i);
    }
  }


  if(img->UseConstrainedIntraPred && (img->type==INTER_IMG_1 || img->type==INTER_IMG_MULT))        // inter frame
  {
    if (!IS_NEWINTRA (currMB) && currMB->b8mode[0]!=IBLOCK) img->intra_block[img->current_mb_nr][0] = 0;
    if (!IS_NEWINTRA (currMB) && currMB->b8mode[1]!=IBLOCK) img->intra_block[img->current_mb_nr][1] = 0;
    if (!IS_NEWINTRA (currMB) && currMB->b8mode[2]!=IBLOCK) img->intra_block[img->current_mb_nr][2] = 0;
    if (!IS_NEWINTRA (currMB) && currMB->b8mode[3]!=IBLOCK) img->intra_block[img->current_mb_nr][3] = 0;
  }
  
  
  //! TO for Error Concelament
  //! If we have an INTRA Macroblock and we lost the partition
  //! which contains the intra coefficients Copy MB would be better 
  //! than just a grey block.
  //! Seems to be a bit at the wrong place to do this right here, but for this case 
  //! up to now there is no other way.
  dP = &(currSlice->partArr[partMap[SE_CBP_INTRA]]);
  if(IS_INTRA (currMB) && dP->bitstream->ei_flag && img->number)
  {
    currMB->mb_type = 0;
    currMB->ei_flag = 1;
    for (i=0;i<4;i++) {currMB->b8mode[i]=currMB->b8pdir[i]=0; }
  }
  if(img->type == B_IMG_1 || img->type == B_IMG_MULT)  dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
  else                                                 dP = &(currSlice->partArr[partMap[currSE.type]]);
  //! End TO


  //--- init macroblock data ---
  if ((img->type==B_IMG_1) || (img->type==B_IMG_MULT))  init_macroblock_Bframe(img);
  else                                                  init_macroblock       (img);



  if (inp->symbol_mode != CABAC && IS_DIRECT (currMB) && img->cod_counter >= 0)
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
    return DECODE_MB;
  }

  if (IS_COPY (currMB)) //keep last macroblock
  {
    return DECODE_COPY_MB;
  }


  // intra prediction modes for a macroblock 4x4 **********************************************
  currSE.type = SE_INTRAPREDMODE;
  if(img->type == B_IMG_1 || img->type == B_IMG_MULT)     dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
  else                                                    dP = &(currSlice->partArr[partMap[currSE.type]]);
  if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag) currSE.mapping = linfo;
  else                                                    currSE.reading = readIntraPredModeFromBuffer_CABAC;
  for(i=0; i<8; i++)
  {
    if (currMB->b8mode[i/2]==IBLOCK)
    {
#if TRACE
      sprintf(currSE.tracestring, "Intra mode ");
#endif
      currSE.context=2*i;
      dP->readSyntaxElement(&currSE,img,inp,dP);

      i1 = img->block_x + 2*((i%4)/2);
      j1 = img->block_y + 2*(i/4) + (i%2);

      if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag)
      {
        dbl_ipred_word = currSE.value1;
        /* find intra prediction mode for two blocks */
        img->ipredmode[i1+1][j1+1] = PRED_IPRED[img->ipredmode[i1+1][j1]+1][img->ipredmode[i1][j1+1]+1][IPRED_ORDER[dbl_ipred_word][0]];
        img->ipredmode[i1+2][j1+1] = PRED_IPRED[img->ipredmode[i1+2][j1]+1][img->ipredmode[i1+1][j1+1]+1][IPRED_ORDER[dbl_ipred_word][1]];
      }
      else
      {
        currMB->intra_pred_modes[2*i  ] = currSE.value1;
        currMB->intra_pred_modes[2*i+1] = currSE.value2;
        img->ipredmode[i1+1][j1+1] = PRED_IPRED[img->ipredmode[i1+1][j1]+1][img->ipredmode[i1][j1+1]+1][currSE.value1];
        img->ipredmode[i1+2][j1+1] = PRED_IPRED[img->ipredmode[i1+2][j1]+1][img->ipredmode[i1+1][j1+1]+1][currSE.value2];
      }
    }
  }


  /* read inter frame vector data *********************************************************/
  if (IS_INTERMV (currMB))
  {
    readMotionInfoFromNAL (img, inp);
  }


  // read CBP and Coeffs  ***************************************************************
  readCBPandCoeffsFromNAL (img,inp);

  return DECODE_MB;
}





/*!
 ************************************************************************
 * \brief
 *    Set motion vector predictor
 ************************************************************************
 */
void
SetMotionVectorPredictor (struct img_par  *img,
                          int             *pmv_x,
                          int             *pmv_y,
                          int             ref_frame,
                          int             **refFrArr,
                          int             ***tmp_mv,
                          int             block_x,
                          int             block_y,
                          int             blockshape_x,
                          int             blockshape_y)
{
  Macroblock* currMB = &img->mb_data[img->current_mb_nr];

  int mb_x                 = 4*block_x;
  int mb_y                 = 4*block_y;
  int pic_block_x          = img->block_x + block_x;
  int pic_block_y          = img->block_y + block_y;
  int mb_nr                = img->current_mb_nr;
  int mb_width             = img->width/16;
  int mb_available_up =      (img->mb_y == 0        ) ? 0 : (currMB->slice_nr == img->mb_data[mb_nr-mb_width].slice_nr);
  int mb_available_left =    (img->mb_x == 0        ) ? 0 : (currMB->slice_nr == img->mb_data[mb_nr-1].slice_nr);
  int mb_available_upleft  = (img->mb_x == 0          ||
                              img->mb_y == 0        ) ? 0 : (currMB->slice_nr == img->mb_data[mb_nr-mb_width-1].slice_nr);
  int mb_available_upright = (img->mb_x >= mb_width-1 ||
                              img->mb_y == 0        ) ? 0 : (currMB->slice_nr == img->mb_data[mb_nr-mb_width+1].slice_nr);
  int block_available_up, block_available_left, block_available_upright, block_available_upleft;
  int mv_a, mv_b, mv_c, mv_d, pred_vec=0;
  int mvPredType, rFrameL, rFrameU, rFrameUR;
  int hv;

  /* D B C */
  /* A X   */

  /* 1 A, B, D are set to 0 if unavailable       */
  /* 2 If C is not available it is replaced by D */
  block_available_up   = mb_available_up   || (mb_y > 0);
  block_available_left = mb_available_left || (mb_x > 0);

  if (mb_y > 0)
  {
    if (mb_x < 8)  // first column of 8x8 blocks
    {
      if (mb_y==8)
      {
        if (blockshape_x == 16)      block_available_upright = 0;
        else                         block_available_upright = 1;
      }
      else
      {
        if (mb_x+blockshape_x != 8)  block_available_upright = 1;
        else                         block_available_upright = 0;
      }
    }
    else
    {
      if (mb_x+blockshape_x != 16)   block_available_upright = 1;
      else                           block_available_upright = 0;
    }
  }
  else if (mb_x+blockshape_x != MB_BLOCK_SIZE)
  {
    block_available_upright = mb_available_up;
  }
  else
  {
    block_available_upright = mb_available_upright;
  }

  if (mb_x > 0)
  {
    block_available_upleft = (mb_y > 0 ? 1 : mb_available_up);
  }
  else if (mb_y > 0)
  {
    block_available_upleft = mb_available_left;
  }
  else
  {
    block_available_upleft = mb_available_upleft;
  }

  mvPredType = MVPRED_MEDIAN;
  rFrameL    = block_available_left    ? refFrArr[pic_block_y]  [pic_block_x-1             ] : -1;
  rFrameU    = block_available_up      ? refFrArr[pic_block_y-1][pic_block_x               ]   : -1;
  rFrameUR   = block_available_upright ? refFrArr[pic_block_y-1][pic_block_x+blockshape_x/4] :
               block_available_upleft  ? refFrArr[pic_block_y-1][pic_block_x-1             ] : -1;


  /* Prediction if only one of the neighbors uses the reference frame
   * we are checking
   */
  if(rFrameL == ref_frame && rFrameU != ref_frame && rFrameUR != ref_frame)       mvPredType = MVPRED_L;
  else if(rFrameL != ref_frame && rFrameU == ref_frame && rFrameUR != ref_frame)  mvPredType = MVPRED_U;
  else if(rFrameL != ref_frame && rFrameU != ref_frame && rFrameUR == ref_frame)  mvPredType = MVPRED_UR;
  // Directional predictions 
  else if(blockshape_x == 8 && blockshape_y == 16)
  {
    if(mb_x == 0)
    {
      if(rFrameL == ref_frame)
        mvPredType = MVPRED_L;
    }
    else
    {
      if(rFrameUR == ref_frame)
        mvPredType = MVPRED_UR;
    }
  }
  else if(blockshape_x == 16 && blockshape_y == 8)
  {
    if(mb_y == 0)
    {
      if(rFrameU == ref_frame)
        mvPredType = MVPRED_U;
    }
    else
    {
      if(rFrameL == ref_frame)
        mvPredType = MVPRED_L;
    }
  }

  for (hv=0; hv < 2; hv++)
  {
    mv_a = block_available_left    ? tmp_mv[4+pic_block_x-1             ][pic_block_y  ][hv] : 0;
    mv_b = block_available_up      ? tmp_mv[4+pic_block_x               ][pic_block_y-1][hv] : 0;
    mv_d = block_available_upleft  ? tmp_mv[4+pic_block_x-1             ][pic_block_y-1][hv] : 0;
    mv_c = block_available_upright ? tmp_mv[4+pic_block_x+blockshape_x/4][pic_block_y-1][hv] : mv_d;

    switch (mvPredType)
    {
    case MVPRED_MEDIAN:
      if(!(block_available_upleft || block_available_up || block_available_upright))
        pred_vec = mv_a;
      else
        pred_vec = mv_a+mv_b+mv_c-min(mv_a,min(mv_b,mv_c))-max(mv_a,max(mv_b,mv_c));
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

    if (hv==0)  *pmv_x = pred_vec;
    else        *pmv_y = pred_vec;
  }
}


/*!
 ************************************************************************
 * \brief
 *    Set context for reference frames
 ************************************************************************
 */
int
BType2CtxRef (int btype)
{
  if (btype<4)  return 0;
  else          return 1;
}

/*!
 ************************************************************************
 * \brief
 *    Read motion info
 ************************************************************************
 */
void readMotionInfoFromNAL (struct img_par *img, struct inp_par *inp)
{
  int i,j,k,l,m;
  int step_h,step_v;
  int curr_mvd;
  int mb_nr           = img->current_mb_nr;
  Macroblock *currMB  = &img->mb_data[mb_nr];
  SyntaxElement currSE;
  Slice *currSlice    = img->currentSlice;
  DataPartition *dP;
  int *partMap        = assignSE2partition[inp->partition_mode];
  int bframe          = (img->type==B_IMG_1 || img->type==B_IMG_MULT);
  int partmode        = (IS_P8x8(currMB)?4:currMB->mb_type);
  int step_h0         = BLOCK_STEP [partmode][0];
  int step_v0         = BLOCK_STEP [partmode][1];


  int mv_mode, i0, j0, refframe;
  int pmv[2];
  int j4, i4, ii,jj;
  int vec;

  //  If multiple ref. frames, read reference frame for the MB *********************************
  if(img->type==INTER_IMG_MULT || img->type == SP_IMG_MULT || img->type == B_IMG_MULT)
  {
    currSE.type = SE_REFFRAME;
    if (bframe)                                               dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
    else                                                      dP = &(currSlice->partArr[partMap[SE_REFFRAME]]);
    if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag)   currSE.mapping = linfo;
    else                                                      currSE.reading = readRefFrameFromBuffer_CABAC;

    for (j0=0; j0<4; j0+=step_v0)
    for (i0=0; i0<4; i0+=step_h0)
    {
      k=2*(j0/2)+(i0/2);
      if ((currMB->b8pdir[k]==0 || currMB->b8pdir[k]==2) && currMB->b8mode[k]!=0)
      {
#if TRACE
        strncpy(currSE.tracestring,  "Reference frame no ", TRACESTRING_SIZE);
#endif
        img->subblock_x = i0;
        img->subblock_y = j0;

        if (!IS_P8x8 (currMB) || bframe || (!bframe && !img->allrefzero))
        {
          currSE.context = BType2CtxRef (currMB->b8mode[k]);
          dP->readSyntaxElement (&currSE,img,inp,dP);
          refframe = currSE.value1;
        }
        else
        {
          refframe = 0;
        }

        if (bframe && refframe>img->buf_cycle) //??? copied from readMotionInfoFrameNAL
        {
          set_ec_flag(SE_REFFRAME);
          refframe = 1;
        }

        if (!bframe)
        {
          for (j=j0; j<j0+step_v0;j++)
          for (i=i0; i<i0+step_h0;i++)
            refFrArr[img->block_y+j][img->block_x+i] = refframe;
        }
        else
        {
          for (j=j0; j<j0+step_v0;j++)
          for (i=i0; i<i0+step_h0;i++)
            img->fw_refFrArr[img->block_y+j][img->block_x+i] = refframe;
        }
      }
    }
  }


  //=====  READ FORWARD MOTION VECTORS =====
  currSE.type = SE_MVD;
  if (bframe)   dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
  else          dP = &(currSlice->partArr[partMap[SE_MVD]]);

  if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag) currSE.mapping = linfo_mvd;
  else if (bframe)                                        currSE.reading = readBiMVD2Buffer_CABAC;
  else                                                    currSE.reading = readMVDFromBuffer_CABAC;

  for (j0=0; j0<4; j0+=step_v0)
  for (i0=0; i0<4; i0+=step_h0)
  {
    k=2*(j0/2)+(i0/2);
    if ((currMB->b8pdir[k]==0 || currMB->b8pdir[k]==2) && (currMB->b8mode[k] !=0))//has forward vector
    {
      mv_mode  = currMB->b8mode[k];
      step_h   = BLOCK_STEP [mv_mode][0];
      step_v   = BLOCK_STEP [mv_mode][1];

      if (!bframe)  refframe = refFrArr        [img->block_y+j0][img->block_x+i0];
      else          refframe = img->fw_refFrArr[img->block_y+j0][img->block_x+i0];

      for (j=j0; j<j0+step_v0; j+=step_v)
      for (i=i0; i<i0+step_h0; i+=step_h)
      {
        j4 = img->block_y+j;
        i4 = img->block_x+i;

        // first make mv-prediction
        if (!bframe)  SetMotionVectorPredictor (img, pmv, pmv+1, refframe, refFrArr,         img->mv,    i, j, 4*step_h, 4*step_v);
        else          SetMotionVectorPredictor (img, pmv, pmv+1, refframe, img->fw_refFrArr, img->fw_mv, i, j, 4*step_h, 4*step_v);

        for (k=0; k < 2; k++) 
        {
#if TRACE
          snprintf(currSE.tracestring, TRACESTRING_SIZE, " MVD");
#endif
          img->subblock_x = i; // position used for context determination
          img->subblock_y = j; // position used for context determination
          currSE.value2 = (!bframe ? k : 2*k); // identifies the component; only used for context determination
          dP->readSyntaxElement(&currSE,img,inp,dP);
          curr_mvd = currSE.value1; 
  
          vec=curr_mvd+pmv[k];           /* find motion vector */

          if (bframe)
          {
            for(ii=0;ii<step_h;ii++)
              for(jj=0;jj<step_v;jj++)
                img->fw_mv[i4+ii+BLOCK_SIZE][j4+jj][k]=vec;
          }
          else
          {
            for(ii=0;ii<step_h;ii++)
              for(jj=0;jj<step_v;jj++)
                img->mv[i4+ii+BLOCK_SIZE][j4+jj][k]=vec;
          }

          /* store (oversampled) mvd */
          for (l=0; l < step_v; l++) 
            for (m=0; m < step_h; m++)  
              currMB->mvd[0][j+l][i+m][k] =  curr_mvd;
        }
      }
    }
  }


  //=====  READ FORWARD MOTION VECTORS =====
  currSE.type = SE_MVD;
  dP          = &(currSlice->partArr[partMap[SE_BFRAME]]);

  if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag) currSE.mapping = linfo_mvd;
  else                                                    currSE.reading = readBiMVD2Buffer_CABAC;

  for (j0=0; j0<4; j0+=step_v0)
  for (i0=0; i0<4; i0+=step_h0)
  {
    k=2*(j0/2)+(i0/2);
    if ((currMB->b8pdir[k]==1 || currMB->b8pdir[k]==2) && (currMB->b8mode[k]!=0))//has forward vector
    {
      mv_mode  = currMB->b8mode[k];
      step_h   = BLOCK_STEP [mv_mode][0];
      step_v   = BLOCK_STEP [mv_mode][1];

      refframe = img->bw_refFrArr[img->block_y+j0][img->block_x+i0]; // always 0

      for (j=j0; j<j0+step_v0; j+=step_v)
      for (i=i0; i<i0+step_h0; i+=step_h)
      {
        j4 = img->block_y+j;
        i4 = img->block_x+i;

        // first make mv-prediction
        SetMotionVectorPredictor (img, pmv, pmv+1, refframe, img->bw_refFrArr, img->bw_mv, i, j, 4*step_h, 4*step_v);

        for (k=0; k < 2; k++) 
        {
#if TRACE
          snprintf(currSE.tracestring, TRACESTRING_SIZE, " MVD");
#endif
          img->subblock_x = i; // position used for context determination
          img->subblock_y = j; // position used for context determination
          currSE.value2   = 2*k+1; // identifies the component; only used for context determination
          dP->readSyntaxElement(&currSE,img,inp,dP);
          curr_mvd = currSE.value1; 
  
          vec=curr_mvd+pmv[k];           /* find motion vector */

          for(ii=0;ii<step_h;ii++)
            for(jj=0;jj<step_v;jj++)
              img->bw_mv[i4+ii+BLOCK_SIZE][j4+jj][k]=vec;

          /* store (oversampled) mvd */
          for (l=0; l < step_v; l++) 
            for (m=0; m < step_h; m++)  
              currMB->mvd[1][j+l][i+m][k] =  curr_mvd;
        }
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
  int coef_ctr, len, i0, j0, b8;
  int ll;
  int scan_loop_ctr;
  int block_x,block_y;
  int start_scan;
  int uv;

  int qp_per;
  int qp_rem;
  int qp_per_uv;
  int qp_rem_uv;


  // read CBP if not new intra mode
  if (!IS_NEWINTRA (currMB))
  {
    if (IS_OLDINTRA (currMB))   currSE.type = SE_CBP_INTRA;
    else                        currSE.type = SE_CBP_INTER;

    if(img->type == B_IMG_1 || img->type == B_IMG_MULT)  dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
    else                                                 dP = &(currSlice->partArr[partMap[currSE.type]]);
    
    if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag)
    {
      if (IS_OLDINTRA (currMB))  currSE.mapping = linfo_cbp_intra;
      else                       currSE.mapping = linfo_cbp_inter;
    }
    else
    {
      currSE.reading = readCBPFromBuffer_CABAC;
    }

#if TRACE
    snprintf(currSE.tracestring, TRACESTRING_SIZE, " CBP ");
#endif
    dP->readSyntaxElement(&currSE,img,inp,dP);
    currMB->cbp = cbp = currSE.value1;
    // Delta quant only if nonzero coeffs
    if (cbp !=0)
    {
      if (IS_INTER (currMB))  currSE.type = SE_DELTA_QUANT_INTER;
      else                    currSE.type = SE_DELTA_QUANT_INTRA;

      if(img->type == B_IMG_1 || img->type == B_IMG_MULT)  dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
      else                                                 dP = &(currSlice->partArr[partMap[currSE.type]]);
      
      if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag)
      {
        currSE.mapping = linfo_dquant;
      }
      else
      {
        if (IS_INTER (currMB))   currSE.reading= readDquant_inter_FromBuffer_CABAC;
        else                     currSE.reading= readDquant_intra_FromBuffer_CABAC;
      } 
#if TRACE
      snprintf(currSE.tracestring, TRACESTRING_SIZE, " Delta quant ");
#endif
      dP->readSyntaxElement(&currSE,img,inp,dP);
      currMB->delta_quant = currSE.value1;
      img->qp= (img->qp-MIN_QP+currMB->delta_quant+(MAX_QP-MIN_QP+1))%(MAX_QP-MIN_QP+1)+MIN_QP;
    }
  }
  else
  {
    cbp = currMB->cbp;
  }

  for (i=0;i<BLOCK_SIZE;i++)
    for (j=0;j<BLOCK_SIZE;j++)
      for(iii=0;iii<BLOCK_SIZE;iii++)
        for(jjj=0;jjj<BLOCK_SIZE;jjj++)
          img->cof[i][j][iii][jjj]=0;// reset luma coeffs


  if (IS_NEWINTRA (currMB)) // read DC coeffs for new intra modes
  {
    currSE.type = SE_DELTA_QUANT_INTRA;

    if(img->type == B_IMG_1 || img->type == B_IMG_MULT)  dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
    else                                                 dP = &(currSlice->partArr[partMap[currSE.type]]);
    
    if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag)
    {
      currSE.mapping = linfo_dquant;
    }
    else
    {
      currSE.reading= readDquant_intra_FromBuffer_CABAC;
    }
#if TRACE
    snprintf(currSE.tracestring, TRACESTRING_SIZE, " Delta quant ");
#endif
    dP->readSyntaxElement(&currSE,img,inp,dP);
    currMB->delta_quant = currSE.value1;
    img->qp= (img->qp-MIN_QP+currMB->delta_quant+(MAX_QP-MIN_QP+1))%(MAX_QP-MIN_QP+1)+MIN_QP;

    for (i=0;i<BLOCK_SIZE;i++)
      for (j=0;j<BLOCK_SIZE;j++)
        img->ipredmode[img->block_x+i+1][img->block_y+j+1]=0;


    currSE.type = SE_LUM_DC_INTRA;
    if(img->type == B_IMG_1 || img->type == B_IMG_MULT)  dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
    else                                                 dP = &(currSlice->partArr[partMap[currSE.type]]);

    if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag)
    {
      currSE.mapping = linfo_levrun_inter;
    }
    else
    {
      currSE.reading = readRunLevelFromBuffer_CABAC;
      currSE.context = 3; // for choosing context model
    }

    coef_ctr=-1;
    level = 1;                            // just to get inside the loop
    for(k=0;(k<17) && (level!=0);k++)
    {
#if TRACE
      snprintf(currSE.tracestring, TRACESTRING_SIZE, "DC luma 16x16 ");
#endif
      dP->readSyntaxElement(&currSE,img,inp,dP);
      level = currSE.value1;
      run   = currSE.value2;
      len   = currSE.len;

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

  qp_per    = (img->qp-MIN_QP)/6;
  qp_rem    = (img->qp-MIN_QP)%6;
  qp_per_uv = ((img->qp<0?img->qp:QP_SCALE_CR[img->qp])-MIN_QP)/6;
  qp_rem_uv = ((img->qp<0?img->qp:QP_SCALE_CR[img->qp])-MIN_QP)%6;
  currMB->qp = img->qp;

  // luma coefficients
  for (block_y=0; block_y < 4; block_y += 2) /* all modes */
  {
    for (block_x=0; block_x < 4; block_x += 2)
    {
      for (j=block_y; j < block_y+2; j++)
      {
        jj=j/2;
        for (i=block_x; i < block_x+2; i++)
        {
          ii = i/2;
          b8 = 2*jj+ii;

          if (IS_NEWINTRA (currMB))   start_scan = 1; /* skip DC coeff */
          else                        start_scan = 0; /* take all coeffs */

          img->subblock_x = i; // position for coeff_count ctx
          img->subblock_y = j; // position for coeff_count ctx
          if (cbp & (1<<b8))  /* are there any coeff in current block at all */
          {
            if (currMB->b8mode[b8]!=IBLOCK || (inp->symbol_mode!=CABAC && img->qp>=24))
            {
              coef_ctr = start_scan-1;
              level    = 1;      
              for(k=start_scan;(k<17) && (level!=0);k++)
              {
                /* 
                 * make distinction between INTRA and INTER coded
                 * luminance coefficients
                 */
                if (k == 0)
                { 
                  if (currMB->b8mode[b8]==IBLOCK || IS_NEWINTRA (currMB))
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
                  if (currMB->b8mode[b8]==IBLOCK)
                  {
                    currSE.context = 2; // for choosing context model
                    currSE.type  = SE_LUM_AC_INTRA;
                  }
                  else if (IS_NEWINTRA (currMB))
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
                sprintf(currSE.tracestring, " Luma sng ");
#endif
                if(img->type == B_IMG_1 || img->type == B_IMG_MULT)  dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
                else                                                 dP = &(currSlice->partArr[partMap[currSE.type]]);

                if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag)  currSE.mapping = linfo_levrun_inter;
                else                                                     currSE.reading = readRunLevelFromBuffer_CABAC;

                dP->readSyntaxElement(&currSE,img,inp,dP);
                level = currSE.value1;
                run   =  currSE.value2;
                len   = currSE.len;
                
                if (level != 0)    /* leave if len=1 */
                {
                  coef_ctr             += run+1;
                  i0=SNGL_SCAN[coef_ctr][0];
                  j0=SNGL_SCAN[coef_ctr][1];
                  currMB->cbp_blk      |= 1 << ((j<<2) + i) ;
                  img->cof[i][j][i0][j0]= level*dequant_coef[qp_rem][i0][j0]<<qp_per;
                }
              }
            }
            else    /* double scan (old intra with QP<24*/
            {
              for(scan_loop_ctr=0;scan_loop_ctr<2;scan_loop_ctr++)
              {
                coef_ctr=start_scan-1;
                level=1;                          /* just to get inside the loop */
                for(k=0; k<9 && level!=0;k++)
                {
                  if (k == 0)  currSE.type  = SE_LUM_DC_INTRA; /* element is of type DC */
                  else         currSE.type  = SE_LUM_AC_INTRA; /* element is of type AC */
#if TRACE
                  sprintf(currSE.tracestring, "Luma dbl(%2d,%2d)  ",scan_loop_ctr,k);
#endif
                  if(img->type == B_IMG_1 || img->type == B_IMG_MULT)  dP = &(currSlice->partArr[partMap[SE_BFRAME]]);
                  else                                                 dP = &(currSlice->partArr[partMap[currSE.type]]);

                  if (inp->symbol_mode == UVLC || dP->bitstream->ei_flag)
                  {
                    currSE.mapping = linfo_levrun_intra;
                  }
                  else
                  {
                    currSE.context = 0; // for choosing context model
                    currSE.reading = readRunLevelFromBuffer_CABAC;
                  }
                  dP->readSyntaxElement(&currSE,img,inp,dP);
                  level = currSE.value1;
                  run   = currSE.value2;
                  len   = currSE.len;

                  if (level != 0)    /* leave if len=1 */
                  {
                    coef_ctr              = coef_ctr+run+1;
                    currMB->cbp_blk      |= 1 << ((j<<2) + i) ;
                    i0=DBL_SCAN[coef_ctr][0][scan_loop_ctr];
                    j0=DBL_SCAN[coef_ctr][1][scan_loop_ctr];
                    img->cof[i][j][i0][j0]= level*dequant_coef[qp_rem][i0][j0]<<qp_per;
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

  m2 =img->mb_x*2;
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
        if (IS_INTRA (currMB))
        {
          currSE.context = 6; // for choosing context model
          currSE.type  = SE_CHR_DC_INTRA;
        }
        else
        {
          currSE.context = 5; // for choosing context model
          currSE.type  = SE_CHR_DC_INTER;
        }
        currSE.k = ll; //coeff_count ctx
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
          img->cofu[coef_ctr]=level*dequant_coef[qp_rem_uv][0][0]<<qp_per_uv;
        }
      }

      if (((img->type==SP_IMG_1 || img->type==SP_IMG_MULT) && IS_INTER (currMB)))
      {
        img->cof[0+ll][4][0][0]=(img->cofu[0]>>qp_per_uv)/dequant_coef[qp_rem_uv][0][0];
        img->cof[1+ll][4][0][0]=(img->cofu[1]>>qp_per_uv)/dequant_coef[qp_rem_uv][0][0];
        img->cof[0+ll][5][0][0]=(img->cofu[2]>>qp_per_uv)/dequant_coef[qp_rem_uv][0][0];
        img->cof[1+ll][5][0][0]=(img->cofu[3]>>qp_per_uv)/dequant_coef[qp_rem_uv][0][0];
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
  uv=-1;
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
            uv++;
            for(k=0;(k<16)&&(level!=0);k++)
            {
              currSE.k = uv;  //coeff_count ctx
              if (IS_INTRA (currMB))
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
                img->cof[i][j][i0][j0]=level*dequant_coef[qp_rem_uv][i0][j0]<<qp_per_uv;
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
  int ref_frame = 0; //currMB->ref_frame;
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
  int tmp_blockbw[BLOCK_SIZE][BLOCK_SIZE];
  int js[2][2];
  int i=0,j=0,k,ii=0,jj=0,i1=0,j1=0,j4=0,i4=0;
  int js0=0,js1=0,js2=0,js3=0,jf=0;
  int uv, hv;
  int vec1_x=0,vec1_y=0,vec2_x=0,vec2_y=0;
  int ioff,joff;

  int bw_pred, fw_pred, ifx;
  int ii0,jj0,ii1,jj1,if1,jf1,if0,jf0;
  int mv_mul,f1,f2,f3,f4;

  const byte decode_block_scan[16] = {0,1,4,5,2,3,6,7,8,9,12,13,10,11,14,15};

  Macroblock *currMB   = &img->mb_data[img->current_mb_nr];
  int refframe, fw_refframe, bw_refframe, mv_mode, pred_dir, intra_prediction; // = currMB->ref_frame;
  int*** mv_array, ***fw_mv_array, ***bw_mv_array;
  int bframe = (img->type==B_IMG_1 || img->type==B_IMG_MULT);
  byte refP_tr, TRb, TRp;

#ifdef _ADAPT_LAST_GROUP_
  extern int *last_P_no;
#endif

  int mb_nr             = img->current_mb_nr;
  int mb_width          = img->width/16;
  int mb_available_up   = (img->mb_y == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width].slice_nr);
  int mb_available_left = (img->mb_x == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-1].slice_nr);

  if(img->UseConstrainedIntraPred)
  {
    if (mb_available_up   && (img->intra_block[mb_nr-mb_width][2]==0 || img->intra_block[mb_nr-mb_width][3]==0))
      mb_available_up   = 0;
    if (mb_available_left && (img->intra_block[mb_nr-       1][1]==0 || img->intra_block[mb_nr       -1][3]==0))
      mb_available_left = 0;
  }


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
  if (IS_NEWINTRA (currMB))
  {
    intrapred_luma_2(img, currMB->i16mode);
  }

  for (k = 0; k < (MB_BLOCK_SIZE/BLOCK_SIZE)*(MB_BLOCK_SIZE/BLOCK_SIZE); k ++)
  {
    i = (decode_block_scan[k] & 3);
    j = ((decode_block_scan[k] >> 2) & 3);

    ioff=i*4;
    i4=img->block_x+i;

    joff=j*4;
    j4=img->block_y+j;

    mv_mode  = currMB->b8mode[2*(j/2)+(i/2)];
    pred_dir = currMB->b8pdir[2*(j/2)+(i/2)];
          
    // PREDICTION
    if (mv_mode==IBLOCK)
    {
      //===== INTRA PREDICTION =====
      if (intrapred(img,ioff,joff,i4,j4)==SEARCH_SYNC)  /* make 4x4 prediction block mpr from given prediction img->mb_mode */
        return SEARCH_SYNC;                   /* bit error */
    }
    else if (!IS_NEWINTRA (currMB))
    {
      if (pred_dir != 2)
      {
        //===== FORWARD/BACKWARD PREDICTION =====
        if (!bframe)
        {
          //refframe = (img->frame_cycle + img->buf_cycle - refFrArr[j4][i4]) % img->buf_cycle;
          refframe = refFrArr[j4][i4];
          mv_array = img->mv;
        }
        else if (!pred_dir)
        {
          //refframe = (img->number - 1 - img->fw_refFrArr[j4][i4] + img->buf_cycle) % img->buf_cycle;
          refframe = img->fw_refFrArr[j4][i4]+1;
          mv_array = img->fw_mv;
        }
        else
        {
          //refframe = (img->frame_cycle + img->buf_cycle) % img->buf_cycle;
          refframe = 0;
          mv_array = img->bw_mv;
        }

        vec1_x = i4*4*mv_mul + mv_array[i4+BLOCK_SIZE][j4][0];
        vec1_y = j4*4*mv_mul + mv_array[i4+BLOCK_SIZE][j4][1];

        get_block (refframe, vec1_x, vec1_y, img, tmp_block);

        for(ii=0;ii<BLOCK_SIZE;ii++)
        for(jj=0;jj<BLOCK_SIZE;jj++)  img->mpr[ii+ioff][jj+joff] = tmp_block[ii][jj];
      }
      else
      {
        if (mv_mode != 0)
        {
          //===== BI-DIRECTIONAL PREDICTION =====
          fw_mv_array = img->fw_mv;
          bw_mv_array = img->bw_mv;
          fw_refframe = img->fw_refFrArr[j4][i4]+1;
          bw_refframe = 0;
        }
        else
        {
          //===== DIRECT PREDICTION =====
          fw_mv_array = img->dfMV;
          bw_mv_array = img->dbMV;
          bw_refframe = 0;

          if(refFrArr[j4][i4]==-1) // next P is intra mode
          {
            for(hv=0; hv<2; hv++)   img->dfMV[i4+BLOCK_SIZE][j4][hv]=img->dbMV[i4+BLOCK_SIZE][j4][hv]=0;
            fw_refframe = 1;
          }
          else // next P is skip or inter mode
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
            //fw_refframe = (img->number - 1 - refFrArr[j4][i4] + img->buf_cycle) % img->buf_cycle;
            fw_refframe = refFrArr[j4][i4]+1;
          }
        }

        vec1_x = i4*4*mv_mul + fw_mv_array[i4+BLOCK_SIZE][j4][0];
        vec1_y = j4*4*mv_mul + fw_mv_array[i4+BLOCK_SIZE][j4][1];
        vec2_x = i4*4*mv_mul + bw_mv_array[i4+BLOCK_SIZE][j4][0];
        vec2_y = j4*4*mv_mul + bw_mv_array[i4+BLOCK_SIZE][j4][1];

        get_block(fw_refframe, vec1_x, vec1_y, img, tmp_block);
        get_block(bw_refframe, vec2_x, vec2_y, img, tmp_blockbw);

        for(ii=0;ii<BLOCK_SIZE;ii++)
        for(jj=0;jj<BLOCK_SIZE;jj++)  img->mpr[ii+ioff][jj+joff] = (tmp_block[ii][jj]+tmp_blockbw[ii][jj]+1)/2;
      }
    }

    if ((img->type==SP_IMG_1 || img->type==SP_IMG_MULT) && (IS_INTER (currMB) && mv_mode!=IBLOCK))
    {
      itrans_sp(img,ioff,joff,i,j);
    }
    else
    {
      itrans   (img,ioff,joff,i,j);      // use DCT transform and make 4x4 block m7 from prediction block mpr
    }

    for(ii=0;ii<BLOCK_SIZE;ii++)
    {
      for(jj=0;jj<BLOCK_SIZE;jj++)
      {
        imgY[j4*BLOCK_SIZE+jj][i4*BLOCK_SIZE+ii]=img->m7[ii][jj]; // contruct picture from 4x4 blocks
      }
    }
  }


  // chroma decoding *******************************************************
  for(uv=0;uv<2;uv++)
  {
    intra_prediction = (IS_NEWINTRA (currMB)        ||
                        currMB->b8mode[0] == IBLOCK ||
                        currMB->b8mode[1] == IBLOCK ||
                        currMB->b8mode[2] == IBLOCK ||
                        currMB->b8mode[3] == IBLOCK);

    if (intra_prediction)
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

        mv_mode  = currMB->b8mode[2*(j-4)+i];
        pred_dir = currMB->b8pdir[2*(j-4)+i];

        // PREDICTION
        if (mv_mode==IBLOCK || IS_NEWINTRA (currMB))
        {
          //--- INTRA PREDICTION ---
          for (ii=0; ii<4; ii++)
          for (jj=0; jj<4; jj++)
          {
            img->mpr[ii+ioff][jj+joff]=js[i][j-4];
          }
        }
        else if (pred_dir != 2)
        {
          //--- FORWARD/BACKWARD PREDICTION ---
          if (!bframe)        mv_array = img->mv;
          else if (!pred_dir) mv_array = img->fw_mv;
          else                mv_array = img->bw_mv;

          for(jj=0;jj<4;jj++)
          {
            jf=(j4+jj)/2;
            for(ii=0;ii<4;ii++)
            {
              if1=(i4+ii)/2;

              if (!bframe)        refframe =           refFrArr[jf][if1];
              else if (!pred_dir) refframe = 1+img->fw_refFrArr[jf][if1];
              else                refframe = 0;

              i1=(img->pix_c_x+ii+ioff)*f1+mv_array[if1+4][jf][0];
              j1=(img->pix_c_y+jj+joff)*f1+mv_array[if1+4][jf][1];

              ii0=max (0, min (i1/f1, img->width_cr-1));
              jj0=max (0, min (j1/f1, img->height_cr-1));
              ii1=max (0, min ((i1+f2)/f1, img->width_cr-1));
              jj1=max (0, min ((j1+f2)/f1, img->height_cr-1));

              if1=(i1 & f2);
              jf1=(j1 & f2);
              if0=f1-if1;
              jf0=f1-jf1;
              img->mpr[ii+ioff][jj+joff]=(if0*jf0*mcef[refframe][uv][jj0][ii0]+
                                          if1*jf0*mcef[refframe][uv][jj0][ii1]+
                                          if0*jf1*mcef[refframe][uv][jj1][ii0]+
                                          if1*jf1*mcef[refframe][uv][jj1][ii1]+f4)/f3;
            }
          }
        }
        else
        {
          if (mv_mode != 0)
          {
            //===== BI-DIRECTIONAL PREDICTION =====
            fw_mv_array = img->fw_mv;
            bw_mv_array = img->bw_mv;
          }
          else
          {
            //===== DIRECT PREDICTION =====
            fw_mv_array = img->dfMV;
            bw_mv_array = img->dbMV;
          }

          for(jj=0;jj<4;jj++)
          {
            jf=(j4+jj)/2;
            for(ii=0;ii<4;ii++)
            {
              ifx=(i4+ii)/2;

              if (mv_mode != 0)
              {
                fw_refframe = 1+img->fw_refFrArr[jf][ifx];
                bw_refframe = 0;
              }
              else
              {
                bw_refframe = 0;
                if(refFrArr[jf][ifx]==-1)  fw_refframe = 1;
                else                       fw_refframe = 1+refFrArr[jf][ifx];
              }

              i1=(img->pix_c_x+ii+ioff)*f1+fw_mv_array[ifx+4][jf][0];
              j1=(img->pix_c_y+jj+joff)*f1+fw_mv_array[ifx+4][jf][1];

              ii0=max (0, min (i1/f1, img->width_cr-1));
              jj0=max (0, min (j1/f1, img->height_cr-1));
              ii1=max (0, min ((i1+f2)/f1, img->width_cr-1));
              jj1=max (0, min ((j1+f2)/f1, img->height_cr-1));

              if1=(i1 & f2);
              jf1=(j1 & f2);
              if0=f1-if1;
              jf0=f1-jf1;

              fw_pred=(if0*jf0*mcef[fw_refframe][uv][jj0][ii0]+
                       if1*jf0*mcef[fw_refframe][uv][jj0][ii1]+
                       if0*jf1*mcef[fw_refframe][uv][jj1][ii0]+
                       if1*jf1*mcef[fw_refframe][uv][jj1][ii1]+f4)/f3;

              i1=(img->pix_c_x+ii+ioff)*f1+bw_mv_array[ifx+4][jf][0];
              j1=(img->pix_c_y+jj+joff)*f1+bw_mv_array[ifx+4][jf][1];

              ii0=max (0, min (i1/f1, img->width_cr-1));
              jj0=max (0, min (j1/f1, img->height_cr-1));
              ii1=max (0, min ((i1+f2)/f1, img->width_cr-1));
              jj1=max (0, min ((j1+f2)/f1, img->height_cr-1));

              if1=(i1 & f2);
              jf1=(j1 & f2);
              if0=f1-if1;
              jf0=f1-jf1;

              bw_pred=(if0*jf0*mcef[bw_refframe][uv][jj0][ii0]+
                       if1*jf0*mcef[bw_refframe][uv][jj0][ii1]+
                       if0*jf1*mcef[bw_refframe][uv][jj1][ii0]+
                       if1*jf1*mcef[bw_refframe][uv][jj1][ii1]+f4)/f3;

              img->mpr[ii+ioff][jj+joff]=(int)((fw_pred+bw_pred)/2.+.5);
            }
          }
        }

        if ((img->type!=SP_IMG_1 && img->type!=SP_IMG_MULT) || IS_INTRA (currMB))
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

    if((img->type==SP_IMG_1 || img->type==SP_IMG_MULT) && IS_INTER (currMB))
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
