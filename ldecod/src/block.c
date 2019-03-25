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
 *  \file
 *      block.c
 *
 *  \brief
 *      Block functions
 *
 *  \author
 *      Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Inge Lille-Langøy          <inge.lille-langoy@telenor.com>
 *      - Rickard Sjoberg            <rickard.sjoberg@era.ericsson.se>
 ***********************************************************************
 */

#include "contributors.h"

#include <stdlib.h>
#include <math.h>

#include "block.h"


#define Q_BITS          15

static const int quant_coef[6][4][4] = {
  {{13107, 8066,13107, 8066},{ 8066, 5243, 8066, 5243},{13107, 8066,13107, 8066},{ 8066, 5243, 8066, 5243}},
  {{11916, 7490,11916, 7490},{ 7490, 4660, 7490, 4660},{11916, 7490,11916, 7490},{ 7490, 4660, 7490, 4660}},
  {{10082, 6554,10082, 6554},{ 6554, 4194, 6554, 4194},{10082, 6554,10082, 6554},{ 6554, 4194, 6554, 4194}},
  {{ 9362, 5825, 9362, 5825},{ 5825, 3647, 5825, 3647},{ 9362, 5825, 9362, 5825},{ 5825, 3647, 5825, 3647}},
  {{ 8192, 5243, 8192, 5243},{ 5243, 3355, 5243, 3355},{ 8192, 5243, 8192, 5243},{ 5243, 3355, 5243, 3355}},
  {{ 7282, 4559, 7282, 4559},{ 4559, 2893, 4559, 2893},{ 7282, 4559, 7282, 4559},{ 4559, 2893, 4559, 2893}}
};
static const int A[4][4] = {
  { 16, 20, 16, 20},
  { 20, 25, 20, 25},
  { 16, 20, 16, 20},
  { 20, 25, 20, 25}
};

// Notation for comments regarding prediction and predictors.
// The pels of the 4x4 block are labelled a..p. The predictor pels above
// are labelled A..H, from the left I..L, and from above left X, as follows:
//
//  X A B C D E F G H
//  I a b c d
//  J e f g h
//  K i j k l
//  L m n o p
//

// Predictor array index definitions
#define P_X (PredPel[0])
#define P_A (PredPel[1])
#define P_B (PredPel[2])
#define P_C (PredPel[3])
#define P_D (PredPel[4])
#define P_E (PredPel[5])
#define P_F (PredPel[6])
#define P_G (PredPel[7])
#define P_H (PredPel[8])
#define P_I (PredPel[9])
#define P_J (PredPel[10])
#define P_K (PredPel[11])
#define P_L (PredPel[12])

/*!
 ***********************************************************************
 * \brief
 *    makes and returns 4x4 blocks with all 5 intra prediction modes
 *
 * \return
 *    DECODING_OK   decoding of intraprediction mode was sucessfull            \n
 *    SEARCH_SYNC   search next sync element as errors while decoding occured
 ***********************************************************************
 */

int intrapred(
  struct img_par *img,  //!< image parameters
  int ioff,             //!< pixel offset X within MB
  int joff,             //!< pixel offset Y within MB
  int img_block_x,      //!< location of block X, multiples of 4
  int img_block_y)      //!< location of block Y, multiples of 4
{
  int i,j;
  int s0;
  int img_y,img_x;
  int PredPel[13];  // array of predictor pels

  int block_available_up;
  int block_available_up_right;
  int block_available_left;

  byte predmode = img->ipredmode[img_block_x+1][img_block_y+1];

  img_x=img_block_x*4;
  img_y=img_block_y*4;
  if (img->mb_field)
  {
  if (img->current_mb_nr%2)
  {
    predmode = img->ipredmode_bot[img_block_x+1][img_block_y+1];

    block_available_up = (img->ipredmode_bot[img_block_x+1][img_block_y] >=0);
    block_available_up_right  = (img->ipredmode_bot[img_x/BLOCK_SIZE+2][img_y/BLOCK_SIZE] >=0); // ???
    block_available_left = (img->ipredmode_bot[img_block_x][img_block_y+1] >=0);
  }
  else
  {
    predmode = img->ipredmode_top[img_block_x+1][img_block_y+1];

    block_available_up = (img->ipredmode_top[img_block_x+1][img_block_y] >=0);
    block_available_up_right  = (img->ipredmode_top[img_x/BLOCK_SIZE+2][img_y/BLOCK_SIZE] >=0); // ???
    block_available_left = (img->ipredmode_top[img_block_x][img_block_y+1] >=0);
  }
  }
  else
  {
    block_available_up = (img->ipredmode[img_block_x+1][img_block_y] >=0);              /// can use frm
    block_available_up_right  = (img->ipredmode[img_x/BLOCK_SIZE+2][img_y/BLOCK_SIZE] >=0); // ???  /// can use frm
    block_available_left = (img->ipredmode[img_block_x][img_block_y+1] >=0);            /// can use frm
  }

  if(img_x%MB_BLOCK_SIZE == 12 && img->mb_frame_field_flag)
    block_available_up_right = 0;

  i = (img_x & 15);
  j = (img_y & 15);
  if (block_available_up_right)
  {
    if ((i == 4  && j == 4) ||
        (i == 12 && j == 4) ||
        (i == 12 && j == 8) ||
        (i == 4  && j == 12) ||
        (i == 12 && j == 12))
    {
      block_available_up_right = 0;
    }
  }

  // form predictor pels
  if (block_available_up)
  {
    P_A = imgY[img_y-1][img_x+0];
    P_B = imgY[img_y-1][img_x+1];
    P_C = imgY[img_y-1][img_x+2];
    P_D = imgY[img_y-1][img_x+3];

    if (block_available_up_right)
    {
      P_E = imgY[img_y-1][img_x+4];
      P_F = imgY[img_y-1][img_x+5];
      P_G = imgY[img_y-1][img_x+6];
      P_H = imgY[img_y-1][img_x+7];
    }
    else
    {
      P_E = P_F = P_G = P_H = P_D;
    }
  }
  else
  {
    P_A = P_B = P_C = P_D = P_E = P_F = P_G = P_H = 128;
  }

  if (block_available_left)
  {
    P_I = imgY[img_y+0][img_x-1];
    P_J = imgY[img_y+1][img_x-1];
    P_K = imgY[img_y+2][img_x-1];
    P_L = imgY[img_y+3][img_x-1];
  }
  else
  {
    P_I = P_J = P_K = P_L = 128;
  }

  if (block_available_up && block_available_left)
  {
    P_X = imgY[img_y-1][img_x-1];
  }
  else
  {
    P_X = 128;
  }

  
  switch (predmode)
  {
  case DC_PRED:                         /* DC prediction */

    s0 = 0;
    if (block_available_up && block_available_left)
    {   
      // no edge
      s0 = (P_A + P_B + P_C + P_D + P_I + P_J + P_K + P_L + 4)/(2*BLOCK_SIZE);
    }
    else if (!block_available_up && block_available_left)
    {
      // upper edge
      s0 = (P_I + P_J + P_K + P_L + 2)/BLOCK_SIZE;             
    }
    else if (block_available_up && !block_available_left)
    {
      // left edge
      s0 = (P_A + P_B + P_C + P_D + 2)/BLOCK_SIZE;             
    }
    else //if (!block_available_up && !block_available_left)
    {
      // top left corner, nothing to predict from
      s0 = 128;                           
    }

    for (j=0; j < BLOCK_SIZE; j++)
    {
      for (i=0; i < BLOCK_SIZE; i++)
      {
        // store DC prediction
        img->mpr[i+ioff][j+joff] = s0;
      }
    }
    break;

  case VERT_PRED:                       /* vertical prediction from block above */
    for(j=0;j<BLOCK_SIZE;j++)
      for(i=0;i<BLOCK_SIZE;i++)
        img->mpr[i+ioff][j+joff]=imgY[img_y-1][img_x+i];/* store predicted 4x4 block */
    break;

  case HOR_PRED:                        /* horisontal prediction from left block */
    for(j=0;j<BLOCK_SIZE;j++)
      for(i=0;i<BLOCK_SIZE;i++)
        img->mpr[i+ioff][j+joff]=imgY[img_y+j][img_x-1]; /* store predicted 4x4 block */
    break;

  case DIAG_DOWN_RIGHT_PRED:
    img->mpr[0+ioff][3+joff] = (P_L + 2*P_K + P_J + 2) / 4; 
    img->mpr[0+ioff][2+joff] =
    img->mpr[1+ioff][3+joff] = (P_K + 2*P_J + P_I + 2) / 4; 
    img->mpr[0+ioff][1+joff] =
    img->mpr[1+ioff][2+joff] = 
    img->mpr[2+ioff][3+joff] = (P_J + 2*P_I + P_X + 2) / 4; 
    img->mpr[0+ioff][0+joff] =
    img->mpr[1+ioff][1+joff] =
    img->mpr[2+ioff][2+joff] =
    img->mpr[3+ioff][3+joff] = (P_I + 2*P_X + P_A + 2) / 4; 
    img->mpr[1+ioff][0+joff] =
    img->mpr[2+ioff][1+joff] =
    img->mpr[3+ioff][2+joff] = (P_X + 2*P_A + P_B + 2) / 4;
    img->mpr[2+ioff][0+joff] =
    img->mpr[3+ioff][1+joff] = (P_A + 2*P_B + P_C + 2) / 4;
    img->mpr[3+ioff][0+joff] = (P_B + 2*P_C + P_D + 2) / 4;
    break;

  case DIAG_DOWN_LEFT_PRED:
    img->mpr[0+ioff][0+joff] = (P_A + P_C + 2*(P_B) + 2) / 4;
    img->mpr[1+ioff][0+joff] = 
    img->mpr[0+ioff][1+joff] = (P_B + P_D + 2*(P_C) + 2) / 4;
    img->mpr[2+ioff][0+joff] =
    img->mpr[1+ioff][1+joff] =
    img->mpr[0+ioff][2+joff] = (P_C + P_E + 2*(P_D) + 2) / 4;
    img->mpr[3+ioff][0+joff] = 
    img->mpr[2+ioff][1+joff] = 
    img->mpr[1+ioff][2+joff] = 
    img->mpr[0+ioff][3+joff] = (P_D + P_F + 2*(P_E) + 2) / 4;
    img->mpr[3+ioff][1+joff] = 
    img->mpr[2+ioff][2+joff] = 
    img->mpr[1+ioff][3+joff] = (P_E + P_G + 2*(P_F) + 2) / 4;
    img->mpr[3+ioff][2+joff] = 
    img->mpr[2+ioff][3+joff] = (P_F + P_H + 2*(P_G) + 2) / 4;
    img->mpr[3+ioff][3+joff] = (P_G + 3*(P_H) + 2) / 4;
    break;

  case  VERT_RIGHT_PRED:/* diagonal prediction -22.5 deg to horizontal plane */
    img->mpr[0+ioff][0+joff] = 
    img->mpr[1+ioff][2+joff] = (P_X + P_A + 1) / 2;
    img->mpr[1+ioff][0+joff] = 
    img->mpr[2+ioff][2+joff] = (P_A + P_B + 1) / 2;
    img->mpr[2+ioff][0+joff] = 
    img->mpr[3+ioff][2+joff] = (P_B + P_C + 1) / 2;
    img->mpr[3+ioff][0+joff] = (P_C + P_D + 1) / 2;
    img->mpr[0+ioff][1+joff] = 
    img->mpr[1+ioff][3+joff] = (P_I + 2*P_X + P_A + 2) / 4;
    img->mpr[1+ioff][1+joff] = 
    img->mpr[2+ioff][3+joff] = (P_X + 2*P_A + P_B + 2) / 4;
    img->mpr[2+ioff][1+joff] = 
    img->mpr[3+ioff][3+joff] = (P_A + 2*P_B + P_C + 2) / 4;
    img->mpr[3+ioff][1+joff] = (P_B + 2*P_C + P_D + 2) / 4;
    img->mpr[0+ioff][2+joff] = (P_X + 2*P_I + P_J + 2) / 4;
    img->mpr[0+ioff][3+joff] = (P_I + 2*P_J + P_K + 2) / 4;
    break;

  case  VERT_LEFT_PRED:/* diagonal prediction -22.5 deg to horizontal plane */
    img->mpr[0+ioff][0+joff] = (P_A + P_B + 1) / 2;
    img->mpr[1+ioff][0+joff] = 
    img->mpr[0+ioff][2+joff] = (P_B + P_C + 1) / 2;
    img->mpr[2+ioff][0+joff] = 
    img->mpr[1+ioff][2+joff] = (P_C + P_D + 1) / 2;
    img->mpr[3+ioff][0+joff] = 
    img->mpr[2+ioff][2+joff] = (P_D + P_E + 1) / 2;
    img->mpr[3+ioff][2+joff] = (P_E + P_F + 1) / 2;
    img->mpr[0+ioff][1+joff] = (P_A + 2*P_B + P_C + 2) / 4;
    img->mpr[1+ioff][1+joff] = 
    img->mpr[0+ioff][3+joff] = (P_B + 2*P_C + P_D + 2) / 4;
    img->mpr[2+ioff][1+joff] = 
    img->mpr[1+ioff][3+joff] = (P_C + 2*P_D + P_E + 2) / 4;
    img->mpr[3+ioff][1+joff] = 
    img->mpr[2+ioff][3+joff] = (P_D + 2*P_E + P_F + 2) / 4;
    img->mpr[3+ioff][3+joff] = (P_E + 2*P_F + P_G + 2) / 4;
    break;

  case  HOR_UP_PRED:/* diagonal prediction -22.5 deg to horizontal plane */
    img->mpr[0+ioff][0+joff] = (P_I + P_J + 1) / 2;
    img->mpr[1+ioff][0+joff] = (P_I + 2*P_J + P_K + 2) / 4;
    img->mpr[2+ioff][0+joff] = 
    img->mpr[0+ioff][1+joff] = (P_J + P_K + 1) / 2;
    img->mpr[3+ioff][0+joff] = 
    img->mpr[1+ioff][1+joff] = (P_J + 2*P_K + P_L + 2) / 4;
    img->mpr[2+ioff][1+joff] = 
    img->mpr[0+ioff][2+joff] = (P_K + P_L + 1) / 2;
    img->mpr[3+ioff][1+joff] = 
    img->mpr[1+ioff][2+joff] = (P_K + 2*P_L + P_L + 2) / 4;
    img->mpr[3+ioff][2+joff] = 
    img->mpr[1+ioff][3+joff] = 
    img->mpr[0+ioff][3+joff] = 
    img->mpr[2+ioff][2+joff] = 
    img->mpr[2+ioff][3+joff] = 
    img->mpr[3+ioff][3+joff] = P_L;
    break;

  case  HOR_DOWN_PRED:/* diagonal prediction -22.5 deg to horizontal plane */
    img->mpr[0+ioff][0+joff] = 
    img->mpr[2+ioff][1+joff] = (P_X + P_I + 1) / 2;
    img->mpr[1+ioff][0+joff] = 
    img->mpr[3+ioff][1+joff] = (P_I + 2*P_X + P_A + 2) / 4;
    img->mpr[2+ioff][0+joff] = (P_X + 2*P_A + P_B + 2) / 4;
    img->mpr[3+ioff][0+joff] = (P_A + 2*P_B + P_C + 2) / 4;
    img->mpr[0+ioff][1+joff] = 
    img->mpr[2+ioff][2+joff] = (P_I + P_J + 1) / 2;
    img->mpr[1+ioff][1+joff] = 
    img->mpr[3+ioff][2+joff] = (P_X + 2*P_I + P_J + 2) / 4;
    img->mpr[0+ioff][2+joff] = 
    img->mpr[2+ioff][3+joff] = (P_J + P_K + 1) / 2;
    img->mpr[1+ioff][2+joff] = 
    img->mpr[3+ioff][3+joff] = (P_I + 2*P_J + P_K + 2) / 4;
    img->mpr[0+ioff][3+joff] = (P_K + P_L + 1) / 2;
    img->mpr[1+ioff][3+joff] = (P_J + 2*P_K + P_L + 2) / 4;
    break;

  default:
    printf("Error: illegal prediction mode input: %d\n",predmode);
    return SEARCH_SYNC;
    break;
  }

  return DECODING_OK;
}


/*!
 ***********************************************************************
 * \return
 *    best SAD
 ***********************************************************************
 */
int intrapred_luma_16x16(struct img_par *img, //!< image parameters
                         int predmode)        //!< prediction mode
{
  int s0=0,s1,s2;

  int i,j;

  int ih,iv;
  int ib,ic,iaa;

  int mb_width = img->width/16;
  //int mb_nr_frame = img->mb_y*mb_width+img->mb_x;
  int mb_nr = img->map_mb_nr;// GB img->current_mb_nr;
  int mb_available_up = (img->mb_y == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width].slice_nr);
  int mb_available_left = (img->mb_x == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-1].slice_nr);
  int field_y = img->pix_y;   // For MB level frame/field coding
  if (img->mb_field)
  {
    field_y /= 2;
    mb_available_up = (img->mb_y/2 == 0) ? 0 : 1;

  if (img->current_mb_nr%2) // bottom field
  {
    field_y -= BLOCK_SIZE*2;
      mb_available_up = ((img->mb_y-1)/2 == 0) ? 0 : 1;
  }
  }

  if(img->constrained_intra_pred_flag)
  {
    if (mb_available_up   && (img->intra_block[mb_nr-mb_width][2]==0 || img->intra_block[mb_nr-mb_width][3]==0))
      mb_available_up   = 0;
    if (mb_available_left && (img->intra_block[mb_nr-       1][1]==0 || img->intra_block[mb_nr       -1][3]==0))
      mb_available_left = 0;
  }

  s1=s2=0;

  switch (predmode)
  {
  case VERT_PRED_16:                       // vertical prediction from block above
    for(j=0;j<MB_BLOCK_SIZE;j++)
      for(i=0;i<MB_BLOCK_SIZE;i++)
        img->mpr[i][j]=imgY[field_y-1][img->pix_x+i];// store predicted 16x16 block
    break;

  case HOR_PRED_16:                        // horisontal prediction from left block
    for(j=0;j<MB_BLOCK_SIZE;j++)
      for(i=0;i<MB_BLOCK_SIZE;i++)
        img->mpr[i][j]=imgY[field_y+j][img->pix_x-1]; // store predicted 16x16 block
    break;

  case DC_PRED_16:                         // DC prediction
    s1=s2=0;
    for (i=0; i < MB_BLOCK_SIZE; i++)
    {
      if (mb_available_up)
        s1 += imgY[field_y-1][img->pix_x+i];    // sum hor pix
      if (mb_available_left)
        s2 += imgY[field_y+i][img->pix_x-1];    // sum vert pix
    }
    if (mb_available_up && mb_available_left)
      s0=(s1+s2+16)/(2*MB_BLOCK_SIZE);       // no edge
    if (!mb_available_up && mb_available_left)
      s0=(s2+8)/MB_BLOCK_SIZE;              // upper edge
    if (mb_available_up && !mb_available_left)
      s0=(s1+8)/MB_BLOCK_SIZE;              // left edge
    if (!mb_available_up && !mb_available_left)
      s0=128;                            // top left corner, nothing to predict from
    for(i=0;i<MB_BLOCK_SIZE;i++)
      for(j=0;j<MB_BLOCK_SIZE;j++)
      {
        img->mpr[i][j]=s0;
      }
    break;
  case PLANE_16:// 16 bit integer plan pred
    ih=0;
    iv=0;
    for (i=1;i<9;i++)
    {
      ih += i*(imgY[field_y-1][img->pix_x+7+i] - imgY[field_y-1][img->pix_x+7-i]);
      iv += i*(imgY[field_y+7+i][img->pix_x-1] - imgY[field_y+7-i][img->pix_x-1]);
    }
    ib=(5*ih+32)>>6;
    ic=(5*iv+32)>>6;

    iaa=16*(imgY[field_y-1][img->pix_x+15]+imgY[field_y+15][img->pix_x-1]);
    for (j=0;j< MB_BLOCK_SIZE;j++)
    {
      for (i=0;i< MB_BLOCK_SIZE;i++)
      {
        img->mpr[i][j]=max(0,min((iaa+(i-7)*ib +(j-7)*ic + 16)/32,255));
      }
    }// store plane prediction
    break;

  default:
    {                                    // indication of fault in bitstream,exit
      printf("Error: illegal prediction mode input: %d\n",predmode);
      return SEARCH_SYNC;
    }
  }

  return DECODING_OK;
}

/*!
 ***********************************************************************
 * \brief
 *    Inverse 4x4 transformation, transforms cof to m7
 ***********************************************************************
 */
void itrans(struct img_par *img, //!< image parameters
            int ioff,            //!< index to 4x4 block
            int joff,            //!<
            int i0,              //!<
            int j0)              //!<
{
  int i,j,i1,j1;
  int m5[4];
  int m6[4];

  // horizontal
  for (j=0;j<BLOCK_SIZE;j++)
  {
    for (i=0;i<BLOCK_SIZE;i++)
    {
      m5[i]=img->cof[i0][j0][i][j];
    }
    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (i=0;i<2;i++)
    {
      i1=3-i;
      img->m7[i][j]=m6[i]+m6[i1];
      img->m7[i1][j]=m6[i]-m6[i1];
    }
  }
  // vertical
  for (i=0;i<BLOCK_SIZE;i++)
  {
    for (j=0;j<BLOCK_SIZE;j++)
      m5[j]=img->m7[i][j];

    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (j=0;j<2;j++)
    {
      j1=3-j;
      img->m7[i][j] =mmax(0,mmin(255,(m6[j]+m6[j1]+(img->mpr[i+ioff][j+joff] <<DQ_BITS)+DQ_ROUND)>>DQ_BITS));
      img->m7[i][j1]=mmax(0,mmin(255,(m6[j]-m6[j1]+(img->mpr[i+ioff][j1+joff]<<DQ_BITS)+DQ_ROUND)>>DQ_BITS));
    }
  }

}


/*!
 ***********************************************************************
 * \brief
 *    invers  transform
 ***********************************************************************
 */
void itrans_2(
   struct img_par *img) //!< image parameters
{
  int i,j,i1,j1;
  int M5[4];
  int M6[4];

  int qp_per = (img->qp-MIN_QP)/6;
  int qp_rem = (img->qp-MIN_QP)%6;

  // horizontal
  for (j=0;j<4;j++)
  {
    for (i=0;i<4;i++)
      M5[i]=img->cof[i][j][0][0];

    M6[0]=M5[0]+M5[2];
    M6[1]=M5[0]-M5[2];
    M6[2]=M5[1]-M5[3];
    M6[3]=M5[1]+M5[3];

    for (i=0;i<2;i++)
    {
      i1=3-i;
      img->cof[i ][j][0][0]= M6[i]+M6[i1];
      img->cof[i1][j][0][0]=M6[i]-M6[i1];
    }
  }

  // vertical
  for (i=0;i<4;i++)
  {
    for (j=0;j<4;j++)
      M5[j]=img->cof[i][j][0][0];

    M6[0]=M5[0]+M5[2];
    M6[1]=M5[0]-M5[2];
    M6[2]=M5[1]-M5[3];
    M6[3]=M5[1]+M5[3];

    for (j=0;j<2;j++)
    {
      j1=3-j;
      img->cof[i][j][0][0] = (((M6[j]+M6[j1])*dequant_coef[qp_rem][0][0]<<qp_per)+2)>>2;
      img->cof[i][j1][0][0]= (((M6[j]-M6[j1])*dequant_coef[qp_rem][0][0]<<qp_per)+2)>>2;
    }
  }
}


void itrans_sp(struct img_par *img,  //!< image parameters
               int ioff,             //!< index to 4x4 block
               int joff,             //!<
               int i0,               //!<
               int j0)               //!<
{
  int i,j,i1,j1;
  int m5[4];
  int m6[4];
  int predicted_block[BLOCK_SIZE][BLOCK_SIZE],ilev;
  
  int qp_per = (img->qp-MIN_QP)/6;
  int qp_rem = (img->qp-MIN_QP)%6;
  int q_bits    = Q_BITS+qp_per;

  int qp_per_sp = (img->qpsp-MIN_QP)/6;
  int qp_rem_sp = (img->qpsp-MIN_QP)%6;
  int q_bits_sp    = Q_BITS+qp_per_sp;
  int qp_const2=(1<<q_bits_sp)/2;  //sp_pred

  if (img->sp_switch || img->type == SI_IMG)
  {
    qp_per = (img->qpsp-MIN_QP)/6;
    qp_rem = (img->qpsp-MIN_QP)%6;
    q_bits = Q_BITS+qp_per;
  }
  for (j=0; j< BLOCK_SIZE; j++)
  for (i=0; i< BLOCK_SIZE; i++)
      predicted_block[i][j]=img->mpr[i+ioff][j+joff];
  for (j=0; j < BLOCK_SIZE; j++)
  {
    for (i=0; i < 2; i++)
    {
      i1=3-i;
      m5[i]=predicted_block[i][j]+predicted_block[i1][j];
      m5[i1]=predicted_block[i][j]-predicted_block[i1][j];
    }
    predicted_block[0][j]=(m5[0]+m5[1]);
    predicted_block[2][j]=(m5[0]-m5[1]);
    predicted_block[1][j]=m5[3]*2+m5[2];
    predicted_block[3][j]=m5[3]-m5[2]*2;
  }

  //  Vertival transform

  for (i=0; i < BLOCK_SIZE; i++)
  {
    for (j=0; j < 2; j++)
    {
      j1=3-j;
      m5[j]=predicted_block[i][j]+predicted_block[i][j1];
      m5[j1]=predicted_block[i][j]-predicted_block[i][j1];
    }
    predicted_block[i][0]=(m5[0]+m5[1]);
    predicted_block[i][2]=(m5[0]-m5[1]);
    predicted_block[i][1]=m5[3]*2+m5[2];
    predicted_block[i][3]=m5[3]-m5[2]*2;
  }

  for (j=0;j<BLOCK_SIZE;j++)
  for (i=0;i<BLOCK_SIZE;i++)
  {
    // recovering coefficient since they are already dequantized earlier
    img->cof[i0][j0][i][j]=(img->cof[i0][j0][i][j] >> qp_per) / dequant_coef[qp_rem][i][j]; 
    ilev=((img->cof[i0][j0][i][j]*dequant_coef[qp_rem][i][j]*A[i][j]<< qp_per) >>6)+predicted_block[i][j] ;
    img->cof[i0][j0][i][j]=sign((abs(ilev) * quant_coef[qp_rem_sp][i][j] + qp_const2) >> q_bits_sp, ilev) * dequant_coef[qp_rem_sp][i][j] << qp_per_sp;
  }
  // horizontal
  for (j=0;j<BLOCK_SIZE;j++)
  {
    for (i=0;i<BLOCK_SIZE;i++)
    {
      m5[i]=img->cof[i0][j0][i][j];
    }
    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (i=0;i<2;i++)
    {
      i1=3-i;
      img->m7[i][j]=m6[i]+m6[i1];
      img->m7[i1][j]=m6[i]-m6[i1];
    }
  }
  // vertical
  for (i=0;i<BLOCK_SIZE;i++)
  {
    for (j=0;j<BLOCK_SIZE;j++)
      m5[j]=img->m7[i][j];

    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (j=0;j<2;j++)
    {
      j1=3-j;
      img->m7[i][j] =mmax(0,mmin(255,(m6[j]+m6[j1]+DQ_ROUND)>>DQ_BITS));
      img->m7[i][j1]=mmax(0,mmin(255,(m6[j]-m6[j1]+DQ_ROUND)>>DQ_BITS));
    }
  }
}

/*!
 ***********************************************************************
 * \brief
 *    The routine performs transform,quantization,inverse transform, adds the diff.
 *    to the prediction and writes the result to the decoded luma frame. Includes the
 *    RD constrained quantization also.
 *
 * \par Input:
 *    block_x,block_y: Block position inside a macro block (0,4,8,12).
 *
 * \par Output:
 *    nonzero: 0 if no levels are nonzero.  1 if there are nonzero levels. \n
 *    coeff_cost: Counter for nonzero coefficients, used to discard expencive levels.
 ************************************************************************
 */
void copyblock_sp(struct img_par *img,int block_x,int block_y)
{
  int sign(int a,int b);

  int i,j,i1,j1,m5[4],m6[4];

  int predicted_block[BLOCK_SIZE][BLOCK_SIZE];
  int qp_per = (img->qpsp-MIN_QP)/6;
  int qp_rem = (img->qpsp-MIN_QP)%6;
  int q_bits    = Q_BITS+qp_per;
  int qp_const2=(1<<q_bits)/2;  //sp_pred


  //  Horizontal transform
  for (j=0; j< BLOCK_SIZE; j++)
  for (i=0; i< BLOCK_SIZE; i++)
    predicted_block[i][j]=img->mpr[i+block_x][j+block_y];

  for (j=0; j < BLOCK_SIZE; j++)
  {
    for (i=0; i < 2; i++)
    {
      i1=3-i;
      m5[i]=predicted_block[i][j]+predicted_block[i1][j];
      m5[i1]=predicted_block[i][j]-predicted_block[i1][j];
    }
    predicted_block[0][j]=(m5[0]+m5[1]);
    predicted_block[2][j]=(m5[0]-m5[1]);
    predicted_block[1][j]=m5[3]*2+m5[2];
    predicted_block[3][j]=m5[3]-m5[2]*2;
  }

  //  Vertival transform

  for (i=0; i < BLOCK_SIZE; i++)
  {
    for (j=0; j < 2; j++)
    {
      j1=3-j;
      m5[j]=predicted_block[i][j]+predicted_block[i][j1];
      m5[j1]=predicted_block[i][j]-predicted_block[i][j1];
    }
    predicted_block[i][0]=(m5[0]+m5[1]);
    predicted_block[i][2]=(m5[0]-m5[1]);
    predicted_block[i][1]=m5[3]*2+m5[2];
    predicted_block[i][3]=m5[3]-m5[2]*2;
  }

  // Quant
  for (j=0;j < BLOCK_SIZE; j++)
  for (i=0; i < BLOCK_SIZE; i++)
    img->m7[i][j]=sign((abs(predicted_block[i][j])* quant_coef[qp_rem][i][j]+qp_const2)>> q_bits,predicted_block[i][j])*dequant_coef[qp_rem][i][j]<<qp_per;

  //     IDCT.
  //     horizontal

  for (j=0;j<BLOCK_SIZE;j++)
  {
    for (i=0;i<BLOCK_SIZE;i++)
    {
      m5[i]=img->m7[i][j];
    }
    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (i=0;i<2;i++)
    {
      i1=3-i;
      img->m7[i][j]=m6[i]+m6[i1];
      img->m7[i1][j]=m6[i]-m6[i1];
    }
  }
  // vertical
  for (i=0;i<BLOCK_SIZE;i++)
  {
    for (j=0;j<BLOCK_SIZE;j++)
      m5[j]=img->m7[i][j];

    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (j=0;j<2;j++)
    {
      j1=3-j;
      img->m7[i][j] =mmax(0,mmin(255,(m6[j]+m6[j1]+DQ_ROUND)>>DQ_BITS));
      img->m7[i][j1]=mmax(0,mmin(255,(m6[j]-m6[j1]+DQ_ROUND)>>DQ_BITS));
    }
  }

  //  Decoded block moved to frame memory

  for (j=0; j < BLOCK_SIZE; j++)
    for (i=0; i < BLOCK_SIZE; i++)
      imgY[img->pix_y+block_y+j][img->pix_x+block_x+i]=img->m7[i][j];

}

void itrans_sp_chroma(struct img_par *img,int ll)
{
  int i,j,i1,j2,ilev,n2,n1,j1,mb_y;
  int m5[BLOCK_SIZE];
  int predicted_chroma_block[MB_BLOCK_SIZE/2][MB_BLOCK_SIZE/2],mp1[BLOCK_SIZE];
  int qp_per,qp_rem,q_bits;
  int qp_per_sp,qp_rem_sp,q_bits_sp,qp_const2;

  qp_per    = ((img->qp<0?img->qp:QP_SCALE_CR[img->qp])-MIN_QP)/6;
  qp_rem    = ((img->qp<0?img->qp:QP_SCALE_CR[img->qp])-MIN_QP)%6;
  q_bits    = Q_BITS+qp_per;

  qp_per_sp    = ((img->qpsp<0?img->qpsp:QP_SCALE_CR[img->qpsp])-MIN_QP)/6;
  qp_rem_sp    = ((img->qpsp<0?img->qpsp:QP_SCALE_CR[img->qpsp])-MIN_QP)%6;
  q_bits_sp    = Q_BITS+qp_per_sp;
  qp_const2=(1<<q_bits_sp)/2;  //sp_pred

  if (img->sp_switch || img->type == SI_IMG)
  {
    qp_per    = ((img->qpsp < 0 ? img->qpsp : QP_SCALE_CR[img->qpsp]) - MIN_QP) / 6;
    qp_rem    = ((img->qpsp < 0 ? img->qpsp : QP_SCALE_CR[img->qpsp]) - MIN_QP) % 6;
    q_bits    = Q_BITS + qp_per;
  }

  for (j=0; j < MB_BLOCK_SIZE/2; j++)
  for (i=0; i < MB_BLOCK_SIZE/2; i++)
  {
    predicted_chroma_block[i][j]=img->mpr[i][j];
    img->mpr[i][j]=0;
  }
  for (n2=0; n2 <= BLOCK_SIZE; n2 += BLOCK_SIZE)
  {
    for (n1=0; n1 <= BLOCK_SIZE; n1 += BLOCK_SIZE)
    {
      //  Horizontal transform.
      for (j=0; j < BLOCK_SIZE; j++)
      {
        mb_y=n2+j;
        for (i=0; i < 2; i++)
        {
          i1=3-i;
          m5[i]=predicted_chroma_block[i+n1][mb_y]+predicted_chroma_block[i1+n1][mb_y];
          m5[i1]=predicted_chroma_block[i+n1][mb_y]-predicted_chroma_block[i1+n1][mb_y];
        }
        predicted_chroma_block[n1][mb_y]  =(m5[0]+m5[1]);
        predicted_chroma_block[n1+2][mb_y]=(m5[0]-m5[1]);
        predicted_chroma_block[n1+1][mb_y]=m5[3]*2+m5[2];
        predicted_chroma_block[n1+3][mb_y]=m5[3]-m5[2]*2;
      }

      //  Vertical transform.

      for (i=0; i < BLOCK_SIZE; i++)
      {
        j1=n1+i;
        for (j=0; j < 2; j++)
        {
          j2=3-j;
          m5[j]=predicted_chroma_block[j1][n2+j]+predicted_chroma_block[j1][n2+j2];
          m5[j2]=predicted_chroma_block[j1][n2+j]-predicted_chroma_block[j1][n2+j2];
        }
        predicted_chroma_block[j1][n2+0]=(m5[0]+m5[1]);
        predicted_chroma_block[j1][n2+2]=(m5[0]-m5[1]);
        predicted_chroma_block[j1][n2+1]=m5[3]*2+m5[2];
        predicted_chroma_block[j1][n2+3]=m5[3]-m5[2]*2;
      }
    }
  }

  //     2X2 transform of DC coeffs.
  mp1[0]=(predicted_chroma_block[0][0]+predicted_chroma_block[4][0]+predicted_chroma_block[0][4]+predicted_chroma_block[4][4]);
  mp1[1]=(predicted_chroma_block[0][0]-predicted_chroma_block[4][0]+predicted_chroma_block[0][4]-predicted_chroma_block[4][4]);
  mp1[2]=(predicted_chroma_block[0][0]+predicted_chroma_block[4][0]-predicted_chroma_block[0][4]-predicted_chroma_block[4][4]);
  mp1[3]=(predicted_chroma_block[0][0]-predicted_chroma_block[4][0]-predicted_chroma_block[0][4]+predicted_chroma_block[4][4]);

  for (n1=0; n1 < 2; n1 ++)
  for (n2=0; n2 < 2; n2 ++)
  {
    ilev=((img->cof[n1+ll][4+n2][0][0]*dequant_coef[qp_rem][0][0]*A[0][0]<< qp_per) >>5)+mp1[n1+n2*2] ;
    mp1[n1+n2*2]=sign((abs(ilev)* quant_coef[qp_rem_sp][0][0]+ 2 * qp_const2)>> (q_bits_sp+1),ilev)*dequant_coef[qp_rem_sp][0][0]<<qp_per_sp;
  }

  for (n2=0; n2 < 2; n2 ++)
  for (n1=0; n1 < 2; n1 ++)
  for (i=0;i< BLOCK_SIZE; i++)
  for (j=0;j< BLOCK_SIZE; j++)
  {
  // recovering coefficient since they are already dequantized earlier
    img->cof[n1+ll][4+n2][i][j] = (img->cof[n1+ll][4+n2][i][j] >> qp_per) / dequant_coef[qp_rem][i][j];
    ilev=((img->cof[n1+ll][4+n2][i][j]*dequant_coef[qp_rem][i][j]*A[i][j]<< qp_per) >>6)+predicted_chroma_block[n1*BLOCK_SIZE+i][n2*BLOCK_SIZE+j] ;
    img->cof[n1+ll][4+n2][i][j] = sign((abs(ilev) * quant_coef[qp_rem_sp][i][j] + qp_const2)>> q_bits_sp,ilev)*dequant_coef[qp_rem_sp][i][j]<<qp_per_sp;
  }
  img->cof[0+ll][4][0][0]=(mp1[0]+mp1[1]+mp1[2]+mp1[3])>>1;
  img->cof[1+ll][4][0][0]=(mp1[0]-mp1[1]+mp1[2]-mp1[3])>>1;
  img->cof[0+ll][5][0][0]=(mp1[0]+mp1[1]-mp1[2]-mp1[3])>>1;
  img->cof[1+ll][5][0][0]=(mp1[0]-mp1[1]-mp1[2]+mp1[3])>>1;
}

int sign(int a , int b)
{
  int x;

  x=abs(a);
  if (b>0)
    return(x);
  else return(-x);
}
