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
 * \file block.c
 *
 * \brief
 *    Process one block
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Inge Lille-Langøy               <inge.lille-langoy@telenor.com>
 *    - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *    - Stephan Wenger                  <stewe@cs.tu-berlin.de>
 *    - Jani Lainema                    <jani.lainema@nokia.com>
 *    - Detlev Marpe                    <marpe@hhi.de>
 *    - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
 *    - Ragip Kurceren                  <ragip.kurceren@nokia.com>
 *************************************************************************************
 */

#include "contributors.h"


#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "block.h"
#include "refbuf.h"


#define Q_BITS          15
#define DQ_BITS         6
#define DQ_ROUND        (1<<(DQ_BITS-1))


static const int quant_coef[6][4][4] = {
  {{13107, 8066,13107, 8066},{ 8066, 5243, 8066, 5243},{13107, 8066,13107, 8066},{ 8066, 5243, 8066, 5243}},
  {{11916, 7490,11916, 7490},{ 7490, 4660, 7490, 4660},{11916, 7490,11916, 7490},{ 7490, 4660, 7490, 4660}},
  {{10082, 6554,10082, 6554},{ 6554, 4194, 6554, 4194},{10082, 6554,10082, 6554},{ 6554, 4194, 6554, 4194}},
  {{ 9362, 5825, 9362, 5825},{ 5825, 3647, 5825, 3647},{ 9362, 5825, 9362, 5825},{ 5825, 3647, 5825, 3647}},
  {{ 8192, 5243, 8192, 5243},{ 5243, 3355, 5243, 3355},{ 8192, 5243, 8192, 5243},{ 5243, 3355, 5243, 3355}},
  {{ 7282, 4559, 7282, 4559},{ 4559, 2893, 4559, 2893},{ 7282, 4559, 7282, 4559},{ 4559, 2893, 4559, 2893}}
};

static const int dequant_coef[6][4][4] = {
  {{10, 13, 10, 13},{ 13, 16, 13, 16},{10, 13, 10, 13},{ 13, 16, 13, 16}},
  {{11, 14, 11, 14},{ 14, 18, 14, 18},{11, 14, 11, 14},{ 14, 18, 14, 18}},
  {{13, 16, 13, 16},{ 16, 20, 16, 20},{13, 16, 13, 16},{ 16, 20, 16, 20}},
  {{14, 18, 14, 18},{ 18, 23, 18, 23},{14, 18, 14, 18},{ 18, 23, 18, 23}},
  {{16, 20, 16, 20},{ 20, 25, 20, 25},{16, 20, 16, 20},{ 20, 25, 20, 25}},
  {{18, 23, 18, 23},{ 23, 29, 23, 29},{18, 23, 18, 23},{ 23, 29, 23, 29}}
};
static const int A[4][4] = {
  { 16, 20, 16, 20},
  { 20, 25, 20, 25},
  { 16, 20, 16, 20},
  { 20, 25, 20, 25}
};


// Notation for comments regarding prediction and predictors.
// The pels of the 4x4 block are labelled a..p. The predictor pels above
// are labelled A..H, from the left I..P, and from above left X, as follows:
//
//  X A B C D E F G H
//  I a b c d
//  J e f g h
//  K i j k l
//  L m n o p
//  M
//  N
//  O
//  P
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
#define P_M (PredPel[13])
#define P_N (PredPel[14])
#define P_O (PredPel[15])
#define P_P (PredPel[16])

/*!
 ************************************************************************
 * \brief
 *    Make intra 4x4 prediction according to all 9 prediction modes.
 *    The routine uses left and upper neighbouring points from
 *    previous coded blocks to do this (if available). Notice that
 *    inaccessible neighbouring points are signalled with a negative
 *    value in the predmode array .
 *
 *  \para Input:
 *     Starting point of current 4x4 block image posision
 *
 *  \para Output:
 *      none
 ************************************************************************
 */
void intrapred_luma(int img_x,int img_y)
{
  int i,j;
  int s0;
  int PredPel[17];  // array of predictor pels
  byte **imgY_pred = imgY;  // For MB level frame/field coding tools -- set default to imgY

  int block_available_up        = (img->ipredmode[img_x/BLOCK_SIZE+1][img_y/BLOCK_SIZE] >=0);
  int block_available_up_right  = (img->ipredmode[img_x/BLOCK_SIZE+2][img_y/BLOCK_SIZE] >=0);
  int block_available_left      = (img->ipredmode[img_x/BLOCK_SIZE][img_y/BLOCK_SIZE+1] >=0);
  int block_available_left_down = (img->ipredmode[img_x/BLOCK_SIZE][img_y/BLOCK_SIZE+2] >=0);

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    imgY_pred     = img->top_field ? imgY_top:imgY_bot;
    if(img->top_field)
    {
      block_available_up        = (img->ipredmode_top[img_x/BLOCK_SIZE+1][img_y/BLOCK_SIZE] >=0);
      block_available_up_right  = (img->ipredmode_top[img_x/BLOCK_SIZE+2][img_y/BLOCK_SIZE] >=0);
      block_available_left      = (img->ipredmode_top[img_x/BLOCK_SIZE][img_y/BLOCK_SIZE+1] >=0);
      block_available_left_down = (img->ipredmode_top[img_x/BLOCK_SIZE][img_y/BLOCK_SIZE+2] >=0);
    }
    else
    {
      block_available_up        = (img->ipredmode_bot[img_x/BLOCK_SIZE+1][img_y/BLOCK_SIZE] >=0);
      block_available_up_right  = (img->ipredmode_bot[img_x/BLOCK_SIZE+2][img_y/BLOCK_SIZE] >=0);
      block_available_left      = (img->ipredmode_bot[img_x/BLOCK_SIZE][img_y/BLOCK_SIZE+1] >=0);
      block_available_left_down = (img->ipredmode_bot[img_x/BLOCK_SIZE][img_y/BLOCK_SIZE+2] >=0);
    }
  }
  
  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
  {
    if(img_x%MB_BLOCK_SIZE == 12)   
      block_available_up_right = 0; // for MB pairs some blocks will not have block available up right  
  }
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

  if (block_available_left_down)
  {
    if (!(i == 0 && j == 0) &&
        !(i == 8 && j == 0) &&
        !(i == 0 && j == 4) &&
        !(i == 0 && j == 8) &&
        !(i == 8 && j == 8))
    {
      block_available_left_down = 0;
    }
  }

  // form predictor pels
  if (block_available_up)
  {
    P_A = imgY_pred[img_y-1][img_x+0];
    P_B = imgY_pred[img_y-1][img_x+1];
    P_C = imgY_pred[img_y-1][img_x+2];
    P_D = imgY_pred[img_y-1][img_x+3];

    if (block_available_up_right)
    {
      P_E = imgY_pred[img_y-1][img_x+4];
      P_F = imgY_pred[img_y-1][img_x+5];
      P_G = imgY_pred[img_y-1][img_x+6];
      P_H = imgY_pred[img_y-1][img_x+7];
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
    P_I = imgY_pred[img_y+0][img_x-1];
    P_J = imgY_pred[img_y+1][img_x-1];
    P_K = imgY_pred[img_y+2][img_x-1];
    P_L = imgY_pred[img_y+3][img_x-1];

    if (block_available_left_down)
    {
      P_M = imgY_pred[img_y+4][img_x-1];
      P_N = imgY_pred[img_y+5][img_x-1];
      P_O = imgY_pred[img_y+6][img_x-1];
      P_P = imgY_pred[img_y+7][img_x-1];
    }
    else
    {
      P_M = P_N = P_O = P_P = P_L;
    }
  }
  else
  {
    P_I = P_J = P_K = P_L = P_M = P_N = P_O = P_P = 128;
  }

  // XXXGJC -- not quite right for this boundary
  if (block_available_up && block_available_left)
  {
    P_X = imgY_pred[img_y-1][img_x-1];
  }
  else
  {
    P_X = 128;
  }

  for(i=0;i<9;i++)
    img->mprr[i][0][0]=-1;

  ///////////////////////////////
  // make DC prediction
  ///////////////////////////////
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
      img->mprr[DC_PRED][i][j] = s0;
    }
  }

  ///////////////////////////////
  // make horiz and vert prediction
  ///////////////////////////////

  for (i=0; i < BLOCK_SIZE; i++)
  {
    img->mprr[VERT_PRED][0][i] = 
    img->mprr[VERT_PRED][1][i] = 
    img->mprr[VERT_PRED][2][i] = 
    img->mprr[VERT_PRED][3][i] = (&P_A)[i];
    img->mprr[HOR_PRED][i][0]  = 
    img->mprr[HOR_PRED][i][1]  = 
    img->mprr[HOR_PRED][i][2]  = 
    img->mprr[HOR_PRED][i][3]  = (&P_I)[i];
  }

  if(!block_available_up)img->mprr[VERT_PRED][0][0]=-1;
  if(!block_available_left)img->mprr[HOR_PRED][0][0]=-1;

  /*  Prediction according to 'diagonal' modes */
  if (block_available_up && block_available_left)
  {
    // Mode DIAG_DOWN_RIGHT_PRED
    img->mprr[DIAG_DOWN_RIGHT_PRED][3][0] = (P_L + 2*P_K + P_J + 2) / 4; 
    img->mprr[DIAG_DOWN_RIGHT_PRED][2][0] =
    img->mprr[DIAG_DOWN_RIGHT_PRED][3][1] = (P_K + 2*P_J + P_I + 2) / 4; 
    img->mprr[DIAG_DOWN_RIGHT_PRED][1][0] =
    img->mprr[DIAG_DOWN_RIGHT_PRED][2][1] = 
    img->mprr[DIAG_DOWN_RIGHT_PRED][3][2] = (P_J + 2*P_I + P_X + 2) / 4; 
    img->mprr[DIAG_DOWN_RIGHT_PRED][0][0] =
    img->mprr[DIAG_DOWN_RIGHT_PRED][1][1] =
    img->mprr[DIAG_DOWN_RIGHT_PRED][2][2] =
    img->mprr[DIAG_DOWN_RIGHT_PRED][3][3] = (P_I + 2*P_X + P_A + 2) / 4; 
    img->mprr[DIAG_DOWN_RIGHT_PRED][0][1] =
    img->mprr[DIAG_DOWN_RIGHT_PRED][1][2] =
    img->mprr[DIAG_DOWN_RIGHT_PRED][2][3] = (P_X + 2*P_A + P_B + 2) / 4;
    img->mprr[DIAG_DOWN_RIGHT_PRED][0][2] =
    img->mprr[DIAG_DOWN_RIGHT_PRED][1][3] = (P_A + 2*P_B + P_C + 2) / 4;
    img->mprr[DIAG_DOWN_RIGHT_PRED][0][3] = (P_B + 2*P_C + P_D + 2) / 4;

    // Mode DIAG_DOWN_LEFT_PRED
    img->mprr[DIAG_DOWN_LEFT_PRED][0][0] = (P_A + P_C + P_I + P_K + 2*(P_B + P_J) + 4) / 8;
    img->mprr[DIAG_DOWN_LEFT_PRED][0][1] = 
    img->mprr[DIAG_DOWN_LEFT_PRED][1][0] = (P_B + P_D + P_J + P_L + 2*(P_C + P_K) + 4) / 8;
    img->mprr[DIAG_DOWN_LEFT_PRED][0][2] =
    img->mprr[DIAG_DOWN_LEFT_PRED][1][1] =
    img->mprr[DIAG_DOWN_LEFT_PRED][2][0] = (P_C + P_E + P_K + P_M + 2*(P_D + P_L) + 4) / 8;
    img->mprr[DIAG_DOWN_LEFT_PRED][0][3] = 
    img->mprr[DIAG_DOWN_LEFT_PRED][1][2] = 
    img->mprr[DIAG_DOWN_LEFT_PRED][2][1] = 
    img->mprr[DIAG_DOWN_LEFT_PRED][3][0] = (P_D + P_F + P_L + P_N + 2*(P_E + P_M) + 4) / 8;
    img->mprr[DIAG_DOWN_LEFT_PRED][1][3] = 
    img->mprr[DIAG_DOWN_LEFT_PRED][2][2] = 
    img->mprr[DIAG_DOWN_LEFT_PRED][3][1] = (P_E + P_G + P_M + P_O + 2*(P_F + P_N) + 4) / 8;
    img->mprr[DIAG_DOWN_LEFT_PRED][2][3] = 
    img->mprr[DIAG_DOWN_LEFT_PRED][3][2] = (P_F + P_H + P_N + P_P + 2*(P_G + P_O) + 4) / 8;
    img->mprr[DIAG_DOWN_LEFT_PRED][3][3] = (P_G + P_O + P_H + P_P + 2) / 4;

    // Mode VERT_RIGHT_PRED
    img->mprr[VERT_RIGHT_PRED][0][0] = 
    img->mprr[VERT_RIGHT_PRED][2][1] = (P_X + P_A + 1) / 2;
    img->mprr[VERT_RIGHT_PRED][0][1] = 
    img->mprr[VERT_RIGHT_PRED][2][2] = (P_A + P_B + 1) / 2;
    img->mprr[VERT_RIGHT_PRED][0][2] = 
    img->mprr[VERT_RIGHT_PRED][2][3] = (P_B + P_C + 1) / 2;
    img->mprr[VERT_RIGHT_PRED][0][3] = (P_C + P_D + 1) / 2;
    img->mprr[VERT_RIGHT_PRED][1][0] = 
    img->mprr[VERT_RIGHT_PRED][3][1] = (P_I + 2*P_X + P_A + 2) / 4;
    img->mprr[VERT_RIGHT_PRED][1][1] = 
    img->mprr[VERT_RIGHT_PRED][3][2] = (P_X + 2*P_A + P_B + 2) / 4;
    img->mprr[VERT_RIGHT_PRED][1][2] = 
    img->mprr[VERT_RIGHT_PRED][3][3] = (P_A + 2*P_B + P_C + 2) / 4;
    img->mprr[VERT_RIGHT_PRED][1][3] = (P_B + 2*P_C + P_D + 2) / 4;
    img->mprr[VERT_RIGHT_PRED][2][0] = (P_X + 2*P_I + P_J + 2) / 4;
    img->mprr[VERT_RIGHT_PRED][3][0] = (P_I + 2*P_J + P_K + 2) / 4;

    // Mode VERT_LEFT_PRED
    img->mprr[VERT_LEFT_PRED][0][0] = (2*(P_A + P_B + P_K) + P_J + P_L + 4) / 8;
    img->mprr[VERT_LEFT_PRED][0][1] = 
    img->mprr[VERT_LEFT_PRED][2][0] = (P_B + P_C + 1) / 2;
    img->mprr[VERT_LEFT_PRED][0][2] = 
    img->mprr[VERT_LEFT_PRED][2][1] = (P_C + P_D + 1) / 2;
    img->mprr[VERT_LEFT_PRED][0][3] = 
    img->mprr[VERT_LEFT_PRED][2][2] = (P_D + P_E + 1) / 2;
    img->mprr[VERT_LEFT_PRED][2][3] = (P_E + P_F + 1) / 2;
    img->mprr[VERT_LEFT_PRED][1][0] = (2*(P_B + P_L) + P_A + P_C + P_K + P_M + 4) / 8;
    img->mprr[VERT_LEFT_PRED][1][1] = 
    img->mprr[VERT_LEFT_PRED][3][0] = (P_B + 2*P_C + P_D + 2) / 4;
    img->mprr[VERT_LEFT_PRED][1][2] = 
    img->mprr[VERT_LEFT_PRED][3][1] = (P_C + 2*P_D + P_E + 2) / 4;
    img->mprr[VERT_LEFT_PRED][1][3] = 
    img->mprr[VERT_LEFT_PRED][3][2] = (P_D + 2*P_E + P_F + 2) / 4;
    img->mprr[VERT_LEFT_PRED][3][3] = (P_E + 2*P_F + P_G + 2) / 4;

    // Mode HOR_UP_PRED
    img->mprr[HOR_UP_PRED][0][0] = (2*(P_C + P_I + P_J) + P_B + P_D + 4) / 8;
    img->mprr[HOR_UP_PRED][0][1] = (2*(P_D + P_J) + P_C + P_E + P_I + P_K + 4) / 8;
    img->mprr[HOR_UP_PRED][0][2] = 
    img->mprr[HOR_UP_PRED][1][0] = (P_J + P_K + 1) / 2;
    img->mprr[HOR_UP_PRED][0][3] = 
    img->mprr[HOR_UP_PRED][1][1] = (P_J + 2*P_K + P_L + 2) / 4;
    img->mprr[HOR_UP_PRED][1][2] = 
    img->mprr[HOR_UP_PRED][2][0] = (P_K + P_L + 1) / 2;
    img->mprr[HOR_UP_PRED][1][3] = 
    img->mprr[HOR_UP_PRED][2][1] = (P_K + 2*P_L + P_M + 2) / 4;
    img->mprr[HOR_UP_PRED][3][0] = 
    img->mprr[HOR_UP_PRED][2][2] = (P_L + P_M + 1) / 2;
    img->mprr[HOR_UP_PRED][2][3] = 
    img->mprr[HOR_UP_PRED][3][1] = (P_L + 2*P_M + P_N + 2) / 4;
    img->mprr[HOR_UP_PRED][3][2] = (P_M + P_N + 1) / 2;
    img->mprr[HOR_UP_PRED][3][3] = (P_M + 2*P_N + P_O + 2) / 4;

    // Mode HOR_DOWN_PRED
    img->mprr[HOR_DOWN_PRED][0][0] = 
    img->mprr[HOR_DOWN_PRED][1][2] = (P_X + P_I + 1) / 2;
    img->mprr[HOR_DOWN_PRED][0][1] = 
    img->mprr[HOR_DOWN_PRED][1][3] = (P_I + 2*P_X + P_A + 2) / 4;
    img->mprr[HOR_DOWN_PRED][0][2] = (P_X + 2*P_A + P_B + 2) / 4;
    img->mprr[HOR_DOWN_PRED][0][3] = (P_A + 2*P_B + P_C + 2) / 4;
    img->mprr[HOR_DOWN_PRED][1][0] = 
    img->mprr[HOR_DOWN_PRED][2][2] = (P_I + P_J + 1) / 2;
    img->mprr[HOR_DOWN_PRED][1][1] = 
    img->mprr[HOR_DOWN_PRED][2][3] = (P_X + 2*P_I + P_J + 2) / 4;
    img->mprr[HOR_DOWN_PRED][2][0] = 
    img->mprr[HOR_DOWN_PRED][3][2] = (P_J + P_K + 1) / 2;
    img->mprr[HOR_DOWN_PRED][2][1] = 
    img->mprr[HOR_DOWN_PRED][3][3] = (P_I + 2*P_J + P_K + 2) / 4;
    img->mprr[HOR_DOWN_PRED][3][0] = (P_K + P_L + 1) / 2;
    img->mprr[HOR_DOWN_PRED][3][1] = (P_J + 2*P_K + P_L + 2) / 4;
  }
}

/*!
 ************************************************************************
 * \brief
 *    16x16 based luma prediction
 *
 * \para Input:
 *    Image parameters
 *
 * \para Output:
 *    none
 ************************************************************************
 */
void intrapred_luma_2()
{
  int s0=0,s1,s2;
  int i,j;
  int s[16][2];

  int ih,iv;
  int ib,ic,iaa;
  byte **imgY_pred = imgY;  // For Mb level field/frame coding tools -- default to frame pred
  int pix_y = img->pix_y; // For MB level field/frame coding tools


  int mb_nr = img->current_mb_nr;
  int mb_width = img->width/16;
  int mb_available_left = (img->mb_x == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-1].slice_nr);
  int mb_available_up = (img->mb_y == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width].slice_nr);

  if(input->UseConstrainedIntraPred)
  {
    if (mb_available_up   && (img->intra_block[mb_nr-mb_width][2]==0 || img->intra_block[mb_nr-mb_width][3]==0))
      mb_available_up   = 0;
    if (mb_available_left && (img->intra_block[mb_nr-       1][1]==0 || img->intra_block[mb_nr       -1][3]==0))
      mb_available_left = 0;
  }

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    if(img->top_field)
    {
      mb_available_up = (img->mb_y/2 == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width].slice_nr);
      pix_y   = img->field_pix_y; // set pix_y to field pix_y
      imgY_pred = imgY_top; // set the prediction image to top field
    }
    else
    {
      mb_available_up = ((img->mb_y-1)/2 == 0) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width].slice_nr);
      imgY_pred = imgY_bot;
      pix_y   = img->field_pix_y; // set pix_y to field pix_y
    }
  }

  
  s1=s2=0;
  // make DC prediction
  for (i=0; i < MB_BLOCK_SIZE; i++)
  {
    if (mb_available_up)
      s1 += imgY_pred[pix_y-1][img->pix_x+i];    // sum hor pix
    if (mb_available_left)
      s2 += imgY_pred[pix_y+i][img->pix_x-1];    // sum vert pix
  }
  if (mb_available_up && mb_available_left)
    s0=(s1+s2+16)/(2*MB_BLOCK_SIZE);             // no edge
  if (!mb_available_up && mb_available_left)
    s0=(s2+8)/MB_BLOCK_SIZE;                     // upper edge
  if (mb_available_up && !mb_available_left)
    s0=(s1+8)/MB_BLOCK_SIZE;                     // left edge
  if (!mb_available_up && !mb_available_left)
    s0=128;                                      // top left corner, nothing to predict from

  for (i=0; i < MB_BLOCK_SIZE; i++)
  {
    // vertical prediction
    if (mb_available_up)
      s[i][0]=imgY_pred[pix_y-1][img->pix_x+i];
    // horizontal prediction
    if (mb_available_left)
      s[i][1]=imgY_pred[pix_y+i][img->pix_x-1];
  }

  for (j=0; j < MB_BLOCK_SIZE; j++)
  {
    for (i=0; i < MB_BLOCK_SIZE; i++)
    {
      img->mprr_2[VERT_PRED_16][j][i]=s[i][0]; // store vertical prediction
      img->mprr_2[HOR_PRED_16 ][j][i]=s[j][1]; // store horizontal prediction
      img->mprr_2[DC_PRED_16  ][j][i]=s0;      // store DC prediction
    }
  }
  if (!mb_available_up || !mb_available_left) // edge
    return;

  // 16 bit integer plan pred

  ih=0;
  iv=0;
  for (i=1;i<9;i++)
  {
    ih += i*(imgY_pred[pix_y-1][img->pix_x+7+i] - imgY_pred[pix_y-1][img->pix_x+7-i]);
    iv += i*(imgY_pred[pix_y+7+i][img->pix_x-1] - imgY_pred[pix_y+7-i][img->pix_x-1]);
  }
  ib=(5*ih+32)>>6;
  ic=(5*iv+32)>>6;

  iaa=16*(imgY_pred[pix_y-1][img->pix_x+15]+imgY_pred[pix_y+15][img->pix_x-1]);
  for (j=0;j< MB_BLOCK_SIZE;j++)
  {
    for (i=0;i< MB_BLOCK_SIZE;i++)
    {
      img->mprr_2[PLANE_16][j][i]=max(0,min(255,(iaa+(i-7)*ib +(j-7)*ic + 16)/32));// store plane prediction
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    For new intra pred routines
 *
 * \para Input:
 *    Image par, 16x16 based intra mode
 *
 * \para Output:
 *    none
 ************************************************************************
 */
int dct_luma2(int new_intra_mode)
{
  int qp_const;
  int i,j;
  int ii,jj;
  int i1,j1;
  int M1[16][16];
  int M4[4][4];
  int M5[4],M6[4];
  int M0[4][4][4][4];
  int run,scan_pos,coeff_ctr,level;
  int qp_per,qp_rem,q_bits;
  int ac_coef = 0;
  int incr=1;
  int offset=0;

  int   b8, b4;
  int*  DCLevel = img->cofDC[0][0];
  int*  DCRun   = img->cofDC[0][1];
  int*  ACLevel;
  int*  ACRun;

  qp_per    = (img->qp-MIN_QP)/6;
  qp_rem    = (img->qp-MIN_QP)%6;
  q_bits    = Q_BITS+qp_per;
  qp_const  = (1<<q_bits)/3;


  if(input->InterlaceCodingOption>=MB_CODING && img->field_mode && mb_adaptive)
  {
    if(img->top_field)
    {
      offset = 0;
      incr   = 2;
    }
    else
    {
      offset = -15;
      incr   = 2;
    }
  }

  for (j=0;j<16;j++)
  {
    for (i=0;i<16;i++)
    {
      M1[i][j]=imgY_org[img->pix_y+(incr*j)+offset][img->pix_x+i]-img->mprr_2[new_intra_mode][j][i];
      M0[i%4][i/4][j%4][j/4]=M1[i][j];
    }
  }

  for (jj=0;jj<4;jj++)
  {
    for (ii=0;ii<4;ii++)
    {
      for (j=0;j<4;j++)
      {
        for (i=0;i<2;i++)
        {
          i1=3-i;
          M5[i]=  M0[i][ii][j][jj]+M0[i1][ii][j][jj];
          M5[i1]= M0[i][ii][j][jj]-M0[i1][ii][j][jj];
        }
        M0[0][ii][j][jj]=M5[0]+M5[1];
        M0[2][ii][j][jj]=M5[0]-M5[1];
        M0[1][ii][j][jj]=M5[3]*2+M5[2];
        M0[3][ii][j][jj]=M5[3]-M5[2]*2;
      }
      // vertical
      for (i=0;i<4;i++)
      {
        for (j=0;j<2;j++)
        {
          j1=3-j;
          M5[j] = M0[i][ii][j][jj]+M0[i][ii][j1][jj];
          M5[j1]= M0[i][ii][j][jj]-M0[i][ii][j1][jj];
        }
        M0[i][ii][0][jj]=M5[0]+M5[1];
        M0[i][ii][2][jj]=M5[0]-M5[1];
        M0[i][ii][1][jj]=M5[3]*2+M5[2];
        M0[i][ii][3][jj]=M5[3] -M5[2]*2;
      }
    }
  }

  // pick out DC coeff

  for (j=0;j<4;j++)
    for (i=0;i<4;i++)
      M4[i][j]= M0[0][i][0][j];

  for (j=0;j<4;j++)
  {
    for (i=0;i<2;i++)
    {
      i1=3-i;
      M5[i]= M4[i][j]+M4[i1][j];
      M5[i1]=M4[i][j]-M4[i1][j];
    }
    M4[0][j]=M5[0]+M5[1];
    M4[2][j]=M5[0]-M5[1];
    M4[1][j]=M5[3]+M5[2];
    M4[3][j]=M5[3]-M5[2];
  }

  // vertical

  for (i=0;i<4;i++)
  {
    for (j=0;j<2;j++)
    {
      j1=3-j;
      M5[j]= M4[i][j]+M4[i][j1];
      M5[j1]=M4[i][j]-M4[i][j1];
    }
    M4[i][0]=(M5[0]+M5[1])>>1;
    M4[i][2]=(M5[0]-M5[1])>>1;
    M4[i][1]=(M5[3]+M5[2])>>1;
    M4[i][3]=(M5[3]-M5[2])>>1;
  }

  // quant

  run=-1;
  scan_pos=0;

  for (coeff_ctr=0;coeff_ctr<16;coeff_ctr++)
  {
    i=SNGL_SCAN[coeff_ctr][0];
    j=SNGL_SCAN[coeff_ctr][1];

    run++;

    level= (abs(M4[i][j]) * quant_coef[qp_rem][0][0] + 2*qp_const)>>(q_bits+1);

    if (level != 0)
    {
      DCLevel[scan_pos] = sign(level,M4[i][j]);
      DCRun  [scan_pos] = run;
      ++scan_pos;
      run=-1;
    }
    M4[i][j]=sign(level,M4[i][j]);
  }
  DCLevel[scan_pos]=0;

  // invers DC transform
  for (j=0;j<4;j++)
  {
    for (i=0;i<4;i++)
      M5[i]=M4[i][j];

    M6[0]=M5[0]+M5[2];
    M6[1]=M5[0]-M5[2];
    M6[2]=M5[1]-M5[3];
    M6[3]=M5[1]+M5[3];

    for (i=0;i<2;i++)
    {
      i1=3-i;
      M4[i][j]= M6[i]+M6[i1];
      M4[i1][j]=M6[i]-M6[i1];
    }
  }

  for (i=0;i<4;i++)
  {
    for (j=0;j<4;j++)
      M5[j]=M4[i][j];

    M6[0]=M5[0]+M5[2];
    M6[1]=M5[0]-M5[2];
    M6[2]=M5[1]-M5[3];
    M6[3]=M5[1]+M5[3];

    for (j=0;j<2;j++)
    {
      j1=3-j;
      M0[0][i][0][j] = (((M6[j]+M6[j1])*dequant_coef[qp_rem][0][0]<<qp_per)+2)>>2;
      M0[0][i][0][j1]= (((M6[j]-M6[j1])*dequant_coef[qp_rem][0][0]<<qp_per)+2)>>2;
    }
  }

  // AC invers trans/quant for MB
  for (jj=0;jj<4;jj++)
  {
    for (ii=0;ii<4;ii++)
    {
      run      = -1;
      scan_pos =  0;
      b8       = 2*(jj/2) + (ii/2);
      b4       = 2*(jj%2) + (ii%2);
      ACLevel  = img->cofAC [b8][b4][0];
      ACRun    = img->cofAC [b8][b4][1];

      for (coeff_ctr=1;coeff_ctr<16;coeff_ctr++) // set in AC coeff
      {
        i=SNGL_SCAN[coeff_ctr][0];
        j=SNGL_SCAN[coeff_ctr][1];
        run++;

        level= ( abs( M0[i][ii][j][jj]) * quant_coef[qp_rem][i][j] + qp_const) >> q_bits;

        if (level != 0)
        {
          ac_coef = 15;
          ACLevel[scan_pos] = sign(level,M0[i][ii][j][jj]);
          ACRun  [scan_pos] = run;
          ++scan_pos;
          run=-1;
        }
        M0[i][ii][j][jj]=sign(level*dequant_coef[qp_rem][i][j]<<qp_per,M0[i][ii][j][jj]);
      }
      ACLevel[scan_pos] = 0;


      // IDCT horizontal

      for (j=0;j<4;j++)
      {
        for (i=0;i<4;i++)
        {
          M5[i]=M0[i][ii][j][jj];
        }

        M6[0]= M5[0]+M5[2];
        M6[1]= M5[0]-M5[2];
        M6[2]=(M5[1]>>1) -M5[3];
        M6[3]= M5[1]+(M5[3]>>1);

        for (i=0;i<2;i++)
        {
          i1=3-i;
          M0[i][ii][j][jj] =M6[i]+M6[i1];
          M0[i1][ii][j][jj]=M6[i]-M6[i1];
        }
      }

      // vert
      for (i=0;i<4;i++)
      {
        for (j=0;j<4;j++)
          M5[j]=M0[i][ii][j][jj];

        M6[0]= M5[0]+M5[2];
        M6[1]= M5[0]-M5[2];
        M6[2]=(M5[1]>>1) -M5[3];
        M6[3]= M5[1]+(M5[3]>>1);

        for (j=0;j<2;j++)
        {
          j1=3-j;
          M0[i][ii][ j][jj]=M6[j]+M6[j1];
          M0[i][ii][j1][jj]=M6[j]-M6[j1];

        }
      }

    }
  }

  for (j=0;j<16;j++)
  {
    for (i=0;i<16;i++)
    {
      M1[i][j]=M0[i%4][i/4][j%4][j/4];
    }
  }

  for (j=0;j<16;j++)
    for (i=0;i<16;i++)
      imgY[img->pix_y+j][img->pix_x+i]=(byte)min(255,max(0,(M1[i][j]+(img->mprr_2[new_intra_mode][j][i]<<DQ_BITS)+DQ_ROUND)>>DQ_BITS));

    return ac_coef;
}


/*!
 ************************************************************************
 * \brief
 *    The routine performs transform,quantization,inverse transform, adds the diff.
 *    to the prediction and writes the result to the decoded luma frame. Includes the
 *    RD constrained quantization also.
 *
 * \para Input:
 *    block_x,block_y: Block position inside a macro block (0,4,8,12).
 *
 * \para Output_
 *    nonzero: 0 if no levels are nonzero.  1 if there are nonzero levels.             \n
 *    coeff_cost: Counter for nonzero coefficients, used to discard expencive levels.
 ************************************************************************
 */
int dct_luma(int block_x,int block_y,int *coeff_cost, int old_intra_mode)
{
  int sign(int a,int b);

  int i,j,i1,j1,ilev,m5[4],m6[4],coeff_ctr;
  int qp_const,level,scan_pos,run;
  int nonzero;
  int qp_per,qp_rem,q_bits;

  int   pos_x   = block_x/BLOCK_SIZE;
  int   pos_y   = block_y/BLOCK_SIZE;
  int   b8      = 2*(pos_y/2) + (pos_x/2);
  int   b4      = 2*(pos_y%2) + (pos_x%2);
  int*  ACLevel = img->cofAC[b8][b4][0];
  int*  ACRun   = img->cofAC[b8][b4][1];

  qp_per    = (img->qp-MIN_QP)/6;
  qp_rem    = (img->qp-MIN_QP)%6;
  q_bits    = Q_BITS+qp_per;

  if (img->type == INTRA_IMG)
    qp_const=(1<<q_bits)/3;    // intra
  else
    qp_const=(1<<q_bits)/6;    // inter

  //  Horizontal transform
  for (j=0; j < BLOCK_SIZE; j++)
  {
    for (i=0; i < 2; i++)
    {
      i1=3-i;
      m5[i]=img->m7[i][j]+img->m7[i1][j];
      m5[i1]=img->m7[i][j]-img->m7[i1][j];
    }
    img->m7[0][j]=(m5[0]+m5[1]);
    img->m7[2][j]=(m5[0]-m5[1]);
    img->m7[1][j]=m5[3]*2+m5[2];
    img->m7[3][j]=m5[3]-m5[2]*2;
  }

  //  Vertival transform
  for (i=0; i < BLOCK_SIZE; i++)
  {
    for (j=0; j < 2; j++)
    {
      j1=3-j;
      m5[j]=img->m7[i][j]+img->m7[i][j1];
      m5[j1]=img->m7[i][j]-img->m7[i][j1];
    }
    img->m7[i][0]=(m5[0]+m5[1]);
    img->m7[i][2]=(m5[0]-m5[1]);
    img->m7[i][1]=m5[3]*2+m5[2];
    img->m7[i][3]=m5[3]-m5[2]*2;
  }

  // Quant

  nonzero=FALSE;

  run=-1;
  scan_pos=0;
  
  for (coeff_ctr=0;coeff_ctr < 16;coeff_ctr++)
  {
    i=SNGL_SCAN[coeff_ctr][0];
    j=SNGL_SCAN[coeff_ctr][1];
    
    run++;
    ilev=0;
    
    level = (abs (img->m7[i][j]) * quant_coef[qp_rem][i][j] + qp_const) >> q_bits;
    
    if (level != 0)
    {
      nonzero=TRUE;
      if (level > 1)
        *coeff_cost += MAX_VALUE;                // set high cost, shall not be discarded
      else
        *coeff_cost += COEFF_COST[run];
      ACLevel[scan_pos] = sign(level,img->m7[i][j]);
      ACRun  [scan_pos] = run;
      ++scan_pos;
      run=-1;                     // reset zero level counter
      ilev=level*dequant_coef[qp_rem][i][j]<<qp_per;
    }
    img->m7[i][j]=sign(ilev,img->m7[i][j]);
  }
  ACLevel[scan_pos] = 0;
  
  
  //     IDCT.
  //     horizontal

  for (j=0; j < BLOCK_SIZE; j++)
  {
    for (i=0; i < BLOCK_SIZE; i++)
    {
      m5[i]=img->m7[i][j];
    }
    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (i=0; i < 2; i++)
    {
      i1=3-i;
      img->m7[i][j]=m6[i]+m6[i1];
      img->m7[i1][j]=m6[i]-m6[i1];
    }
  }

  //  vertical

  for (i=0; i < BLOCK_SIZE; i++)
  {
    for (j=0; j < BLOCK_SIZE; j++)
    {
      m5[j]=img->m7[i][j];
    }
    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (j=0; j < 2; j++)
    {
      j1=3-j;
      img->m7[i][j] =min(255,max(0,(m6[j]+m6[j1]+(img->mpr[i+block_x][j+block_y] <<DQ_BITS)+DQ_ROUND)>>DQ_BITS));
      img->m7[i][j1]=min(255,max(0,(m6[j]-m6[j1]+(img->mpr[i+block_x][j1+block_y]<<DQ_BITS)+DQ_ROUND)>>DQ_BITS));
    }
  }

  //  Decoded block moved to frame memory

  for (j=0; j < BLOCK_SIZE; j++)
    for (i=0; i < BLOCK_SIZE; i++)
      imgY[img->pix_y+block_y+j][img->pix_x+block_x+i]=img->m7[i][j];


  return nonzero;
}



/*!
 ************************************************************************
 * \brief
 *    Transform,quantization,inverse transform for chroma.
 *    The main reason why this is done in a separate routine is the
 *    additional 2x2 transform of DC-coeffs. This routine is called
 *    ones for each of the chroma components.
 *
 * \para Input:
 *    uv    : Make difference between the U and V chroma component  \n
 *    cr_cbp: chroma coded block pattern
 *
 * \para Output:
 *    cr_cbp: Updated chroma coded block pattern.
 ************************************************************************
 */
int dct_chroma(int uv,int cr_cbp)
{
  int i,j,i1,j2,ilev,n2,n1,j1,mb_y,coeff_ctr,qp_const,level ,scan_pos,run;
  int m1[BLOCK_SIZE],m5[BLOCK_SIZE],m6[BLOCK_SIZE];
  int coeff_cost;
  int cr_cbp_tmp;
  int nn0,nn1;
  int DCcoded=0 ;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  int qp_per,qp_rem,q_bits;

  int   b4;
  int*  DCLevel = img->cofDC[uv+1][0];
  int*  DCRun   = img->cofDC[uv+1][1];
  int*  ACLevel;
  int*  ACRun;

  qp_per    = QP_SCALE_CR[img->qp-MIN_QP]/6;
  qp_rem    = QP_SCALE_CR[img->qp-MIN_QP]%6;
  q_bits    = Q_BITS+qp_per;

  if (img->type == INTRA_IMG)
    qp_const=(1<<q_bits)/3;    // intra
  else
    qp_const=(1<<q_bits)/6;    // inter

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
          m5[i]=img->m7[i+n1][mb_y]+img->m7[i1+n1][mb_y];
          m5[i1]=img->m7[i+n1][mb_y]-img->m7[i1+n1][mb_y];
        }
        img->m7[n1][mb_y]  =(m5[0]+m5[1]);
        img->m7[n1+2][mb_y]=(m5[0]-m5[1]);
        img->m7[n1+1][mb_y]=m5[3]*2+m5[2];
        img->m7[n1+3][mb_y]=m5[3]-m5[2]*2;
      }

      //  Vertical transform.

      for (i=0; i < BLOCK_SIZE; i++)
      {
        j1=n1+i;
        for (j=0; j < 2; j++)
        {
          j2=3-j;
          m5[j]=img->m7[j1][n2+j]+img->m7[j1][n2+j2];
          m5[j2]=img->m7[j1][n2+j]-img->m7[j1][n2+j2];
        }
        img->m7[j1][n2+0]=(m5[0]+m5[1]);
        img->m7[j1][n2+2]=(m5[0]-m5[1]);
        img->m7[j1][n2+1]=m5[3]*2+m5[2];
        img->m7[j1][n2+3]=m5[3]-m5[2]*2;
      }
    }
  }

  //     2X2 transform of DC coeffs.
  m1[0]=(img->m7[0][0]+img->m7[4][0]+img->m7[0][4]+img->m7[4][4]);
  m1[1]=(img->m7[0][0]-img->m7[4][0]+img->m7[0][4]-img->m7[4][4]);
  m1[2]=(img->m7[0][0]+img->m7[4][0]-img->m7[0][4]-img->m7[4][4]);
  m1[3]=(img->m7[0][0]-img->m7[4][0]-img->m7[0][4]+img->m7[4][4]);

  //     Quant of chroma 2X2 coeffs.
  run=-1;
  scan_pos=0;

  for (coeff_ctr=0; coeff_ctr < 4; coeff_ctr++)
  {
    run++;
    ilev=0;

    level =(abs(m1[coeff_ctr]) * quant_coef[qp_rem][0][0] + 2*qp_const) >> (q_bits+1);

    if (level  != 0)
    {
      currMB->cbp_blk |= 0xf0000 << (uv << 2) ;    // if one of the 2x2-DC levels is != 0 set the
      cr_cbp=max(1,cr_cbp);                     // coded-bit all 4 4x4 blocks (bit 16-19 or 20-23)
      DCcoded = 1 ;
      DCLevel[scan_pos] = sign(level ,m1[coeff_ctr]);
      DCRun  [scan_pos] = run;
      scan_pos++;
      run=-1;
      ilev=level*dequant_coef[qp_rem][0][0]<<qp_per;
    }
    m1[coeff_ctr]=sign(ilev,m1[coeff_ctr]);
  }
  DCLevel[scan_pos] = 0;

  //  Invers transform of 2x2 DC levels

  img->m7[0][0]=(m1[0]+m1[1]+m1[2]+m1[3])>>1;
  img->m7[4][0]=(m1[0]-m1[1]+m1[2]-m1[3])>>1;
  img->m7[0][4]=(m1[0]+m1[1]-m1[2]-m1[3])>>1;
  img->m7[4][4]=(m1[0]-m1[1]-m1[2]+m1[3])>>1;

  //     Quant of chroma AC-coeffs.
  coeff_cost=0;
  cr_cbp_tmp=0;

  for (n2=0; n2 <= BLOCK_SIZE; n2 += BLOCK_SIZE)
  {
    for (n1=0; n1 <= BLOCK_SIZE; n1 += BLOCK_SIZE)
    {
      b4      = 2*(n2/4) + (n1/4);
      ACLevel = img->cofAC[uv+4][b4][0];
      ACRun   = img->cofAC[uv+4][b4][1];
      run=-1;
      scan_pos=0;

      for (coeff_ctr=1; coeff_ctr < 16; coeff_ctr++)// start change rd_quant
      {
        i = SNGL_SCAN[coeff_ctr][0];
        j = SNGL_SCAN[coeff_ctr][1];
        ++run;
        ilev=0;

        level=(abs(img->m7[n1+i][n2+j])*quant_coef[qp_rem][i][j]+qp_const)>>q_bits;

        if (level  != 0)
        {
          currMB->cbp_blk |= 1 << (16 + (uv << 2) + ((n2 >> 1) + (n1 >> 2))) ;
          if (level > 1)
            coeff_cost += MAX_VALUE;                // set high cost, shall not be discarded
          else
            coeff_cost += COEFF_COST[run];

          cr_cbp_tmp=2;
          ACLevel[scan_pos] = sign(level,img->m7[n1+i][n2+j]);
          ACRun  [scan_pos] = run;
          ++scan_pos;
          run=-1;
          ilev=level*dequant_coef[qp_rem][i][j]<<qp_per;
        }
        img->m7[n1+i][n2+j]=sign(ilev,img->m7[n1+i][n2+j]); // for use in IDCT
      }
      ACLevel[scan_pos] = 0;
    }
  }

  // * reset chroma coeffs
  if(coeff_cost<7)
  {
    cr_cbp_tmp = 0 ;
    for (n2=0; n2 <= BLOCK_SIZE; n2 += BLOCK_SIZE)
    {
      for (n1=0; n1 <= BLOCK_SIZE; n1 += BLOCK_SIZE)
      {
        b4      = 2*(n2/4) + (n1/4);
        ACLevel = img->cofAC[uv+4][b4][0];
        ACRun   = img->cofAC[uv+4][b4][1];
        if( DCcoded == 0) currMB->cbp_blk &= ~(0xf0000 << (uv << 2));  // if no chroma DC's: then reset coded-bits of this chroma subblock
        nn0 = (n1>>2) + (uv<<1);
        nn1 = 4 + (n2>>2) ;
        ACLevel[0] = 0;
        for (coeff_ctr=1; coeff_ctr < 16; coeff_ctr++)// ac coeff
        {
          i=SNGL_SCAN[coeff_ctr][0];
          j=SNGL_SCAN[coeff_ctr][1];
          img->m7[n1+i][n2+j]=0;
          ACLevel[coeff_ctr] = 0;
        }
      }
    }
  }
  if(cr_cbp_tmp==2)
    cr_cbp = 2;
  //     IDCT.

      //     Horizontal.
  for (n2=0; n2 <= BLOCK_SIZE; n2 += BLOCK_SIZE)
  {
    for (n1=0; n1 <= BLOCK_SIZE; n1 += BLOCK_SIZE)
    {
      for (j=0; j < BLOCK_SIZE; j++)
      {
        for (i=0; i < BLOCK_SIZE; i++)
        {
          m5[i]=img->m7[n1+i][n2+j];
        }
        m6[0]=(m5[0]+m5[2]);
        m6[1]=(m5[0]-m5[2]);
        m6[2]=(m5[1]>>1)-m5[3];
        m6[3]=m5[1]+(m5[3]>>1);

        for (i=0; i < 2; i++)
        {
          i1=3-i;
          img->m7[n1+i][n2+j]=m6[i]+m6[i1];
          img->m7[n1+i1][n2+j]=m6[i]-m6[i1];
        }
      }

      //     Vertical.
      for (i=0; i < BLOCK_SIZE; i++)
      {
        for (j=0; j < BLOCK_SIZE; j++)
        {
          m5[j]=img->m7[n1+i][n2+j];
        }
        m6[0]=(m5[0]+m5[2]);
        m6[1]=(m5[0]-m5[2]);
        m6[2]=(m5[1]>>1)-m5[3];
        m6[3]=m5[1]+(m5[3]>>1);

        for (j=0; j < 2; j++)
        {
          j2=3-j;
          img->m7[n1+i][n2+j] =min(255,max(0,(m6[j]+m6[j2]+(img->mpr[n1+i][n2+j] <<DQ_BITS)+DQ_ROUND)>>DQ_BITS));
          img->m7[n1+i][n2+j2]=min(255,max(0,(m6[j]-m6[j2]+(img->mpr[n1+i][n2+j2]<<DQ_BITS)+DQ_ROUND)>>DQ_BITS));
        }
      }
    }
  }

  //  Decoded block moved to memory
  for (j=0; j < BLOCK_SIZE*2; j++)
    for (i=0; i < BLOCK_SIZE*2; i++)
      imgUV[uv][img->pix_c_y+j][img->pix_c_x+i]= img->m7[i][j];

  return cr_cbp;
}


/*!
 ************************************************************************
 * \brief
 *    The routine performs transform,quantization,inverse transform, adds the diff.
 *    to the prediction and writes the result to the decoded luma frame. Includes the
 *    RD constrained quantization also.
 *
 * \para Input:
 *    block_x,block_y: Block position inside a macro block (0,4,8,12).
 *
 * \para Output:
 *    nonzero: 0 if no levels are nonzero.  1 if there are nonzero levels.              \n
 *    coeff_cost: Counter for nonzero coefficients, used to discard expencive levels.
 *
 *
 ************************************************************************
 */
int dct_luma_sp(int block_x,int block_y,int *coeff_cost)
{
  int sign(int a,int b);

  int i,j,i1,j1,ilev,m5[4],m6[4],coeff_ctr;
  int qp_const,level,scan_pos,run;
  int nonzero;

  int predicted_block[BLOCK_SIZE][BLOCK_SIZE],c_err,qp_const2;
  int qp_per,qp_rem,q_bits;
  int qp_per_sp,qp_rem_sp,q_bits_sp;

  int   pos_x   = block_x/BLOCK_SIZE;
  int   pos_y   = block_y/BLOCK_SIZE;
  int   b8      = 2*(pos_y/2) + (pos_x/2);
  int   b4      = 2*(pos_y%2) + (pos_x%2);
  int*  ACLevel = img->cofAC[b8][b4][0];
  int*  ACRun   = img->cofAC[b8][b4][1];

  // For encoding optimization
  int c_err1, c_err2, level1, level2;
  double D_dis1, D_dis2;
  int len, info;
  double lambda_mode   = 0.85 * pow (2, img->qp/3.0) * 4; 

  qp_per    = (img->qp-MIN_QP)/6;
  qp_rem    = (img->qp-MIN_QP)%6;
  q_bits    = Q_BITS+qp_per;
  qp_per_sp    = (img->qpsp-MIN_QP)/6;
  qp_rem_sp    = (img->qpsp-MIN_QP)%6;
  q_bits_sp    = Q_BITS+qp_per_sp;

  qp_const=(1<<q_bits)/6;    // inter
  qp_const2=(1<<q_bits_sp)/2;  //sp_pred

  //  Horizontal transform
  for (j=0; j< BLOCK_SIZE; j++)
    for (i=0; i< BLOCK_SIZE; i++)
    {
      img->m7[i][j]+=img->mpr[i+block_x][j+block_y];
      predicted_block[i][j]=img->mpr[i+block_x][j+block_y];
    }

  for (j=0; j < BLOCK_SIZE; j++)
  {
    for (i=0; i < 2; i++)
    {
      i1=3-i;
      m5[i]=img->m7[i][j]+img->m7[i1][j];
      m5[i1]=img->m7[i][j]-img->m7[i1][j];
    }
    img->m7[0][j]=(m5[0]+m5[1]);
    img->m7[2][j]=(m5[0]-m5[1]);
    img->m7[1][j]=m5[3]*2+m5[2];
    img->m7[3][j]=m5[3]-m5[2]*2;
  }

  //  Vertival transform

  for (i=0; i < BLOCK_SIZE; i++)
  {
    for (j=0; j < 2; j++)
    {
      j1=3-j;
      m5[j]=img->m7[i][j]+img->m7[i][j1];
      m5[j1]=img->m7[i][j]-img->m7[i][j1];
    }
    img->m7[i][0]=(m5[0]+m5[1]);
    img->m7[i][2]=(m5[0]-m5[1]);
    img->m7[i][1]=m5[3]*2+m5[2];
    img->m7[i][3]=m5[3]-m5[2]*2;
  }

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
  nonzero=FALSE;

  run=-1;
  scan_pos=0;
  
  for (coeff_ctr=0;coeff_ctr < 16;coeff_ctr++)     // 8 times if double scan, 16 normal scan
  {
    i=SNGL_SCAN[coeff_ctr][0];
    j=SNGL_SCAN[coeff_ctr][1];
    
    run++;
    ilev=0;
    
    // decide prediction
    
    // case 1
    level1 = (abs (predicted_block[i][j]) * quant_coef[qp_rem_sp][i][j] + qp_const2) >> q_bits_sp; 
    level1 = (level1 << q_bits_sp) / quant_coef[qp_rem_sp][i][j];                 
    c_err1 = img->m7[i][j]-sign(level1, predicted_block[i][j]);                   
    level1 = (abs (c_err1) * quant_coef[qp_rem][i][j] + qp_const) >> q_bits;
    
    // case 2
    c_err2=img->m7[i][j]-predicted_block[i][j];
    level2 = (abs (c_err2) * quant_coef[qp_rem][i][j] + qp_const) >> q_bits;
    
    // select prediction
    if ((level1 != level2) && (level1 != 0) && (level2 != 0))
    {
      D_dis1 = img->m7[i][j] - ((sign(level1,c_err1)*dequant_coef[qp_rem][i][j]*A[i][j]<< qp_per) >>6) - predicted_block[i][j]; 
      levrun_linfo_inter(level1, run, &len, &info);
      D_dis1 = D_dis1*D_dis1 + lambda_mode * len;
      
      D_dis2 = img->m7[i][j] - ((sign(level2,c_err2)*dequant_coef[qp_rem][i][j]*A[i][j]<< qp_per) >>6) - predicted_block[i][j]; 
      levrun_linfo_inter(level2, run, &len, &info);
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
      nonzero=TRUE;
      if (level > 1)
        *coeff_cost += MAX_VALUE;                // set high cost, shall not be discarded
      else
        *coeff_cost += COEFF_COST[run];
      ACLevel[scan_pos] = sign(level,c_err);
      ACRun  [scan_pos] = run;
      ++scan_pos;
      run=-1;                     // reset zero level counter
      ilev=((sign(level,c_err)*dequant_coef[qp_rem][i][j]*A[i][j]<< qp_per) >>6);
    }
    ilev+=predicted_block[i][j] ; 
    img->m7[i][j] = sign((abs(ilev) * quant_coef[qp_rem_sp][i][j] + qp_const2)>> q_bits_sp, ilev) * dequant_coef[qp_rem_sp][i][j] << qp_per_sp;
  }
  ACLevel[scan_pos] = 0;
  
    
  //     IDCT.
  //     horizontal

  for (j=0; j < BLOCK_SIZE; j++)
  {
    for (i=0; i < BLOCK_SIZE; i++)
    {
      m5[i]=img->m7[i][j];
    }
    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (i=0; i < 2; i++)
    {
      i1=3-i;
      img->m7[i][j]=m6[i]+m6[i1];
      img->m7[i1][j]=m6[i]-m6[i1];
    }
  }

  //  vertical

  for (i=0; i < BLOCK_SIZE; i++)
  {
    for (j=0; j < BLOCK_SIZE; j++)
    {
      m5[j]=img->m7[i][j];
    }
    m6[0]=(m5[0]+m5[2]);
    m6[1]=(m5[0]-m5[2]);
    m6[2]=(m5[1]>>1)-m5[3];
    m6[3]=m5[1]+(m5[3]>>1);

    for (j=0; j < 2; j++)
    {
      j1=3-j;
      img->m7[i][j] =min(255,max(0,(m6[j]+m6[j1]+DQ_ROUND)>>DQ_BITS));
      img->m7[i][j1]=min(255,max(0,(m6[j]-m6[j1]+DQ_ROUND)>>DQ_BITS));
    }
  }

  //  Decoded block moved to frame memory

  for (j=0; j < BLOCK_SIZE; j++)
  for (i=0; i < BLOCK_SIZE; i++)
    imgY[img->pix_y+block_y+j][img->pix_x+block_x+i]=img->m7[i][j];

  return nonzero;
}

/*!
 ************************************************************************
 * \brief
 *    Transform,quantization,inverse transform for chroma.
 *    The main reason why this is done in a separate routine is the
 *    additional 2x2 transform of DC-coeffs. This routine is called
 *    ones for each of the chroma components.
 *
 * \para Input:
 *    uv    : Make difference between the U and V chroma component               \n
 *    cr_cbp: chroma coded block pattern
 *
 * \para Output:
 *    cr_cbp: Updated chroma coded block pattern.
 ************************************************************************
 */
int dct_chroma_sp(int uv,int cr_cbp)
{
  int i,j,i1,j2,ilev,n2,n1,j1,mb_y,coeff_ctr,qp_const,c_err,level ,scan_pos,run;
  int m1[BLOCK_SIZE],m5[BLOCK_SIZE],m6[BLOCK_SIZE];
  int coeff_cost;
  int cr_cbp_tmp;
  int predicted_chroma_block[MB_BLOCK_SIZE/2][MB_BLOCK_SIZE/2],qp_const2,mp1[BLOCK_SIZE];
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  int qp_per,qp_rem,q_bits;
  int qp_per_sp,qp_rem_sp,q_bits_sp;

  int   b4;
  int*  DCLevel = img->cofDC[uv+1][0];
  int*  DCRun   = img->cofDC[uv+1][1];
  int*  ACLevel;
  int*  ACRun;

  int c_err1, c_err2, level1, level2;
  int len, info;
  double D_dis1, D_dis2;
  double lambda_mode   = 0.85 * pow (2, img->qp/3.0) * 4; 

  qp_per    = ((img->qp<0?img->qp:QP_SCALE_CR[img->qp])-MIN_QP)/6;
  qp_rem    = ((img->qp<0?img->qp:QP_SCALE_CR[img->qp])-MIN_QP)%6;
  q_bits    = Q_BITS+qp_per;
  qp_const=(1<<q_bits)/6;    // inter
  qp_per_sp    = ((img->qpsp<0?img->qpsp:QP_SCALE_CR[img->qpsp])-MIN_QP)/6;
  qp_rem_sp    = ((img->qpsp<0?img->qpsp:QP_SCALE_CR[img->qpsp])-MIN_QP)%6;
  q_bits_sp    = Q_BITS+qp_per_sp;
  qp_const2=(1<<q_bits_sp)/2;  //sp_pred

  for (j=0; j < MB_BLOCK_SIZE/2; j++)
    for (i=0; i < MB_BLOCK_SIZE/2; i++)
    {
      img->m7[i][j]+=img->mpr[i][j];
      predicted_chroma_block[i][j]=img->mpr[i][j];
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
          m5[i]=img->m7[i+n1][mb_y]+img->m7[i1+n1][mb_y];
          m5[i1]=img->m7[i+n1][mb_y]-img->m7[i1+n1][mb_y];
        }
        img->m7[n1][mb_y]  =(m5[0]+m5[1]);
        img->m7[n1+2][mb_y]=(m5[0]-m5[1]);
        img->m7[n1+1][mb_y]=m5[3]*2+m5[2];
        img->m7[n1+3][mb_y]=m5[3]-m5[2]*2;
      }

      //  Vertical transform.

      for (i=0; i < BLOCK_SIZE; i++)
      {
        j1=n1+i;
        for (j=0; j < 2; j++)
        {
          j2=3-j;
          m5[j]=img->m7[j1][n2+j]+img->m7[j1][n2+j2];
          m5[j2]=img->m7[j1][n2+j]-img->m7[j1][n2+j2];
        }
        img->m7[j1][n2+0]=(m5[0]+m5[1]);
        img->m7[j1][n2+2]=(m5[0]-m5[1]);
        img->m7[j1][n2+1]=m5[3]*2+m5[2];
        img->m7[j1][n2+3]=m5[3]-m5[2]*2;
      }
    }
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
  m1[0]=(img->m7[0][0]+img->m7[4][0]+img->m7[0][4]+img->m7[4][4]);
  m1[1]=(img->m7[0][0]-img->m7[4][0]+img->m7[0][4]-img->m7[4][4]);
  m1[2]=(img->m7[0][0]+img->m7[4][0]-img->m7[0][4]-img->m7[4][4]);
  m1[3]=(img->m7[0][0]-img->m7[4][0]-img->m7[0][4]+img->m7[4][4]);

  //     2X2 transform of DC coeffs.
  mp1[0]=(predicted_chroma_block[0][0]+predicted_chroma_block[4][0]+predicted_chroma_block[0][4]+predicted_chroma_block[4][4]);
  mp1[1]=(predicted_chroma_block[0][0]-predicted_chroma_block[4][0]+predicted_chroma_block[0][4]-predicted_chroma_block[4][4]);
  mp1[2]=(predicted_chroma_block[0][0]+predicted_chroma_block[4][0]-predicted_chroma_block[0][4]-predicted_chroma_block[4][4]);
  mp1[3]=(predicted_chroma_block[0][0]-predicted_chroma_block[4][0]-predicted_chroma_block[0][4]+predicted_chroma_block[4][4]);

  run=-1;
  scan_pos=0;

  for (coeff_ctr=0; coeff_ctr < 4; coeff_ctr++)
  {
    run++;
    ilev=0;

  // case 1
    c_err1 = (abs (mp1[coeff_ctr]) * quant_coef[qp_rem_sp][0][0] + 2 * qp_const2) >> (q_bits_sp + 1);
    c_err1 = (c_err1 << (q_bits_sp + 1)) / quant_coef[qp_rem_sp][0][0];
    c_err1 = m1[coeff_ctr] - sign(c_err1, mp1[coeff_ctr]);
    level1 = (abs(c_err1) * quant_coef[qp_rem][0][0] + 2 * qp_const) >> (q_bits+1);

  // case 2
    c_err2 = m1[coeff_ctr] - mp1[coeff_ctr];
    level2 = (abs(c_err2) * quant_coef[qp_rem][0][0] + 2 * qp_const) >> (q_bits+1);

  if (level1 != level2 && level1 != 0 && level2 != 0)
  {
    D_dis1 = m1[coeff_ctr] - ((sign(level1,c_err1)*dequant_coef[qp_rem][0][0]*A[0][0]<< qp_per) >>5)- mp1[coeff_ctr];
    levrun_linfo_c2x2(level1, run, &len, &info);
    D_dis1 = D_dis1 * D_dis1 + lambda_mode * len;

    D_dis2 = m1[coeff_ctr] - ((sign(level2,c_err2)*dequant_coef[qp_rem][0][0]*A[0][0]<< qp_per) >>5)- mp1[coeff_ctr];
    levrun_linfo_c2x2(level2, run, &len, &info);
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

    if (level  != 0)
    {
      currMB->cbp_blk |= 0xf0000 << (uv << 2) ;  // if one of the 2x2-DC levels is != 0 the coded-bit
      cr_cbp=max(1,cr_cbp);
      DCLevel[scan_pos] = sign(level ,c_err);
      DCRun  [scan_pos] = run;
      scan_pos++;
      run=-1;
      ilev=((sign(level,c_err)*dequant_coef[qp_rem][0][0]*A[0][0]<< qp_per) >>5);
    }
    ilev+= mp1[coeff_ctr];
    m1[coeff_ctr]=sign((abs(ilev)  * quant_coef[qp_rem_sp][0][0] + 2 * qp_const2) >> (q_bits_sp+1), ilev) * dequant_coef[qp_rem_sp][0][0] << qp_per_sp;
  }
  DCLevel[scan_pos] = 0;

  //  Invers transform of 2x2 DC levels

  img->m7[0][0]=(m1[0]+m1[1]+m1[2]+m1[3])/2;
  img->m7[4][0]=(m1[0]-m1[1]+m1[2]-m1[3])/2;
  img->m7[0][4]=(m1[0]+m1[1]-m1[2]-m1[3])/2;
  img->m7[4][4]=(m1[0]-m1[1]-m1[2]+m1[3])/2;

  //     Quant of chroma AC-coeffs.
  coeff_cost=0;
  cr_cbp_tmp=0;

  for (n2=0; n2 <= BLOCK_SIZE; n2 += BLOCK_SIZE)
  {
    for (n1=0; n1 <= BLOCK_SIZE; n1 += BLOCK_SIZE)
    {
      b4      = 2*(n2/4) + (n1/4);
      ACLevel = img->cofAC[uv+4][b4][0];
      ACRun   = img->cofAC[uv+4][b4][1];

      run      = -1;
      scan_pos =  0;

      for (coeff_ctr=1; coeff_ctr < 16; coeff_ctr++)// start change rd_quant
      {
        i=SNGL_SCAN[coeff_ctr][0];
        j=SNGL_SCAN[coeff_ctr][1];
        ++run;
        ilev=0;

    // quantization on prediction
    c_err1 = (abs(predicted_chroma_block[n1+i][n2+j]) * quant_coef[qp_rem_sp][i][j] + qp_const2) >> q_bits_sp;
    c_err1 = (c_err1 << q_bits_sp) / quant_coef[qp_rem_sp][i][j];
    c_err1 = img->m7[n1+i][n2+j] - sign(c_err1, predicted_chroma_block[n1+i][n2+j]);
    level1 = (abs(c_err1) * quant_coef[qp_rem][i][j] + qp_const) >> q_bits;

    // no quantization on prediction
    c_err2 = img->m7[n1+i][n2+j] - predicted_chroma_block[n1+i][n2+j];
    level2 = (abs(c_err2) * quant_coef[qp_rem][i][j] + qp_const) >> q_bits;

    if (level1 != level2 && level1 != 0 && level2 != 0)
    {
      D_dis1 = img->m7[n1+i][n2+j] - ((sign(level1,c_err1)*dequant_coef[qp_rem][i][j]*A[i][j]<< qp_per) >>6) - predicted_chroma_block[n1+i][n2+j]; 

      levrun_linfo_inter(level1, run, &len, &info);
      D_dis1 = D_dis1 * D_dis1 + lambda_mode * len;

      D_dis2 = img->m7[n1+i][n2+j] - ((sign(level2,c_err2)*dequant_coef[qp_rem][i][j]*A[i][j]<< qp_per) >>6) - predicted_chroma_block[n1+i][n2+j]; 
      levrun_linfo_inter(level2, run, &len, &info);
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

        if (level  != 0)
        {
          currMB->cbp_blk |=  1 << (16 + (uv << 2) + ((n2 >> 1) + (n1 >> 2))) ;
          if (level > 1)
            coeff_cost += MAX_VALUE;                // set high cost, shall not be discarded
          else
            coeff_cost += COEFF_COST[run];

          cr_cbp_tmp=2;
          ACLevel[scan_pos] = sign(level,c_err);
          ACRun  [scan_pos] = run;
          ++scan_pos;
          run=-1;
          ilev=((sign(level,c_err)*dequant_coef[qp_rem_sp][i][j]*A[i][j]<< qp_per) >>6);
        }
        ilev+=predicted_chroma_block[n1+i][n2+j];
        img->m7[n1+i][n2+j] = sign((abs(ilev) * quant_coef[qp_rem_sp][i][j] + qp_const2) >> q_bits_sp,ilev) * dequant_coef[qp_rem_sp][i][j] << qp_per_sp;
      }
      ACLevel[scan_pos] = 0;
    }
  }

  // * reset chroma coeffs

  if(cr_cbp_tmp==2)
      cr_cbp=2;
  //     IDCT.

      //     Horizontal.
  for (n2=0; n2 <= BLOCK_SIZE; n2 += BLOCK_SIZE)
  {
    for (n1=0; n1 <= BLOCK_SIZE; n1 += BLOCK_SIZE)
    {
      for (j=0; j < BLOCK_SIZE; j++)
      {
        for (i=0; i < BLOCK_SIZE; i++)
        {
          m5[i]=img->m7[n1+i][n2+j];
        }
        m6[0]=(m5[0]+m5[2]);
        m6[1]=(m5[0]-m5[2]);
        m6[2]=(m5[1]>>1)-m5[3];
        m6[3]=m5[1]+(m5[3]>>1);

        for (i=0; i < 2; i++)
        {
          i1=3-i;
          img->m7[n1+i][n2+j]=m6[i]+m6[i1];
          img->m7[n1+i1][n2+j]=m6[i]-m6[i1];
        }
      }

      //     Vertical.
      for (i=0; i < BLOCK_SIZE; i++)
      {
        for (j=0; j < BLOCK_SIZE; j++)
        {
          m5[j]=img->m7[n1+i][n2+j];
        }
        m6[0]=(m5[0]+m5[2]);
        m6[1]=(m5[0]-m5[2]);
        m6[2]=(m5[1]>>1)-m5[3];
        m6[3]=m5[1]+(m5[3]>>1);

        for (j=0; j < 2; j++)
        {
          j2=3-j;
          img->m7[n1+i][n2+j] =min(255,max(0,(m6[j]+m6[j2]+DQ_ROUND)>>DQ_BITS));
          img->m7[n1+i][n2+j2]=min(255,max(0,(m6[j]-m6[j2]+DQ_ROUND)>>DQ_BITS));
        }
      }
    }
  }

  //  Decoded block moved to memory
  for (j=0; j < BLOCK_SIZE*2; j++)
    for (i=0; i < BLOCK_SIZE*2; i++)
    {
      imgUV[uv][img->pix_c_y+j][img->pix_c_x+i]= img->m7[i][j];
    }

  return cr_cbp;
}

/*!
 ************************************************************************
 * \brief
 *    The routine performs transform,quantization,inverse transform, adds the diff.
 *    to the prediction and writes the result to the decoded luma frame. Includes the
 *    RD constrained quantization also.
 *
 * \para Input:
 *    block_x,block_y: Block position inside a macro block (0,4,8,12).
 *
 * \para Output:
 *    nonzero: 0 if no levels are nonzero.  1 if there are nonzero levels.            \n
 *    coeff_cost: Counter for nonzero coefficients, used to discard expencive levels.
 ************************************************************************
 */
void copyblock_sp(int block_x,int block_y)
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
    {
      predicted_block[i][j]=img->mpr[i+block_x][j+block_y];
    }

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
      img->m7[i][j] =min(255,max(0,(m6[j]+m6[j1]+DQ_ROUND)>>DQ_BITS));
      img->m7[i][j1]=min(255,max(0,(m6[j]-m6[j1]+DQ_ROUND)>>DQ_BITS));
    }
  }

  //  Decoded block moved to frame memory

  for (j=0; j < BLOCK_SIZE; j++)
    for (i=0; i < BLOCK_SIZE; i++)
      imgY[img->pix_y+block_y+j][img->pix_x+block_x+i]=img->m7[i][j];
}
