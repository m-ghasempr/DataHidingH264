/*!
 *************************************************************************************
 * \file intra8x8.c
 *
 * \brief
 *    Intra 8x8 mode functions
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Alexis Michael Tourapis         <alexismt@ieee.org>
 *************************************************************************************
 */

#include "global.h"
#include "image.h"
#include "mb_access.h"
#include "intra8x8.h"


// Notation for comments regarding prediction and predictors.
// The pels of the 4x4 block are labelled a..p. The predictor pels above
// are labelled A..H, from the left I..P, and from above left X, as follows:
//
//  Z  A  B  C  D  E  F  G  H  I  J  K  L  M   N  O  P
//  Q  a1 b1 c1 d1 e1 f1 g1 h1
//  R  a2 b2 c2 d2 e2 f2 g2 h2
//  S  a3 b3 c3 d3 e3 f3 g3 h3
//  T  a4 b4 c4 d4 e4 f4 g4 h4
//  U  a5 b5 c5 d5 e5 f5 g5 h5
//  V  a6 b6 c6 d6 e6 f6 g6 h6
//  W  a7 b7 c7 d7 e7 f7 g7 h7
//  X  a8 b8 c8 d8 e8 f8 g8 h8

// Predictor array index definitions
#define P_Z (PredPel[0])
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
#define P_Q (PredPel[17])
#define P_R (PredPel[18])
#define P_S (PredPel[19])
#define P_T (PredPel[20])
#define P_U (PredPel[21])
#define P_V (PredPel[22])
#define P_W (PredPel[23])
#define P_X (PredPel[24])

/*!
 *************************************************************************************
 * \brief
 *    Prefiltering for Intra8x8 prediction
 *************************************************************************************
 */
void LowPassForIntra8x8Pred(imgpel *PredPel, int block_up_left, int block_up, int block_left)
{
  int i;
  imgpel LoopArray[25];

  memcpy(LoopArray,PredPel, 25 * sizeof(imgpel));

  if(block_up)
  {
    if(block_up_left)
    {
      LoopArray[1] = (imgpel) ((PredPel[0] + (PredPel[1]<<1) + PredPel[2] + 2) >> 2);
    }
    else
      LoopArray[1] = (imgpel) ((PredPel[1] + (PredPel[1]<<1) + PredPel[2] + 2) >> 2);


    for(i = 2; i <16; i++)
    {
      LoopArray[i] = (imgpel) ((PredPel[i-1] + (PredPel[i]<<1) + PredPel[i+1] + 2) >> 2);
    }
    LoopArray[16] = (imgpel) ((PredPel[15] + PredPel[16] + (PredPel[16]<<1) + 2) >> 2);
  }

  if(block_up_left)
  {
    if(block_up && block_left)
    {
      LoopArray[0] = (imgpel) (((PredPel[0] << 1) + PredPel[1] + PredPel[17] + 2) >> 2);
    }
    else
    {
      if(block_up)
        LoopArray[0] = (imgpel) ((PredPel[0] + (PredPel[0] << 1) + PredPel[1] + 2) >> 2);
      else
        if(block_left)
          LoopArray[0] = (imgpel) ((PredPel[0] + (PredPel[0] << 1) + PredPel[17] + 2) >> 2);
    }
  }

  if(block_left)
  {
    if(block_up_left)
      LoopArray[17] = (imgpel) ((PredPel[0] + (PredPel[17] << 1) + PredPel[18] + 2) >> 2);
    else
      LoopArray[17] = (imgpel) ((PredPel[17] + (PredPel[17] << 1) + PredPel[18] + 2) >> 2);

    for(i = 18; i <24; i++)
    {
      LoopArray[i] = (imgpel) ((PredPel[i-1] + (PredPel[i]<<1) + PredPel[i+1] + 2) >> 2);
    }
    LoopArray[24] = (imgpel) ((PredPel[23] + (PredPel[24] << 1) + PredPel[24] + 2) >> 2);
  }

  memcpy(PredPel, LoopArray, 25 * sizeof(imgpel));
}

/*!
 ************************************************************************
 * \brief
 *    Vertical 8x8 Prediction
 ************************************************************************
 */
static inline void get_i8x8_vertical(imgpel **cur_pred, imgpel *PredPel)
{
  int j;
  for (j=0; j < BLOCK_SIZE_8x8; j++)
    memcpy(cur_pred[j], &P_A, BLOCK_SIZE_8x8 * sizeof(imgpel));
}


/*!
 ************************************************************************
 * \brief
 *    Horizontal 8x8 Prediction
 ************************************************************************
 */
static inline void get_i8x8_horizontal(imgpel **cur_pred, imgpel *PredPel)
{
  int i;
  for (i=0; i < BLOCK_SIZE_8x8; i++)
  {
    cur_pred[i][0]  =
    cur_pred[i][1]  =
    cur_pred[i][2]  =
    cur_pred[i][3]  =
    cur_pred[i][4]  =
    cur_pred[i][5]  =
    cur_pred[i][6]  =
    cur_pred[i][7]  = (imgpel) (&P_Q)[i];
  }
}

/*!
 ************************************************************************
 * \brief
 *    DC 8x8 Prediction
 ************************************************************************
 */
static inline void get_i8x8_dc(imgpel **cur_pred, imgpel *PredPel, int left_available, int up_available)
{
  int i, j, s0 = 0;
  if (up_available && left_available)
  {
    // no edge
    s0 = rshift_rnd_sf((P_A + P_B + P_C + P_D + P_E + P_F + P_G + P_H + P_Q + P_R + P_S + P_T + P_U + P_V + P_W + P_X), 4);
  }
  else if (!up_available && left_available)
  {
    // upper edge
    s0 = rshift_rnd_sf((P_Q + P_R + P_S + P_T + P_U + P_V + P_W + P_X), 3);
  }
  else if (up_available && !left_available)
  {
    // left edge
    s0 = rshift_rnd_sf((P_A + P_B + P_C + P_D + P_E + P_F + P_G + P_H), 3);
  }
  else //if (!up_available && !left_available)
  {
    // top left corner, nothing to predict from
    s0 = P_A; //p_Vid->dc_pred_value;
  }

  // First line
  for (i=0; i < BLOCK_SIZE_8x8; i++)
  {
    cur_pred[0][i] = (imgpel) s0;
  }

  for (j=0; j < BLOCK_SIZE_8x8 - 1; j++)
  {
    memcpy (cur_pred[j + 1], cur_pred[j], BLOCK_SIZE_8x8 * sizeof(imgpel));
  }
}

/*!
 ************************************************************************
 * \brief
 *    Diagonal Down Left 8x8 Prediction
 ************************************************************************
 */
static inline void get_i8x8_downleft(imgpel **cur_pred, imgpel *PredPel)
{
  cur_pred[0][0] = (imgpel) ((P_A + P_C + ((P_B) << 1) + 2) >> 2);
  cur_pred[0][1] =
  cur_pred[1][0] = (imgpel) ((P_B + P_D + ((P_C) << 1) + 2) >> 2);
  cur_pred[0][2] =
  cur_pred[1][1] =
  cur_pred[2][0] = (imgpel) ((P_C + P_E + ((P_D) << 1) + 2) >> 2);
  cur_pred[0][3] =
  cur_pred[1][2] =
  cur_pred[2][1] =
  cur_pred[3][0] = (imgpel) ((P_D + P_F + ((P_E) << 1) + 2) >> 2);
  cur_pred[0][4] =
  cur_pred[1][3] =
  cur_pred[2][2] =
  cur_pred[3][1] =
  cur_pred[4][0] = (imgpel) ((P_E + P_G + ((P_F) << 1) + 2) >> 2);
  cur_pred[0][5] =
  cur_pred[1][4] =
  cur_pred[2][3] =
  cur_pred[3][2] =
  cur_pred[4][1] =
  cur_pred[5][0] = (imgpel) ((P_F + P_H + ((P_G) << 1) + 2) >> 2);
  cur_pred[0][6] =
  cur_pred[1][5] =
  cur_pred[2][4] =
  cur_pred[3][3] =
  cur_pred[4][2] =
  cur_pred[5][1] =
  cur_pred[6][0] = (imgpel) ((P_G + P_I + ((P_H) << 1) + 2) >> 2);
  cur_pred[0][7] =
  cur_pred[1][6] =
  cur_pred[2][5] =
  cur_pred[3][4] =
  cur_pred[4][3] =
  cur_pred[5][2] =
  cur_pred[6][1] =
  cur_pred[7][0] = (imgpel) ((P_H + P_J + ((P_I) << 1) + 2) >> 2);
  cur_pred[1][7] =
  cur_pred[2][6] =
  cur_pred[3][5] =
  cur_pred[4][4] =
  cur_pred[5][3] =
  cur_pred[6][2] =
  cur_pred[7][1] = (imgpel) ((P_I + P_K + ((P_J) << 1) + 2) >> 2);
  cur_pred[2][7] =
  cur_pred[3][6] =
  cur_pred[4][5] =
  cur_pred[5][4] =
  cur_pred[6][3] =
  cur_pred[7][2] = (imgpel) ((P_J + P_L + (P_K << 1) + 2) >> 2);
  cur_pred[3][7] =
  cur_pred[4][6] =
  cur_pred[5][5] =
  cur_pred[6][4] =
  cur_pred[7][3] = (imgpel) ((P_K + P_M + (P_L << 1) + 2) >> 2);
  cur_pred[4][7] =
  cur_pred[5][6] =
  cur_pred[6][5] =
  cur_pred[7][4] = (imgpel) ((P_L + P_N + (P_M << 1) + 2) >> 2);
  cur_pred[5][7] =
  cur_pred[6][6] =
  cur_pred[7][5] = (imgpel) ((P_M + P_O + (P_N << 1) + 2) >> 2);
  cur_pred[6][7] =
  cur_pred[7][6] = (imgpel) ((P_N + P_P + (P_O << 1) + 2) >> 2);
  cur_pred[7][7] = (imgpel) ((P_O + P_P + (P_P << 1) + 2) >> 2);
}

/*!
 ************************************************************************
 * \brief
 *    Diagonal Down Right 8x8 Prediction
 ************************************************************************
 */
static inline void get_i8x8_downright(imgpel **cur_pred, imgpel *PredPel)
{
    cur_pred[7][0] = (imgpel) ((P_X + P_V + (P_W << 1) + 2) >> 2);
    cur_pred[6][0] =
    cur_pred[7][1] = (imgpel) ((P_W + P_U + (P_V << 1) + 2) >> 2);
    cur_pred[5][0] =
    cur_pred[6][1] =
    cur_pred[7][2] = (imgpel) ((P_V + P_T + (P_U << 1) + 2) >> 2);
    cur_pred[4][0] =
    cur_pred[5][1] =
    cur_pred[6][2] =
    cur_pred[7][3] = (imgpel) ((P_U + P_S + (P_T << 1) + 2) >> 2);
    cur_pred[3][0] =
    cur_pred[4][1] =
    cur_pred[5][2] =
    cur_pred[6][3] =
    cur_pred[7][4] = (imgpel) ((P_T + P_R + (P_S << 1) + 2) >> 2);
    cur_pred[2][0] =
    cur_pred[3][1] =
    cur_pred[4][2] =
    cur_pred[5][3] =
    cur_pred[6][4] =
    cur_pred[7][5] = (imgpel) ((P_S + P_Q + (P_R << 1) + 2) >> 2);
    cur_pred[1][0] =
    cur_pred[2][1] =
    cur_pred[3][2] =
    cur_pred[4][3] =
    cur_pred[5][4] =
    cur_pred[6][5] =
    cur_pred[7][6] = (imgpel) ((P_R + P_Z + (P_Q << 1) + 2) >> 2);
    cur_pred[0][0] =
    cur_pred[1][1] =
    cur_pred[2][2] =
    cur_pred[3][3] =
    cur_pred[4][4] =
    cur_pred[5][5] =
    cur_pred[6][6] =
    cur_pred[7][7] = (imgpel) ((P_Q + P_A + (P_Z << 1) + 2) >> 2);
    cur_pred[0][1] =
    cur_pred[1][2] =
    cur_pred[2][3] =
    cur_pred[3][4] =
    cur_pred[4][5] =
    cur_pred[5][6] =
    cur_pred[6][7] = (imgpel) ((P_Z + P_B + (P_A << 1) + 2) >> 2);
    cur_pred[0][2] =
    cur_pred[1][3] =
    cur_pred[2][4] =
    cur_pred[3][5] =
    cur_pred[4][6] =
    cur_pred[5][7] = (imgpel) ((P_A + P_C + (P_B << 1) + 2) >> 2);
    cur_pred[0][3] =
    cur_pred[1][4] =
    cur_pred[2][5] =
    cur_pred[3][6] =
    cur_pred[4][7] = (imgpel) ((P_B + P_D + (P_C << 1) + 2) >> 2);
    cur_pred[0][4] =
    cur_pred[1][5] =
    cur_pred[2][6] =
    cur_pred[3][7] = (imgpel) ((P_C + P_E + (P_D << 1) + 2) >> 2);
    cur_pred[0][5] =
    cur_pred[1][6] =
    cur_pred[2][7] = (imgpel) ((P_D + P_F + (P_E << 1) + 2) >> 2);
    cur_pred[0][6] =
    cur_pred[1][7] = (imgpel) ((P_E + P_G + (P_F << 1) + 2) >> 2);
    cur_pred[0][7] = (imgpel) ((P_F + P_H + (P_G << 1) + 2) >> 2);
}


/*!
 ************************************************************************
 * \brief
 *    Vertical Left 8x8 Prediction
 ************************************************************************
 */
static inline void get_i8x8_vertleft(imgpel **cur_pred, imgpel *PredPel)
{
    cur_pred[0][0] = (imgpel) ((P_A + P_B + 1) >> 1);
    cur_pred[0][1] =
    cur_pred[2][0] = (imgpel) ((P_B + P_C + 1) >> 1);
    cur_pred[0][2] =
    cur_pred[2][1] =
    cur_pred[4][0] = (imgpel) ((P_C + P_D + 1) >> 1);
    cur_pred[0][3] =
    cur_pred[2][2] =
    cur_pred[4][1] =
    cur_pred[6][0] = (imgpel) ((P_D + P_E + 1) >> 1);
    cur_pred[0][4] =
    cur_pred[2][3] =
    cur_pred[4][2] =
    cur_pred[6][1] = (imgpel) ((P_E + P_F + 1) >> 1);
    cur_pred[0][5] =
    cur_pred[2][4] =
    cur_pred[4][3] =
    cur_pred[6][2] = (imgpel) ((P_F + P_G + 1) >> 1);
    cur_pred[0][6] =
    cur_pred[2][5] =
    cur_pred[4][4] =
    cur_pred[6][3] = (imgpel) ((P_G + P_H + 1) >> 1);
    cur_pred[0][7] =
    cur_pred[2][6] =
    cur_pred[4][5] =
    cur_pred[6][4] = (imgpel) ((P_H + P_I + 1) >> 1);
    cur_pred[2][7] =
    cur_pred[4][6] =
    cur_pred[6][5] = (imgpel) ((P_I + P_J + 1) >> 1);
    cur_pred[4][7] =
    cur_pred[6][6] = (imgpel) ((P_J + P_K + 1) >> 1);
    cur_pred[6][7] = (imgpel) ((P_K + P_L + 1) >> 1);
    cur_pred[1][0] = (imgpel) ((P_A + P_C + 2*P_B + 2) >> 2);
    cur_pred[1][1] =
    cur_pred[3][0] = (imgpel) ((P_B + P_D + 2*P_C + 2) >> 2);
    cur_pred[1][2] =
    cur_pred[3][1] =
    cur_pred[5][0] = (imgpel) ((P_C + P_E + 2*P_D + 2) >> 2);
    cur_pred[1][3] =
    cur_pred[3][2] =
    cur_pred[5][1] =
    cur_pred[7][0] = (imgpel) ((P_D + P_F + 2*P_E + 2) >> 2);
    cur_pred[1][4] =
    cur_pred[3][3] =
    cur_pred[5][2] =
    cur_pred[7][1] = (imgpel) ((P_E + P_G + 2*P_F + 2) >> 2);
    cur_pred[1][5] =
    cur_pred[3][4] =
    cur_pred[5][3] =
    cur_pred[7][2] = (imgpel) ((P_F + P_H + 2*P_G + 2) >> 2);
    cur_pred[1][6] =
    cur_pred[3][5] =
    cur_pred[5][4] =
    cur_pred[7][3] = (imgpel) ((P_G + P_I + 2*P_H + 2) >> 2);
    cur_pred[1][7] =
    cur_pred[3][6] =
    cur_pred[5][5] =
    cur_pred[7][4] = (imgpel) ((P_H + P_J + 2*P_I + 2) >> 2);
    cur_pred[3][7] =
    cur_pred[5][6] =
    cur_pred[7][5] = (imgpel) ((P_I + P_K + 2*P_J + 2) >> 2);
    cur_pred[5][7] =
    cur_pred[7][6] = (imgpel) ((P_J + P_L + 2*P_K + 2) >> 2);
    cur_pred[7][7] = (imgpel) ((P_K + P_M + 2*P_L + 2) >> 2);
}

/*!
 ************************************************************************
 * \brief
 *    Vertical Right 8x8 Prediction
 ************************************************************************
 */
static inline void get_i8x8_vertright(imgpel **cur_pred, imgpel *PredPel)
{
    cur_pred[0][0] =
    cur_pred[2][1] =
    cur_pred[4][2] =
    cur_pred[6][3] = (imgpel) ((P_Z + P_A + 1) >> 1);
    cur_pred[0][1] =
    cur_pred[2][2] =
    cur_pred[4][3] =
    cur_pred[6][4] = (imgpel) ((P_A + P_B + 1) >> 1);
    cur_pred[0][2] =
    cur_pred[2][3] =
    cur_pred[4][4] =
    cur_pred[6][5] = (imgpel) ((P_B + P_C + 1) >> 1);
    cur_pred[0][3] =
    cur_pred[2][4] =
    cur_pred[4][5] =
    cur_pred[6][6] = (imgpel) ((P_C + P_D + 1) >> 1);
    cur_pred[0][4] =
    cur_pred[2][5] =
    cur_pred[4][6] =
    cur_pred[6][7] = (imgpel) ((P_D + P_E + 1) >> 1);
    cur_pred[0][5] =
    cur_pred[2][6] =
    cur_pred[4][7] = (imgpel) ((P_E + P_F + 1) >> 1);
    cur_pred[0][6] =
    cur_pred[2][7] = (imgpel) ((P_F + P_G + 1) >> 1);
    cur_pred[0][7] = (imgpel) ((P_G + P_H + 1) >> 1);
    cur_pred[1][0] =
    cur_pred[3][1] =
    cur_pred[5][2] =
    cur_pred[7][3] = (imgpel) ((P_Q + P_A + (P_Z << 1) + 2) >> 2);
    cur_pred[1][1] =
    cur_pred[3][2] =
    cur_pred[5][3] =
    cur_pred[7][4] = (imgpel) ((P_Z + P_B + (P_A << 1) + 2) >> 2);
    cur_pred[1][2] =
    cur_pred[3][3] =
    cur_pred[5][4] =
    cur_pred[7][5] = (imgpel) ((P_A + P_C + (P_B << 1) + 2) >> 2);
    cur_pred[1][3] =
    cur_pred[3][4] =
    cur_pred[5][5] =
    cur_pred[7][6] = (imgpel) ((P_B + P_D + (P_C << 1) + 2) >> 2);
    cur_pred[1][4] =
    cur_pred[3][5] =
    cur_pred[5][6] =
    cur_pred[7][7] = (imgpel) ((P_C + P_E + (P_D << 1) + 2) >> 2);
    cur_pred[1][5] =
    cur_pred[3][6] =
    cur_pred[5][7] = (imgpel) ((P_D + P_F + (P_E << 1) + 2) >> 2);
    cur_pred[1][6] =
    cur_pred[3][7] = (imgpel) ((P_E + P_G + (P_F << 1) + 2) >> 2);
    cur_pred[1][7] = (imgpel) ((P_F + P_H + (P_G << 1) + 2) >> 2);
    cur_pred[2][0] =
    cur_pred[4][1] =
    cur_pred[6][2] = (imgpel) ((P_R + P_Z + (P_Q << 1) + 2) >> 2);
    cur_pred[3][0] =
    cur_pred[5][1] =
    cur_pred[7][2] = (imgpel) ((P_S + P_Q + (P_R << 1) + 2) >> 2);
    cur_pred[4][0] =
    cur_pred[6][1] = (imgpel) ((P_T + P_R + (P_S << 1) + 2) >> 2);
    cur_pred[5][0] =
    cur_pred[7][1] = (imgpel) ((P_U + P_S + (P_T << 1) + 2) >> 2);
    cur_pred[6][0] = (imgpel) ((P_V + P_T + (P_U << 1) + 2) >> 2);
    cur_pred[7][0] = (imgpel) ((P_W + P_U + (P_V << 1) + 2) >> 2);
}

/*!
 ************************************************************************
 * \brief
 *    Horizontal Down 8x8 Prediction
 ************************************************************************
 */
static inline void get_i8x8_hordown(imgpel **cur_pred, imgpel *PredPel)
{
    cur_pred[0][0] =
    cur_pred[1][2] =
    cur_pred[2][4] =
    cur_pred[3][6] = (imgpel) ((P_Q + P_Z + 1) >> 1);
    cur_pred[1][0] =
    cur_pred[2][2] =
    cur_pred[3][4] =
    cur_pred[4][6] = (imgpel) ((P_R + P_Q + 1) >> 1);
    cur_pred[2][0] =
    cur_pred[3][2] =
    cur_pred[4][4] =
    cur_pred[5][6] = (imgpel) ((P_S + P_R + 1) >> 1);
    cur_pred[3][0] =
    cur_pred[4][2] =
    cur_pred[5][4] =
    cur_pred[6][6] = (imgpel) ((P_T + P_S + 1) >> 1);
    cur_pred[4][0] =
    cur_pred[5][2] =
    cur_pred[6][4] =
    cur_pred[7][6] = (imgpel) ((P_U + P_T + 1) >> 1);
    cur_pred[5][0] =
    cur_pred[6][2] =
    cur_pred[7][4] = (imgpel) ((P_V + P_U + 1) >> 1);
    cur_pred[6][0] =
    cur_pred[7][2] = (imgpel) ((P_W + P_V + 1) >> 1);
    cur_pred[7][0] = (imgpel) ((P_X + P_W + 1) >> 1);
    cur_pred[0][1] =
    cur_pred[1][3] =
    cur_pred[2][5] =
    cur_pred[3][7] = (imgpel) ((P_Q + P_A + (P_Z << 1) + 2) >> 2);
    cur_pred[1][1] =
    cur_pred[2][3] =
    cur_pred[3][5] =
    cur_pred[4][7] = (imgpel) ((P_Z + P_R + (P_Q << 1) + 2) >> 2);
    cur_pred[2][1] =
    cur_pred[3][3] =
    cur_pred[4][5] =
    cur_pred[5][7] = (imgpel) ((P_Q + P_S + (P_R << 1) + 2) >> 2);
    cur_pred[3][1] =
    cur_pred[4][3] =
    cur_pred[5][5] =
    cur_pred[6][7] = (imgpel) ((P_R + P_T + (P_S << 1) + 2) >> 2);
    cur_pred[4][1] =
    cur_pred[5][3] =
    cur_pred[6][5] =
    cur_pred[7][7] = (imgpel) ((P_S + P_U + (P_T << 1) + 2) >> 2);
    cur_pred[5][1] =
    cur_pred[6][3] =
    cur_pred[7][5] = (imgpel) ((P_T + P_V + (P_U << 1) + 2) >> 2);
    cur_pred[6][1] =
    cur_pred[7][3] = (imgpel) ((P_U + P_W + (P_V << 1) + 2) >> 2);
    cur_pred[7][1] = (imgpel) ((P_V + P_X + (P_W << 1) + 2) >> 2);
    cur_pred[0][2] =
    cur_pred[1][4] =
    cur_pred[2][6] = (imgpel) ((P_Z + P_B + (P_A << 1) + 2) >> 2);
    cur_pred[0][3] =
    cur_pred[1][5] =
    cur_pred[2][7] = (imgpel) ((P_A + P_C + (P_B << 1) + 2) >> 2);
    cur_pred[0][4] =
    cur_pred[1][6] = (imgpel) ((P_B + P_D + (P_C << 1) + 2) >> 2);
    cur_pred[0][5] =
    cur_pred[1][7] = (imgpel) ((P_C + P_E + (P_D << 1) + 2) >> 2);
    cur_pred[0][6] = (imgpel) ((P_D + P_F + (P_E << 1) + 2) >> 2);
    cur_pred[0][7] = (imgpel) ((P_E + P_G + (P_F << 1) + 2) >> 2);
}


/*!
 ************************************************************************
 * \brief
 *    Horizontal Up 8x8 Prediction
 ************************************************************************
 */
static inline void get_i8x8_horup(imgpel **cur_pred, imgpel *PredPel)
{
  cur_pred[0][0] = (imgpel) ((P_Q + P_R + 1) >> 1);
  cur_pred[1][0] =
  cur_pred[0][2] = (imgpel) ((P_R + P_S + 1) >> 1);
  cur_pred[2][0] =
  cur_pred[1][2] =
  cur_pred[0][4] = (imgpel) ((P_S + P_T + 1) >> 1);
  cur_pred[3][0] =
  cur_pred[2][2] =
  cur_pred[1][4] =
  cur_pred[0][6] = (imgpel) ((P_T + P_U + 1) >> 1);
  cur_pred[4][0] =
  cur_pred[3][2] =
  cur_pred[2][4] =
  cur_pred[1][6] = (imgpel) ((P_U + P_V + 1) >> 1);
  cur_pred[5][0] =
  cur_pred[4][2] =
  cur_pred[3][4] =
  cur_pred[2][6] = (imgpel) ((P_V + P_W + 1) >> 1);
  cur_pred[6][0] =
  cur_pred[5][2] =
  cur_pred[4][4] =
  cur_pred[3][6] = (imgpel) ((P_W + P_X + 1) >> 1);
  cur_pred[4][6] =
  cur_pred[4][7] =
  cur_pred[5][4] =
  cur_pred[5][5] =
  cur_pred[5][6] =
  cur_pred[5][7] =
  cur_pred[6][2] =
  cur_pred[6][3] =
  cur_pred[6][4] =
  cur_pred[6][5] =
  cur_pred[6][6] =
  cur_pred[6][7] =
  cur_pred[7][0] =
  cur_pred[7][1] =
  cur_pred[7][2] =
  cur_pred[7][3] =
  cur_pred[7][4] =
  cur_pred[7][5] =
  cur_pred[7][6] =
  cur_pred[7][7] = (imgpel) P_X;
  cur_pred[6][1] =
  cur_pred[5][3] =
  cur_pred[4][5] =
  cur_pred[3][7] = (imgpel) ((P_W + P_X + (P_X << 1) + 2) >> 2);
  cur_pred[5][1] =
  cur_pred[4][3] =
  cur_pred[3][5] =
  cur_pred[2][7] = (imgpel) ((P_X + P_V + (P_W << 1) + 2) >> 2);
  cur_pred[4][1] =
  cur_pred[3][3] =
  cur_pred[2][5] =
  cur_pred[1][7] = (imgpel) ((P_W + P_U + (P_V << 1) + 2) >> 2);
  cur_pred[3][1] =
  cur_pred[2][3] =
  cur_pred[1][5] =
  cur_pred[0][7] = (imgpel) ((P_V + P_T + (P_U << 1) + 2) >> 2);
  cur_pred[2][1] =
  cur_pred[1][3] =
  cur_pred[0][5] = (imgpel) ((P_U + P_S + (P_T << 1) + 2) >> 2);
  cur_pred[1][1] =
  cur_pred[0][3] = (imgpel) ((P_T + P_R + (P_S << 1) + 2) >> 2);
  cur_pred[0][1] = (imgpel) ((P_S + P_Q + (P_R << 1) + 2) >> 2);
}


/*!
 ************************************************************************
 * \brief
 *    Make intra 8x8 prediction according to all 9 prediction modes.
 *    The routine uses left and upper neighbouring points from
 *    previous coded blocks to do this (if available). Notice that
 *    inaccessible neighbouring points are signalled with a negative
 *    value in the predmode array .
 *
 *  \par Input:
 *     Starting point of current 8x8 block image posision
 *
 *  \par Output:
 *      none
 ************************************************************************
 */
void set_intrapred_8x8(Macroblock *currMB, ColorPlane pl, int img_x,int img_y, int *left_available, int *up_available, int *all_available)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;

  int i;
  imgpel *PredPel = currMB->intra8x8_pred[pl];  // array of predictor pels
  imgpel **img_enc = p_Vid->enc_picture->p_curr_img;
  imgpel *img_pel;
  int *mb_size = p_Vid->mb_size[IS_LUMA];

  int ioff = (img_x & 15);
  int joff = (img_y & 15);

  PixelPos pix_a[8];
  PixelPos pix_b, pix_c, pix_d;

  int block_available_up;
  int block_available_left;
  int block_available_up_left;
  int block_available_up_right;

  for (i=0;i<8;i++)
  {
    p_Vid->getNeighbour(currMB, ioff - 1, joff + i , mb_size, &pix_a[i]);
  }

  p_Vid->getNeighbour(currMB, ioff    , joff - 1, mb_size, &pix_b);
  p_Vid->getNeighbour(currMB, ioff + 8, joff - 1, mb_size, &pix_c);
  p_Vid->getNeighbour(currMB, ioff - 1, joff - 1, mb_size, &pix_d);

  pix_c.available = pix_c.available &&!(ioff == 8 && joff == 8);

  if (p_Inp->UseConstrainedIntraPred)
  {
    for (i=0, block_available_left=1; i<8;i++)
      block_available_left  &= pix_a[i].available ? p_Vid->intra_block[pix_a[i].mb_addr]: 0;

    block_available_up       = pix_b.available ? p_Vid->intra_block [pix_b.mb_addr] : 0;
    block_available_up_right = pix_c.available ? p_Vid->intra_block [pix_c.mb_addr] : 0;
    block_available_up_left  = pix_d.available ? p_Vid->intra_block [pix_d.mb_addr] : 0;
  }
  else
  {
    block_available_left     = pix_a[0].available;
    block_available_up       = pix_b.available;
    block_available_up_right = pix_c.available;
    block_available_up_left  = pix_d.available;
  }

  *left_available    = block_available_left;
  *up_available      = block_available_up;
  *all_available  = block_available_up && block_available_left && block_available_up_left;

  // form predictor pels
  if (block_available_up)
  {
    img_pel = &img_enc[pix_b.pos_y][pix_b.pos_x];
    P_A = *(img_pel++);
    P_B = *(img_pel++);
    P_C = *(img_pel++);
    P_D = *(img_pel++);
    P_E = *(img_pel++);
    P_F = *(img_pel++);
    P_G = *(img_pel++);
    P_H = *(img_pel);
  }
  else
  {
    P_A = P_B = P_C = P_D = P_E = P_F = P_G = P_H = p_Vid->dc_pred_value;
  }

  if (block_available_up_right)
  {
    img_pel = &img_enc[pix_c.pos_y][pix_c.pos_x];
    P_I = *(img_pel++);
    P_J = *(img_pel++);
    P_K = *(img_pel++);
    P_L = *(img_pel++);
    P_M = *(img_pel++);
    P_N = *(img_pel++);
    P_O = *(img_pel++);
    P_P = *(img_pel);

  }
  else
  {
    P_I = P_J = P_K = P_L = P_M = P_N = P_O = P_P = P_H;
  }

  if (block_available_left)
  {
    P_Q = img_enc[pix_a[0].pos_y][pix_a[0].pos_x];
    P_R = img_enc[pix_a[1].pos_y][pix_a[1].pos_x];
    P_S = img_enc[pix_a[2].pos_y][pix_a[2].pos_x];
    P_T = img_enc[pix_a[3].pos_y][pix_a[3].pos_x];
    P_U = img_enc[pix_a[4].pos_y][pix_a[4].pos_x];
    P_V = img_enc[pix_a[5].pos_y][pix_a[5].pos_x];
    P_W = img_enc[pix_a[6].pos_y][pix_a[6].pos_x];
    P_X = img_enc[pix_a[7].pos_y][pix_a[7].pos_x];
  }
  else
  {
    P_Q = P_R = P_S = P_T = P_U = P_V = P_W = P_X = p_Vid->dc_pred_value;
  }

  if (block_available_up_left)
  {
    P_Z = img_enc[pix_d.pos_y][pix_d.pos_x];
  }
  else
  {
    P_Z = p_Vid->dc_pred_value;
  }

  LowPassForIntra8x8Pred(&(P_Z), block_available_up_left, block_available_up, block_available_left);
}

/*!
 ************************************************************************
 * \brief
 *    Generate 8x8 intra prediction block
 *
 *  \par Input:
 *     Starting point of current 8x8 block image posision
 *
 *  \par Output:
 *      none
 ************************************************************************
 */
void get_intrapred_8x8(Macroblock *currMB, ColorPlane pl, int i8x8_mode, int left_available, int up_available)
{
  imgpel *PredPel = currMB->intra8x8_pred[pl];  // array of predictor pels
  Slice *currSlice = currMB->p_slice;
  imgpel ***curr_mpr_8x8  = currSlice->mpr_8x8[pl];

  switch (i8x8_mode)
  {
  case VERT_PRED :    
    get_i8x8_vertical(curr_mpr_8x8[VERT_PRED], PredPel);
    break;
  case HOR_PRED :
    get_i8x8_horizontal(curr_mpr_8x8[HOR_PRED], PredPel);
    break;
  case DC_PRED :
    get_i8x8_dc(curr_mpr_8x8[DC_PRED], PredPel, left_available, up_available);
    break;
  case DIAG_DOWN_LEFT_PRED :
    get_i8x8_downleft(curr_mpr_8x8[DIAG_DOWN_LEFT_PRED], PredPel);
    break;
  case DIAG_DOWN_RIGHT_PRED :
    get_i8x8_downright(curr_mpr_8x8[DIAG_DOWN_RIGHT_PRED], PredPel);
    break;
  case VERT_RIGHT_PRED :
    get_i8x8_vertright(curr_mpr_8x8[VERT_RIGHT_PRED], PredPel);
    break;
  case HOR_DOWN_PRED :
    get_i8x8_hordown(curr_mpr_8x8[HOR_DOWN_PRED], PredPel);
    break;
  case VERT_LEFT_PRED :
    get_i8x8_vertleft(curr_mpr_8x8[VERT_LEFT_PRED], PredPel);
    break;
  case HOR_UP_PRED :
    get_i8x8_horup(curr_mpr_8x8[HOR_UP_PRED], PredPel);
    break;
  default:
    printf("invalid prediction mode \n");
    break;
  }
}
