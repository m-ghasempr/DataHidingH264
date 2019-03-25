/*!
 *************************************************************************************
 * \file intra4x4_pred.c
 *
 * \brief
 *    Functions for intra 4x4 prediction
 *
 * \author
 *      Main contributors (see contributors.h for copyright, 
 *                         address and affiliation details)
 *      - Alexis Michael Tourapis  <alexismt@ieee.org>
 *
 *************************************************************************************
 */
#include "global.h"
#include "intra4x4_pred.h"
#include "mb_access.h"
#include "image.h"

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
 *    makes and returns 4x4 DC prediction mode
 *
 * \param currMB
 *    current MB structure
 * \param pl
 *    color plane
 * \param ioff
 *    pixel offset X within MB
 * \param joff
 *    pixel offset Y within MB
 * \return
 *    DECODING_OK   decoding of intra prediction mode was successful            \n
 *
 ***********************************************************************
 */
static inline int intra4x4_dc_pred(Macroblock *currMB, 
                                   ColorPlane pl,               
                                   int ioff,
                                   int joff)
{
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;  

  int i,j;
  int s0 = 0;  
  imgpel **imgY = (pl) ? p_Vid->dec_picture->imgUV[pl - 1] : p_Vid->dec_picture->imgY;

  PixelPos pix_a[4], pix_b;

  int block_available_up;
  int block_available_left;  

  imgpel **mb_pred = currSlice->mb_pred[pl];    

  for (i=0;i<4;++i)
  {
    p_Vid->getNeighbour(currMB, ioff - 1, joff + i, p_Vid->mb_size[IS_LUMA], &pix_a[i]);
  }
  p_Vid->getNeighbour(currMB, ioff    , joff -1 , p_Vid->mb_size[IS_LUMA], &pix_b);

  if (p_Vid->active_pps->constrained_intra_pred_flag)
  {
    for (i=0, block_available_left=1; i<4;++i)
      block_available_left  &= pix_a[i].available ? p_Vid->intra_block[pix_a[i].mb_addr]: 0;
    block_available_up       = pix_b.available ? p_Vid->intra_block [pix_b.mb_addr] : 0;
  }
  else
  {
    block_available_left     = pix_a[0].available;
    block_available_up       = pix_b.available;
  }

  // form predictor pels
  if (block_available_up)
  {
    s0 += imgY[pix_b.pos_y][pix_b.pos_x + 0];
    s0 += imgY[pix_b.pos_y][pix_b.pos_x + 1];
    s0 += imgY[pix_b.pos_y][pix_b.pos_x + 2];
    s0 += imgY[pix_b.pos_y][pix_b.pos_x + 3];
  }

  if (block_available_left)
  {
    s0 += imgY[pix_a[0].pos_y][pix_a[0].pos_x];
    s0 += imgY[pix_a[1].pos_y][pix_a[1].pos_x];
    s0 += imgY[pix_a[2].pos_y][pix_a[2].pos_x];
    s0 += imgY[pix_a[3].pos_y][pix_a[3].pos_x];
  }

  if (block_available_up && block_available_left)
  {
    // no edge
    s0 = (s0 + 4)>>3;
  }
  else if (!block_available_up && block_available_left)
  {
    // upper edge
    s0 = (s0 + 2)>>2;
  }
  else if (block_available_up && !block_available_left)
  {
    // left edge
    s0 = (s0 + 2)>>2;
  }
  else //if (!block_available_up && !block_available_left)
  {
    // top left corner, nothing to predict from
    s0 = p_Vid->dc_pred_value_comp[pl];
  }

  for (j=joff; j < joff + BLOCK_SIZE; ++j)
  {
    for (i=ioff; i < ioff + BLOCK_SIZE; ++i)
    {
      // store DC prediction
      mb_pred[j][i] = (imgpel) s0;
    }
  }
  return DECODING_OK;
}

/*!
 ***********************************************************************
 * \brief
 *    makes and returns 4x4 vertical prediction mode
 *
 * \return
 *    DECODING_OK   decoding of intraprediction mode was sucessfull            \n
 *
 ***********************************************************************
 */
static inline int intra4x4_vert_pred(Macroblock *currMB,    //!< current macroblock
                                     ColorPlane pl,         //!< current image plane
                                     int ioff,              //!< pixel offset X within MB
                                     int joff)              //!< pixel offset Y within MB
{
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  
  int j;
  imgpel **imgY = (pl) ? p_Vid->dec_picture->imgUV[pl - 1] : p_Vid->dec_picture->imgY;
  
  int block_available_up;

  PixelPos pix_b;

  imgpel **mb_pred = currSlice->mb_pred[pl];

  p_Vid->getNeighbour(currMB, ioff, joff - 1 , p_Vid->mb_size[IS_LUMA], &pix_b);

  if (p_Vid->active_pps->constrained_intra_pred_flag)
  {
    block_available_up = pix_b.available ? p_Vid->intra_block [pix_b.mb_addr] : 0;
  }
  else
  {
    block_available_up = pix_b.available;
  }

  if (!block_available_up)
    printf ("warning: Intra_4x4_Vertical prediction mode not allowed at mb %d\n", (int) p_Vid->current_mb_nr);

  for(j = joff; j < joff + BLOCK_SIZE; ++j) /* store predicted 4x4 block */
    memcpy(&(mb_pred[j][ioff]), &(imgY[pix_b.pos_y][pix_b.pos_x]), BLOCK_SIZE * sizeof(imgpel));

  return DECODING_OK;
}

/*!
 ***********************************************************************
 * \brief
 *    makes and returns 4x4 horizontal prediction mode
 *
 * \param currMB
 *    current MB structure
 * \param pl
 *    color plane
 * \param ioff
 *    pixel offset X within MB
 * \param joff
 *    pixel offset Y within MB
 *
 * \return
 *    DECODING_OK   decoding of intra prediction mode was successful
 *
 ***********************************************************************
 */
static inline int intra4x4_hor_pred(Macroblock *currMB, 
                                    ColorPlane pl,               
                                    int ioff,
                                    int joff)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  Slice *currSlice = currMB->p_Slice;

  int i,j;
  imgpel **imgY = (pl) ? p_Vid->dec_picture->imgUV[pl - 1] : p_Vid->dec_picture->imgY;

  PixelPos pix_a[4];

  int block_available_left;

  imgpel *predrow, prediction, **mb_pred = currSlice->mb_pred[pl];    

  for (i=0;i<4;++i)
  {
    p_Vid->getNeighbour(currMB, ioff -1 , joff +i , p_Vid->mb_size[IS_LUMA], &pix_a[i]);
  }

  if (p_Vid->active_pps->constrained_intra_pred_flag)
  {
    for (i=0, block_available_left=1; i<4;++i)
      block_available_left  &= pix_a[i].available ? p_Vid->intra_block[pix_a[i].mb_addr]: 0;
  }
  else
  {
    block_available_left     = pix_a[0].available;
  }

  if (!block_available_left)
    printf ("warning: Intra_4x4_Horizontal prediction mode not allowed at mb %d\n",(int) p_Vid->current_mb_nr);

  for(j=0;j<BLOCK_SIZE;++j)
  {
    predrow = mb_pred[j+joff];
    prediction = imgY[pix_a[j].pos_y][pix_a[j].pos_x];
    for(i = ioff;i < ioff + BLOCK_SIZE;++i)
      predrow[i]= prediction; /* store predicted 4x4 block */
  }

  return DECODING_OK;
}

/*!
 ***********************************************************************
 * \brief
 *    makes and returns 4x4 diagonal down right prediction mode
 *
 * \return
 *    DECODING_OK   decoding of intraprediction mode was sucessfull            \n
 *
 ***********************************************************************
 */
static inline int intra4x4_diag_down_right_pred(Macroblock *currMB,    //!< current macroblock
                                                ColorPlane pl,         //!< current image plane
                                                int ioff,              //!< pixel offset X within MB
                                                int joff)              //!< pixel offset Y within MB
{
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  
  int i;
  imgpel PredPel[13];  // array of predictor pels
  imgpel **imgY = (pl) ? p_Vid->dec_picture->imgUV[pl - 1] : p_Vid->dec_picture->imgY;

  PixelPos pix_a[4];
  PixelPos pix_b, pix_d;

  int block_available_up;
  int block_available_left;
  int block_available_up_left;

  int jpos0 = joff, jpos1 = joff + 1, jpos2 = joff + 2, jpos3 = joff + 3;
  int ipos0 = ioff, ipos1 = ioff + 1, ipos2 = ioff + 2, ipos3 = ioff + 3;
  imgpel **mb_pred = currSlice->mb_pred[pl];    

  for (i=0;i<4;++i)
  {
    p_Vid->getNeighbour(currMB, ioff -1 , joff +i , p_Vid->mb_size[IS_LUMA], &pix_a[i]);
  }

  p_Vid->getNeighbour(currMB, ioff    , joff -1 , p_Vid->mb_size[IS_LUMA], &pix_b);
  p_Vid->getNeighbour(currMB, ioff -1 , joff -1 , p_Vid->mb_size[IS_LUMA], &pix_d);

  if (p_Vid->active_pps->constrained_intra_pred_flag)
  {
    for (i=0, block_available_left=1; i<4;++i)
      block_available_left  &= pix_a[i].available ? p_Vid->intra_block[pix_a[i].mb_addr]: 0;
    block_available_up       = pix_b.available ? p_Vid->intra_block [pix_b.mb_addr] : 0;
    block_available_up_left  = pix_d.available ? p_Vid->intra_block [pix_d.mb_addr] : 0;
  }
  else
  {
    block_available_left     = pix_a[0].available;
    block_available_up       = pix_b.available;
    block_available_up_left  = pix_d.available;
  }

  if ((!block_available_up)||(!block_available_left)||(!block_available_up_left))
    printf ("warning: Intra_4x4_Diagonal_Down_Right prediction mode not allowed at mb %d\n",(int) p_Vid->current_mb_nr);

  // form predictor pels
  P_A = imgY[pix_b.pos_y][pix_b.pos_x + 0];
  P_B = imgY[pix_b.pos_y][pix_b.pos_x + 1];
  P_C = imgY[pix_b.pos_y][pix_b.pos_x + 2];
  P_D = imgY[pix_b.pos_y][pix_b.pos_x + 3];

  P_I = imgY[pix_a[0].pos_y][pix_a[0].pos_x];
  P_J = imgY[pix_a[1].pos_y][pix_a[1].pos_x];
  P_K = imgY[pix_a[2].pos_y][pix_a[2].pos_x];
  P_L = imgY[pix_a[3].pos_y][pix_a[3].pos_x];

  P_X = imgY[pix_d.pos_y][pix_d.pos_x];

  mb_pred[jpos3][ipos0] = (imgpel) ((P_L + 2*P_K + P_J + 2) >> 2);
  mb_pred[jpos2][ipos0] =
  mb_pred[jpos3][ipos1] = (imgpel) ((P_K + 2*P_J + P_I + 2) >> 2);
  mb_pred[jpos1][ipos0] =
  mb_pred[jpos2][ipos1] =
  mb_pred[jpos3][ipos2] = (imgpel) ((P_J + 2*P_I + P_X + 2) >> 2);
  mb_pred[jpos0][ipos0] =
  mb_pred[jpos1][ipos1] =
  mb_pred[jpos2][ipos2] =
  mb_pred[jpos3][ipos3] = (imgpel) ((P_I + 2*P_X + P_A + 2) >> 2);
  mb_pred[jpos0][ipos1] =
  mb_pred[jpos1][ipos2] =
  mb_pred[jpos2][ipos3] = (imgpel) ((P_X + 2*P_A + P_B + 2) >> 2);
  mb_pred[jpos0][ipos2] =
  mb_pred[jpos1][ipos3] = (imgpel) ((P_A + 2*P_B + P_C + 2) >> 2);
  mb_pred[jpos0][ipos3] = (imgpel) ((P_B + 2*P_C + P_D + 2) >> 2);

  return DECODING_OK;
}

/*!
 ***********************************************************************
 * \brief
 *    makes and returns 4x4 diagonal down left prediction mode
 *
 * \return
 *    DECODING_OK   decoding of intraprediction mode was sucessfull            \n
 *
 ***********************************************************************
 */
static inline int intra4x4_diag_down_left_pred(Macroblock *currMB,    //!< current macroblock
                                               ColorPlane pl,         //!< current image plane
                                               int ioff,              //!< pixel offset X within MB
                                               int joff)              //!< pixel offset Y within MB
{
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  
  imgpel PredPel[13];  // array of predictor pels
  imgpel **imgY = (pl) ? p_Vid->dec_picture->imgUV[pl - 1] : p_Vid->dec_picture->imgY;

  PixelPos pix_b, pix_c;

  int block_available_up;
  int block_available_up_right;

  int jpos0 = joff, jpos1 = joff + 1, jpos2 = joff + 2, jpos3 = joff + 3;
  int ipos0 = ioff, ipos1 = ioff + 1, ipos2 = ioff + 2, ipos3 = ioff + 3;
  imgpel **mb_pred = currSlice->mb_pred[pl];    

  p_Vid->getNeighbour(currMB, ioff    , joff - 1, p_Vid->mb_size[IS_LUMA], &pix_b);
  p_Vid->getNeighbour(currMB, ioff + 4, joff - 1, p_Vid->mb_size[IS_LUMA], &pix_c);

  pix_c.available = pix_c.available && !((ioff==4) && ((joff==4)||(joff==12)));

  if (p_Vid->active_pps->constrained_intra_pred_flag)
  {
    block_available_up       = pix_b.available ? p_Vid->intra_block [pix_b.mb_addr] : 0;
    block_available_up_right = pix_c.available ? p_Vid->intra_block [pix_c.mb_addr] : 0;
  }
  else
  {
    block_available_up       = pix_b.available;
    block_available_up_right = pix_c.available;
  }

  if (!block_available_up)
    printf ("warning: Intra_4x4_Diagonal_Down_Left prediction mode not allowed at mb %d\n", (int) p_Vid->current_mb_nr);

  // form predictor pels
  P_A = imgY[pix_b.pos_y][pix_b.pos_x + 0];
  P_B = imgY[pix_b.pos_y][pix_b.pos_x + 1];
  P_C = imgY[pix_b.pos_y][pix_b.pos_x + 2];
  P_D = imgY[pix_b.pos_y][pix_b.pos_x + 3];

  if (block_available_up_right)
  {
    P_E = imgY[pix_c.pos_y][pix_c.pos_x + 0];
    P_F = imgY[pix_c.pos_y][pix_c.pos_x + 1];
    P_G = imgY[pix_c.pos_y][pix_c.pos_x + 2];
    P_H = imgY[pix_c.pos_y][pix_c.pos_x + 3];
  }
  else
  {
    P_E = P_F = P_G = P_H = P_D;
  }

  mb_pred[jpos0][ipos0] = (imgpel) ((P_A + P_C + 2*(P_B) + 2) >> 2);
  mb_pred[jpos0][ipos1] =
  mb_pred[jpos1][ipos0] = (imgpel) ((P_B + P_D + 2*(P_C) + 2) >> 2);
  mb_pred[jpos0][ipos2] =
  mb_pred[jpos1][ipos1] =
  mb_pred[jpos2][ipos0] = (imgpel) ((P_C + P_E + 2*(P_D) + 2) >> 2);
  mb_pred[jpos0][ipos3] =
  mb_pred[jpos1][ipos2] =
  mb_pred[jpos2][ipos1] =
  mb_pred[jpos3][ipos0] = (imgpel) ((P_D + P_F + 2*(P_E) + 2) >> 2);
  mb_pred[jpos1][ipos3] =
  mb_pred[jpos2][ipos2] =
  mb_pred[jpos3][ipos1] = (imgpel) ((P_E + P_G + 2*(P_F) + 2) >> 2);
  mb_pred[jpos2][ipos3] =
  mb_pred[jpos3][ipos2] = (imgpel) ((P_F + P_H + 2*(P_G) + 2) >> 2);
  mb_pred[jpos3][ipos3] = (imgpel) ((P_G + 3*(P_H) + 2) >> 2);

  return DECODING_OK;
}

/*!
 ***********************************************************************
 * \brief
 *    makes and returns 4x4 vertical right prediction mode
 *
 * \return
 *    DECODING_OK   decoding of intraprediction mode was sucessfull            \n
 *
 ***********************************************************************
 */
static inline int intra4x4_vert_right_pred(Macroblock *currMB,    //!< current macroblock
                                           ColorPlane pl,         //!< current image plane
                                           int ioff,              //!< pixel offset X within MB
                                           int joff)              //!< pixel offset Y within MB
{
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  
  int i;
  imgpel PredPel[13];  // array of predictor pels
  imgpel **imgY = (pl) ? p_Vid->dec_picture->imgUV[pl - 1] : p_Vid->dec_picture->imgY;

  PixelPos pix_a[4];
  PixelPos pix_b, pix_d;

  int block_available_up;
  int block_available_left;
  int block_available_up_left;

  int jpos0 = joff, jpos1 = joff + 1, jpos2 = joff + 2, jpos3 = joff + 3;
  int ipos0 = ioff, ipos1 = ioff + 1, ipos2 = ioff + 2, ipos3 = ioff + 3;
  imgpel **mb_pred = currSlice->mb_pred[pl];    

  for (i=0;i<4;++i)
  {
    p_Vid->getNeighbour(currMB, ioff -1 , joff +i , p_Vid->mb_size[IS_LUMA], &pix_a[i]);
  }

  p_Vid->getNeighbour(currMB, ioff    , joff -1 , p_Vid->mb_size[IS_LUMA], &pix_b);
  p_Vid->getNeighbour(currMB, ioff -1 , joff -1 , p_Vid->mb_size[IS_LUMA], &pix_d);

  if (p_Vid->active_pps->constrained_intra_pred_flag)
  {
    for (i=0, block_available_left=1; i<4;++i)
      block_available_left  &= pix_a[i].available ? p_Vid->intra_block[pix_a[i].mb_addr]: 0;
    block_available_up       = pix_b.available ? p_Vid->intra_block [pix_b.mb_addr] : 0;
    block_available_up_left  = pix_d.available ? p_Vid->intra_block [pix_d.mb_addr] : 0;
  }
  else
  {
    block_available_left     = pix_a[0].available;
    block_available_up       = pix_b.available;
    block_available_up_left  = pix_d.available;
  }

  if ((!block_available_up)||(!block_available_left)||(!block_available_up_left))
    printf ("warning: Intra_4x4_Vertical_Right prediction mode not allowed at mb %d\n", (int) p_Vid->current_mb_nr);

  // form predictor pels
  P_A = imgY[pix_b.pos_y][pix_b.pos_x + 0];
  P_B = imgY[pix_b.pos_y][pix_b.pos_x + 1];
  P_C = imgY[pix_b.pos_y][pix_b.pos_x + 2];
  P_D = imgY[pix_b.pos_y][pix_b.pos_x + 3];

  P_I = imgY[pix_a[0].pos_y][pix_a[0].pos_x];
  P_J = imgY[pix_a[1].pos_y][pix_a[1].pos_x];
  P_K = imgY[pix_a[2].pos_y][pix_a[2].pos_x];
  P_L = imgY[pix_a[3].pos_y][pix_a[3].pos_x];

  P_X = imgY[pix_d.pos_y][pix_d.pos_x];

  mb_pred[jpos0][ipos0] =
  mb_pred[jpos2][ipos1] = (imgpel) ((P_X + P_A + 1) >> 1);
  mb_pred[jpos0][ipos1] =
  mb_pred[jpos2][ipos2] = (imgpel) ((P_A + P_B + 1) >> 1);
  mb_pred[jpos0][ipos2] =
  mb_pred[jpos2][ipos3] = (imgpel) ((P_B + P_C + 1) >> 1);
  mb_pred[jpos0][ipos3] = (imgpel) ((P_C + P_D + 1) >> 1);
  mb_pred[jpos1][ipos0] =
  mb_pred[jpos3][ipos1] = (imgpel) ((P_I + 2*P_X + P_A + 2) >> 2);
  mb_pred[jpos1][ipos1] =
  mb_pred[jpos3][ipos2] = (imgpel) ((P_X + 2*P_A + P_B + 2) >> 2);
  mb_pred[jpos1][ipos2] =
  mb_pred[jpos3][ipos3] = (imgpel) ((P_A + 2*P_B + P_C + 2) >> 2);
  mb_pred[jpos1][ipos3] = (imgpel) ((P_B + 2*P_C + P_D + 2) >> 2);
  mb_pred[jpos2][ipos0] = (imgpel) ((P_X + 2*P_I + P_J + 2) >> 2);
  mb_pred[jpos3][ipos0] = (imgpel) ((P_I + 2*P_J + P_K + 2) >> 2);

  return DECODING_OK;
}


/*!
 ***********************************************************************
 * \brief
 *    makes and returns 4x4 vertical left prediction mode
 *
 * \return
 *    DECODING_OK   decoding of intraprediction mode was sucessfull            \n
 *
 ***********************************************************************
 */
static inline int intra4x4_vert_left_pred(Macroblock *currMB,    //!< current macroblock
                                          ColorPlane pl,         //!< current image plane
                                          int ioff,              //!< pixel offset X within MB
                                          int joff)              //!< pixel offset Y within MB
{
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  
  imgpel PredPel[13];  // array of predictor pels  
  imgpel **imgY = (pl) ? p_Vid->dec_picture->imgUV[pl - 1] : p_Vid->dec_picture->imgY;

  PixelPos pix_b, pix_c;

  int block_available_up;
  int block_available_up_right;

  int jpos0 = joff, jpos1 = joff + 1, jpos2 = joff + 2, jpos3 = joff + 3;
  int ipos0 = ioff, ipos1 = ioff + 1, ipos2 = ioff + 2, ipos3 = ioff + 3;
  imgpel **mb_pred = currSlice->mb_pred[pl];    

  p_Vid->getNeighbour(currMB, ioff    , joff -1 , p_Vid->mb_size[IS_LUMA], &pix_b);
  p_Vid->getNeighbour(currMB, ioff +4 , joff -1 , p_Vid->mb_size[IS_LUMA], &pix_c);

  pix_c.available = pix_c.available && !((ioff==4) && ((joff==4)||(joff==12)));
  
  if (p_Vid->active_pps->constrained_intra_pred_flag)
  {
    block_available_up       = pix_b.available ? p_Vid->intra_block [pix_b.mb_addr] : 0;
    block_available_up_right = pix_c.available ? p_Vid->intra_block [pix_c.mb_addr] : 0;
  }
  else
  {
    block_available_up       = pix_b.available;
    block_available_up_right = pix_c.available;
  }


  if (!block_available_up)
    printf ("warning: Intra_4x4_Vertical_Left prediction mode not allowed at mb %d\n", (int) p_Vid->current_mb_nr);

  // form predictor pels
  P_A = imgY[pix_b.pos_y][pix_b.pos_x + 0];
  P_B = imgY[pix_b.pos_y][pix_b.pos_x + 1];
  P_C = imgY[pix_b.pos_y][pix_b.pos_x + 2];
  P_D = imgY[pix_b.pos_y][pix_b.pos_x + 3];

  if (block_available_up_right)
  {
    P_E = imgY[pix_c.pos_y][pix_c.pos_x + 0];
    P_F = imgY[pix_c.pos_y][pix_c.pos_x + 1];
    P_G = imgY[pix_c.pos_y][pix_c.pos_x + 2];
    P_H = imgY[pix_c.pos_y][pix_c.pos_x + 3];
  }
  else
  {
    P_E = P_F = P_G = P_H = P_D;
  }

  mb_pred[jpos0][ipos0] = (imgpel) ((P_A + P_B + 1) >> 1);
  mb_pred[jpos0][ipos1] =
  mb_pred[jpos2][ipos0] = (imgpel) ((P_B + P_C + 1) >> 1);
  mb_pred[jpos0][ipos2] =
  mb_pred[jpos2][ipos1] = (imgpel) ((P_C + P_D + 1) >> 1);
  mb_pred[jpos0][ipos3] =
  mb_pred[jpos2][ipos2] = (imgpel) ((P_D + P_E + 1) >> 1);
  mb_pred[jpos2][ipos3] = (imgpel) ((P_E + P_F + 1) >> 1);
  mb_pred[jpos1][ipos0] = (imgpel) ((P_A + 2*P_B + P_C + 2) >> 2);
  mb_pred[jpos1][ipos1] =
  mb_pred[jpos3][ipos0] = (imgpel) ((P_B + 2*P_C + P_D + 2) >> 2);
  mb_pred[jpos1][ipos2] =
  mb_pred[jpos3][ipos1] = (imgpel) ((P_C + 2*P_D + P_E + 2) >> 2);
  mb_pred[jpos1][ipos3] =
  mb_pred[jpos3][ipos2] = (imgpel) ((P_D + 2*P_E + P_F + 2) >> 2);
  mb_pred[jpos3][ipos3] = (imgpel) ((P_E + 2*P_F + P_G + 2) >> 2);

  return DECODING_OK;
}

/*!
 ***********************************************************************
 * \brief
 *    makes and returns 4x4 horizontal up prediction mode
 *
 * \return
 *    DECODING_OK   decoding of intraprediction mode was sucessfull            \n
 *
 ***********************************************************************
 */
static inline int intra4x4_hor_up_pred(Macroblock *currMB,    //!< current macroblock
                                       ColorPlane pl,         //!< current image plane
                                       int ioff,              //!< pixel offset X within MB
                                       int joff)              //!< pixel offset Y within MB
{
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  
  int i;
  imgpel PredPel[13];  // array of predictor pels
  imgpel **imgY = (pl) ? p_Vid->dec_picture->imgUV[pl - 1] : p_Vid->dec_picture->imgY;

  PixelPos pix_a[4];

  int block_available_left;

  int jpos0 = joff, jpos1 = joff + 1, jpos2 = joff + 2, jpos3 = joff + 3;
  int ipos0 = ioff, ipos1 = ioff + 1, ipos2 = ioff + 2, ipos3 = ioff + 3;
  imgpel **mb_pred = currSlice->mb_pred[pl];    

  for (i=0;i<4;++i)
  {
    p_Vid->getNeighbour(currMB, ioff -1 , joff +i , p_Vid->mb_size[IS_LUMA], &pix_a[i]);
  }

  if (p_Vid->active_pps->constrained_intra_pred_flag)
  {
    for (i=0, block_available_left=1; i<4;++i)
      block_available_left  &= pix_a[i].available ? p_Vid->intra_block[pix_a[i].mb_addr]: 0;
  }
  else
  {
    block_available_left     = pix_a[0].available;
  }

  if (!block_available_left)
    printf ("warning: Intra_4x4_Horizontal_Up prediction mode not allowed at mb %d\n",(int) p_Vid->current_mb_nr);

  // form predictor pels
  P_I = imgY[pix_a[0].pos_y][pix_a[0].pos_x];
  P_J = imgY[pix_a[1].pos_y][pix_a[1].pos_x];
  P_K = imgY[pix_a[2].pos_y][pix_a[2].pos_x];
  P_L = imgY[pix_a[3].pos_y][pix_a[3].pos_x];

  mb_pred[jpos0][ipos0] = (imgpel) ((P_I + P_J + 1) >> 1);
  mb_pred[jpos0][ipos1] = (imgpel) ((P_I + 2*P_J + P_K + 2) >> 2);
  mb_pred[jpos0][ipos2] =
  mb_pred[jpos1][ipos0] = (imgpel) ((P_J + P_K + 1) >> 1);
  mb_pred[jpos0][ipos3] =
  mb_pred[jpos1][ipos1] = (imgpel) ((P_J + 2*P_K + P_L + 2) >> 2);
  mb_pred[jpos1][ipos2] =
  mb_pred[jpos2][ipos0] = (imgpel) ((P_K + P_L + 1) >> 1);
  mb_pred[jpos1][ipos3] =
  mb_pred[jpos2][ipos1] = (imgpel) ((P_K + 2*P_L + P_L + 2) >> 2);
  mb_pred[jpos2][ipos3] =
  mb_pred[jpos3][ipos1] =
  mb_pred[jpos3][ipos0] =
  mb_pred[jpos2][ipos2] =
  mb_pred[jpos3][ipos2] =
  mb_pred[jpos3][ipos3] = (imgpel) P_L;

  return DECODING_OK;
}

/*!
 ***********************************************************************
 * \brief
 *    makes and returns 4x4 horizontal down prediction mode
 *
 * \return
 *    DECODING_OK   decoding of intraprediction mode was sucessfull            \n
 *
 ***********************************************************************
 */
static inline int intra4x4_hor_down_pred(Macroblock *currMB,    //!< current macroblock
                                         ColorPlane pl,         //!< current image plane
                                         int ioff,              //!< pixel offset X within MB
                                         int joff)              //!< pixel offset Y within MB
{
  Slice *currSlice = currMB->p_Slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  
  int i;
  imgpel PredPel[13];  // array of predictor pels
  imgpel **imgY = (pl) ? p_Vid->dec_picture->imgUV[pl - 1] : p_Vid->dec_picture->imgY;

  PixelPos pix_a[4];
  PixelPos pix_b, pix_d;

  int block_available_up;
  int block_available_left;
  int block_available_up_left;

  int jpos0 = joff, jpos1 = joff + 1, jpos2 = joff + 2, jpos3 = joff + 3;
  int ipos0 = ioff, ipos1 = ioff + 1, ipos2 = ioff + 2, ipos3 = ioff + 3;
  imgpel **mb_pred = currSlice->mb_pred[pl];    
  
  for (i=0;i<4;++i)
  {
    p_Vid->getNeighbour(currMB, ioff -1 , joff +i , p_Vid->mb_size[IS_LUMA], &pix_a[i]);
  }

  p_Vid->getNeighbour(currMB, ioff    , joff -1 , p_Vid->mb_size[IS_LUMA], &pix_b);
  p_Vid->getNeighbour(currMB, ioff -1 , joff -1 , p_Vid->mb_size[IS_LUMA], &pix_d);

  if (p_Vid->active_pps->constrained_intra_pred_flag)
  {
    for (i=0, block_available_left=1; i<4;++i)
      block_available_left  &= pix_a[i].available ? p_Vid->intra_block[pix_a[i].mb_addr]: 0;
    block_available_up       = pix_b.available ? p_Vid->intra_block [pix_b.mb_addr] : 0;
    block_available_up_left  = pix_d.available ? p_Vid->intra_block [pix_d.mb_addr] : 0;
  }
  else
  {
    block_available_left     = pix_a[0].available;
    block_available_up       = pix_b.available;
    block_available_up_left  = pix_d.available;
  }

  if ((!block_available_up)||(!block_available_left)||(!block_available_up_left))
    printf ("warning: Intra_4x4_Horizontal_Down prediction mode not allowed at mb %d\n", (int) p_Vid->current_mb_nr);

  // form predictor pels
  P_A = imgY[pix_b.pos_y][pix_b.pos_x + 0];
  P_B = imgY[pix_b.pos_y][pix_b.pos_x + 1];
  P_C = imgY[pix_b.pos_y][pix_b.pos_x + 2];
  P_D = imgY[pix_b.pos_y][pix_b.pos_x + 3];

  P_I = imgY[pix_a[0].pos_y][pix_a[0].pos_x];
  P_J = imgY[pix_a[1].pos_y][pix_a[1].pos_x];
  P_K = imgY[pix_a[2].pos_y][pix_a[2].pos_x];
  P_L = imgY[pix_a[3].pos_y][pix_a[3].pos_x];

  P_X = imgY[pix_d.pos_y][pix_d.pos_x];

  mb_pred[jpos0][ipos0] =
  mb_pred[jpos1][ipos2] = (imgpel) ((P_X + P_I + 1) >> 1);
  mb_pred[jpos0][ipos1] =
  mb_pred[jpos1][ipos3] = (imgpel) ((P_I + 2*P_X + P_A + 2) >> 2);
  mb_pred[jpos0][ipos2] = (imgpel) ((P_X + 2*P_A + P_B + 2) >> 2);
  mb_pred[jpos0][ipos3] = (imgpel) ((P_A + 2*P_B + P_C + 2) >> 2);
  mb_pred[jpos1][ipos0] =
  mb_pred[jpos2][ipos2] = (imgpel) ((P_I + P_J + 1) >> 1);
  mb_pred[jpos1][ipos1] =
  mb_pred[jpos2][ipos3] = (imgpel) ((P_X + 2*P_I + P_J + 2) >> 2);
  mb_pred[jpos2][ipos0] =
  mb_pred[jpos3][ipos2] = (imgpel) ((P_J + P_K + 1) >> 1);
  mb_pred[jpos2][ipos1] =
  mb_pred[jpos3][ipos3] = (imgpel) ((P_I + 2*P_J + P_K + 2) >> 2);
  mb_pred[jpos3][ipos0] = (imgpel) ((P_K + P_L + 1) >> 1);
  mb_pred[jpos3][ipos1] = (imgpel) ((P_J + 2*P_K + P_L + 2) >> 2);

  return DECODING_OK;
}


/*!
 ***********************************************************************
 * \brief
 *    makes and returns 4x4 intra prediction blocks 
 *
 * \return
 *    DECODING_OK   decoding of intraprediction mode was sucessfull            \n
 *    SEARCH_SYNC   search next sync element as errors while decoding occured
 ***********************************************************************
 */
int intrapred(Macroblock *currMB,    //!< current macroblock
              ColorPlane pl,         //!< current image plane
              int ioff,              //!< pixel offset X within MB
              int joff,              //!< pixel offset Y within MB
              int img_block_x,       //!< location of block X, multiples of 4
              int img_block_y)       //!< location of block Y, multiples of 4
{
  VideoParameters *p_Vid = currMB->p_Vid;
  byte predmode = p_Vid->ipredmode[img_block_y][img_block_x];
  currMB->ipmode_DPCM = predmode; //For residual DPCM

  switch (predmode)
  {
  case DC_PRED:
    return (intra4x4_dc_pred(currMB, pl, ioff, joff));
    break;
  case VERT_PRED:
    return (intra4x4_vert_pred(currMB, pl, ioff, joff));
    break;
  case HOR_PRED:
    return (intra4x4_hor_pred(currMB, pl, ioff, joff));
    break;
  case DIAG_DOWN_RIGHT_PRED:
    return (intra4x4_diag_down_right_pred(currMB, pl, ioff, joff));
    break;
  case DIAG_DOWN_LEFT_PRED:
    return (intra4x4_diag_down_left_pred(currMB, pl, ioff, joff));
    break;
  case VERT_RIGHT_PRED:
    return (intra4x4_vert_right_pred(currMB, pl, ioff, joff));
    break;
  case VERT_LEFT_PRED:
    return (intra4x4_vert_left_pred(currMB, pl, ioff, joff));
    break;
  case HOR_UP_PRED:
    return (intra4x4_hor_up_pred(currMB, pl, ioff, joff));
    break;
  case HOR_DOWN_PRED:  
    return (intra4x4_hor_down_pred(currMB, pl, ioff, joff));
  default:
    printf("Error: illegal intra_4x4 prediction mode: %d\n", (int) predmode);
    return SEARCH_SYNC;
    break;
  }
}
