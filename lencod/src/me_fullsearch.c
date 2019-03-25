
/*!
*************************************************************************************
* \file me_fullsearch.c
*
* \brief
*    Motion Estimation using Fullsearch
*
* \author
*    Main contributors (see contributors.h for copyright, address and affiliation details)
*      - Alexis Michael Tourapis <alexismt@ieee.org>
*      - Athanasios Leontaris    <aleon@dolby.com>
*
*************************************************************************************
*/

// Includes
#include "contributors.h"
#include <limits.h>

#include "global.h"
#include "image.h"
#include "memalloc.h"
#include "mb_access.h"
#include "refbuf.h"

#include "me_distortion.h"
#include "me_fullsearch.h"
#include "mv-search.h"

// Define Global Parameters
extern short*  spiral_search_x;
extern short*  spiral_search_y;
extern short*  spiral_hpel_search_x;
extern short*  spiral_hpel_search_y;

// Functions
/*!
 ***********************************************************************
 * \brief
 *    Full pixel block motion search
 ***********************************************************************
 */
int                                                //  ==> minimum motion cost after search
FullPelBlockMotionSearch (Macroblock *currMB ,     // <--  current Macroblock
                          imgpel*   orig_pic,      // <--  original pixel values for the AxB block
                          short     ref,           // <--  reference frame (0... or -1 (backward))
                          int       list,          // <--  current list
                          char   ***refPic,        // <--  reference array
                          short ****tmp_mv,        // <--  mv array
                          int       pic_pix_x,     // <--  absolute x-coordinate of regarded AxB block
                          int       pic_pix_y,     // <--  absolute y-coordinate of regarded AxB block
                          int       blocktype,     // <--  block type (1-16x16 ... 7-4x4)
                          MotionVector *pred,    // <--  motion vector predictor (x) in sub-pel units
                          MotionVector *mv,         // <--> in: search center (x) / out: motion vector (x) - in pel units
                          int       search_range,  // <--  1-d search range in pel units
                          int       min_mcost,     // <--  minimum motion cost (cost for center or huge value)
                          int       lambda_factor, // <--  lagrangian parameter for determining motion cost
                          int       apply_weights) // <--  use weight based ME
{
  int   pos, cand_x, cand_y, mcost;

  StorablePicture *ref_picture = listX[list+currMB->list_offset][ref];

  int   best_pos      = 0;                                        // position with minimum motion cost
  int   max_pos       = (2*search_range+1)*(2*search_range+1);    // number of search positions
  int   blocksize_y   = params->blc_size[blocktype][1];            // vertical block size
  int   blocksize_x   = params->blc_size[blocktype][0];            // horizontal block size

  int   pred_x        = (pic_pix_x << 2) + pred->mv_x;       // predicted position x (in sub-pel units)
  int   pred_y        = (pic_pix_y << 2) + pred->mv_y;       // predicted position y (in sub-pel units)
  int   center_x      = pic_pix_x + mv->mv_x;                        // center position x (in pel units)
  int   center_y      = pic_pix_y + mv->mv_y;                        // center position y (in pel units)
  int   check_for_00  = (blocktype==1 && !params->rdopt && img->type!=B_SLICE && ref==0);
  int   dist_method = F_PEL + 3 * apply_weights;

  ref_pic_sub.luma = ref_picture->p_curr_img_sub;

  img_width  = ref_picture->size_x;
  img_height = ref_picture->size_y;
  width_pad  = ref_picture->size_x_pad;
  height_pad = ref_picture->size_y_pad;

  if (ChromaMEEnable)
  {
    ref_pic_sub.crcb[0] = ref_picture->imgUV_sub[0];
    ref_pic_sub.crcb[1] = ref_picture->imgUV_sub[1];
    width_pad_cr  = ref_picture->size_x_cr_pad;
    height_pad_cr = ref_picture->size_y_cr_pad;
  }

  //===== loop over all search positions =====
  for (pos=0; pos<max_pos; pos++)
  {
    //--- set candidate position (absolute position in pel units) ---
    cand_x = (center_x + spiral_search_x[pos])<<2;
    cand_y = (center_y + spiral_search_y[pos])<<2;

    //--- initialize motion cost (cost for motion vector) and check ---
    mcost = MV_COST_SMP (lambda_factor, cand_x, cand_y, pred_x, pred_y);
    if (check_for_00 && cand_x==pic_pix_x && cand_y==pic_pix_y)
    {
      mcost -= WEIGHTED_COST (lambda_factor, 16);
    }
    if (mcost >= min_mcost)   continue;

    //--- add residual cost to motion cost ---
    mcost += computeUniPred[dist_method](orig_pic, blocksize_y, blocksize_x,
      min_mcost - mcost, cand_x + IMG_PAD_SIZE_TIMES4, cand_y + IMG_PAD_SIZE_TIMES4);

    //--- check if motion cost is less than minimum cost ---
    if (mcost < min_mcost)
    {
      best_pos  = pos;
      min_mcost = mcost;
    }
  }


  //===== set best motion vector and return minimum motion cost =====
  if (best_pos)
  {
    mv->mv_x += spiral_search_x[best_pos];
    mv->mv_y += spiral_search_y[best_pos];
  }
  return min_mcost;
}

/*!
 ***********************************************************************
 * \brief
 *    Full pixel block motion search
 ***********************************************************************
 */
int                                                //  ==> minimum motion cost after search
FullPelBlockMotionBiPred (Macroblock *currMB,      // <--  current Macroblock
                          imgpel*   orig_pic,      // <--  original pixel values for the AxB block
                          short     ref,           // <--  reference frame (0... or -1 (backward))
                          int       list,          // <--  reference list
                          char   ***refPic,        // <--  reference array
                          short ****tmp_mv,        // <--  mv array
                          int       pic_pix_x,     // <--  absolute x-coordinate of regarded AxB block
                          int       pic_pix_y,     // <--  absolute y-coordinate of regarded AxB block
                          int       blocktype,     // <--  block type (1-16x16 ... 7-4x4)
                          MotionVector *pred_mv1,   // <--  motion vector predictor from first list (x|y) in sub-pel units
                          MotionVector *pred_mv2,   // <--  motion vector predictor from second list (x|y) in sub-pel units
                          MotionVector *mv1,         // <--> in: search center (x|y) / out: motion vector (x|y) - in pel units
                          MotionVector *mv2,       // <--> in: search center (x|y) 
                          int       search_range,  // <--  1-d search range in pel units
                          int       min_mcost,     // <--  minimum motion cost (cost for center or huge value)
                          int       iteration_no,  // <--  bi pred iteration number
                          int       lambda_factor, // <--  lagrangian parameter for determining motion cost
                          int       apply_weights  // <--  perform weight based ME
                          ) 
{
  StorablePicture *ref_picture1 = listX[list       + currMB->list_offset][ref];
  StorablePicture *ref_picture2 = listX[(list ^ 1) + currMB->list_offset][0];
  short blocksize_y   = params->blc_size[blocktype][1];        // vertical block size
  short blocksize_x   = params->blc_size[blocktype][0];        // horizontal block size
  
  int   pos, mcost;

  int   best_pos      = 0;                                     // position with minimum motion cost
  int   max_pos       = (2*search_range+1)*(2*search_range+1); // number of search positions

  static MotionVector center1, center2, cand, pred1, pred2;
  pred1.mv_x   = (pic_pix_x << 2) + pred_mv1->mv_x; // predicted position x (in sub-pel units)
  pred1.mv_y   = (pic_pix_y << 2) + pred_mv1->mv_y; // predicted position y (in sub-pel units)
  pred2.mv_x   = (pic_pix_x << 2) + pred_mv2->mv_x; // predicted position x (in sub-pel units)
  pred2.mv_y   = (pic_pix_y << 2) + pred_mv2->mv_y; // predicted position y (in sub-pel units)
  center1.mv_x = pic_pix_x + mv1->mv_x;              // center position x (in pel units)
  center1.mv_y = pic_pix_y + mv1->mv_y;              // center position y (in pel units)
  center2.mv_x = pic_pix_x + mv2->mv_x;       // center position x of static mv (in pel units)
  center2.mv_y = pic_pix_y + mv2->mv_y;       // center position y of static mv (in pel units)

  ref_pic1_sub.luma = ref_picture1->p_curr_img_sub;
  ref_pic2_sub.luma = ref_picture2->p_curr_img_sub;

  img_width  = ref_picture1->size_x;
  img_height = ref_picture1->size_y;
  width_pad  = ref_picture1->size_x_pad;
  height_pad = ref_picture1->size_y_pad;

  if (apply_weights)
  {
    computeBiPred = computeBiPred2[F_PEL];
  }
  else
  {
    computeBiPred = computeBiPred1[F_PEL];
  }

  if ( ChromaMEEnable )
  {
    ref_pic1_sub.crcb[0] = ref_picture1->imgUV_sub[0];
    ref_pic1_sub.crcb[1] = ref_picture1->imgUV_sub[1];
    ref_pic2_sub.crcb[0] = ref_picture2->imgUV_sub[0];
    ref_pic2_sub.crcb[1] = ref_picture2->imgUV_sub[1];
    width_pad_cr  = ref_picture1->size_x_cr_pad;
    height_pad_cr = ref_picture1->size_y_cr_pad;
  }


  //===== loop over all search positions =====
  for (pos=0; pos<max_pos; pos++)
  {
    //--- set candidate position (absolute position in pel units) ---
    cand.mv_x = (center1.mv_x + spiral_search_x[pos])<<2;
    cand.mv_y = (center1.mv_y + spiral_search_y[pos])<<2;

    //--- initialize motion cost (cost for motion vector) and check ---
    mcost  = MV_COST_SMP (lambda_factor, cand.mv_x, cand.mv_y, pred1.mv_x, pred1.mv_y);
    mcost += MV_COST_SMP (lambda_factor, (center2.mv_x << 2), (center2.mv_y<<2), pred2.mv_x, pred2.mv_y);

    if (mcost >= min_mcost)   continue;

    //--- add residual cost to motion cost ---
    mcost += computeBiPred(orig_pic,
      blocksize_y, blocksize_x, min_mcost - mcost,
      cand.mv_x + IMG_PAD_SIZE_TIMES4, cand.mv_y + IMG_PAD_SIZE_TIMES4,
      (center2.mv_x << 2) + IMG_PAD_SIZE_TIMES4,
      (center2.mv_y << 2) + IMG_PAD_SIZE_TIMES4);

    //--- check if motion cost is less than minimum cost ---
    if (mcost < min_mcost)
    {
      best_pos  = pos;
      min_mcost = mcost;
    }
  }

  //===== set best motion vector and return minimum motion cost =====
  if (best_pos)
  {
    mv1->mv_x += spiral_search_x[best_pos];
    mv1->mv_y += spiral_search_y[best_pos];
  }
  return min_mcost;
}

/*!
 ***********************************************************************
 * \brief
 *    Sub pixel block motion search
 ***********************************************************************
 */
int                                               //  ==> minimum motion cost after search
SubPelBlockMotionSearch (imgpel*   orig_pic,      // <--  original pixel values for the AxB block
                         short     ref,           // <--  reference frame (0... or -1 (backward))
                         int       list,          // <--  reference picture list
                         int       list_offset,   // <--  MBAFF list offset
                         int       pic_pix_x,     // <--  absolute x-coordinate of regarded AxB block
                         int       pic_pix_y,     // <--  absolute y-coordinate of regarded AxB block
                         int       blocktype,     // <--  block type (1-16x16 ... 7-4x4)
                         MotionVector *pred,    // <--  motion vector predictor in sub-pel units
                         MotionVector *mv,         // <--> in: search center / out: motion vector - in pel units
                         int       search_pos2,   // <--  search positions for    half-pel search  (default: 9)
                         int       search_pos4,   // <--  search positions for quarter-pel search  (default: 9)
                         int       min_mcost,     // <--  minimum motion cost (cost for center or huge value)
                         int*      lambda,        // <--  lagrangian parameter for determining motion cost
                         int       apply_weights  // <--  use weight based ME
                         )
{
  int   pos, best_pos, mcost;

  int   cand_mv_x, cand_mv_y;

  int   check_position0 = (!params->rdopt && img->type!=B_SLICE && ref==0 && blocktype==1 && mv->mv_x == 0 && mv->mv_y ==0);
  int   blocksize_x     = params->blc_size[blocktype][0];
  int   blocksize_y     = params->blc_size[blocktype][1];
  int   pic4_pix_x      = ((pic_pix_x + IMG_PAD_SIZE)<< 2);
  int   pic4_pix_y      = ((pic_pix_y + IMG_PAD_SIZE)<< 2);
  int   max_pos2        = ( !start_me_refinement_hp ? imax(1,search_pos2) : search_pos2);  
  int   cmv_x, cmv_y;
  int dist_method = H_PEL + 3 * apply_weights;
  StorablePicture *ref_picture = listX[list+list_offset][ref];

  int lambda_factor = lambda[H_PEL];

  ref_pic_sub.luma = ref_picture->p_curr_img_sub;
  width_pad  = ref_picture->size_x_pad;
  height_pad = ref_picture->size_y_pad;

  if (ChromaMEEnable)
  {
    ref_pic_sub.crcb[0] = ref_picture->imgUV_sub[0];
    ref_pic_sub.crcb[1] = ref_picture->imgUV_sub[1];
    width_pad_cr  = ref_picture->size_x_cr_pad;
    height_pad_cr = ref_picture->size_y_cr_pad;
  }

  /*********************************
   *****                       *****
   *****  HALF-PEL REFINEMENT  *****
   *****                       *****
   *********************************/

  //===== loop over search positions =====
  for (best_pos = 0, pos = start_me_refinement_hp; pos < max_pos2; pos++)
  {
    cand_mv_x = mv->mv_x + (spiral_hpel_search_x[pos]);    // quarter-pel units
    cand_mv_y = mv->mv_y + (spiral_hpel_search_y[pos]);    // quarter-pel units

    //----- set motion vector cost -----
    mcost = MV_COST_SMP (lambda_factor, cand_mv_x, cand_mv_y, pred->mv_x, pred->mv_y);


    if (mcost >= min_mcost) continue;

    cmv_x = cand_mv_x + pic4_pix_x;
    cmv_y = cand_mv_y + pic4_pix_y;

    mcost += computeUniPred[dist_method]( orig_pic, blocksize_y, blocksize_x, min_mcost - mcost, cmv_x, cmv_y);

    if (pos==0 && check_position0)
    {
      mcost -= WEIGHTED_COST (lambda_factor, 16);
    }

    if (mcost < min_mcost)
    {
      min_mcost = mcost;
      best_pos  = pos;
    }
  }
  if (best_pos)
  {
    mv->mv_x += (spiral_hpel_search_x [best_pos]);
    mv->mv_y += (spiral_hpel_search_y [best_pos]);
  }

  if ( !start_me_refinement_qp )
    min_mcost = INT_MAX;

  /************************************
  *****                          *****
  *****  QUARTER-PEL REFINEMENT  *****
  *****                          *****
  ************************************/

  dist_method = Q_PEL + 3 * apply_weights;
  lambda_factor = lambda[Q_PEL];

  //===== loop over search positions =====
  for (best_pos = 0, pos = start_me_refinement_qp; pos < search_pos4; pos++)
  {
    cand_mv_x = mv->mv_x + spiral_search_x[pos];    // quarter-pel units
    cand_mv_y = mv->mv_y + spiral_search_y[pos];    // quarter-pel units

    //----- set motion vector cost -----
    mcost = MV_COST_SMP (lambda_factor, cand_mv_x, cand_mv_y, pred->mv_x, pred->mv_y);

    if (mcost >= min_mcost) continue;

    cmv_x = cand_mv_x + pic4_pix_x;
    cmv_y = cand_mv_y + pic4_pix_y;

    mcost += computeUniPred[dist_method]( orig_pic, blocksize_y, blocksize_x, min_mcost - mcost, cmv_x, cmv_y);

    if (mcost < min_mcost)
    {
      min_mcost = mcost;
      best_pos  = pos;
    }
  }
  if (best_pos)
  {
    mv->mv_x += spiral_search_x [best_pos];
    mv->mv_y += spiral_search_y [best_pos];
  }

  //===== return minimum motion cost =====
  return min_mcost;
}

/*!
***********************************************************************
* \brief
*    Bipred Sub pixel block motion search
***********************************************************************
*/
int                                               //  ==> minimum motion cost after search
SubPelBlockSearchBiPred (imgpel*   orig_pic,      // <--  original pixel values for the AxB block
                         short     ref,           // <--  reference frame (0... or -1 (backward))
                         int       list,          // <--  reference picture list
                         int       pic_pix_x,     // <--  absolute x-coordinate of regarded AxB block
                         int       pic_pix_y,     // <--  absolute y-coordinate of regarded AxB block
                         int       blocktype,     // <--  block type (1-16x16 ... 7-4x4)
                         MotionVector *pred_mv1,  // <--  motion vector predictor (x) in sub-pel units
                         MotionVector *pred_mv2,  // <--  motion vector predictor (x) in sub-pel units
                         MotionVector *mv1,       // <--> in: search center (x) / out: motion vector (x) - in pel units
                         MotionVector *mv2,       // <--> in: search center (x) / out: motion vector (x) - in pel units
                         int       search_pos2,   // <--  search positions for    half-pel search  (default: 9)
                         int       search_pos4,   // <--  search positions for quarter-pel search  (default: 9)
                         int       min_mcost,     // <--  minimum motion cost (cost for center or huge value)
                         int*      lambda,        // <--  lagrangian parameter for determining motion cost
                         int       apply_weights  // <--  perform weight based ME
                         )
{
  int   list_offset   = img->mb_data[img->current_mb_nr].list_offset;

  int   pos, best_pos, mcost;
  int   cand_mv_x, cand_mv_y;

  int   blocksize_x     = params->blc_size[blocktype][0];
  int   blocksize_y     = params->blc_size[blocktype][1];

  int   pic4_pix_x      = ((pic_pix_x + IMG_PAD_SIZE)<< 2);
  int   pic4_pix_y      = ((pic_pix_y + IMG_PAD_SIZE)<< 2);

  int   max_pos2        = ( !start_me_refinement_hp ? imax(1,search_pos2) : search_pos2);
  int   cmv_x, cmv_y;
  int   smv_x = mv2->mv_x + pic4_pix_x;
  int   smv_y = mv2->mv_y + pic4_pix_y;

  StorablePicture *ref_picture1 = listX[list       + list_offset][ref];
  StorablePicture *ref_picture2 = listX[(list ^ 1) + list_offset][0];

  int lambda_factor = lambda[H_PEL];

  ref_pic1_sub.luma = ref_picture1->p_curr_img_sub;
  ref_pic2_sub.luma = ref_picture2->p_curr_img_sub;
  img_width    = ref_picture1->size_x;
  img_height   = ref_picture1->size_y;
  width_pad    = ref_picture1->size_x_pad;
  height_pad   = ref_picture1->size_y_pad;

  if (apply_weights)
  {
    computeBiPred = computeBiPred2[H_PEL];
  }
  else
  {
    computeBiPred = computeBiPred1[H_PEL];
  }


  if ( ChromaMEEnable )
  {
    ref_pic1_sub.crcb[0] = ref_picture1->imgUV_sub[0];
    ref_pic1_sub.crcb[1] = ref_picture1->imgUV_sub[1];
    ref_pic2_sub.crcb[0] = ref_picture2->imgUV_sub[0];
    ref_pic2_sub.crcb[1] = ref_picture2->imgUV_sub[1];
    width_pad_cr  = ref_picture1->size_x_cr_pad;
    height_pad_cr = ref_picture1->size_y_cr_pad;
  }

  /*********************************
   *****                       *****
   *****  HALF-PEL REFINEMENT  *****
   *****                       *****
   *********************************/

  //===== loop over search positions =====
  for (best_pos = 0, pos = start_me_refinement_hp; pos < max_pos2; pos++)
  {
    cand_mv_x = mv1->mv_x + (spiral_hpel_search_x[pos]);    // quarter-pel units
    cand_mv_y = mv1->mv_y + (spiral_hpel_search_y[pos]);    // quarter-pel units

    //----- set motion vector cost -----
    mcost  = MV_COST_SMP (lambda_factor, cand_mv_x, cand_mv_y, pred_mv1->mv_x, pred_mv1->mv_y);
    mcost += MV_COST_SMP (lambda_factor, mv2->mv_x, mv2->mv_y, pred_mv2->mv_x, pred_mv2->mv_y);

    if (mcost >= min_mcost) continue;

    cmv_x = cand_mv_x + pic4_pix_x;
    cmv_y = cand_mv_y + pic4_pix_y;

    mcost += computeBiPred(orig_pic, blocksize_y, blocksize_x,
      min_mcost - mcost, cmv_x, cmv_y, smv_x, smv_y);

    if (mcost < min_mcost)
    {
      min_mcost = mcost;
      best_pos  = pos;
    }
  }

  if (best_pos)
  {
    mv1->mv_x += (spiral_hpel_search_x [best_pos]);
    mv1->mv_y += (spiral_hpel_search_y [best_pos]);
  }

  computeBiPred = apply_weights? computeBiPred2[Q_PEL] : computeBiPred1[Q_PEL];

  /************************************
  *****                          *****
  *****  QUARTER-PEL REFINEMENT  *****
  *****                          *****
  ************************************/

  if ( !start_me_refinement_qp )
    min_mcost = INT_MAX;

  lambda_factor = lambda[Q_PEL];

  //===== loop over search positions =====
  for (best_pos = 0, pos = start_me_refinement_qp; pos < search_pos4; pos++)
  {
    cand_mv_x = mv1->mv_x + spiral_search_x[pos];    // quarter-pel units
    cand_mv_y = mv1->mv_y + spiral_search_y[pos];    // quarter-pel units

    //----- set motion vector cost -----
    mcost  = MV_COST_SMP (lambda_factor, cand_mv_x, cand_mv_y, pred_mv1->mv_x, pred_mv1->mv_y);
    mcost += MV_COST_SMP (lambda_factor, mv2->mv_x, mv2->mv_y, pred_mv2->mv_x, pred_mv2->mv_y);

    if (mcost >= min_mcost) continue;
    cmv_x = cand_mv_x + pic4_pix_x;
    cmv_y = cand_mv_y + pic4_pix_y;

    mcost += computeBiPred(orig_pic, blocksize_y, blocksize_x,
      min_mcost - mcost, cmv_x, cmv_y, smv_x, smv_y);

    if (mcost < min_mcost)
    {
      min_mcost = mcost;
      best_pos  = pos;
    }

  }

  if (best_pos)
  {
    mv1->mv_x += spiral_search_x [best_pos];
    mv1->mv_y += spiral_search_y [best_pos];
  }

  //===== return minimum motion cost =====
  return min_mcost;
}

