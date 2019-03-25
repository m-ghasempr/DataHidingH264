/*!
 ***************************************************************************
 * \file rd_intra_jm_low.c
 *
 * \brief
 *    Rate-Distortion optimized mode decision
 *
 * \author
 *    - Heiko Schwarz              <hschwarz@hhi.de>
 *    - Valeri George              <george@hhi.de>
 *    - Lowell Winger              <lwinger@lsil.com>
 *    - Alexis Michael Tourapis    <alexismt@ieee.org>
 * \date
 *    12. April 2001
 **************************************************************************
 */

#include <limits.h>

#include "global.h"

#include "image.h"
#include "macroblock.h"
#include "mb_access.h"
#include "rdopt_coding_state.h"
#include "mode_decision.h"
#include "rdopt.h"
#include "rd_intra_jm.h"
#include "q_around.h"
#include "intra4x4.h"
#include "intra8x8.h"

/*!
 *************************************************************************************
 * \brief
 *    Mode Decision for an 4x4 Intra block
 *************************************************************************************
 */
int Mode_Decision_for_4x4IntraBlocks_JM_Low (Macroblock *currMB, int  b8,  int  b4,  int  lambda,  distblk*  min_cost)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;
  Slice *currSlice = currMB->p_slice;

  int     ipmode, best_ipmode = 0, dummy;
  int     nonzero = 0;
  distblk cost;
  int  block_x     = ((b8 & 0x01) << 3) + ((b4 & 0x01) << 2);
  int  block_y     = ((b8 >> 1) << 3)  + ((b4 >> 1) << 2);
  int  pic_pix_x   = currMB->pix_x  + block_x;
  int  pic_pix_y   = currMB->pix_y  + block_y;
  int  pic_opix_x  = currMB->pix_x + block_x;
  int  pic_opix_y  = currMB->opix_y + block_y;
  int  pic_block_x = pic_pix_x >> 2;
  int  pic_block_y = pic_pix_y >> 2;

  int left_available, up_available, all_available;

  char   upMode, leftMode;
  int    mostProbableMode;

  PixelPos left_block, top_block;

  distblk  fixedcost = weighted_cost(lambda, 4); //(int) floor(4 * lambda );
  int  *mb_size  = p_Vid->mb_size[IS_LUMA];


#ifdef BEST_NZ_COEFF
  int best_nz_coeff = 0;
  int best_coded_block_flag = 0;
  int bit_pos = 1 + ((((b8>>1)<<1)+(b4>>1))<<2) + (((b8&1)<<1)+(b4&1));
  int64 cbp_bits;

  if (b8==0 && b4==0)
    cbp_bits = 0;
#endif

  get4x4Neighbour(currMB, block_x - 1, block_y    , mb_size, &left_block);
  get4x4Neighbour(currMB, block_x,     block_y - 1, mb_size, &top_block );

  // constrained intra pred
  if (p_Inp->UseConstrainedIntraPred)
  {
    left_block.available = left_block.available ? p_Vid->intra_block[left_block.mb_addr] : 0;
    top_block.available  = top_block.available  ? p_Vid->intra_block[top_block.mb_addr]  : 0;
  }

  upMode            =  top_block.available ? p_Vid->ipredmode[top_block.pos_y ][top_block.pos_x ] : -1;
  leftMode          = left_block.available ? p_Vid->ipredmode[left_block.pos_y][left_block.pos_x] : -1;

  mostProbableMode  = (upMode < 0 || leftMode < 0) ? DC_PRED : upMode < leftMode ? upMode : leftMode;
  *min_cost = DISTBLK_MAX;
  currMB->ipmode_DPCM = NO_INTRA_PMODE; //For residual DPCM

  //===== INTRA PREDICTION FOR 4x4 BLOCK =====
  // set intra prediction values for 4x4 intra prediction
  set_intrapred_4x4(currMB, PLANE_Y, pic_pix_x, pic_pix_y, &left_available, &up_available, &all_available);  

  //===== LOOP OVER ALL 4x4 INTRA PREDICTION MODES =====
  for (ipmode = 0; ipmode < NO_INTRA_PMODE; ipmode++)
  {
    int available_mode =  (all_available) || (ipmode==DC_PRED) ||
      (up_available && (ipmode==VERT_PRED||ipmode==VERT_LEFT_PRED||ipmode==DIAG_DOWN_LEFT_PRED)) ||
      (left_available && (ipmode==HOR_PRED||ipmode==HOR_UP_PRED));

    if (valid_intra_mode(currSlice, ipmode) == 0)
      continue;

    if( available_mode)
    {
      // generate intra 4x4 prediction block given availability
      // Note that some checks may not be necessary internally since the "available_mode" test has already tested
      // mode availability.
      get_intrapred_4x4(currMB, PLANE_Y, ipmode, block_x, block_y, left_available, up_available);
      cost  = (ipmode == mostProbableMode) ? 0 : fixedcost;
      currSlice->compute_cost4x4(p_Vid, &p_Vid->pCurImg[pic_opix_y], currSlice->mpr_4x4[0][ipmode], pic_opix_x, &cost, *min_cost);
      if (cost < *min_cost)
      {
        best_ipmode = ipmode;
        *min_cost   = cost;
      }
    }
  }

#ifdef BEST_NZ_COEFF
  p_Vid->nz_coeff [p_Vid->current_mb_nr][block_x4][block_y4] = best_nz_coeff;
  cbp_bits &= (~(int64)(1<<bit_pos));
  cbp_bits |= (int64)(best_coded_block_flag<<bit_pos);
#endif
  //===== set intra mode prediction =====
  p_Vid->ipredmode[pic_block_y][pic_block_x] = (char) best_ipmode;
  currMB->intra_pred_modes[4*b8+b4] =
    (char) (mostProbableMode == best_ipmode ? -1 : (best_ipmode < mostProbableMode ? best_ipmode : best_ipmode-1));

  // get prediction and prediction error
  generate_pred_error_4x4(&p_Vid->pCurImg[pic_opix_y], currSlice->mpr_4x4[0][best_ipmode], &currSlice->mb_pred[0][block_y], &currSlice->mb_ores[0][block_y], pic_opix_x, block_x);

  currMB->ipmode_DPCM = (short) best_ipmode;  

  select_dct(currMB);

  if (currMB->mb_type == I4MB)
	  nonzero = currMB->cr_cbp[0] = dct_4x4 (currMB, PLANE_Y, block_x, block_y, &dummy, 1);
  else
	  nonzero = currMB->cr_cbp[0] = currMB->trans_4x4 (currMB, PLANE_Y, block_x, block_y, &dummy, 1);

  return nonzero;
}

/*!
 *************************************************************************************
 * \brief
 *    8x8 Intra mode decision for a macroblock - Low complexity
 *************************************************************************************
 */
int Mode_Decision_for_8x8IntraBlocks_JM_Low (Macroblock *currMB, int b8, int lambda, distblk *min_cost)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;
  Slice *currSlice = currMB->p_slice;

  int     ipmode, best_ipmode = 0, i, j, dummy;
  distblk     cost;
  int     nonzero = 0;
  int     block_x     = (b8 & 0x01) << 3;
  int     block_y     = (b8 >> 1) << 3;
  int     pic_pix_x   = currMB->pix_x + block_x;
  int     pic_pix_y   = currMB->pix_y + block_y;
  int     pic_opix_x   = currMB->pix_x + block_x;
  int     pic_opix_y   = currMB->opix_y + block_y;
  int     pic_block_x = pic_pix_x >> 2;
  int     pic_block_y = pic_pix_y >> 2;
  int     mb_block_y  = (currMB->block_y) + ((b8 >> 1) << 1);
  int     mb_block_x  = (currMB->block_x) + ((b8 & 0x01) << 1);

  imgpel    *img_org, *img_prd;
  int       *residual;
  int left_available, up_available, all_available;
  int    **mb_ores = currSlice->mb_ores[0]; 
  imgpel **mb_pred = currSlice->mb_pred[0];
  int *mb_size = p_Vid->mb_size[IS_LUMA];

  char   upMode;
  char   leftMode;
  int    mostProbableMode;

  PixelPos left_block;
  PixelPos top_block;

  get4x4Neighbour(currMB, block_x - 1, block_y    , mb_size, &left_block);
  get4x4Neighbour(currMB, block_x,     block_y - 1, mb_size, &top_block );

  if (p_Inp->UseConstrainedIntraPred)
  {
    top_block.available  = top_block.available ? p_Vid->intra_block [top_block.mb_addr] : 0;
    left_block.available = left_block.available ? p_Vid->intra_block [left_block.mb_addr] : 0;
  }

  if(b8 >> 1)
    upMode    =  top_block.available ? p_Vid->ipredmode8x8[top_block.pos_y ][top_block.pos_x ] : -1;
  else
    upMode    =  top_block.available ? p_Vid->ipredmode   [top_block.pos_y ][top_block.pos_x ] : -1;
  if(b8 & 0x01)
    leftMode  = left_block.available ? p_Vid->ipredmode8x8[left_block.pos_y][left_block.pos_x] : -1;
  else
    leftMode  = left_block.available ? p_Vid->ipredmode[left_block.pos_y][left_block.pos_x] : -1;

  mostProbableMode  = (upMode < 0 || leftMode < 0) ? DC_PRED : upMode < leftMode ? upMode : leftMode;
  *min_cost = DISTBLK_MAX;
  currMB->ipmode_DPCM = NO_INTRA_PMODE; //For residual DPCM

  //===== INTRA PREDICTION FOR 8x8 BLOCK =====
  set_intrapred_8x8(currMB, PLANE_Y, pic_pix_x, pic_pix_y, &left_available, &up_available, &all_available);

  //===== LOOP OVER ALL 8x8 INTRA PREDICTION MODES =====
  for (ipmode = 0; ipmode < NO_INTRA_PMODE; ipmode++)
  {
    if( (ipmode==DC_PRED) ||
      ((ipmode==VERT_PRED||ipmode==VERT_LEFT_PRED||ipmode==DIAG_DOWN_LEFT_PRED) && up_available ) ||
      ((ipmode==HOR_PRED||ipmode==HOR_UP_PRED) && left_available ) ||
      (all_available) )
    {
      get_intrapred_8x8(currMB, PLANE_Y, ipmode, left_available, up_available);
      cost  = (ipmode == mostProbableMode) ? 0 : (weighted_cost(lambda, 4)); //(int)floor(4 * lambda );
      currSlice->compute_cost8x8(p_Vid, &p_Vid->pImgOrg[0][pic_opix_y], currSlice->mpr_8x8[0][ipmode], pic_opix_x, &cost, *min_cost);
      if (cost < *min_cost)
      {
        best_ipmode = ipmode;
        *min_cost   = cost;
      }
    }
  }

  //===== set intra mode prediction =====
  p_Vid->ipredmode8x8[pic_block_y][pic_block_x] = (char) best_ipmode;
  currMB->ipmode_DPCM = (short) best_ipmode; //For residual DPCM

  currMB->intra_pred_modes8x8[4*b8] = (char)((mostProbableMode == best_ipmode)
    ? -1
    : (best_ipmode < mostProbableMode ? best_ipmode : best_ipmode-1));

  for(j = mb_block_y; j < mb_block_y + 2; j++)   //loop 4x4s in the subblock for 8x8 prediction setting
    memset(&p_Vid->ipredmode8x8[j][mb_block_x], best_ipmode, 2 * sizeof(char));

  // get prediction and prediction error
  for (j = block_y; j < block_y + 8; j++)
  {
    memcpy(&mb_pred[j][block_x],currSlice->mpr_8x8[0][best_ipmode][j - block_y], 8 * sizeof(imgpel));
    img_org  = &p_Vid->pCurImg[currMB->opix_y + j][pic_opix_x];
    img_prd  = &mb_pred[j][block_x];
    residual = &mb_ores[j][block_x];
    for (i=0; i<8; i++)
    {
      *residual++ = *img_org++ - *img_prd++;
    }
  }

  currMB->ipmode_DPCM = (short) best_ipmode;

  nonzero = currMB->trans_8x8 (currMB, PLANE_Y, b8, &dummy, 1);    
  return nonzero;
}

