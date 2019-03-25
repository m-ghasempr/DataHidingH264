/*!
 ***************************************************************************
 * \file rd_intra_jm444.c
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

/*!
 *************************************************************************************
 * \brief
 *    Mode Decision for an 4x4 Intra block
 *************************************************************************************
 */
int Mode_Decision_for_4x4IntraBlocks_JM_High444 (Macroblock *currMB, int  b8,  int  b4,  int  lambda,  distblk*  min_cost)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;
  Slice *currSlice = currMB->p_slice;
  RDOPTStructure  *p_RDO = currSlice->p_RDO;

  int    ipmode, best_ipmode = 0, i, j, y, dummy;
  int    c_nz, nonzero = 0;
  int*   ACLevel = currSlice->cofAC[b8][b4][0];
  int*   ACRun   = currSlice->cofAC[b8][b4][1];
  int    uv;
  distblk rdcost = 0;
  distblk min_rdcost  = DISTBLK_MAX;
  int    block_x     = ((b8 & 0x01) << 3) + ((b4 & 0x01) << 2);
  int    block_y     = ((b8 >> 1) << 3)  + ((b4 >> 1) << 2);
  int    pic_pix_x   = currMB->pix_x  + block_x;
  int    pic_pix_y   = currMB->pix_y  + block_y;
  int    pic_opix_x  = currMB->pix_x + block_x;
  int    pic_opix_y  = currMB->opix_y + block_y;
  int    pic_block_x = pic_pix_x >> 2;
  int    pic_block_y = pic_pix_y >> 2;

  int left_available, up_available, all_available;
  int *mb_size = p_Vid->mb_size[IS_LUMA];

  char   upMode, leftMode;
  int    mostProbableMode;

  PixelPos left_block, top_block;

  int  lrec4x4[4][4];
  int best_nz_coeff = 0;
  int block_x4 = block_x>>2;
  int block_y4 = block_y>>2;

#ifdef BEST_NZ_COEFF
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
  currMB->ipmode_DPCM = NO_INTRA_PMODE; ////For residual DPCM

  //===== INTRA PREDICTION FOR 4x4 BLOCK =====
  // set intra prediction values for 4x4 intra prediction
  set_intrapred_4x4(currMB, PLANE_Y, pic_pix_x, pic_pix_y, &left_available, &up_available, &all_available);  

  if (currSlice->P444_joined)
  {
    select_plane(p_Vid, PLANE_U);
    set_intrapred_4x4(currMB, PLANE_U, pic_pix_x, pic_pix_y, &left_available, &up_available, &all_available);  
    select_plane(p_Vid, PLANE_V);
    set_intrapred_4x4(currMB, PLANE_V, pic_pix_x, pic_pix_y, &left_available, &up_available, &all_available);  
    select_plane(p_Vid, PLANE_Y);
  }

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
      get_intrapred_4x4(currMB, PLANE_Y, ipmode, block_x, block_y, left_available, up_available);

      // get prediction and prediction error
      generate_pred_error_4x4(&p_Vid->pCurImg[pic_opix_y], currSlice->mpr_4x4[0][ipmode], &currSlice->mb_pred[0][block_y], &currSlice->mb_ores[0][block_y], pic_opix_x, block_x);     

      if (p_Vid->yuv_format == YUV444)
      {
        currMB->ipmode_DPCM = (short) ipmode;
        if (!IS_INDEPENDENT(p_Inp)) 
        {
          // generate intra 4x4 prediction block given availability
          get_intrapred_4x4(currMB, PLANE_U, ipmode, block_x, block_y, left_available, up_available);
          generate_pred_error_4x4(&p_Vid->pImgOrg[1][pic_opix_y], currSlice->mpr_4x4[1][ipmode], &currSlice->mb_pred[1][block_y], &currSlice->mb_ores[1][block_y], pic_opix_x, block_x);
          // generate intra 4x4 prediction block given availability
          get_intrapred_4x4(currMB, PLANE_V, ipmode, block_x, block_y, left_available, up_available);
          generate_pred_error_4x4(&p_Vid->pImgOrg[2][pic_opix_y], currSlice->mpr_4x4[2][ipmode], &currSlice->mb_pred[2][block_y], &currSlice->mb_ores[2][block_y], pic_opix_x, block_x);     
        }
      }

      // get and check rate-distortion cost
#ifdef BEST_NZ_COEFF
      currMB->cbp_bits[0] = cbp_bits;
#endif      

      rdcost = currSlice->rdcost_for_4x4_intra_blocks (currMB, &c_nz, b8, b4, ipmode, lambda, mostProbableMode, min_rdcost);
      if ((rdcost < min_rdcost) || (rdcost == min_rdcost && ipmode == mostProbableMode))
      {
        //--- set coefficients ---
        memcpy(p_RDO->cofAC4x4[0], ACLevel, 18 * sizeof(int));
        memcpy(p_RDO->cofAC4x4[1], ACRun,   18 * sizeof(int));

        //--- set reconstruction ---
        copy_4x4block(p_RDO->rec4x4[PLANE_Y], &p_Vid->enc_picture->imgY[pic_pix_y], 0, pic_pix_x);

        // SP/SI reconstruction
        if(currSlice->slice_type == SP_SLICE && !currSlice->sp2_frame_indicator)
        {
          for (y=0; y<4; y++)
          {
            memcpy(lrec4x4[y],&p_Vid->lrec[pic_pix_y+y][pic_pix_x], BLOCK_SIZE * sizeof(int));// stores the mode coefficients
          }
        }

        if(currSlice->P444_joined) 
        { 
          //--- set coefficients ---
          for (uv=0; uv < 2; uv++)
          {
            memcpy(p_RDO->cofAC4x4CbCr[uv][0],currSlice->cofAC[b8+4+uv*4][b4][0], 18 * sizeof(int));
            memcpy(p_RDO->cofAC4x4CbCr[uv][1],currSlice->cofAC[b8+4+uv*4][b4][1], 18 * sizeof(int));
            currMB->cr_cbp[uv + 1] = currMB->c_nzCbCr[uv + 1];

            //--- set reconstruction ---
            copy_4x4block(p_RDO->rec4x4[uv + 1], &p_Vid->enc_picture->imgUV[uv][pic_pix_y], 0, pic_pix_x);
          }
        }
        //--- flag if dct-coefficients must be coded ---
        nonzero = c_nz;

        //--- set best mode update minimum cost ---
        *min_cost     = rdcost;
        min_rdcost    = rdcost;
        best_ipmode   = ipmode;

        best_nz_coeff = p_Vid->nz_coeff [p_Vid->current_mb_nr][block_x4][block_y4];
#ifdef BEST_NZ_COEFF
        best_coded_block_flag = (int)((currMB->cbp_bits[0] >> bit_pos)&(int64)(1));
#endif
        if (p_Vid->AdaptiveRounding)
        {
          store_adaptive_rounding_4x4 (p_Vid, p_Vid->ARCofAdj4x4, I4MB, block_y, block_x);
        }
      }
    }
  }
#if INTRA_RDCOSTCALC_NNZ
  p_Vid->nz_coeff [p_Vid->current_mb_nr][block_x4][block_y4] = best_nz_coeff;
#endif
#ifdef BEST_NZ_COEFF
  cbp_bits &= (~(int64)(1<<bit_pos));
  cbp_bits |= (int64)(best_coded_block_flag<<bit_pos);
#endif

  //===== set intra mode prediction =====
  p_Vid->ipredmode[pic_block_y][pic_block_x] = (char) best_ipmode;
  currMB->intra_pred_modes[4*b8+b4] =
    (char) (mostProbableMode == best_ipmode ? -1 : (best_ipmode < mostProbableMode ? best_ipmode : best_ipmode-1));

  if(currSlice->P444_joined)
  {
    ColorPlane k;
    for (k = PLANE_U; k <= PLANE_V; k++)
    {
      select_plane(p_Vid, k);

      copy_4x4block(&currSlice->mb_pred[k][block_y], currSlice->mpr_4x4[k][best_ipmode], block_x, 0);
      for (j=0; j<4; j++)
      {
        for (i=0; i<4; i++)
        {
          currSlice->mb_ores[k][block_y+j][block_x+i]   = p_Vid->pImgOrg[k][currMB->pix_y+block_y+j][currMB->pix_x+block_x+i] - currSlice->mpr_4x4[k][best_ipmode][j][i];
        }
      }
      currMB->cr_cbp[k] = currMB->trans_4x4(currMB, k, block_x,block_y,&dummy,1);
    }
    select_plane(p_Vid, PLANE_Y);
  }

  //===== restore coefficients =====
  memcpy (ACLevel, p_RDO->cofAC4x4[0], 18 * sizeof(int));
  memcpy (ACRun,   p_RDO->cofAC4x4[1], 18 * sizeof(int));

  //===== restore reconstruction and prediction (needed if single coeffs are removed) =====
  copy_4x4block(&p_Vid->enc_picture->imgY[pic_pix_y], p_RDO->rec4x4[PLANE_Y], pic_pix_x, 0);
  copy_4x4block(&currSlice->mb_pred[0][block_y], currSlice->mpr_4x4[0][best_ipmode], block_x, 0);

  // SP/SI reconstuction
  if(currSlice->slice_type == SP_SLICE && !p_Vid->sp2_frame_indicator)
  {
    for (y=0; y<BLOCK_SIZE; y++)
    {
      memcpy (&p_Vid->lrec[pic_pix_y+y][pic_pix_x], lrec4x4[y], BLOCK_SIZE * sizeof(int));//restore coefficients when encoding primary SP frame
    }
  }
  if (currSlice->P444_joined) 
  {
    for (uv=0; uv < 2; uv++ )
    {
      //===== restore coefficients =====
      memcpy(currSlice->cofAC[b8+4+uv*4][b4][0], p_RDO->cofAC4x4CbCr[uv][0], 18 * sizeof(int));
      memcpy(currSlice->cofAC[b8+4+uv*4][b4][1], p_RDO->cofAC4x4CbCr[uv][1], 18 * sizeof(int));
      //===== restore reconstruction and prediction (needed if single coeffs are removed) =====
      copy_4x4block(&p_Vid->enc_picture->imgUV[uv][pic_pix_y], p_RDO->rec4x4[uv + 1], pic_pix_x, 0);
      copy_4x4block(&currSlice->mb_pred[uv + 1][block_y], currSlice->mpr_4x4[uv + 1][best_ipmode], block_x, 0);
    }
  }

  if (p_Vid->AdaptiveRounding)
  {
    update_adaptive_rounding_4x4 (p_Vid,p_Vid->ARCofAdj4x4, I4MB, block_y, block_x);
  }

  return nonzero;
}

/*!
 *************************************************************************************
 * \brief
 *    Mode Decision for an 4x4 Intra block
 *************************************************************************************
 */
int Mode_Decision_for_4x4IntraBlocks_JM_Low444 (Macroblock *currMB, int  b8,  int  b4,  int  lambda,  distblk*  min_cost)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;
  Slice *currSlice = currMB->p_slice;

  int     ipmode, best_ipmode = 0, i, j, dummy;
  distblk cost;
  int     nonzero = 0;

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

  PixelPos left_block;
  PixelPos top_block;

  distblk  fixedcost = weighted_cost(lambda, 4);   //(int) floor(4 * lambda );
  int  *mb_size  = p_Vid->mb_size[IS_LUMA];
  int best_nz_coeff = 0;
  int block_x4 = block_x>>2;
  int block_y4 = block_y>>2;

#ifdef BEST_NZ_COEFF
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
  currMB->ipmode_DPCM = NO_INTRA_PMODE; ////For residual DPCM

  //===== INTRA PREDICTION FOR 4x4 BLOCK =====
  // set intra prediction values for 4x4 intra prediction
  set_intrapred_4x4(currMB, PLANE_Y, pic_pix_x, pic_pix_y, &left_available, &up_available, &all_available);  

  if (currSlice->P444_joined)
  {
    select_plane(p_Vid, PLANE_U);
    set_intrapred_4x4(currMB, PLANE_U, pic_pix_x, pic_pix_y, &left_available, &up_available, &all_available);  
    select_plane(p_Vid, PLANE_V);
    set_intrapred_4x4(currMB, PLANE_V, pic_pix_x, pic_pix_y, &left_available, &up_available, &all_available);  
    select_plane(p_Vid, PLANE_Y);
  }

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
      if (currSlice->P444_joined)
      {
        get_intrapred_4x4(currMB, PLANE_U, ipmode, block_x, block_y, left_available, up_available);
        currSlice->compute_cost4x4(p_Vid, &p_Vid->pImgOrg[1][pic_opix_y], currSlice->mpr_4x4[1][ipmode], pic_opix_x, &cost, (int) *min_cost);
        get_intrapred_4x4(currMB, PLANE_V, ipmode, block_x, block_y, left_available, up_available);
        currSlice->compute_cost4x4(p_Vid, &p_Vid->pImgOrg[2][pic_opix_y], currSlice->mpr_4x4[2][ipmode], pic_opix_x, &cost, (int) *min_cost);
      }
      
      if (cost < *min_cost)
      {
        best_ipmode = ipmode;
        *min_cost   = cost;
      }
    }
  }

#if INTRA_RDCOSTCALC_NNZ
  p_Vid->nz_coeff [p_Vid->current_mb_nr][block_x4][block_y4] = best_nz_coeff;
#endif
#ifdef BEST_NZ_COEFF  
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
  nonzero = currMB->cr_cbp[0] = currMB->trans_4x4 (currMB, PLANE_Y, block_x, block_y, &dummy, 1);

  if (currSlice->P444_joined)
  {
    ColorPlane k;
    for (k = PLANE_U; k <= PLANE_V; k++)
    {
      select_plane(p_Vid, k);
      for (j=0; j<4; j++)
      {
        for (i=0; i<4; i++)
        {
          currSlice->mb_pred[k][block_y + j][block_x+i] = currSlice->mpr_4x4[k][best_ipmode][j][i];    
          currSlice->mb_ores[k][block_y + j][block_x+i] = p_Vid->pImgOrg[k][pic_opix_y+j][pic_opix_x+i] - currSlice->mpr_4x4[k][best_ipmode][j][i]; 
        }
      }

      currMB->cr_cbp[k] = currMB->trans_4x4 (currMB, k, block_x, block_y, &dummy, 1);
    }
    select_plane(p_Vid, PLANE_Y);
  }

  return nonzero;
}

/*!
*************************************************************************************
* \brief
*    Intra 16x16 mode decision
*************************************************************************************
*/
void Intra16x16_Mode_Decision444 (Macroblock* currMB)
{
  Slice *currSlice = currMB->p_slice;
  if (!currSlice->P444_joined)
  {
    /* generate intra prediction samples for all 4 16x16 modes */
    intrapred_16x16 (currMB, PLANE_Y);
    currSlice->find_sad_16x16 (currMB);   /* get best new intra mode */
    currMB->cbp = currMB->trans_16x16 (currMB, PLANE_Y);    
  }
  else
  {
    VideoParameters *p_Vid = currSlice->p_Vid;

    /* generate intra prediction samples for all 4 16x16 modes */
    intrapred_16x16 (currMB, PLANE_Y);
    select_plane(p_Vid, PLANE_U);
    intrapred_16x16 (currMB, PLANE_U);
    select_plane(p_Vid, PLANE_V);
    intrapred_16x16 (currMB, PLANE_V);
    select_plane(p_Vid, PLANE_Y);

    currSlice->find_sad_16x16 (currMB);   /* get best new intra mode */

    currMB->cbp = currMB->trans_16x16 (currMB, PLANE_Y);
    select_plane(p_Vid, PLANE_U);
    currSlice->cmp_cbp[1] = currMB->trans_16x16 (currMB, PLANE_U);
    select_plane(p_Vid, PLANE_V);
    currSlice->cmp_cbp[2] = currMB->trans_16x16 (currMB, PLANE_V);
    select_plane(p_Vid, PLANE_Y);

    currMB->cbp |= (currSlice->cmp_cbp[1] | currSlice->cmp_cbp[2]);
    currSlice->cmp_cbp[1] = currMB->cbp;
    currSlice->cmp_cbp[2] = currMB->cbp;
  }
}

