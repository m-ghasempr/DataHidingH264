
/*!
 ***************************************************************************
 * \file mode_decision.c
 *
 * \brief
 *    Main macroblock mode decision functions and helpers
 *
 **************************************************************************
 */

#include <math.h>
#include <limits.h>
#include <float.h>

#include "global.h"
#include "rdopt_coding_state.h"
#include "mb_access.h"
#include "intrarefresh.h"
#include "image.h"
#include "transform8x8.h"
#include "ratectl.h"
#include "mode_decision.h"
#include "fmo.h"
#include "me_umhex.h"
#include "me_umhexsmp.h"
#include "macroblock.h"
#include "rdoq.h"
#include "errdo.h"
#include "q_around.h"
#include "slice.h"
#include "md_common.h"
#include "conformance.h"
#include "me_umhex.h"
#include "rdopt.h"


/*!
*************************************************************************************
* \brief
*    Reset Valid Modes
*************************************************************************************
*/
void reset_valid_modes(RD_PARAMS *enc_mb)
{
  memset(enc_mb->valid, 0, MAXMODE * sizeof(short));
}


/*!
*************************************************************************************
* \brief
*    Checks whether a primary SP slice macroblock was encoded as I16
*************************************************************************************
*/
static inline int check_for_SI16(int **lrec, int pix_x, int pix_y)
{
  int i,j;
  for (i = pix_y; i < pix_y + MB_BLOCK_SIZE; i++)
  {
    for (j = pix_x;j < pix_x + MB_BLOCK_SIZE; j++)
      if (lrec[i][j] != -16)
        return 0;
  }
  return 1;
}

/*!
*************************************************************************************
* \brief
*    Update parameters after encoding a macroblock
*************************************************************************************
*/
void end_encode_one_macroblock(Macroblock *currMB)
{
  Slice *currSlice = currMB->p_slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;

  int bslice = (currSlice->slice_type == B_SLICE);

  update_qp_cbp(currMB);

  if ( (currSlice->mb_aff_frame_flag)
    && (currMB->mbAddrX & 0x01)
    && (currMB->mb_type ? 0:((bslice) ? !currMB->cbp:1))  // bottom is skip
    && ((currMB->PrevMB)->mb_type ? 0:((bslice) ? !(currMB->PrevMB)->cbp:1))
    && !(field_flag_inference(currMB) == currMB->mb_field)) // top is skip
  {
    currSlice->rddata->min_rdcost = 1e30;  // don't allow coding of a MB pair as skip if wrong inference
  }
  else
    currSlice->rddata->min_rdcost = (double)currMB->min_rdcost;

  currSlice->rddata->min_dcost  = (double)currMB->min_dcost;
  currSlice->rddata->min_rate   = (double)currMB->min_rate;

  if(p_Inp->SearchMode == UM_HEX)
  {
    UMHEX_skip_intrabk_SAD(currMB, p_Vid->listXsize[currMB->list_offset]);
  }
  else if(p_Inp->SearchMode == UM_HEX_SIMPLE)
  {
    smpUMHEX_skip_intrabk_SAD(currMB);
  }

  //--- constrain intra prediction ---
  if(p_Inp->UseConstrainedIntraPred && (currSlice->slice_type == P_SLICE || currSlice->slice_type == B_SLICE))
  {
    p_Vid->intra_block[currMB->mbAddrX] = IS_INTRA(currMB);
  }
}

/*!
*************************************************************************************
* \brief
*    Initialize Encoding parameters for Macroblock
*************************************************************************************
*/
void init_enc_mb_params(Macroblock* currMB, RD_PARAMS *enc_mb, int intra)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;
  Slice *currSlice = currMB->p_slice;
  int bslice = (currSlice->slice_type == B_SLICE);

  int l,k;

  enc_mb->curr_mb_field = (short) ((currSlice->mb_aff_frame_flag)&&(currMB->mb_field));

  // Set valid modes  
  enc_mb->valid[I8MB]  = (short) ((!p_Inp->DisableIntraInInter || intra )?   p_Inp->Transform8x8Mode : 0);
  enc_mb->valid[I4MB]  = (short) ((!p_Inp->DisableIntraInInter || intra )? ((p_Inp->Transform8x8Mode == 2) ? 0 : 1) : 0);
  enc_mb->valid[I4MB]  = (short) ((!p_Inp->DisableIntra4x4  ) ? enc_mb->valid[I4MB] : 0);
  enc_mb->valid[I16MB] = (short) ((!p_Inp->DisableIntraInInter || intra )? 1 : 0);
  enc_mb->valid[I16MB] = (short) ((!p_Inp->DisableIntra16x16) ? enc_mb->valid[I16MB] : 0);
  enc_mb->valid[IPCM]  = (short) ((!p_Inp->DisableIntraInInter || intra )? p_Inp->EnableIPCM : 0);
  enc_mb->valid[SI4MB] = 0;

  enc_mb->valid[0]     = (short) (!intra && p_Inp->InterSearch[bslice][0]);
  enc_mb->valid[1]     = (short) (!intra && p_Inp->InterSearch[bslice][1]);
  enc_mb->valid[2]     = (short) (!intra && p_Inp->InterSearch[bslice][2]);
  enc_mb->valid[3]     = (short) (!intra && p_Inp->InterSearch[bslice][3]);
  enc_mb->valid[4]     = (short) (!intra && p_Inp->InterSearch[bslice][4]);
  enc_mb->valid[5]     = (short) (!intra && p_Inp->InterSearch[bslice][5] && !(p_Inp->Transform8x8Mode==2));
  enc_mb->valid[6]     = (short) (!intra && p_Inp->InterSearch[bslice][6] && !(p_Inp->Transform8x8Mode==2));
  enc_mb->valid[7]     = (short) (!intra && p_Inp->InterSearch[bslice][7] && !(p_Inp->Transform8x8Mode==2));
  enc_mb->valid[P8x8]  = (short) (enc_mb->valid[4] || enc_mb->valid[5] || enc_mb->valid[6] || enc_mb->valid[7]);


  if (currSlice->UseRDOQuant && p_Inp->RDOQ_CP_Mode && (p_Vid->qp != p_Vid->masterQP) )  
    RDOQ_update_mode(currSlice, enc_mb);
  
  if(currSlice->slice_type == SP_SLICE || currSlice->slice_type == SI_SLICE)
  {
    if(currSlice->slice_type == SI_SLICE)
    {
      reset_valid_modes(enc_mb);
      if(check_for_SI16(p_Vid->lrec, currMB->pix_x, currMB->pix_y))
      {
        enc_mb->valid[I16MB] = 1;
      }
      else
      {
        enc_mb->valid[I4MB]  = 1;
      }
    }

    if(p_Vid->sp2_frame_indicator)
    {
      if(check_for_SI16(p_Vid->lrec, currMB->pix_x, currMB->pix_y))
      {
        reset_valid_modes(enc_mb);
        enc_mb->valid[I16MB] = 1;
      }
      else
      {
        enc_mb->valid[I8MB]  = 0;
        enc_mb->valid[IPCM]  = 0;
        enc_mb->valid[0]     = 0;
        //enc_mb->valid[I16MB] = 0;
      }
    }
  }

  //===== SET LAGRANGE PARAMETERS =====
  // Note that these are now computed at the slice level to reduce
  // computations and cleanup code.
  if (bslice && p_Vid->nal_reference_idc)
  {
    enc_mb->lambda_md = p_Vid->lambda_md[5][p_Vid->masterQP];
//#if JCOST_CALC_SCALEUP
    enc_mb->lambda_mdfp = LAMBDA_FACTOR(enc_mb->lambda_md);
//#endif
    enc_mb->lambda_me[F_PEL] = p_Vid->lambda_me[5][p_Vid->masterQP][F_PEL];
    enc_mb->lambda_me[H_PEL] = p_Vid->lambda_me[5][p_Vid->masterQP][H_PEL];
    enc_mb->lambda_me[Q_PEL] = p_Vid->lambda_me[5][p_Vid->masterQP][Q_PEL];

    enc_mb->lambda_mf[F_PEL] = p_Vid->lambda_mf[5][p_Vid->masterQP][F_PEL];
    enc_mb->lambda_mf[H_PEL] = p_Vid->lambda_mf[5][p_Vid->masterQP][H_PEL];
    enc_mb->lambda_mf[Q_PEL] = p_Vid->lambda_mf[5][p_Vid->masterQP][Q_PEL];
  }
  else
  {
    enc_mb->lambda_md = p_Vid->lambda_md[currSlice->slice_type][p_Vid->masterQP];

//#if JCOST_CALC_SCALEUP
    enc_mb->lambda_mdfp = LAMBDA_FACTOR(enc_mb->lambda_md);
//#endif

    enc_mb->lambda_me[F_PEL] = p_Vid->lambda_me[currSlice->slice_type][p_Vid->masterQP][F_PEL];
    enc_mb->lambda_me[H_PEL] = p_Vid->lambda_me[currSlice->slice_type][p_Vid->masterQP][H_PEL];
    enc_mb->lambda_me[Q_PEL] = p_Vid->lambda_me[currSlice->slice_type][p_Vid->masterQP][Q_PEL];

    enc_mb->lambda_mf[F_PEL] = p_Vid->lambda_mf[currSlice->slice_type][p_Vid->masterQP][F_PEL];
    enc_mb->lambda_mf[H_PEL] = p_Vid->lambda_mf[currSlice->slice_type][p_Vid->masterQP][H_PEL];
    enc_mb->lambda_mf[Q_PEL] = p_Vid->lambda_mf[currSlice->slice_type][p_Vid->masterQP][Q_PEL];
  }

  if (!currSlice->mb_aff_frame_flag)
  {
    for (l = LIST_0; l < BI_PRED; l++)
    {
      for(k = 0; k < p_Vid->listXsize[l]; k++)
      {
        if(currSlice->structure != p_Vid->listX[l][k]->structure)
        {
          if (currSlice->structure == TOP_FIELD)
            p_Vid->listX[l][k]->chroma_vector_adjustment = -2;
          else if (currSlice->structure == BOTTOM_FIELD)
            p_Vid->listX[l][k]->chroma_vector_adjustment = 2;
          else
            p_Vid->listX[l][k]->chroma_vector_adjustment= 0;
        }
        else
          p_Vid->listX[l][k]->chroma_vector_adjustment= 0;
      }
    }
  }
  else
  {
    if (enc_mb->curr_mb_field)
    {
      for (l = currMB->list_offset; l <= currMB->list_offset + LIST_1; l++)
      {
        for(k = 0; k < p_Vid->listXsize[l]; k++)
        {
          p_Vid->listX[l][k]->chroma_vector_adjustment= 0;
          if((currMB->mbAddrX & 0x01) == 0 && p_Vid->listX[l][k]->structure == BOTTOM_FIELD)
            p_Vid->listX[l][k]->chroma_vector_adjustment = -2;
          if((currMB->mbAddrX & 0x01) == 1 && p_Vid->listX[l][k]->structure == TOP_FIELD)
            p_Vid->listX[l][k]->chroma_vector_adjustment = 2;
        }
      }
    }
    else
    {
      for (l = currMB->list_offset; l <= currMB->list_offset + LIST_1; l++)
      {
        for(k = 0; k < p_Vid->listXsize[l]; k++)
          p_Vid->listX[l][k]->chroma_vector_adjustment = 0;
      }
    }
  }

  // reset chroma intra predictor to default
  currMB->c_ipred_mode = DC_PRED_8;

  if(p_Inp->SearchMode == UM_HEX)
  {
    UMHEX_decide_intrabk_SAD(currMB);
    UMHEX_DefineThresholdMB(p_Vid, p_Inp);
  }
  else if (p_Inp->SearchMode == UM_HEX_SIMPLE)
  {
    smpUMHEX_decide_intrabk_SAD(currMB);
  }
}

/*!
*************************************************************************************
* \brief
*    computation of prediction list (including biprediction) cost
*************************************************************************************
*/
void list_prediction_cost(Macroblock *currMB, int list, int block, int mode, RD_PARAMS *enc_mb, distblk bmcost[5], char best_ref[2])
{
  Slice *currSlice = currMB->p_slice;
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;

  short ref;
  distblk mcost;
  int cur_list = list < BI_PRED ? currMB->list_offset + list : currMB->list_offset;
  int ref_lambda = (p_Inp->rdopt) ? enc_mb->lambda_mf[Q_PEL] :  enc_mb->lambda_mf[Q_PEL] >> 2;

  //--- get cost and reference frame for forward prediction ---

  if (list < BI_PRED)
  {
    for (ref=0; ref < p_Vid->listXsize[cur_list]; ref++)
    {
      if (!p_Vid->checkref || list || ref==0 || (p_Inp->RestrictRef && CheckReliabilityOfRef (currMB, block, list, ref, mode)))
      {
        // limit the number of reference frames to 1 when switching SP frames are used
        if((!p_Inp->sp2_frame_indicator && !p_Inp->sp_output_indicator)||
          ((p_Inp->sp2_frame_indicator || p_Inp->sp_output_indicator) && (currSlice->slice_type != P_SLICE && currSlice->slice_type != SP_SLICE))||
          ((p_Inp->sp2_frame_indicator || p_Inp->sp_output_indicator) && ((currSlice->slice_type == P_SLICE || currSlice->slice_type == SP_SLICE) &&(ref==0))))
        {
          mcost  = ref_cost(p_Vid, ref_lambda, (short) ref, cur_list);

          mcost += p_Vid->motion_cost[mode][list][ref][block];

          if (mcost < bmcost[list])
          {
            bmcost[list]   = mcost;
            best_ref[list] = (char)ref;
          }
        }
      }
    }
  }
  else if (list == BI_PRED)
  {
    if (p_Vid->active_pps->weighted_bipred_idc == 1)
    {
      int weight_sum = currSlice->wbp_weight[0][(int) best_ref[LIST_0]][(int) best_ref[LIST_1]][0] + currSlice->wbp_weight[1][(int) best_ref[LIST_0]][(int) best_ref[LIST_1]][0];

      if (weight_sum < -128 ||  weight_sum > 127)
      {
        bmcost[list] = DISTBLK_MAX;
      }
      else
      {
        bmcost[list]  = ref_cost(p_Vid, ref_lambda, (short) best_ref[LIST_0], cur_list) +
          ref_cost(p_Vid, ref_lambda, (short) best_ref[LIST_1], cur_list + LIST_1);
        bmcost[list] += BIDPartitionCost (currMB, mode, block, best_ref, enc_mb->lambda_mf[Q_PEL]);
      }
    }
    else
    {
      bmcost[list]  = ref_cost(p_Vid, ref_lambda, (short) best_ref[LIST_0], cur_list) +
        ref_cost(p_Vid, ref_lambda, (short) best_ref[LIST_1], cur_list + LIST_1);
      bmcost[list] += BIDPartitionCost (currMB, mode, block, best_ref, enc_mb->lambda_mf[Q_PEL]);
    }
  }
  else
  {
    bmcost[list]  = ref_cost(p_Vid, ref_lambda, 0, cur_list) + ref_cost(p_Vid, ref_lambda, 0, cur_list + LIST_1);
    bmcost[list] += BPredPartitionCost(currMB, mode, block, 0, 0, enc_mb->lambda_mf[Q_PEL], !(list&1));
  }
}

static inline distblk compute_ref_cost(VideoParameters *p_Vid, RD_PARAMS *enc_mb, int ref, int list)
{
  return weighted_cost(enc_mb->lambda_mf[Q_PEL],((p_Vid->listXsize[list] <= 1)? 0 : p_Vid->refbits[ref]));
}

/*!
*************************************************************************************
* \brief
*    Determination of prediction list based on simple distortion computation
*************************************************************************************
*/
void determine_prediction_list( distblk bmcost[5], Info8x8 *best, distblk *cost)
{
  int bestlist;  
  *cost += distblkminarray ( bmcost, 5, &bestlist);
  if (bestlist <= BI_PRED)  //LIST_0, LIST_1 & BI_DIR
  {
    best->pdir = (char) bestlist; 
    best->bipred= 0;
  }
  else                      //BI_PRED_L0 & BI_PRED_L1
  {
    best->pdir = 2;    
    best->bipred = (char) (bestlist - 2);
    best->ref[LIST_0] = 0;
    best->ref[LIST_1] = 0;
  }
}

/*!
*************************************************************************************
* \brief
*    RD decision process
*************************************************************************************
*/
void compute_mode_RD_cost(Macroblock *currMB,
                          RD_PARAMS *enc_mb,
                          short mode,
                          short *inter_skip)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;
  int terminate_16x16 = 0, terminate_trans = 0, ctr16x16 = 0;
  Slice *currSlice = currMB->p_slice;
  RDOPTStructure  *p_RDO = currSlice->p_RDO;
  int bslice = (currSlice->slice_type == B_SLICE);

  //--- transform size ---
  currMB->luma_transform_size_8x8_flag = (byte) (p_Inp->Transform8x8Mode==2
    ?  (mode >= 1 && mode <= 3)
    || (mode == 0 && bslice && p_Vid->active_sps->direct_8x8_inference_flag)
    || ((mode == P8x8) && (enc_mb->valid[4]))
    :  0);

  SetModesAndRefframeForBlocks (currMB, (short) mode);
  memset( currSlice->cofAC[0][0][0], 0, 2080 * sizeof(int)); // 4 * 4 * 2 * 65

  // Encode with coefficients
  currSlice->NoResidueDirect = 0;

  if ((p_Inp->FastCrIntraDecision ) || (currMB->c_ipred_mode == DC_PRED_8 || (IS_INTRA(currMB) )))
  {
    do
    {
      // This seems to have a problem since we are not properly copying the right b8x8info (transform based)
      terminate_16x16 = bslice_16x16_termination_control(p_Inp, p_Vid->b8x8info, &ctr16x16, mode, bslice);

      do
      {

        if (mode == P8x8)
        {
          int i;
          if (currMB->luma_transform_size_8x8_flag)
          {
            for (i = 0; i < 4; i++)
              p_Vid->b8x8info->best[P8x8][i] = p_RDO->tr8x8->part[i];
          }
          else
          {
            for (i = 0; i < 4; i++)
              p_Vid->b8x8info->best[P8x8][i] = p_RDO->tr4x4->part[i];
          }          
        }
        // check if prediction parameters are in valid range.
        if (CheckPredictionParams(currMB, p_Vid->b8x8info, mode) == TRUE)
        {
          if (RDCost_for_macroblocks (currMB, enc_mb->lambda_mdfp, mode))
          {
            //Rate control
            if (p_Inp->RCEnable)
            {
              if(mode == P8x8)
              {
                rc_store_diff(currSlice->diffy, &p_Vid->pCurImg[currMB->opix_y], currMB->pix_x, currMB->luma_transform_size_8x8_flag == TRUE ? p_RDO->tr8x8->mpr8x8 : p_RDO->tr4x4->mpr8x8);
              }
              else
                rc_store_diff(currSlice->diffy, &p_Vid->pCurImg[currMB->opix_y], currMB->pix_x, p_RDO->pred);
            }

            store_macroblock_parameters (currMB, mode);

            if(p_Inp->rdopt == 2 && mode == 0 && p_Inp->EarlySkipEnable)
            {
              // check transform quantized coeff.
              if(currMB->cbp == 0)
                *inter_skip = 1;
            }
          }        
        }
        // This code needs to be fixed - BUG? ATOUR

        terminate_trans = transform_termination_control(currMB, mode);

      }while (!terminate_trans);
    }while (!terminate_16x16);

    // Encode with no coefficients. Currently only for direct. This could be extended to all other modes as in example.
    //if (mode < P8x8 && (*inter_skip == 0) && enc_mb->valid[mode] && currMB->cbp && (currMB->cbp&15) != 15 && !p_Inp->nobskip)
    if (bslice && mode == 0 && (*inter_skip == 0) && enc_mb->valid[mode] 
    && currMB->cbp && (currMB->cbp&15) != 15 && !p_Inp->nobskip
      && !(currMB->qp_scaled[0] == 0 && p_Vid->lossless_qpprime_flag==1) )
    {
      currSlice->NoResidueDirect = 1;

      if (CheckPredictionParams(currMB, p_Vid->b8x8info, mode) == TRUE)
      {
        if (RDCost_for_macroblocks (currMB, enc_mb->lambda_mdfp, mode))
        {
          //Rate control
          if (p_Inp->RCEnable)
            rc_store_diff(currSlice->diffy, &p_Vid->pCurImg[currMB->opix_y], currMB->pix_x, p_RDO->pred);

          if (p_Vid->AdaptiveRounding)
            reset_adaptive_rounding_direct(p_Vid);

          store_macroblock_parameters (currMB, mode);
        }
      }
    }

    //modes 0 and 1 of a B frame 
    if (p_Vid->AdaptiveRounding && bslice && mode <= 1)
    { 
      if (currMB->temp_transform_size_8x8_flag)
        update_adaptive_rounding_16x16( p_Vid, p_Vid->ARCofAdj8x8, mode);
      else
        update_adaptive_rounding_16x16( p_Vid, p_Vid->ARCofAdj4x4, mode);
    }
  }
}


/*!
*************************************************************************************
* \brief
*    Mode Decision for an 8x8 sub-macroblock
*************************************************************************************
*/
void submacroblock_mode_decision(Macroblock *currMB,
                                 RD_PARAMS *enc_mb,
                                 RD_8x8DATA *dataTr,                                 
                                 int ****cofACtr,
                                 int block,                                 
                                 distblk *cost)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;
  Slice *currSlice = currMB->p_slice;
  PicMotionParams *motion = &p_Vid->enc_picture->motion;
  int transform8x8 = currMB->luma_transform_size_8x8_flag;
  int64 curr_cbp_blk;
  int j, k;
  int index;
  int mode;
  distblk min_rdcost, rdcost = 0;
  distblk min_cost8x8;
  distblk bmcost[5] = {DISTBLK_MAX};
  int cnt_nonz = 0;
  int best_cnt_nonz = 0;
  int maxindex =  (transform8x8) ? 2 : 5;
  int block_x, block_y;
  int lambda_mf[3];

  int ****fadjust = transform8x8? p_Vid->ARCofAdj8x8 : p_Vid->ARCofAdj4x4;

  //--- set coordinates ---
  int j0 = ((block>>1)<<3);
  int j1 = (j0>>2);
  int i0 = ((block&0x01)<<3);
  int i1 = (i0>>2);

  Boolean valid_8x8 = FALSE;
  Boolean stored_state_8x8 = FALSE;

#ifdef BEST_NZ_COEFF
  int best_nz_coeff[2][2];
#endif

  Info8x8 *partition = &(dataTr->part[block]);
  Info8x8 best;
  // Init best (need to create simple function)
  best.mode = 0;
  best.pdir = 0;
  best.bipred = 0;
  best.ref[LIST_0] = 0;
  best.ref[LIST_1] = -1;
  *partition = best;

#ifdef BEST_NZ_COEFF
  for(j = 0; j <= 1; j++)
  {
    for(i = 0; i <= 1; i++)
      best_nz_coeff[i][j] = p_Vid->nz_coeff[currMB->mbAddrX][i1 + i][j1 + j] = 0;
  }
#endif

  if(p_Inp->subMBCodingState == 2)
    currSlice->store_coding_state(currMB, currSlice->p_RDO->cs_tmp);

  //=====  LOOP OVER POSSIBLE CODING MODES FOR 8x8 SUB-PARTITION  =====
  for (min_cost8x8 = DISTBLK_MAX, min_rdcost = DISTBLK_MAX, index = (currSlice->slice_type == B_SLICE ? 0 : 1); index < maxindex; index++)
  {
    mode = b8_mode_table[index];
    best.mode = (char) mode;
    *cost = 0;

    if (enc_mb->valid[mode] && (!(transform8x8 == 1 && mode > 4)) && (transform8x8 == 0 || mode != 0 || (mode == 0 && p_Vid->active_sps->direct_8x8_inference_flag)))
    {
      if (transform8x8)
      {
        currMB->valid_8x8 = TRUE;
      }

      valid_8x8 = TRUE;
      curr_cbp_blk = 0;

      if (mode==0)  //--- Direct8x8 Mode ---
      {        
        block_x = currMB->block_x + (block & 0x01)*2;
        block_y = currMB->block_y + (block & 0x02);
        best.ref[LIST_0] = currSlice->direct_ref_idx[block_y][block_x][LIST_0];
        best.ref[LIST_1] = currSlice->direct_ref_idx[block_y][block_x][LIST_1];
        best.pdir        = currSlice->direct_pdir[block_y][block_x];
      } // if (mode==0)
      else
      {
        int64 ref_pic_num;
        char b_ref;

        //======= motion estimation for all reference frames ========
        //-----------------------------------------------------------
        memcpy(lambda_mf, enc_mb->lambda_mf, 3 * sizeof(int));

        if (p_Inp->CtxAdptLagrangeMult == 1)
        {
          RDOPTStructure *p_RDO = currSlice->p_RDO;
          lambda_mf[F_PEL] = (int)(lambda_mf[F_PEL] * p_RDO->lambda_mf_factor);
          lambda_mf[H_PEL] = (int)(lambda_mf[H_PEL] * p_RDO->lambda_mf_factor);
          lambda_mf[Q_PEL] = (int)(lambda_mf[Q_PEL] * p_RDO->lambda_mf_factor);
        }

        SubPartitionMotionSearch (currMB, mode, block, lambda_mf);

        //--- get cost and reference frame for LIST 0 prediction ---
        bmcost[LIST_0] = DISTBLK_MAX;
        list_prediction_cost(currMB, LIST_0, block, mode, enc_mb, bmcost, best.ref);

        //store LIST 0 reference index for every block
        block_x = currMB->block_x + (block & 0x01)*2;
        block_y = currMB->block_y + (block & 0x02);
        b_ref = best.ref[LIST_0];
        ref_pic_num = p_Vid->enc_picture->ref_pic_num[currMB->list_offset][(short) b_ref];

        memset(&motion->ref_idx [LIST_0][block_y    ][block_x], b_ref, 2 * sizeof(char));
        memset(&motion->ref_idx [LIST_0][block_y + 1][block_x], b_ref, 2 * sizeof(char));

        motion->ref_pic_id[LIST_0][block_y    ][block_x    ] = ref_pic_num;
        motion->ref_pic_id[LIST_0][block_y    ][block_x + 1] = ref_pic_num;
        motion->ref_pic_id[LIST_0][block_y + 1][block_x    ] = ref_pic_num;
        motion->ref_pic_id[LIST_0][block_y + 1][block_x + 1] = ref_pic_num;

        memcpy(motion->mv[LIST_0][block_y   ][block_x], currSlice->all_mv[LIST_0][(short) b_ref][mode][j1    ][i1], 4 * sizeof(short));
        memcpy(motion->mv[LIST_0][block_y + 1][block_x], currSlice->all_mv[LIST_0][(short) b_ref][mode][j1 + 1][i1], 4 * sizeof(short));

        if (currSlice->slice_type == B_SLICE)
        {
          //--- get cost and reference frame for LIST 1 prediction ---
          bmcost[LIST_1] = DISTBLK_MAX;
          bmcost[BI_PRED] = DISTBLK_MAX;
          list_prediction_cost(currMB, LIST_1, block, mode, enc_mb, bmcost, best.ref);

          // Compute bipredictive cost between best list 0 and best list 1 references
          list_prediction_cost(currMB, BI_PRED, block, mode, enc_mb, bmcost, best.ref);

          // currently Bi prediction ME is only supported for modes 1, 2, 3 and only for ref 0 and only for ref 0
          if (is_bipred_enabled(p_Inp, mode))
          {
            get_bipred_cost(currMB, mode, block, i1, j1, &best, enc_mb, bmcost);
          }
          else
          {
            bmcost[BI_PRED_L0] = DISTBLK_MAX;
            bmcost[BI_PRED_L1] = DISTBLK_MAX;
          }

          //--- get prediction direction ----
          determine_prediction_list(bmcost, &best, cost);

          //store backward reference index for every block
          for (k = LIST_0; k <= LIST_1; k++)
          {
            memset(&motion->ref_idx[k][block_y    ][block_x], best.ref[k], 2 * sizeof(char));
            memset(&motion->ref_idx[k][block_y + 1][block_x], best.ref[k], 2 * sizeof(char));

            if (best.bipred)
            {              
              memcpy(motion->mv[k][block_y    ][block_x], currSlice->bipred_mv[best.bipred - 1][k][(short) best.ref[k]][mode][j1    ][i1], 4 * sizeof(short));
              memcpy(motion->mv[k][block_y + 1][block_x], currSlice->bipred_mv[best.bipred - 1][k][(short) best.ref[k]][mode][j1 + 1][i1], 4 * sizeof(short));
            }
            else
            {
              memcpy(motion->mv[k][block_y    ][block_x], currSlice->all_mv[k][(short) best.ref[k]][mode][j1    ][i1], 4 * sizeof(short));
              memcpy(motion->mv[k][block_y + 1][block_x], currSlice->all_mv[k][(short) best.ref[k]][mode][j1 + 1][i1], 4 * sizeof(short));
            }
          }
        } // if (currSlice->slice_type == B_SLICE)
        else
        {
          best.pdir = 0;
          *cost     = bmcost[LIST_0];
        }
      } // if (mode!=0)

      //--- get and check rate-distortion cost ---
      rdcost = RDCost_for_8x8blocks (currMB, dataTr, &cnt_nonz, &curr_cbp_blk, enc_mb->lambda_mdfp, block, (short) mode, &best, min_rdcost);
      //--- set variables if best mode has changed ---
      if (rdcost < min_rdcost)
      {
        min_cost8x8              = *cost;
        min_rdcost               = rdcost;        
        *partition               = best;
        partition->mode          = (char) mode;
        currMB->b8x8[block].mode = (char) mode;

#ifdef BEST_NZ_COEFF
        if (cnt_nonz)
        {
          for(i = 0; i <= 1; i++)
          {
            best_nz_coeff[i][0]= p_Vid->nz_coeff[currMB->mbAddrX][i1 + i][j1    ];
            best_nz_coeff[i][1]= p_Vid->nz_coeff[currMB->mbAddrX][i1 + i][j1 + 1];
          }
        }
        else
        {
          for(i = 0; i <= 1; i++)
          {
            best_nz_coeff[i][0]= 0;
            best_nz_coeff[i][1]= 0;
          }
        }
#endif

        //--- store number of nonzero coefficients ---
        best_cnt_nonz  = cnt_nonz;

        //--- store block cbp ---
        dataTr->cbp_blk8x8 &= (~(0x33 << (((block>>1)<<3)+((block & 0x01)<<1)))); // delete bits for block
        dataTr->cbp_blk8x8 |= curr_cbp_blk;

        //--- store coefficients ---
        memcpy(&cofACtr[0][0][0][0],&currSlice->cofAC[block][0][0][0], 4 * 2 * 65 * sizeof(int));

        if( currSlice->P444_joined ) 
        {
          //--- store coefficients ---
          memcpy(&cofACtr[1][0][0][0],&currSlice->cofAC[block + 4][0][0][0], 4 * 2 * 65 * sizeof(int));
          memcpy(&cofACtr[2][0][0][0],&currSlice->cofAC[block + 8][0][0][0], 4 * 2 * 65 * sizeof(int));
        }

        //--- store reconstruction and prediction ---
        copy_image_data_8x8(&dataTr->rec_mbY8x8[j0], &p_Vid->enc_picture->imgY[currMB->pix_y + j0], i0, currMB->pix_x + i0);
        copy_image_data_8x8(&dataTr->mpr8x8[j0], &currSlice->mb_pred[0][j0], i0, i0);

        if (p_Inp->rdopt == 3)
        {
          errdo_store_best_block(p_Inp, p_Vid->p_decs->dec_mbY_best8x8[transform8x8], p_Vid->enc_picture->p_dec_img[0], i0, j0, currMB->pix_x + i0, currMB->pix_y, BLOCK_SIZE_8x8);
          errdo_store_best_block(p_Inp, p_Vid->p_decs->dec_mb_pred_best8x8[transform8x8], p_Vid->p_decs->dec_mb_pred, i0, j0, i0, 0, BLOCK_SIZE_8x8); 
        }

        if(currSlice->slice_type == SP_SLICE)
        {
          for (j = j0; j < j0 + BLOCK_SIZE_8x8; j++)
          {
            memcpy(&dataTr->lrec[j][i0],&p_Vid->lrec[currMB->pix_y + j][currMB->pix_x + i0], BLOCK_SIZE_8x8 * sizeof(int));
          }
        }

        if(currSlice->P444_joined) 
        {
          copy_image_data_8x8(&dataTr->rec_mb8x8_cr[0][j0], &p_Vid->enc_picture->imgUV[0][currMB->pix_y + j0], i0, currMB->pix_x + i0);
          copy_image_data_8x8(&dataTr->mpr8x8CbCr[0][j0], &currSlice->mb_pred[1][j0], i0, i0);

          copy_image_data_8x8(&dataTr->rec_mb8x8_cr[1][j0], &p_Vid->enc_picture->imgUV[1][currMB->pix_y + j0], i0, currMB->pix_x + i0);
          copy_image_data_8x8(&dataTr->mpr8x8CbCr[1][j0], &currSlice->mb_pred[2][j0], i0, i0);
        }


        //--- store best 8x8 coding state ---
        if (block < 3)
        {
          currSlice->store_coding_state (currMB, currSlice->p_RDO->cs_b8);
          stored_state_8x8 = TRUE;
        }
      } // if (rdcost <= min_rdcost)

      //--- re-set coding state as it was before coding with current mode was performed ---
      if (index != maxindex - 1)
      {
        if(p_Inp->subMBCodingState == 1)
          currSlice->reset_coding_state (currMB, currSlice->p_RDO->cs_cm);
        else if(p_Inp->subMBCodingState == 2)
          currSlice->reset_coding_state (currMB, currSlice->p_RDO->cs_tmp);
      }
    } // if ((enc_mb->valid[mode] && (transform8x8 == 0 || mode != 0 || (mode == 0 && active_sps->direct_8x8_inference_flag)))
  } // for (min_rdcost=1e30, index=(currSlice->slice_type == B_SLICE ? 0 : 1); index<6; index++)

  if (valid_8x8 == TRUE)
  {
#ifdef BEST_NZ_COEFF
    for(i = 0; i <= 1; i++)  
    {
      for(j = 0; j <= 1; j++)
        p_Vid->nz_coeff[currMB->mbAddrX][i1 + i][j1 + j] = best_nz_coeff[i][j];
    }
#endif

    if (!transform8x8)
    {
      if (min_cost8x8 != LLONG_MAX)
        dataTr->mb_p8x8_cost += min_cost8x8;
      else
        dataTr->mb_p8x8_cost = DISTBLK_MAX;
    }

    //----- set cbp and count of nonzero coefficients ---
    if (best_cnt_nonz)
    {
      dataTr->cbp8x8       |= (1 << block);
      dataTr->cnt_nonz_8x8 += best_cnt_nonz;
    }

    if (!transform8x8)
    {
      if (block < 3)
      {
        //===== re-set reconstructed block =====
        j0   = 8*(block >> 1);
        i0   = 8*(block & 0x01);

        // need to double check code since original did not use i0
        copy_image_data_8x8(&p_Vid->enc_picture->imgY[currMB->pix_y + j0], &dataTr->rec_mbY8x8[j0], currMB->pix_x + i0, i0);

        if (p_Inp->rdopt == 3)
        {
          errdo_get_best_block(currMB, p_Vid->enc_picture->p_dec_img[0], p_Vid->p_decs->dec_mbY_best8x8[transform8x8], j0, BLOCK_SIZE_8x8);
        }

        if(currSlice->slice_type == SP_SLICE)
        {
          for (j = j0; j < j0 + BLOCK_SIZE_8x8; j++)
          {
            memcpy(&p_Vid->lrec[currMB->pix_y + j][currMB->pix_x], dataTr->lrec[j], 2 * BLOCK_SIZE * sizeof(int)); // reset the coefficients for SP slice
          }
        }

        if(currSlice->P444_joined) 
        {

          for (k=0; k<2; k++)
          {
            copy_image_data_8x8(&p_Vid->enc_picture->imgUV[k][currMB->pix_y + j0], &dataTr->rec_mb8x8_cr[k][j0], currMB->pix_x + i0, i0);
          }
        }
      } // if (block<3)
    }
    else
    {
      //======= save motion data for 8x8 partition for transform size 8x8 ========    
      currSlice->store_8x8_motion_vectors(currSlice, 0, block, partition);
    }

    //===== set motion vectors and reference frames (prediction) =====    
    currSlice->set_ref_and_motion_vectors (currMB, motion, partition, block);

    //===== set the coding state after current block =====
    //if (transform8x8 == 0 || block < 3)
    if (stored_state_8x8 == TRUE)
      currSlice->reset_coding_state (currMB, currSlice->p_RDO->cs_b8);
    else
    {
      currSlice->reset_coding_state (currMB, currSlice->p_RDO->cs_cm);
      update_adaptive_rounding_8x8(p_Vid, p_Inp, dataTr, fadjust);
    }
  }
}

/*!
*************************************************************************************
* \brief
*    Low Complexity Mode Decision for an 8x8 sub-macroblock
*************************************************************************************
*/
void submacroblock_mode_decision_low(Macroblock *currMB,
                                     RD_PARAMS *enc_mb,
                                     RD_8x8DATA *dataTr,                                     
                                     int ****cofACtr,
                                     int *have_direct,
                                     int block,
                                     distblk *cost_direct,
                                     distblk *cost,
                                     distblk *cost8x8_direct,
                                     int transform8x8)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;
  Slice *currSlice = currMB->p_slice;
  PicMotionParams *motion = &p_Vid->enc_picture->motion;

  int64 curr_cbp_blk;
  double min_rdcost, rdcost = 0.0;
  int j0, i0, j1, i1;
  int i,j, k;
  int index;
  int mode;
  distblk min_cost8x8;
  distblk direct4x4_tmp, direct8x8_tmp;
  distblk bmcost[5] = {DISTBLK_MAX};
  int cnt_nonz = 0;
  int dummy;
  int best_cnt_nonz = 0;
  int maxindex =  (transform8x8) ? 2 : 5;
  int block_x, block_y;
  int lambda_mf[3];

  int ****fadjust = transform8x8? p_Vid->ARCofAdj8x8 : p_Vid->ARCofAdj4x4;
  short pdir;

  Boolean valid_8x8 = FALSE;
  Boolean stored_state_8x8 = FALSE;

#ifdef BEST_NZ_COEFF
  int best_nz_coeff[2][2];
#endif

  Info8x8 best;
  // Init best (need to create simple function)

  best.pdir = 0;
  best.bipred = 0;
  best.ref[LIST_0] = 0;
  best.ref[LIST_1] = -1;
  
  dataTr->part[block].mode = 0;

  //--- set coordinates ---
  j0 = ((block>>1)<<3);
  j1 = (j0>>2);
  i0 = ((block&0x01)<<3);
  i1 = (i0>>2);

#ifdef BEST_NZ_COEFF
  for(j = 0; j <= 1; j++)
  {
    for(i = 0; i <= 1; i++)
      best_nz_coeff[i][j] = p_Vid->nz_coeff[currMB->mbAddrX][i1 + i][j1 + j] = 0;
  }
#endif

  if (transform8x8)
    currMB->luma_transform_size_8x8_flag = TRUE; //switch to transform size 8x8

  //--- store coding state before coding ---
  currSlice->store_coding_state (currMB, currSlice->p_RDO->cs_cm);

  //=====  LOOP OVER POSSIBLE CODING MODES FOR 8x8 SUB-PARTITION  =====
  for (min_cost8x8 = DISTBLK_MAX, min_rdcost = 1e20, index = (currSlice->slice_type == B_SLICE ? 0 : 1); index < maxindex; index++)
  {
    mode = b8_mode_table[index];
    best.mode = (char) mode;
    *cost = 0;

    if (enc_mb->valid[mode] && (!(transform8x8 == 1 && mode > 4)) && (transform8x8 == 0 || mode != 0 || (mode == 0 && p_Vid->active_sps->direct_8x8_inference_flag)))
    {
      if (transform8x8)
      {
        currMB->valid_8x8 = TRUE;
      }

      valid_8x8 = TRUE;
      curr_cbp_blk = 0;

      if (mode==0)
      {
        //--- Direct Mode ---       
        direct4x4_tmp = 0;
        direct8x8_tmp = 0;
        direct4x4_tmp = GetDirectCost8x8 ( currMB, block, &direct8x8_tmp);
        if ((direct4x4_tmp==DISTBLK_MAX)||(*cost_direct==DISTBLK_MAX))
        {
          *cost_direct = DISTBLK_MAX;
          if (transform8x8)
            *cost8x8_direct = DISTBLK_MAX;
        }
        else
        {
          *cost_direct += direct4x4_tmp;
          if (transform8x8)
            *cost8x8_direct += direct8x8_tmp;
        }
        (*have_direct) ++;

        if (transform8x8)
        {
          switch(p_Inp->Transform8x8Mode)
          {
          case 1: // Mixture of 8x8 & 4x4 transform
            if((direct8x8_tmp < direct4x4_tmp) || !(enc_mb->valid[5] && enc_mb->valid[6] && enc_mb->valid[7]))
              *cost = direct8x8_tmp;
            else
              *cost = direct4x4_tmp;
            break;
          case 2: // 8x8 Transform only
            *cost = direct8x8_tmp;
            break;
          default: // 4x4 Transform only
            *cost = direct4x4_tmp;
            break;
          }
          if (p_Inp->Transform8x8Mode==2)
            *cost = DISTBLK_MAX;
        }
        else
        {
          *cost = direct4x4_tmp;
        }

        block_x = currMB->block_x + (block & 0x01)*2;
        block_y = currMB->block_y + (block & 0x02);
        best.ref[LIST_0] = currSlice->direct_ref_idx[block_y][block_x][LIST_0];
        best.ref[LIST_1] = currSlice->direct_ref_idx[block_y][block_x][LIST_1];
        best.pdir        = currSlice->direct_pdir[block_y][block_x];
      } // if (mode==0)
      else
      {
        int64 ref_pic_num;
        char b_ref;

        //======= motion estimation for all reference frames ========
        //-----------------------------------------------------------
        memcpy(lambda_mf, enc_mb->lambda_mf, 3 * sizeof(int));
        if (p_Inp->CtxAdptLagrangeMult == 1)
        {
          RDOPTStructure *p_RDO = currSlice->p_RDO;
          lambda_mf[F_PEL] = (int)(lambda_mf[F_PEL] * p_RDO->lambda_mf_factor);
          lambda_mf[H_PEL] = (int)(lambda_mf[H_PEL] * p_RDO->lambda_mf_factor);
          lambda_mf[Q_PEL] = (int)(lambda_mf[Q_PEL] * p_RDO->lambda_mf_factor);
        }

        SubPartitionMotionSearch (currMB, mode, block, lambda_mf);

        //--- get cost and reference frame for LIST 0 prediction ---
        bmcost[LIST_0] = DISTBLK_MAX;
        list_prediction_cost(currMB, LIST_0, block, mode, enc_mb, bmcost, best.ref);

        //store LIST 0 reference index for every block
        block_x = currMB->block_x + (block & 0x01)*2;
        block_y = currMB->block_y + (block & 0x02);
        b_ref = best.ref[LIST_0];
        ref_pic_num = p_Vid->enc_picture->ref_pic_num[currMB->list_offset][(short) b_ref];

        for (j = block_y; j< block_y + 2; j++)
        {
          memset(&motion->ref_idx [LIST_0][j][block_x], b_ref, 2 * sizeof(char));
        }

        for (j = block_y; j< block_y + 2; j++)
        {
          for (i = block_x; i < block_x + 2; i++)
          {
            motion->ref_pic_id[LIST_0][j][i] = ref_pic_num;
          }
        }

        if (currSlice->slice_type == B_SLICE)
        {
          //--- get cost and reference frame for LIST 1 prediction ---
          bmcost[LIST_1] = DISTBLK_MAX;
          bmcost[BI_PRED] = DISTBLK_MAX;
          list_prediction_cost(currMB, LIST_1, block, mode, enc_mb, bmcost, best.ref);

          // Compute bipredictive cost between best list 0 and best list 1 references
          list_prediction_cost(currMB, BI_PRED, block, mode, enc_mb, bmcost, best.ref);

          // currently Bi prediction ME is only supported for modes 1, 2, 3 and only for ref 0 and only for ref 0
          if (is_bipred_enabled(p_Inp, mode))
          {
            list_prediction_cost(currMB, BI_PRED_L0, block, mode, enc_mb, bmcost, 0);
            list_prediction_cost(currMB, BI_PRED_L1, block, mode, enc_mb, bmcost, 0);
          }
          else
          {
            bmcost[BI_PRED_L0] = DISTBLK_MAX;
            bmcost[BI_PRED_L1] = DISTBLK_MAX;
          }

          //--- get prediction direction ----
          determine_prediction_list(bmcost, &best, cost);

          //store backward reference index for every block
          for (k = LIST_0; k <= LIST_1; k++)
          {
            for (j = block_y; j< block_y + 2; j++)
            {
              memset(&motion->ref_idx[k][j][block_x], best.ref[k], 2 * sizeof(char));
            }
          }
        } // if (currSlice->slice_type == B_SLICE)
        else
        {
          best.pdir = 0;
          *cost     = bmcost[LIST_0];
        }
      } // if (mode!=0)
      if (*cost != DISTBLK_MAX)
        *cost += (ref_cost(p_Vid, enc_mb->lambda_mf[Q_PEL], (short) B8Mode2Value (currSlice, (short) mode, best.pdir), 
        (best.pdir < 1 ? currMB->list_offset : currMB->list_offset + LIST_1)) - 1);

      //--- set variables if best mode has changed ---
      if (*cost < min_cost8x8)
      {
        min_cost8x8             = *cost;
        min_rdcost               = rdcost;
        
        dataTr->part[block] = best;
        currMB->b8x8[block].mode = (char) mode;

#ifdef BEST_NZ_COEFF
        if (cnt_nonz)
        {
          best_nz_coeff[0][0]= p_Vid->nz_coeff[currMB->mbAddrX][i1    ][j1    ];
          best_nz_coeff[0][1]= p_Vid->nz_coeff[currMB->mbAddrX][i1    ][j1 + 1];
          best_nz_coeff[1][0]= p_Vid->nz_coeff[currMB->mbAddrX][i1 + 1][j1    ];
          best_nz_coeff[1][1]= p_Vid->nz_coeff[currMB->mbAddrX][i1 + 1][j1 + 1];
        }
        else
        {
          best_nz_coeff[0][0]= 0;
          best_nz_coeff[0][1]= 0;
          best_nz_coeff[1][0]= 0;
          best_nz_coeff[1][1]= 0;
        }
#endif

        //--- store number of nonzero coefficients ---
        best_cnt_nonz  = cnt_nonz;                

        //--- store best 8x8 coding state ---
        if (block < 3)
        {
          currSlice->store_coding_state (currMB, currSlice->p_RDO->cs_b8);
          stored_state_8x8 = TRUE;
        }
      } // if (rdcost <= min_rdcost)

      //--- re-set coding state as it was before coding with current mode was performed ---
      currSlice->reset_coding_state (currMB, currSlice->p_RDO->cs_cm);
    } // if ((enc_mb->valid[mode] && (transform8x8 == 0 || mode != 0 || (mode == 0 && active_sps->direct_8x8_inference_flag)))
  } // for (min_rdcost=1e30, index=(currSlice->slice_type == B_SLICE ? 0 : 1); index<6; index++)

  if (valid_8x8 == TRUE)
  {
#ifdef BEST_NZ_COEFF
    for(i = 0; i <= 1; i++)  
    {
      for(j = 0; j <= 1; j++)
        p_Vid->nz_coeff[currMB->mbAddrX][i1 + i][j1 + j] = best_nz_coeff[i][j];
    }
#endif

    if (!transform8x8)
    {
      if (min_cost8x8 != DISTBLK_MAX)
        dataTr->mb_p8x8_cost += min_cost8x8;
      else
        dataTr->mb_p8x8_cost = DISTBLK_MAX;
    }


    if (transform8x8)
    {
      if (min_cost8x8 != DISTBLK_MAX)
        dataTr->mb_p8x8_cost += min_cost8x8;
      else
        dataTr->mb_p8x8_cost = DISTBLK_MAX;

      mode = dataTr->part[block].mode;
      pdir = dataTr->part[block].pdir;
    }
    else
    {
      mode = dataTr->part[block].mode;
      pdir = dataTr->part[block].pdir;
    }

    curr_cbp_blk  = 0;
    currMB->b8x8[block].bipred = dataTr->part[block].bipred;
    currMB->ar_mode = (short) ((mode != 0)? mode: P8x8);

    if (currSlice->P444_joined)
    {
      int list_mode[2];
      list_mode[0] = (pdir == 0 || pdir == 2 ? mode : 0);
      list_mode[1] = (pdir == 1 || pdir == 2 ? mode : 0);

      best_cnt_nonz = luma_residual_coding_p444_8x8 (currMB, &dummy, &curr_cbp_blk, block, pdir, list_mode, dataTr->part[block].ref);
      best_cnt_nonz += currSlice->coeff_cost_cr[1] + currSlice->coeff_cost_cr[2];
    }
    else
    {
      int list_mode[2];
      list_mode[0] = (pdir == 0 || pdir == 2 ? mode : 0);
      list_mode[1] = (pdir == 1 || pdir == 2 ? mode : 0);
      best_cnt_nonz = luma_residual_coding_8x8 (currMB, &dummy, &curr_cbp_blk, block, pdir, list_mode, dataTr->part[block].ref);
    }


    dataTr->cbp_blk8x8   &= (~(0x33 << (((block>>1)<<3)+((block & 0x01)<<1)))); // delete bits for block
    dataTr->cbp_blk8x8   |= curr_cbp_blk;

    //--- store coefficients ---
    memcpy(cofACtr[0][0][0],currSlice->cofAC[block][0][0], 4 * 2 * 65 * sizeof(int));

    if(currSlice->P444_joined) 
    {
      //--- store coefficients ---
      memcpy(cofACtr[1][0][0],currSlice->cofAC[block + 4][0][0], 4 * 2 * 65 * sizeof(int));
      memcpy(cofACtr[2][0][0],currSlice->cofAC[block + 8][0][0], 4 * 2 * 65 * sizeof(int));
    }


    //--- store reconstruction and prediction ---
    copy_image_data_8x8(&dataTr->rec_mbY8x8[j0], &p_Vid->enc_picture->imgY[currMB->pix_y + j0], i0, currMB->pix_x + i0);
    copy_image_data_8x8(&dataTr->mpr8x8[j0], &currSlice->mb_pred[0][j0], i0, i0);


    //--- store reconstruction and prediction ---
    if(currSlice->slice_type == SP_SLICE)
    {
      for (j=j0; j < j0 + BLOCK_SIZE_8x8; j++)
      {
        memcpy(&dataTr->lrec[j][i0],&p_Vid->lrec[currMB->pix_y+j][currMB->pix_x+i0],BLOCK_SIZE_8x8 * sizeof(int)); // store coefficients for primary SP slice
      }
    }
    if(currSlice->P444_joined) 
    {
      copy_image_data_8x8(&dataTr->rec_mb8x8_cr[0][j0], &p_Vid->enc_picture->imgUV[0][currMB->pix_y + j0], i0, currMB->pix_x + i0);
      copy_image_data_8x8(&dataTr->mpr8x8CbCr[0][j0], &currSlice->mb_pred[1][j0], i0, i0);
      copy_image_data_8x8(&dataTr->rec_mb8x8_cr[1][j0], &p_Vid->enc_picture->imgUV[1][currMB->pix_y + j0], i0, currMB->pix_x + i0);
      copy_image_data_8x8(&dataTr->mpr8x8CbCr[1][j0], &currSlice->mb_pred[2][j0], i0, i0);
    }   


    //----- set cbp and count of nonzero coefficients ---
    if (best_cnt_nonz)
    {
      dataTr->cbp8x8       |= (1 << block);
      dataTr->cnt_nonz_8x8 += best_cnt_nonz;
    }

    if (!transform8x8)
    {
      if (block < 3)
      {
        //===== re-set reconstructed block =====
        j0   = 8*(block >> 1);
        i0   = 8*(block & 0x01);

        copy_image_data_8x8(&p_Vid->enc_picture->imgY[currMB->pix_y + j0], &dataTr->rec_mbY8x8[j0], currMB->pix_x + i0, i0);

        if(currSlice->slice_type == SP_SLICE)
        {
          for (j = j0; j < j0 + BLOCK_SIZE_8x8; j++)
          {
            memcpy(&p_Vid->lrec[currMB->pix_y + j][currMB->pix_x], dataTr->lrec[j], 2 * BLOCK_SIZE * sizeof(int)); // reset the coefficients for SP slice
          }
        }

        if(currSlice->P444_joined) 
        {
          copy_image_data_8x8(&p_Vid->enc_picture->imgUV[0][currMB->pix_y + j0], &dataTr->rec_mb8x8_cr[0][j0], currMB->pix_x + i0, i0);
          copy_image_data_8x8(&p_Vid->enc_picture->imgUV[1][currMB->pix_y + j0], &dataTr->rec_mb8x8_cr[1][j0], currMB->pix_x + i0, i0);
        }
      } // if (block<3)
    }
    else
    {
      //======= save motion data for 8x8 partition for transform size 8x8 ========    
      currSlice->store_8x8_motion_vectors(currSlice, 0, block, &dataTr->part[block]);
    }

    //===== set motion vectors and reference frames (prediction) =====
    currSlice->set_ref_and_motion_vectors (currMB, motion, &dataTr->part[block], block);

    //===== set the coding state after current block =====
    //if (transform8x8 == 0 || block < 3)
    if (stored_state_8x8 == TRUE)
      currSlice->reset_coding_state (currMB, currSlice->p_RDO->cs_b8);
    else
    {
      update_adaptive_rounding_8x8(p_Vid, p_Inp, dataTr, fadjust);
    }
  }
}


void get_initial_mb16x16_cost(Macroblock* currMB)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  RDOPTStructure *p_RDO = currMB->p_slice->p_RDO;
  if (currMB->mb_left && currMB->mb_up)
  {
    p_Vid->mb16x16_cost = (p_Vid->mb16x16_cost_frame[currMB->mbAddrX - 1] +
      p_Vid->mb16x16_cost_frame[currMB->mbAddrX - (p_Vid->width>>4)] + 1)/2.0;
  }
  else if (currMB->mb_left)
  {
    p_Vid->mb16x16_cost = p_Vid->mb16x16_cost_frame[currMB->mbAddrX - 1];
  }
  else if (currMB->mb_up)
  {
    p_Vid->mb16x16_cost = p_Vid->mb16x16_cost_frame[currMB->mbAddrX - (p_Vid->width>>4)];
  }
  else
  {
    p_Vid->mb16x16_cost = CALM_MF_FACTOR_THRESHOLD;
  }

  p_RDO->lambda_mf_factor = p_Vid->mb16x16_cost < CALM_MF_FACTOR_THRESHOLD ? 1.0 : sqrt(p_Vid->mb16x16_cost / (CALM_MF_FACTOR_THRESHOLD * p_Vid->lambda_mf_factor[p_Vid->type][p_Vid->qp]));
}

void adjust_mb16x16_cost(Macroblock *currMB, distblk cost)
{
  VideoParameters *p_Vid = currMB->p_Vid;
  RDOPTStructure  *p_RDO = currMB->p_slice->p_RDO;

  p_Vid->mb16x16_cost = (double) cost;
  p_Vid->mb16x16_cost_frame[p_Vid->current_mb_nr] = p_Vid->mb16x16_cost;
#if JCOST_CALC_SCALEUP
  p_RDO->lambda_mf_factor = (p_Vid->mb16x16_cost < CALM_MF_FACTOR_THRESHOLD*(1<<LAMBDA_ACCURACY_BITS))
  ? 1.0
  : sqrt(p_Vid->mb16x16_cost / (CALM_MF_FACTOR_THRESHOLD*(1<<LAMBDA_ACCURACY_BITS) * p_Vid->lambda_mf_factor[p_Vid->type][p_Vid->qp]));
#else
  p_RDO->lambda_mf_factor = (p_Vid->mb16x16_cost < CALM_MF_FACTOR_THRESHOLD)
  ? 1.0
  : sqrt(p_Vid->mb16x16_cost / (CALM_MF_FACTOR_THRESHOLD * p_Vid->lambda_mf_factor[p_Vid->type][p_Vid->qp]));
#endif
}

void update_lambda_costs(Macroblock *currMB, RD_PARAMS *enc_mb, int lambda_mf[3])
{
  InputParameters *p_Inp = currMB->p_Inp;
  RDOPTStructure  *p_RDO = currMB->p_slice->p_RDO;

  int MEPos;
  for (MEPos = 0; MEPos < 3; MEPos ++)
  {
    lambda_mf[MEPos] = p_Inp->CtxAdptLagrangeMult == 0 ? enc_mb->lambda_mf[MEPos] : (int)(enc_mb->lambda_mf[MEPos] * sqrt(p_RDO->lambda_mf_factor));
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Return array's minimum and its index
 *************************************************************************************
 */
int iminarray ( int arr[], int size, int *minind )
{
  int i; 
  int mincand = arr[0];
  *minind = 0;
  for ( i = 1; i < size; i++ )
  {
    if (arr[i] < mincand)
    {
      mincand = arr[i];
      *minind = i;
    }
  }
  return mincand;
} 

distblk distblkminarray ( distblk arr[], int size, int *minind )
{
  int i; 
  distblk mincand = arr[0];
  *minind = 0;
  for ( i = 1; i < size; i++ )
  {
    if (arr[i] < mincand)
    {
      mincand = arr[i];
      *minind = i;
    }
  }
  return mincand;
} 

/*!
 *************************************************************************************
 * \brief
 *    Determines whether bi prediction is enabaled for current mode
 *************************************************************************************
 */
int is_bipred_enabled(InputParameters *p_Inp, int mode) 
{
  int enabled = 0;
  mode = (mode == P8x8) ? 4: mode;

  if (p_Inp->BiPredMotionEstimation)
  {
    if (mode > 0 && mode < 5)
    {
      enabled = (p_Inp->BiPredSearch[mode - 1]) ? 1: 0;
    }    
    else
    {
      enabled = 0;
    }
  }
  else
  {
    enabled = 0;
  }
  return enabled;
}


/*!
 *************************************************************************************
 * \brief
 *    Decides whether to perform tranform 8x8 for this mode
 *************************************************************************************
 */
int transform_termination_control(Macroblock* currMB, int mode) 
{  
  VideoParameters *p_Vid = currMB->p_Vid;
  InputParameters *p_Inp = currMB->p_Inp;

  // Go through transform modes.
  // Note that if currMB->cbp is 0 one could choose to skip 8x8 mode
  // although this could be due to deadzoning decisions.
  //if (p_Inp->Transform8x8Mode==1 && currMB->cbp!=0)
  if (p_Inp->Transform8x8Mode == 1)
  {
    int bslice = currMB->p_slice->slice_type == B_SLICE;
    //=========== try the 8x8 transform with mb_types 16x16,16x8, 8x16, 8x8, and DIRECT 16x16 ===========
    if (currMB->luma_transform_size_8x8_flag == FALSE && 
      ((mode >= 1 && mode <= 3) || (bslice && mode == 0 && p_Vid->active_sps->direct_8x8_inference_flag) || (mode == P8x8)))
        //if (currMB->luma_transform_size_8x8_flag == FALSE && 
      //((mode >= 1 && mode <= 3) || (bslice && mode == 0 && active_sps->direct_8x8_inference_flag)))
    {
      //try with 8x8 transform size
      currMB->luma_transform_size_8x8_flag = TRUE;
      return 0;
    }
    else
    {
      currMB->luma_transform_size_8x8_flag = FALSE;
      return 1;
    }
  }
  else
  {
    return 1;
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Update prediction direction for mode P16x16 to check all prediction directions
 *************************************************************************************
 */
int bslice_16x16_termination_control(InputParameters *p_Inp, Block8x8Info *b8x8info, int *ctr16x16, int mode, int bslice)
{
  int lastcheck = 1;
  //--- for INTER16x16 in BSLICEs check all prediction directions ---
  if (mode == 1 && bslice)
  {
    char pdir = 0;
    short i;
    char bipred_me = 0;
 
    switch (*ctr16x16)
    {
    case 0:
      pdir = 0;
      lastcheck = 0;
      break;
    case 1:
      pdir = 1;
      lastcheck = 0;
      break;
    case 2:
      pdir = 2;
      if (p_Inp->BiPredMotionEstimation)
      {
        lastcheck = 0;
      }
      break;
    case 3:
      pdir = 2;
      bipred_me = 1;
      lastcheck = 0;
      break;
    case 4:
      pdir = 2;
      bipred_me = 2;
      break;
    default:
      error("invalid 'ctr16x16' value", -1);
      break;
    }
    for (i = 0; i< 4; i++)
    {
      b8x8info->best[1][i].bipred = bipred_me;
      b8x8info->best[1][i].pdir = pdir;
    }
    (*ctr16x16)++;
  }
  return lastcheck;
}

/*!
**************************************************************************************
* \brief
*     Compute bipred costs
**************************************************************************************
*/
void get_bipred_cost(Macroblock *currMB, int mode, int block, int i, int j, Info8x8  *best, RD_PARAMS *enc_mb, distblk bmcost[5])
{
  Slice *currSlice = currMB->p_slice;
  short   *bi0_mv_l0    = currSlice->bipred_mv[0][LIST_0][0][mode][j][i];
  short   *bi0_mv_l1    = currSlice->bipred_mv[0][LIST_1][0][mode][j][i];
  short   *bi1_mv_l0    = currSlice->bipred_mv[1][LIST_0][0][mode][j][i];
  short   *bi1_mv_l1    = currSlice->bipred_mv[1][LIST_1][0][mode][j][i];
  distblk  MaxVal = DISTBLK_MAX;
  if (best->ref[0] == 0 && best->ref[1] == 0)
  {
    short   *single_mv_l0 = currSlice->all_mv [LIST_0][0][mode][j][i];
    short   *single_mv_l1 = currSlice->all_mv [LIST_1][0][mode][j][i];

    if ((single_mv_l0[0] != bi0_mv_l0[0]) || (single_mv_l0[1] != bi0_mv_l0[1]) ||
      (single_mv_l1[0] != bi0_mv_l1[0]) || (single_mv_l1[1] != bi0_mv_l1[1]))
      list_prediction_cost(currMB, BI_PRED_L0, block, mode, enc_mb, bmcost, 0);
    else
      bmcost[BI_PRED_L0] = MaxVal;

    if ((single_mv_l0[0] != bi1_mv_l0[0]) || (single_mv_l0[1] != bi1_mv_l0[1]) ||
      (single_mv_l1[0] != bi1_mv_l1[0]) || (single_mv_l1[1] != bi1_mv_l1[1]))
    {
      if ((bi0_mv_l0[0] != bi1_mv_l0[0]) || (bi0_mv_l0[1] != bi1_mv_l0[1]) ||
        (bi0_mv_l1[0] != bi1_mv_l1[0]) || (bi0_mv_l1[1] != bi1_mv_l1[1]))
        list_prediction_cost(currMB, BI_PRED_L1, block, mode, enc_mb, bmcost, 0);
      else
        bmcost[BI_PRED_L1] = MaxVal;
    }
    else
      bmcost[BI_PRED_L1] = MaxVal;
  }
  else
  {
    list_prediction_cost(currMB, BI_PRED_L0, block, mode, enc_mb, bmcost, 0);
    if ((bi0_mv_l0[0] != bi1_mv_l0[0]) || (bi0_mv_l0[1] != bi1_mv_l0[1]) ||
      (bi0_mv_l1[0] != bi1_mv_l1[0]) || (bi0_mv_l1[1] != bi1_mv_l1[1]))
      list_prediction_cost(currMB, BI_PRED_L1, block, mode, enc_mb, bmcost, 0);
    else
      bmcost[BI_PRED_L1] = MaxVal;
  }
}




