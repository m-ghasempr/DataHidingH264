/*!
*************************************************************************************
* \file rd_pic_coding.c
*
* \brief
*    Rate Distortion Picture level multiple pass encoding functions
*
* \author
*    Main contributors (see contributors.h for copyright, address and affiliation details)
*     - Alexis Michael Tourapis         <alexismt@ieee.org>
*************************************************************************************
*/

#include <math.h>

#include "global.h"
#include "image.h"
#include "rc_quadratic.h"
#include "wp.h"
#include "pred_struct.h"

/*!
 ************************************************************************
 * \brief
 *    performs multi-pass encoding of same picture using different
 *    coding conditions
 ************************************************************************
 */

void rd_picture_coding(VideoParameters *p_Vid, InputParameters *p_Inp)
{
  int   second_qp = p_Vid->p_curr_frm_struct->qp, rd_qp = p_Vid->p_curr_frm_struct->qp;
  int   previntras = p_Vid->intras;
  int   prevtype = p_Vid->type;
  int   skip_encode = 0;
  pic_parameter_set_rbsp_t *sec_pps;
  int   tmpFrameQP = p_Vid->SumFrameQP;
  int   num_ref_idx_l0 = p_Vid->num_ref_idx_l0_active;
  int   num_ref_idx_l1 = p_Vid->num_ref_idx_l1_active;
  float rateRatio = 1.0F;

  if ( p_Inp->RCEnable )
    rc_save_state(p_Vid, p_Inp);

  if ( p_Inp->WPMCPrecision )
    p_Vid->pWPX->curr_wp_rd_pass = p_Vid->pWPX->wp_rd_passes + 1;

  if (p_Vid->type!=I_SLICE && p_Inp->GenerateMultiplePPS)
  {
    if (p_Vid->type==P_SLICE)
    {
      if ((p_Inp->RDPSliceWeightOnly != 2) && (p_Vid->TestWPPSlice(p_Vid, 0) == 1))
      {
        p_Vid->active_pps = p_Vid->PicParSet[1];
        if ( p_Inp->WPMCPrecision )
        {
          p_Vid->pWPX->curr_wp_rd_pass->algorithm = WP_REGULAR;
        }
      }
      else if ( p_Inp->WPMCPrecision )
      {
        p_Vid->active_pps = p_Vid->PicParSet[1];
      }
      else
      {
        skip_encode = p_Inp->RDPSliceWeightOnly;
        p_Vid->active_pps = p_Vid->PicParSet[0];

        if (!p_Vid->AdaptiveRounding)
        {
          p_Vid->qp-=1;
          if ( p_Inp->RCEnable )
            rateRatio = 1.15F;
        }
      }
    }
    else
    {
      p_Vid->active_pps = p_Vid->PicParSet[2];
      if ( p_Inp->WPMCPrecision )
      {
        p_Vid->pWPX->curr_wp_rd_pass->algorithm = WP_REGULAR;
      }
    }
  }
  else
  {
    if (!p_Vid->AdaptiveRounding)
    {
      p_Vid->qp-=1;
      if ( p_Inp->RCEnable )
        rateRatio = 1.15F;
    }
  }

  sec_pps = p_Vid->active_pps;
  second_qp = p_Vid->qp;

  p_Vid->write_macroblock = FALSE;

  if (skip_encode)
  {
    p_Vid->rd_pass = 0;
    p_Vid->enc_frame_picture[1] = NULL;
  }
  else
  {
    if(p_Inp->RCEnable)
      rc_init_frame_rdpic( p_Vid, p_Inp, rateRatio );

    p_Vid->qp = iClip3( p_Vid->RCMinQP, p_Vid->RCMaxQP, p_Vid->qp );
    p_Vid->p_curr_frm_struct->qp = p_Vid->qp;
    frame_picture (p_Vid, p_Vid->frame_pic[1], &p_Vid->imgData, 1);
    p_Vid->rd_pass = picture_coding_decision(p_Vid, p_Vid->frame_pic[0], p_Vid->frame_pic[1], rd_qp);
  }

  if (p_Vid->rd_pass==0)
  {
    p_Vid->enc_picture=p_Vid->enc_frame_picture[0];
    if (p_Vid->type!=I_SLICE && p_Inp->GenerateMultiplePPS)
    {
      p_Vid->qp = rd_qp;
      p_Vid->active_pps = p_Vid->PicParSet[0];
    }
    else
    {
      p_Vid->qp = rd_qp;
    }
    p_Vid->intras = previntras;
    p_Vid->p_frame_pic = p_Vid->frame_pic[0];
  }
  else
  {
    previntras  = p_Vid->intras;
    p_Vid->p_frame_pic = p_Vid->frame_pic[1];
    tmpFrameQP  = p_Vid->SumFrameQP;
    num_ref_idx_l0 = p_Vid->num_ref_idx_l0_active;
    num_ref_idx_l1 = p_Vid->num_ref_idx_l1_active;


    if(p_Inp->RCEnable)
      rc_save_state(p_Vid, p_Inp);
  }
  // Final Encoding pass - note that we should
  // make this more flexible in a later version.
  if ( p_Inp->RCEnable )
    rateRatio = 1.0F;

  if ( p_Inp->WPMCPrecision )
    p_Vid->pWPX->curr_wp_rd_pass = p_Vid->pWPX->wp_rd_passes + 2;

  if (p_Vid->type != I_SLICE )
  {
    skip_encode = 0;
    p_Vid->qp    = rd_qp;

    if (p_Vid->type == P_SLICE && (p_Vid->intras * 100 )/p_Vid->FrameSizeInMbs >=75)
    {
      set_slice_type(p_Vid, p_Inp, I_SLICE );
      populate_frame_slice_type( p_Inp, p_Vid->p_curr_frm_struct, I_SLICE, p_Vid->p_pred->max_num_slices );
      p_Vid->active_pps = p_Vid->PicParSet[0];
    }
    else if (p_Vid->type==P_SLICE)
    {
      if (p_Inp->GenerateMultiplePPS)
      {
        if ((p_Inp->RDPSliceWeightOnly != 2) && (p_Vid->TestWPPSlice(p_Vid, 1) == 1))
        {
          p_Vid->active_pps = p_Vid->PicParSet[1];
          if ( p_Inp->WPMCPrecision )
            p_Vid->pWPX->curr_wp_rd_pass->algorithm = WP_REGULAR;
        }
        else if ( p_Inp->WPMCPrecision == 2 )
          p_Vid->active_pps = p_Vid->PicParSet[1];
        else if (p_Inp->RDPSliceBTest && p_Vid->active_sps->profile_idc != BASELINE)
        {
          set_slice_type(p_Vid, p_Inp, B_SLICE );
          populate_frame_slice_type( p_Inp, p_Vid->p_curr_frm_struct, B_SLICE, p_Vid->p_pred->max_num_slices );
          p_Vid->active_pps = p_Vid->PicParSet[0];
        }
        else
        {
          skip_encode = p_Inp->RDPSliceWeightOnly;
          p_Vid->active_pps = p_Vid->PicParSet[0];
          if (!p_Vid->AdaptiveRounding)
          {
            p_Vid->qp+=1;
            if ( p_Inp->RCEnable )
              rateRatio = 0.85F;
          }
        }
      }
    }
    else
    {
      if (p_Inp->GenerateMultiplePPS && (p_Inp->RDBSliceWeightOnly != 2) && p_Vid->TestWPBSlice(p_Vid, 0) == 1)
      {
        p_Vid->active_pps = p_Vid->PicParSet[1];
        if ( p_Inp->WPMCPrecision )
          p_Vid->pWPX->curr_wp_rd_pass->algorithm = WP_REGULAR;
      }
      else if ( p_Inp->WPMCPrecision == 2 && (p_Inp->WPMCPrecBSlice == 2 || (p_Inp->WPMCPrecBSlice == 1 && p_Vid->nal_reference_idc) ) )
        p_Vid->active_pps = p_Vid->PicParSet[1];
      else
      {
        skip_encode = (p_Inp->RDBSliceWeightOnly == 1);
        p_Vid->qp = rd_qp + (p_Vid->nal_reference_idc ? - 1 : 1);
        if ( p_Inp->RCEnable )
          rateRatio = p_Vid->nal_reference_idc ? 1.15F : 0.85F;
        if ( p_Inp->WPMCPrecision )
          p_Vid->pWPX->curr_wp_rd_pass->algorithm = WP_REGULAR;
      }
    }
  }
  else
  {
    p_Vid->active_pps = p_Vid->PicParSet[0];
    if (!p_Vid->AdaptiveRounding)
      p_Vid->qp    = (rd_qp + 1);
  }

  p_Vid->write_macroblock = FALSE;

  if (skip_encode)
  {
    p_Vid->enc_frame_picture[2] = NULL;
    p_Vid->qp = rd_qp;
  }
  else
  {
    if(p_Inp->RCEnable)
      rc_init_frame_rdpic( p_Vid, p_Inp, rateRatio );

    p_Vid->qp = iClip3( p_Vid->RCMinQP, p_Vid->RCMaxQP, p_Vid->qp );
    p_Vid->p_curr_frm_struct->qp = p_Vid->qp;
    frame_picture (p_Vid, p_Vid->frame_pic[2], &p_Vid->imgData, 2);

    if (p_Vid->rd_pass==0)
      p_Vid->rd_pass  = 2 * picture_coding_decision(p_Vid, p_Vid->frame_pic[0], p_Vid->frame_pic[2], rd_qp);
    else
      p_Vid->rd_pass +=     picture_coding_decision(p_Vid, p_Vid->frame_pic[1], p_Vid->frame_pic[2], rd_qp);

    if ( p_Inp->RCEnable && p_Vid->rd_pass == 2 )
      rc_save_state(p_Vid, p_Inp);

    if ( p_Vid->rd_pass == 2 )
    {
      tmpFrameQP = p_Vid->SumFrameQP;
      num_ref_idx_l0 = p_Vid->num_ref_idx_l0_active;
      num_ref_idx_l1 = p_Vid->num_ref_idx_l1_active;
    }
  }

  if (p_Vid->rd_pass==0)
  {
    p_Vid->enc_picture = p_Vid->enc_frame_picture[0];
    set_slice_type( p_Vid, p_Inp, prevtype );
    populate_frame_slice_type( p_Inp, p_Vid->p_curr_frm_struct, prevtype, p_Vid->p_pred->max_num_slices );
    p_Vid->active_pps  = p_Vid->PicParSet[0];
    p_Vid->qp     = rd_qp;
    p_Vid->intras = previntras;
  }
  else if (p_Vid->rd_pass==1)
  {
    p_Vid->enc_picture = p_Vid->enc_frame_picture[1];
    set_slice_type( p_Vid, p_Inp, prevtype );
    populate_frame_slice_type( p_Inp, p_Vid->p_curr_frm_struct, prevtype, p_Vid->p_pred->max_num_slices );
    p_Vid->active_pps  = sec_pps;
    p_Vid->qp     = second_qp;
    p_Vid->intras = previntras;
  }

  if ( p_Inp->RCEnable )
    rc_restore_state(p_Vid, p_Inp);

  p_Vid->SumFrameQP = tmpFrameQP;
  p_Vid->num_ref_idx_l0_active = num_ref_idx_l0;
  p_Vid->num_ref_idx_l1_active = num_ref_idx_l1;
  p_Vid->p_curr_frm_struct->qp = p_Vid->qp;
}

