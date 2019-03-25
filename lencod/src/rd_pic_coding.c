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

/*!
 ************************************************************************
 * \brief
 *    performs multi-pass encoding of same picture using different
 *    coding conditions
 ************************************************************************
 */

void rd_picture_coding(ImageParameters *p_Img, InputParameters *p_Inp)
{
  int   second_qp = p_Img->qp, rd_qp = p_Img->qp;
  int   previntras = p_Img->intras;
  int   prevtype = p_Img->type;
  int   skip_encode = 0;
  pic_parameter_set_rbsp_t *sec_pps;
  int   tmpFrameQP = p_Img->SumFrameQP;
  int   num_ref_idx_l0 = p_Img->num_ref_idx_l0_active;
  int   num_ref_idx_l1 = p_Img->num_ref_idx_l1_active;
  float rateRatio = 1.0F;

  if ( p_Inp->RCEnable )
    rc_save_state(p_Img, p_Inp);

  if ( p_Inp->WPMCPrecision )
    p_Img->pWPX->curr_wp_rd_pass = p_Img->pWPX->wp_rd_passes + 1;

  if (p_Img->type!=I_SLICE && p_Inp->GenerateMultiplePPS)
  {
    if ( p_Inp->WPMCPrecision )
    {
      p_Img->pWPX->curr_wp_rd_pass->algorithm = WP_REGULAR;
      if (p_Img->type==P_SLICE)
      {
        p_Img->active_pps = p_Img->PicParSet[1];
      }
      else
      {
        p_Img->active_pps = p_Img->PicParSet[2];
      }
    }
    else
    {
      if (p_Img->type==P_SLICE)
      {
        if ((p_Inp->RDPSliceWeightOnly != 2) && (p_Img->TestWPPSlice(p_Img, p_Inp, 0) == 1))
        {
          p_Img->active_pps = p_Img->PicParSet[1];
        }
        else
        {
          skip_encode = p_Inp->RDPSliceWeightOnly;
          p_Img->active_pps = p_Img->PicParSet[0];

          if (!p_Img->AdaptiveRounding)
          {
            p_Img->qp-=1;
            if ( p_Inp->RCEnable )
              rateRatio = 1.15F;
          }
        }
      }
      else
      {
        p_Img->active_pps = p_Img->PicParSet[2];
      }
    }
  }
  else
  {
    if (!p_Img->AdaptiveRounding)
    {
      p_Img->qp-=1;
      if ( p_Inp->RCEnable )
        rateRatio = 1.15F;
    }
  }

  sec_pps = p_Img->active_pps;
  second_qp = p_Img->qp;

  p_Img->write_macroblock = FALSE;

  if (skip_encode)
  {
    p_Img->rd_pass = 0;
    p_Img->enc_frame_picture[1] = NULL;
  }
  else
  {
    if(p_Inp->RCEnable)
      rc_init_frame_rdpic( p_Img, p_Inp, rateRatio );

    p_Img->qp = iClip3( p_Img->RCMinQP, p_Img->RCMaxQP, p_Img->qp );
    frame_picture (p_Img, p_Inp, p_Img->frame_pic[1], &p_Img->imgData, 1);
    p_Img->rd_pass = picture_coding_decision(p_Img, p_Inp, p_Img->frame_pic[0], p_Img->frame_pic[1], rd_qp);
  }

  if (p_Img->rd_pass==0)
  {
    p_Img->enc_picture=p_Img->enc_frame_picture[0];
    if (p_Img->type!=I_SLICE && p_Inp->GenerateMultiplePPS)
    {
      p_Img->qp = rd_qp;
      p_Img->active_pps = p_Img->PicParSet[0];
    }
    else
    {
      p_Img->qp = rd_qp;
    }
    p_Img->intras = previntras;
    p_Img->p_frame_pic = p_Img->frame_pic[0];
  }
  else
  {
    previntras  = p_Img->intras;
    p_Img->p_frame_pic = p_Img->frame_pic[1];
    tmpFrameQP  = p_Img->SumFrameQP;
    num_ref_idx_l0 = p_Img->num_ref_idx_l0_active;
    num_ref_idx_l1 = p_Img->num_ref_idx_l1_active;


    if(p_Inp->RCEnable)
      rc_save_state(p_Img, p_Inp);
  }
  // Final Encoding pass - note that we should
  // make this more flexible in a later version.
  if ( p_Inp->RCEnable )
    rateRatio = 1.0F;

  if ( p_Inp->WPMCPrecision )
    p_Img->pWPX->curr_wp_rd_pass = p_Img->pWPX->wp_rd_passes + 2;

  if (p_Img->type != I_SLICE )
  {
    skip_encode = 0;
    p_Img->qp    = rd_qp;

    if (p_Img->type == P_SLICE && (p_Img->intras * 100 )/p_Img->FrameSizeInMbs >=75)
    {
      set_slice_type(p_Img, p_Inp, I_SLICE );
      p_Img->active_pps = p_Img->PicParSet[0];
    }
    else if (p_Img->type==P_SLICE)
    {
      if (p_Inp->GenerateMultiplePPS)
      {
        if ((p_Inp->RDPSliceWeightOnly != 2) && (p_Img->TestWPPSlice(p_Img, p_Inp, 1) == 1))
        {
          p_Img->active_pps = p_Img->PicParSet[1];
          if ( p_Inp->WPMCPrecision )
            p_Img->pWPX->curr_wp_rd_pass->algorithm = WP_REGULAR;
        }
        else if ( p_Inp->WPMCPrecision == 2 )
          p_Img->active_pps = p_Img->PicParSet[1];
        else if (p_Inp->RDPSliceBTest && p_Img->active_sps->profile_idc != BASELINE)
        {
          set_slice_type(p_Img,  p_Inp, B_SLICE );
          p_Img->active_pps = p_Img->PicParSet[0];
        }
        else
        {
          skip_encode = p_Inp->RDPSliceWeightOnly;
          p_Img->active_pps = p_Img->PicParSet[0];
          if (!p_Img->AdaptiveRounding)
          {
            p_Img->qp+=1;
            if ( p_Inp->RCEnable )
              rateRatio = 0.85F;
          }
        }
      }
    }
    else
    {
      if (p_Inp->GenerateMultiplePPS && (p_Inp->RDBSliceWeightOnly != 2) && p_Img->TestWPBSlice(p_Img, p_Inp, 0) == 1)
      {
        p_Img->active_pps = p_Img->PicParSet[1];
        if ( p_Inp->WPMCPrecision )
          p_Img->pWPX->curr_wp_rd_pass->algorithm = WP_REGULAR;
      }
      else if ( p_Inp->WPMCPrecision == 2 && (p_Inp->WPMCPrecBSlice == 2 || (p_Inp->WPMCPrecBSlice == 1 && p_Img->nal_reference_idc) ) )
        p_Img->active_pps = p_Img->PicParSet[1];
      else
      {
        skip_encode = (p_Inp->RDBSliceWeightOnly == 1);
        p_Img->qp = rd_qp + (p_Img->nal_reference_idc ? - 1 : 1);
        if ( p_Inp->RCEnable )
          rateRatio = p_Img->nal_reference_idc ? 1.15F : 0.85F;
        if ( p_Inp->WPMCPrecision )
          p_Img->pWPX->curr_wp_rd_pass->algorithm = WP_REGULAR;
      }
    }
  }
  else
  {
    p_Img->active_pps = p_Img->PicParSet[0];
    if (!p_Img->AdaptiveRounding)
      p_Img->qp    = (rd_qp + 1);
  }

  p_Img->write_macroblock = FALSE;

  if (skip_encode)
  {
    p_Img->enc_frame_picture[2] = NULL;
    p_Img->qp = rd_qp;
  }
  else
  {
    if(p_Inp->RCEnable)
      rc_init_frame_rdpic( p_Img, p_Inp, rateRatio );

    p_Img->qp = iClip3( p_Img->RCMinQP, p_Img->RCMaxQP, p_Img->qp );
    frame_picture (p_Img, p_Inp, p_Img->frame_pic[2], &p_Img->imgData, 2);

    if (p_Img->rd_pass==0)
      p_Img->rd_pass  = 2 * picture_coding_decision(p_Img, p_Inp, p_Img->frame_pic[0], p_Img->frame_pic[2], rd_qp);
    else
      p_Img->rd_pass +=     picture_coding_decision(p_Img, p_Inp, p_Img->frame_pic[1], p_Img->frame_pic[2], rd_qp);

    if ( p_Inp->RCEnable && p_Img->rd_pass == 2 )
      rc_save_state(p_Img, p_Inp);

    if ( p_Img->rd_pass == 2 )
    {
      tmpFrameQP = p_Img->SumFrameQP;
      num_ref_idx_l0 = p_Img->num_ref_idx_l0_active;
      num_ref_idx_l1 = p_Img->num_ref_idx_l1_active;
    }
  }

  if (p_Img->rd_pass==0)
  {
    p_Img->enc_picture = p_Img->enc_frame_picture[0];
    set_slice_type( p_Img, p_Inp, prevtype );
    p_Img->active_pps  = p_Img->PicParSet[0];
    p_Img->qp     = rd_qp;
    p_Img->intras = previntras;
  }
  else if (p_Img->rd_pass==1)
  {
    p_Img->enc_picture = p_Img->enc_frame_picture[1];
    set_slice_type( p_Img, p_Inp, prevtype );
    p_Img->active_pps  = sec_pps;
    p_Img->qp     = second_qp;
    p_Img->intras = previntras;
  }

  if ( p_Inp->RCEnable )
    rc_restore_state(p_Img, p_Inp);

  p_Img->SumFrameQP = tmpFrameQP;
  p_Img->num_ref_idx_l0_active = num_ref_idx_l0;
  p_Img->num_ref_idx_l1_active = num_ref_idx_l1;
}

