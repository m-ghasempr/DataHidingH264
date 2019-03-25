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
 * \file header.c
 *
 * \brief
 *    H.26L Slice headers
 *
 *************************************************************************************
 */

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "global.h"
#include "elements.h"
#include "defines.h"
#include "fmo.h"
#include "vlc.h"
#include "mbuffer.h"


#if TRACE
#define SYMTRACESTRING(s) strncpy(sym.tracestring,s,TRACESTRING_SIZE)
#else
#define SYMTRACESTRING(s) // to nothing
#endif

extern int UsedBits;

static void ref_pic_list_reordering();
static void pred_weight_table();
static void dec_ref_pic_marking();


/*!
 ************************************************************************
 * \brief
 *    read the first part of the header (only the pic_parameter_set_id)
 * \return
 *    Length of the first part of the slice header (in bits)
 ************************************************************************
 */
int FirstPartOfSliceHeader()
{
  Slice *currSlice = img->currentSlice;
  int dP_nr = assignSE2partition[currSlice->dp_mode][SE_HEADER];
  DataPartition *partition = &(currSlice->partArr[dP_nr]);

  UsedBits= partition->bitstream->frame_bitoffset; // was hardcoded to 31 for previous start-code. This is better.

  // Get first_mb_in_slice
  currSlice->start_mb_nr = ue_v ("SH: first_mb_in_slice", partition);

  currSlice->picture_type = ue_v ("SH: slice_type", partition);

  currSlice->pic_parameter_set_id = ue_v ("SH: pic_parameter_set_id", partition);
  
  return UsedBits;
}

/*!
 ************************************************************************
 * \brief
 *    read the scond part of the header (without the pic_parameter_set_id 
 * \return
 *    Length of the second part of the Slice header in bits
 ************************************************************************
 */
int RestOfSliceHeader()
{
  Slice *currSlice = img->currentSlice;
  int dP_nr = assignSE2partition[currSlice->dp_mode][SE_HEADER];
  DataPartition *partition = &(currSlice->partArr[dP_nr]);
//  int UsedBits= partition->bitstream->frame_bitoffset; // was hardcoded to 31 for previous start-code. This is better.
  static int last_imgtr_frm=0,modulo_ctr_frm=0,last_imgtr_fld=0,modulo_ctr_fld=0;
  static int last_imgtr_frm_b=0,modulo_ctr_frm_b=0,last_imgtr_fld_b=0,modulo_ctr_fld_b=0;
  static int FirstCall = 1;   // used for FMO initialization
  static unsigned int AbsFrameNum, ExpectedPicOrderCnt, PicOrderCntCycleCnt, FrameNumInPicOrderCntCycle;
  static unsigned int PreviousFrameNum, FrameNumOffset, ExpectedDeltaPerPicOrderCntCycle;
  static unsigned int Previousfield_pic_flag,Previousbottom_field_flag,Previousnal_reference_idc;
  static unsigned int Previousdelta_pic_order_cnt[2], PreviousPOC=11111, ThisPOC, FirstFieldType;
  int i,val;

  int multpred=active_sps->num_ref_frames>1;
  
  // set the img->type to the old style defines
  switch (currSlice->picture_type)
  {
  case 0:
  case 5:
    if (multpred) currSlice->picture_type = INTER_IMG_MULT;
    else currSlice->picture_type = INTER_IMG_1;
    break;
  case 1:
  case 6:
    if (multpred) currSlice->picture_type = B_IMG_MULT;
    else currSlice->picture_type = B_IMG_1;
    break;
  case 2:
  case 7:
    currSlice->picture_type = INTRA_IMG;
    break;
  case 3:
  case 8:
    if (multpred) currSlice->picture_type = SP_IMG_MULT;
    else currSlice->picture_type = SP_IMG_1;
    break;
  default:
    error("Invalid Slice Type", 500);
    break;
  }

  img->type = currSlice->picture_type;

  img->frame_num = u_v (active_sps->log2_max_frame_num_minus4 + 4, "SH: frame_num", partition);

  if (active_sps->frame_mbs_only_flag)
  {
    img->structure = FRAME;
    img->field_pic_flag=0;
  }
  else
  {
    // field_pic_flag   u(1)
    img->field_pic_flag = u_1("SH: field_pic_flag", partition);
    if (img->field_pic_flag)
    {
      // bottom_field_flag  u(1)
      img->bottom_field_flag = u_1("SH: bottom_field_flag", partition);

      img->structure = img->bottom_field_flag ? BOTTOM_FIELD : TOP_FIELD;
    }
    else
    {
      img->structure = FRAME;
      img->bottom_field_flag=0;
    }
  }

  currSlice->structure = img->structure;

  img->mb_frame_field_flag=(active_sps->mb_adaptive_frame_field_flag && (img->field_pic_flag==0));

  if (img->structure == 0) assert (img->field_pic_flag == 0);
  if (img->structure == 1) assert (img->field_pic_flag == 1 && img->bottom_field_flag == 0);
  if (img->structure == 2) assert (img->field_pic_flag == 1 && img->bottom_field_flag == 1);

  if (img->idr_flag)
  {
    img->idr_pic_id = ue_v("SH: idr_pic_id", partition);
  }

  if(!img->delta_pic_order_always_zero_flag)
  {
    img->delta_pic_order_cnt[0] = se_v("SH: delta_pic_order_cnt[0]", partition);
  }
  else 
    img->delta_pic_order_cnt[0] = 0;
                        
  if (img->pic_order_present_flag  && !img->delta_pic_order_always_zero_flag)
  {
    img->delta_pic_order_cnt[1] = se_v("SH: delta_pic_order_cnt[1]", partition);
  }
  else 
    img->delta_pic_order_cnt[1] = 0;
  
  //! redundant_pic_cnt is missing here
  if (active_pps->redundant_pic_cnt_present_flag)
  {
    img->redundant_pic_cnt = u_1 ("SH: redundant_pic_cnt", partition);
  }

  if(img->type==B_IMG_1 || img->type==B_IMG_MULT)
  {
    img->direct_type = u_1 ("SH: direct_spatial_mv_pred_flag", partition);
  }

  img->num_ref_pic_active_fwd = 0;
  img->num_ref_pic_active_bwd = 0;

  val = u_1 ("SH: num_ref_idx_override_flag", partition);
  if (val)
  {
    if(img->type==INTER_IMG_1 || img->type==INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT || img->type==B_IMG_1 || img->type==B_IMG_MULT)
    {
      img->num_ref_pic_active_fwd = 1 + ue_v ("SH: num_ref_pic_active_fwd_minus1", partition);
      
      if(img->type==B_IMG_1 || img->type==B_IMG_MULT)
      {
        img->num_ref_pic_active_bwd = 1 + ue_v ("SH: num_ref_pic_active_bwd_minus1", partition);
      }
    }
  }

  ref_pic_list_reordering();

  img->apply_weights = ((img->weighted_pred_flag && (currSlice->picture_type == INTER_IMG_1 || currSlice->picture_type == INTER_IMG_MULT || currSlice->picture_type == SP_IMG_1 || currSlice->picture_type == SP_IMG_MULT) )
          || ((img->weighted_bipred_explicit_flag  || img->weighted_bipred_implicit_flag ) && (currSlice->picture_type == B_IMG_1 || currSlice->picture_type == B_IMG_MULT )));

  if ((active_pps->weighted_pred_flag&&(img->type==INTER_IMG_1 || img->type==INTER_IMG_MULT || img->type == SP_IMG_1 || img->type == SP_IMG_MULT))||
      (active_pps->weighted_bipred_idc==1 && (img->type==B_IMG_1 || img->type==B_IMG_MULT)))
  {
    pred_weight_table();
  }

  dec_ref_pic_marking();

  if (active_pps->entropy_coding_mode && img->type!=INTRA_IMG && img->type!=SP_IMG_1 && img->type!=SP_IMG_MULT)
  {
    currSlice->cabac_init_idc = ue_v("SH: cabac_init_idc", partition);
  }
  else 
  {
    currSlice->cabac_init_idc=0;
  }

  val = se_v("SH: slice_qp_delta", partition);
  currSlice->qp = img->qp = 26 + active_pps->pic_init_qp_minus26 + val;

  if(img->type==SP_IMG_1 || img->type==SP_IMG_MULT || img->type == SI_IMG) 
  {
    if(img->type==SP_IMG_1 || img->type==SP_IMG_MULT)
    {
      img->sp_switch = u_1 ("SH: sp_for_switch_flag", partition);
    }
    val = se_v("SH: slice_qs_delta", partition);
    img->qpsp = 26 + active_pps->pic_init_qs_minus26 + val;
  }

  if (active_pps->deblocking_filter_parameters_present_flag)
  {
    currSlice->LFDisableIdc = ue_v ("SH: disable_deblocking_filter_idc", partition);

    if (currSlice->LFDisableIdc!=1)
    {
      currSlice->LFAlphaC0Offset = 2 * se_v("SH: slice_alpha_c0_offset_div2", partition);

      currSlice->LFBetaOffset = 2 * se_v("SH: slice_beta_offset_div2", partition);
    }
  }

  if (active_pps->num_slice_groups_minus1>0 && active_pps->mb_slice_group_map_type>=3 &&
      active_pps->mb_slice_group_map_type<=5)
  {
    assert ("FMO not yet completely supported");
    assert (0==1);
  }
  // 5. Finally, read Reference Picture ID (same as TR here).  Note that this is an
  // optional field that is not present if the input parameters do not indicate
  // multiframe prediction ??

  // Ok, the above comment is nonsense.  There is no way how a decoder could
  // know that we use multiple reference frames (except probably through a
  // sequence header).  Hence, it's now an if (1) -- PHRefPicID is always present

  // Of course, the decoder can know. It is indicated by the previously decoded
  // parameter "PHPictureType". So, I changed the if-statement again to be
  // compatible with the encoder.

  // WYK: Oct. 16, 2001. Now I use this for the reference frame ID (non-B frame ID). 
  // Thus, we can know how many  non-B frames are lost, and then we can adjust 
  // the reference frame buffers correctly.
/*  if (1)
  {
    sym.type = SE_HEADER;                 // This will be true for all symbols generated here
    // refPicID, variable length
    sym.mapping = linfo;               // change to unsigned integer
    SYMTRACESTRING("SH RefPicID");
    readSyntaxElement_UVLC (&sym,img,input,partition);
    if (img->refPicID != sym.value1)
    {
      img->refPicID_old = img->refPicID;
      img->refPicID = sym.value1;
    }
    UsedBits += sym.len;
  }
*/

  img->buf_cycle = input->buf_cycle+1;
  img->pn=(((img->structure==BOTTOM_FIELD) ? (img->number/2):img->number)%img->buf_cycle);

  img->max_mb_nr = 2*(img->width * img->height) / (MB_BLOCK_SIZE * MB_BLOCK_SIZE);

  //calculate pocs
  if(img->idr_flag)
  {
    FrameNumOffset=0;     //  first pix of IDRGOP, 
    FirstFieldType = img->bottom_field_flag;              //save type of first field of frame
                                                          //NB may not work with mixed field & frame coding
    img->delta_pic_order_cnt[0]=0;                        //ignore first delta
    assert (img->frame_num == 0);
    if(img->frame_num)error("frame_num != 0 in idr pix", -1020);

  }
  else if (img->frame_num<PreviousFrameNum)
  {             //not first pix of IDRGOP
    FrameNumOffset += img->MaxFrameNum;
  }

  if(img->num_ref_frames_in_pic_order_cnt_cycle) 
    AbsFrameNum = FrameNumOffset+img->frame_num;
  else 
    AbsFrameNum=0;

  if(img->disposable_flag && AbsFrameNum)
    AbsFrameNum--;

  //! calculate ExpectedDeltaPerPicOrderCntCycle (used to be done during the reading
  //! of the SPS symbols (which usd to live in the slice header)

  ExpectedDeltaPerPicOrderCntCycle=0;
  if(active_sps->num_ref_frames_in_pic_order_cnt_cycle)
    for(i=0;i<(int) active_sps->num_ref_frames_in_pic_order_cnt_cycle;i++)
    {
      ExpectedDeltaPerPicOrderCntCycle += active_sps->offset_for_ref_frame[i];
    }


  if(AbsFrameNum)
  {
    PicOrderCntCycleCnt = (AbsFrameNum-1)/img->num_ref_frames_in_pic_order_cnt_cycle;
    FrameNumInPicOrderCntCycle = (AbsFrameNum-1)%img->num_ref_frames_in_pic_order_cnt_cycle;
    ExpectedPicOrderCnt = PicOrderCntCycleCnt*ExpectedDeltaPerPicOrderCntCycle;
    for(i=0;i<=(int)FrameNumInPicOrderCntCycle;i++)
      ExpectedPicOrderCnt += img->offset_for_ref_frame[i];
  }
  else 
    ExpectedPicOrderCnt=0;

  if(img->disposable_flag)
    ExpectedPicOrderCnt += img->offset_for_non_ref_pic;


  if(img->field_pic_flag==0)
  {           //frame pix
    ThisPOC = img->toppoc = ExpectedPicOrderCnt + img->delta_pic_order_cnt[0];
    img->bottompoc = img->toppoc + img->offset_for_top_to_bottom_field + img->delta_pic_order_cnt[1];
    if(PreviousPOC!=ThisPOC)
    {               //new frame detected
      if(img->disposable_flag)
        push_poc(img->toppoc,img->bottompoc,NONREFFRAME);
      else                                    
        push_poc(img->toppoc,img->bottompoc,REFFRAME);
    }
  }
  else if (img->bottom_field_flag==0)
  {  //top field 
    ThisPOC = img->toppoc = ExpectedPicOrderCnt + img->delta_pic_order_cnt[0];
    if(PreviousPOC!=ThisPOC  &&  FirstFieldType==img->bottom_field_flag)
    {           //new frame detected
      if(img->disposable_flag)
        push_poc(img->toppoc,0,NONREFFRAME);
      else                                    
        push_poc(img->toppoc,0,REFFRAME);
    }
    else 
      toprefpoc[0] = img->toppoc;                //2nd field of same frame
  } 
  else
  {                                                   //bottom field
    ThisPOC = img->bottompoc = ExpectedPicOrderCnt + img->offset_for_top_to_bottom_field + img->delta_pic_order_cnt[0];
    if(PreviousPOC!=ThisPOC  &&  FirstFieldType==img->bottom_field_flag)
    {           //new frame detected
      if(img->disposable_flag)
        push_poc(0,img->bottompoc,NONREFFRAME);
      else 
        push_poc(0,img->bottompoc,REFFRAME);
    }
    else 
      bottomrefpoc[0] = img->bottompoc;          //2nd field of same frame
  }

                                                                //temp stuff to track tr
  if(img->field_pic_flag)
  {
    img->tr_fld = ThisPOC;
    if(img->bottom_field_flag)
    {
      currSlice->picture_id = img->tr = img->bottompoc%256;
    }
    else
    {   //top field
      currSlice->picture_id = img->tr = img->toppoc%256;
    }
  }
  else
  {           //frame pix  -  use toppoc/2
    img->tr_frm = ThisPOC/2;
    currSlice->picture_id = img->tr = (img->toppoc/2)%256;
  }
            //update "Previous" stuff for next slice
  PreviousFrameNum=img->frame_num;
  Previousfield_pic_flag=img->field_pic_flag;
  Previousbottom_field_flag=img->bottom_field_flag;
  Previousnal_reference_idc=img->idr_flag;               
  Previousdelta_pic_order_cnt[0]=img->delta_pic_order_cnt[0];
  Previousdelta_pic_order_cnt[1]=img->delta_pic_order_cnt[1];
  PreviousPOC=ThisPOC;
            //moved from above for stuff that still uses img->tr
            //soon to be obsolete
  if(!img->current_slice_nr)
  { 
    if((img->type != B_IMG_MULT && img->type != B_IMG_1) || !img->disposable_flag) 
    {
      img->pstruct_next_P = img->structure;
      if(img->structure == TOP_FIELD)
      {
        img->imgtr_last_P = img->imgtr_next_P;
        img->imgtr_next_P = img->tr_fld;
      }
      else if(img->structure == FRAME)
      {
        img->imgtr_last_P = img->imgtr_next_P;
        img->imgtr_next_P = 2*img->tr_frm;
      }
    }
    
  }
  //note  UsedBits is probably inaccurate
//  dumppoc (img);
  return UsedBits;
}


/*!
 ************************************************************************
 * \brief
 *    read the reference picture reordering information
 ************************************************************************
 */
static void ref_pic_list_reordering()
{
  Slice *currSlice = img->currentSlice;
  int dP_nr = assignSE2partition[currSlice->dp_mode][SE_HEADER];
  DataPartition *partition = &(currSlice->partArr[dP_nr]);
  int i, val;

  alloc_ref_pic_list_reordering_buffer(currSlice);
  
  if (img->type!=INTRA_IMG && img->type!=SI_IMG)
  {
    val = currSlice->ref_pic_list_reordering_flag_l0 = u_1 ("SH: ref_pic_list_reordering_flag_l0", partition);
    
    if (val)
    {
      i=0;
      do
      {
        val = currSlice->remapping_of_pic_nums_idc_l0[i] = ue_v("SH: remapping_of_pic_nums_idc_l0", partition);
        if (val==0 || val==1)
        {
          currSlice->abs_diff_pic_num_minus1_l0[i] = ue_v("SH: abs_diff_pic_num_minus1_l0", partition);
        }
        else
        {
          if (val==2)
          {
            currSlice->long_term_pic_idx_l0[i] = ue_v("SH: long_term_pic_idx_l0", partition);
          }
        }
        i++;
        // assert (i>img->num_ref_pic_active_fwd);
      } while (val != 3);
    }
  }

  if (img->type==B_IMG_1 || img->type==B_IMG_MULT)
  {
    val = currSlice->ref_pic_list_reordering_flag_l1 = u_1 ("SH: ref_pic_list_reordering_flag_l1", partition);
    
    if (val)
    {
      i=0;
      do
      {
        val = currSlice->remapping_of_pic_nums_idc_l1[i] = ue_v("SH: remapping_of_pic_nums_idc_l1", partition);
        if (val==0 || val==1)
        {
          currSlice->abs_diff_pic_num_minus1_l1[i] = ue_v("SH: abs_diff_pic_num_minus1_l1", partition);
        }
        else
        {
          if (val==2)
          {
            currSlice->long_term_pic_idx_l1[i] = ue_v("SH: long_term_pic_idx_l1", partition);
          }
        }
        i++;
        // assert (i>img->num_ref_pic_active_bwd);
      } while (val != 3);
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    read the weighted prediction tables
 ************************************************************************
 */
static void pred_weight_table()
{
  Slice *currSlice = img->currentSlice;
  int dP_nr = assignSE2partition[currSlice->dp_mode][SE_HEADER];
  DataPartition *partition = &(currSlice->partArr[dP_nr]);
  int luma_weight_flag_l0, luma_weight_flag_l1, chroma_weight_flag_l0, chroma_weight_flag_l1;
  int i,j;

  img->luma_log_weight_denom = ue_v ("SH: luma_log_weight_denom", partition);
  img->wp_round_luma = 1<<(img->luma_log_weight_denom - 1);
  
  img->chroma_log_weight_denom = ue_v ("SH: chroma_log_weight_denom", partition);
  img->wp_round_chroma = 1<<(img->chroma_log_weight_denom - 1);

  reset_wp_params(img);

  for (i=0; i<img->num_ref_pic_active_fwd; i++)
  {
    luma_weight_flag_l0 = u_1("SH: luma_weight_flag_l0", partition);
    
    if (luma_weight_flag_l0)
    {
      img->wp_weight[0][i][0] = se_v ("SH: luma_weight_l0", partition);
      img->wp_offset[0][i][0] = se_v ("SH: luma_offset_l0", partition);
    }
    else
    {
      img->wp_weight[0][i][0] = 1<<img->luma_log_weight_denom;
      img->wp_offset[0][i][0] = 0;
    }
    
    chroma_weight_flag_l0 = u_1 ("SH: chroma_weight_flag_l0", partition);
    
    for (j=1; j<3; j++)
    {
      if (chroma_weight_flag_l0)
      {
        img->wp_weight[0][i][j] = se_v("SH: chroma_weight_l0", partition);
        img->wp_offset[0][i][j] = se_v("SH: chroma_offset_l0", partition);
      }
      else
      {
        img->wp_weight[0][i][j] = 1<<img->luma_log_weight_denom;
        img->wp_offset[0][i][j] = 0;
      }
    }
  }
  if ((img->type == B_IMG_1 || img->type == B_IMG_MULT) && img->weighted_bipred_explicit_flag == 1)
  {
    for (i=0; i<img->num_ref_pic_active_bwd; i++)
    {
      luma_weight_flag_l1 = u_1("SH: luma_weight_flag_l1", partition);
      
      if (luma_weight_flag_l1)
      {
        img->wp_weight[1][i][0] = se_v ("SH: luma_weight_l1", partition);
        img->wp_offset[1][i][0] = se_v ("SH: luma_offset_l1", partition);
      }
      else
      {
        img->wp_weight[1][i][0] = 1<<img->luma_log_weight_denom;
        img->wp_offset[1][i][0] = 0;
      }
      
      chroma_weight_flag_l1 = u_1 ("SH: chroma_weight_flag_l1", partition);
      
      for (j=1; j<3; j++)
      {
        if (chroma_weight_flag_l1)
        {
          img->wp_weight[1][i][j] = se_v("SH: chroma_weight_l1", partition);
          img->wp_offset[1][i][j] = se_v("SH: chroma_offset_l1", partition);
        }
        else
        {
          img->wp_weight[1][i][j] = 1<<img->luma_log_weight_denom;
          img->wp_offset[1][i][j] = 0;
        }
      }
    }
  }    
}


/*!
 ************************************************************************
 * \brief
 *    read the memory control operations
 ************************************************************************
 */
static void dec_ref_pic_marking()
{
  Slice *currSlice = img->currentSlice;
  int dP_nr = assignSE2partition[currSlice->dp_mode][SE_HEADER];
  DataPartition *partition = &(currSlice->partArr[dP_nr]);
  int val;

  MMCObuffer_t *tmp_mmco,*tmp_mmco2;

  // free old buffer content
  while (img->mmco_buffer)
  { 
    tmp_mmco=img->mmco_buffer;

    img->mmco_buffer=tmp_mmco->Next;
    free (tmp_mmco);
  } 

  if (img->idr_flag)
  {
    img->no_output_of_prior_pics_flag = u_1("SH: no_output_of_prior_pics_flag", partition);
    img->long_term_reference_flag = u_1("SH: long_term_reference_flag", partition);
  }
  else
  {
    img->adaptive_ref_pic_buffering_flag = u_1("SH: adaptive_ref_pic_buffering_flag", partition);
    if (img->adaptive_ref_pic_buffering_flag)
    {
      // read Memory Management Control Operation 
      do
      {
        tmp_mmco=(MMCObuffer_t*)calloc (1,sizeof (MMCObuffer_t));
        tmp_mmco->Next=NULL;
        
        val = tmp_mmco->MMCO = ue_v("SH: memory_management_control_operation", partition);

        if ((val==1)||(val==3)) 
        {
          tmp_mmco->DPN = 1 + ue_v("SH: difference_of_pic_nums_minus1", partition);
        }
        if ((val==2)||(val==3)||(val==6))
        {
          tmp_mmco->LPIN = ue_v("SH: long_term_pic_idx", partition);
        }
        if (val==4)
        {
          tmp_mmco->MLIP1 = ue_v("SH: max_long_term_pic_idx_plus1", partition);
        }
        
        // add MMCO to list
        if (img->mmco_buffer==NULL) 
        {
          img->mmco_buffer=tmp_mmco;
        }
        else
        {
          tmp_mmco2=img->mmco_buffer;
          while (tmp_mmco2->Next!=NULL) tmp_mmco2=tmp_mmco2->Next;
          tmp_mmco2->Next=tmp_mmco;
        }
        
      }while (val != 0);
      
    }
  }
}



/*!
 ************************************************************************
 * \brief
 *    A little helper for the debugging of POC code
 * \return
 *    none
 ************************************************************************
 */
int dumppoc(struct img_par *img) {
    printf ("\nPOC locals...\n");
    printf ("toppoc                                %d\n", img->toppoc);
    printf ("bottompoc                             %d\n", img->bottompoc);
    printf ("frame_num                             %d\n", img->frame_num);
    printf ("field_pic_flag                        %d\n", img->field_pic_flag);
    printf ("bottom_field_flag                     %d\n", img->bottom_field_flag);
    printf ("POC SPS\n");
    printf ("log2_max_frame_num_minus4             %d\n", 4);
    printf ("pic_order_cnt_type                    %d\n", img->pic_order_cnt_type);
    printf ("num_ref_frames_in_pic_order_cnt_cycle %d\n", img->num_ref_frames_in_pic_order_cnt_cycle);
    printf ("delta_pic_order_always_zero_flag      %d\n", img->delta_pic_order_always_zero_flag);
    printf ("offset_for_non_ref_pic                %d\n", img->offset_for_non_ref_pic);
    printf ("offset_for_top_to_bottom_field        %d\n", img->offset_for_top_to_bottom_field);
    printf ("offset_for_ref_frame[0]               %d\n", img->offset_for_ref_frame[0]);
    printf ("offset_for_ref_frame[1]               %d\n", img->offset_for_ref_frame[1]);
    printf ("POC in SLice Header\n");
    printf ("pic_order_present_flag                %d\n", img->pic_order_present_flag);
    printf ("delta_pic_order_cnt[0]                %d\n", img->delta_pic_order_cnt[0]);
    printf ("delta_pic_order_cnt[1]                %d\n", img->delta_pic_order_cnt[1]);
    printf ("delta_pic_order_cnt[2]                %d\n", img->delta_pic_order_cnt[2]);
    printf ("idr_flag                              %d\n", img->idr_flag);
    printf ("MaxFrameNum                           %d\n", img->MaxFrameNum);
    return 0;
}
