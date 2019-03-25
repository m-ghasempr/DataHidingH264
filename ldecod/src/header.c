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
#include <stdio.h>

#include "global.h"
#include "elements.h"
#include "defines.h"
#include "fmo.h"
#include "vlc.h"
#include "mbuffer.h"
#include "header.h"

#include "ctx_tables.h"

#if TRACE
#define SYMTRACESTRING(s) strncpy(sym.tracestring,s,TRACESTRING_SIZE)
#else
#define SYMTRACESTRING(s) // to nothing
#endif

extern int UsedBits;

static void ref_pic_list_reordering();
static void pred_weight_table();


/*!
 ************************************************************************
 * \brief
 *    calculate Ceil(Log2(uiVal))
 ************************************************************************
 */
unsigned CeilLog2( unsigned uiVal)
{
  unsigned uiTmp = uiVal;
  unsigned uiRet = 0;

  while( uiTmp != 0 )
  {
    uiTmp >>= 1;
    uiRet++;
  }
  return uiRet;
}


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
  Bitstream *currStream = partition->bitstream;
  int tmp;

  UsedBits= partition->bitstream->frame_bitoffset; // was hardcoded to 31 for previous start-code. This is better.

  // Get first_mb_in_slice
  currSlice->start_mb_nr = ue_v ("SH: first_mb_in_slice", currStream);

  tmp = ue_v ("SH: slice_type", currStream);
  
  if (tmp>4) tmp -=5;

  img->type = currSlice->picture_type = (SliceType) tmp;

  currSlice->pic_parameter_set_id = ue_v ("SH: pic_parameter_set_id", currStream);
  
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
  Bitstream *currStream = partition->bitstream;

  int val, len;

  img->frame_num = u_v (active_sps->log2_max_frame_num_minus4 + 4, "SH: frame_num", currStream);

  if (active_sps->frame_mbs_only_flag)
  {
    img->structure = FRAME;
    img->field_pic_flag=0;
  }
  else
  {
    // field_pic_flag   u(1)
    img->field_pic_flag = u_1("SH: field_pic_flag", currStream);
    if (img->field_pic_flag)
    {
      // bottom_field_flag  u(1)
      img->bottom_field_flag = u_1("SH: bottom_field_flag", currStream);

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
    img->idr_pic_id = ue_v("SH: idr_pic_id", currStream);
  }

  // POC200301
  if (img->pic_order_cnt_type == 0)
  {
    img->pic_order_cnt_lsb = u_v(active_sps->log2_max_pic_order_cnt_lsb_minus4 + 4, "SH: pic_order_cnt_lsb", currStream);
    if( img->pic_order_present_flag  ==  1 &&  !img->field_pic_flag )
      img->delta_pic_order_cnt_bottom = se_v("SH: delta_pic_order_cnt_bottom", currStream);
  }
  if( img->pic_order_cnt_type == 1 && !img->delta_pic_order_always_zero_flag ) 
  {
    img->delta_pic_order_cnt[ 0 ] = se_v("SH: delta_pic_order_cnt[0]", currStream);
    if( img->pic_order_present_flag  ==  1  &&  !img->field_pic_flag )
      img->delta_pic_order_cnt[ 1 ] = se_v("SH: delta_pic_order_cnt[1]", currStream);
  }
  
  //! redundant_pic_cnt is missing here
  if (active_pps->redundant_pic_cnt_present_flag)
  {
    img->redundant_pic_cnt = u_1 ("SH: redundant_pic_cnt", currStream);
  }

  if(img->type==B_SLICE)
  {
    img->direct_type = u_1 ("SH: direct_spatial_mv_pred_flag", currStream);
  }

  img->num_ref_pic_active_fwd = active_pps->num_ref_idx_l0_active_minus1 + 1;
  img->num_ref_pic_active_bwd = active_pps->num_ref_idx_l1_active_minus1 + 1;

  if(img->type==P_SLICE || img->type == SP_SLICE || img->type==B_SLICE)
  {
    val = u_1 ("SH: num_ref_idx_override_flag", currStream);
    if (val)
    {
      img->num_ref_pic_active_fwd = 1 + ue_v ("SH: num_ref_pic_active_fwd_minus1", currStream);
      
      if(img->type==B_SLICE)
      {
        img->num_ref_pic_active_bwd = 1 + ue_v ("SH: num_ref_pic_active_bwd_minus1", currStream);
      }
    }
  }

  ref_pic_list_reordering();

  img->apply_weights = ((img->weighted_pred_flag && (currSlice->picture_type == P_SLICE || currSlice->picture_type == SP_SLICE) )
          || ((img->weighted_bipred_idc > 0 ) && (currSlice->picture_type == B_SLICE)));

  if ((active_pps->weighted_pred_flag&&(img->type==P_SLICE|| img->type == SP_SLICE))||
      (active_pps->weighted_bipred_idc==1 && (img->type==B_SLICE)))
  {
    pred_weight_table();
  }

  dec_ref_pic_marking(currStream);

  if (active_pps->entropy_coding_mode && img->type!=I_SLICE && img->type!=SI_SLICE)
  {
    img->model_number = ue_v("SH: cabac_init_idc", currStream);
  }
  else 
  {
    img->model_number = 0;
  }

  val = se_v("SH: slice_qp_delta", currStream);
  currSlice->qp = img->qp = 26 + active_pps->pic_init_qp_minus26 + val;

  if(img->type==SP_SLICE || img->type == SI_SLICE) 
  {
    if(img->type==SP_SLICE)
    {
      img->sp_switch = u_1 ("SH: sp_for_switch_flag", currStream);
    }
    val = se_v("SH: slice_qs_delta", currStream);
    img->qpsp = 26 + active_pps->pic_init_qs_minus26 + val;
  }

  if (active_pps->deblocking_filter_parameters_present_flag)
  {
    currSlice->LFDisableIdc = ue_v ("SH: disable_deblocking_filter_idc", currStream);

    if (currSlice->LFDisableIdc!=1)
    {
      currSlice->LFAlphaC0Offset = 2 * se_v("SH: slice_alpha_c0_offset_div2", currStream);

      currSlice->LFBetaOffset = 2 * se_v("SH: slice_beta_offset_div2", currStream);
    }
  }

  if (active_pps->num_slice_groups_minus1>0 && active_pps->slice_group_map_type>=3 &&
      active_pps->slice_group_map_type<=5)
  {
    len = (active_sps->pic_height_in_map_units_minus1+1)*(active_sps->pic_width_in_mbs_minus1+1)/ 
          (active_pps->slice_group_change_rate_minus1+1);
    if (((active_sps->pic_height_in_map_units_minus1+1)*(active_sps->pic_width_in_mbs_minus1+1))% 
          (active_pps->slice_group_change_rate_minus1+1))
          len +=1;
    len = CeilLog2(len);

    img->slice_group_change_cycle = u_v (len, "SH: slice_group_change_cycle", currStream);
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
  img->PicHeightInMbs = img->FrameHeightInMbs / ( 1 + img->field_pic_flag );
  img->PicSizeInMbs   = img->PicWidthInMbs * img->PicHeightInMbs;
  img->FrameSizeInMbs = img->PicWidthInMbs * img->FrameHeightInMbs;

  img->buf_cycle = input->buf_cycle+1;
  img->pn=(((img->structure==BOTTOM_FIELD) ? (img->number/2):img->number)%img->buf_cycle);

  img->max_mb_nr = 2*(img->width * img->height) / (MB_BLOCK_SIZE * MB_BLOCK_SIZE);

  //calculate pocs  POC200301
  decoding_poc(img);
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
  Bitstream *currStream = partition->bitstream;
  int i, val;

  alloc_ref_pic_list_reordering_buffer(currSlice);
  
  if (img->type!=I_SLICE && img->type!=SI_SLICE)
  {
    val = currSlice->ref_pic_list_reordering_flag_l0 = u_1 ("SH: ref_pic_list_reordering_flag_l0", currStream);
    
    if (val)
    {
      i=0;
      do
      {
        val = currSlice->remapping_of_pic_nums_idc_l0[i] = ue_v("SH: remapping_of_pic_nums_idc_l0", currStream);
        if (val==0 || val==1)
        {
          currSlice->abs_diff_pic_num_minus1_l0[i] = ue_v("SH: abs_diff_pic_num_minus1_l0", currStream);
        }
        else
        {
          if (val==2)
          {
            currSlice->long_term_pic_idx_l0[i] = ue_v("SH: long_term_pic_idx_l0", currStream);
          }
        }
        i++;
        // assert (i>img->num_ref_pic_active_fwd);
      } while (val != 3);
    }
  }

  if (img->type==B_SLICE)
  {
    val = currSlice->ref_pic_list_reordering_flag_l1 = u_1 ("SH: ref_pic_list_reordering_flag_l1", currStream);
    
    if (val)
    {
      i=0;
      do
      {
        val = currSlice->remapping_of_pic_nums_idc_l1[i] = ue_v("SH: remapping_of_pic_nums_idc_l1", currStream);
        if (val==0 || val==1)
        {
          currSlice->abs_diff_pic_num_minus1_l1[i] = ue_v("SH: abs_diff_pic_num_minus1_l1", currStream);
        }
        else
        {
          if (val==2)
          {
            currSlice->long_term_pic_idx_l1[i] = ue_v("SH: long_term_pic_idx_l1", currStream);
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
  Bitstream *currStream = partition->bitstream;
  int luma_weight_flag_l0, luma_weight_flag_l1, chroma_weight_flag_l0, chroma_weight_flag_l1;
  int i,j;

  img->luma_log_weight_denom = ue_v ("SH: luma_log_weight_denom", currStream);
  img->wp_round_luma = 1<<(img->luma_log_weight_denom - 1);
  
  img->chroma_log_weight_denom = ue_v ("SH: chroma_log_weight_denom", currStream);
  img->wp_round_chroma = 1<<(img->chroma_log_weight_denom - 1);

  reset_wp_params(img);

  for (i=0; i<img->num_ref_pic_active_fwd; i++)
  {
    luma_weight_flag_l0 = u_1("SH: luma_weight_flag_l0", currStream);
    
    if (luma_weight_flag_l0)
    {
      img->wp_weight[0][i][0] = se_v ("SH: luma_weight_l0", currStream);
      img->wp_offset[0][i][0] = se_v ("SH: luma_offset_l0", currStream);
    }
    else
    {
      img->wp_weight[0][i][0] = 1<<img->luma_log_weight_denom;
      img->wp_offset[0][i][0] = 0;
    }
    
    chroma_weight_flag_l0 = u_1 ("SH: chroma_weight_flag_l0", currStream);
    
    for (j=1; j<3; j++)
    {
      if (chroma_weight_flag_l0)
      {
        img->wp_weight[0][i][j] = se_v("SH: chroma_weight_l0", currStream);
        img->wp_offset[0][i][j] = se_v("SH: chroma_offset_l0", currStream);
      }
      else
      {
        img->wp_weight[0][i][j] = 1<<img->luma_log_weight_denom;
        img->wp_offset[0][i][j] = 0;
      }
    }
  }
  if ((img->type == B_SLICE) && img->weighted_bipred_idc == 1)
  {
    for (i=0; i<img->num_ref_pic_active_bwd; i++)
    {
      luma_weight_flag_l1 = u_1("SH: luma_weight_flag_l1", currStream);
      
      if (luma_weight_flag_l1)
      {
        img->wp_weight[1][i][0] = se_v ("SH: luma_weight_l1", currStream);
        img->wp_offset[1][i][0] = se_v ("SH: luma_offset_l1", currStream);
      }
      else
      {
        img->wp_weight[1][i][0] = 1<<img->luma_log_weight_denom;
        img->wp_offset[1][i][0] = 0;
      }
      
      chroma_weight_flag_l1 = u_1 ("SH: chroma_weight_flag_l1", currStream);
      
      for (j=1; j<3; j++)
      {
        if (chroma_weight_flag_l1)
        {
          img->wp_weight[1][i][j] = se_v("SH: chroma_weight_l1", currStream);
          img->wp_offset[1][i][j] = se_v("SH: chroma_offset_l1", currStream);
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
void dec_ref_pic_marking(Bitstream *currStream)
{
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
    img->no_output_of_prior_pics_flag = u_1("SH: no_output_of_prior_pics_flag", currStream);
    img->long_term_reference_flag = u_1("SH: long_term_reference_flag", currStream);
  }
  else
  {
    img->adaptive_ref_pic_buffering_flag = u_1("SH: adaptive_ref_pic_buffering_flag", currStream);
    if (img->adaptive_ref_pic_buffering_flag)
    {
      // read Memory Management Control Operation 
      do
      {
        tmp_mmco=(MMCObuffer_t*)calloc (1,sizeof (MMCObuffer_t));
        tmp_mmco->Next=NULL;
        
        val = tmp_mmco->MMCO = ue_v("SH: memory_management_control_operation", currStream);

        if ((val==1)||(val==3)) 
        {
          tmp_mmco->DPN = 1 + ue_v("SH: difference_of_pic_nums_minus1", currStream);
        }
        if ((val==2)||(val==3)||(val==6))
        {
          tmp_mmco->LPIN = ue_v("SH: long_term_pic_idx", currStream);
        }
        if (val==4)
        {
          tmp_mmco->MLIP1 = ue_v("SH: max_long_term_pic_idx_plus1", currStream);
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
 *    To calculate the poc values
 *        based upon JVT-F100d2
 *  POC200301: Until Jan 2003, this function will calculate the correct POC
 *    values, but the management of POCs in buffered pictures may need more work.
 * \return
 *    none
 ************************************************************************
 */
void decoding_poc(struct img_par *img)
{
  int i;
  static int flag = 0;
  Slice *currSlice = img->currentSlice;
  // for POC mode 0:
  unsigned int        MaxPicOrderCntLsb = (1<<(active_sps->log2_max_pic_order_cnt_lsb_minus4+4));
  flag ++;
  // for POC mode 1:

  switch ( img->pic_order_cnt_type )
  {
  case 0: // POC MODE 0
    // 1st
    if(img->idr_flag)
      img->PicOrderCntMsb = 0;

    // Calculate the MSBs of current picture
    if( img->pic_order_cnt_lsb  <  img->PrevPicOrderCntLsb  &&  
      ( img->PrevPicOrderCntLsb - img->pic_order_cnt_lsb )  >=  ( MaxPicOrderCntLsb / 2 ) )
	    img->CurrPicOrderCntMsb = img->PicOrderCntMsb + MaxPicOrderCntLsb;
    else if ( img->pic_order_cnt_lsb  >  img->PrevPicOrderCntLsb  &&
      ( img->pic_order_cnt_lsb - img->PrevPicOrderCntLsb )  >  ( MaxPicOrderCntLsb / 2 ) )
	    img->CurrPicOrderCntMsb = img->PicOrderCntMsb - MaxPicOrderCntLsb;
    else
	    img->CurrPicOrderCntMsb = img->PicOrderCntMsb;
    
    // 2nd

    img->toppoc = img->CurrPicOrderCntMsb + img->pic_order_cnt_lsb;
    if( img->bottom_field_flag ) 
	    img->bottompoc = img->CurrPicOrderCntMsb + img->pic_order_cnt_lsb;
    else
	    img->bottompoc = img->toppoc + img->delta_pic_order_cnt_bottom;

    // last: some post-processing. 
    if ( img->toppoc <= img->bottompoc )
      img->ThisPOC = img->framepoc = img->toppoc;
    else
      img->ThisPOC = img->framepoc = img->bottompoc;
    break;


  case 1: // POC MODE 1
    // 1st
    if(img->idr_flag)
    {
      img->FrameNumOffset=0;     //  first pix of IDRGOP, 
      // why the following two lines?????
      img->FirstFieldType = img->bottom_field_flag;              //save type of first field of frame
                                                            //NB may not work with mixed field & frame coding
      img->delta_pic_order_cnt[0]=0;                        //ignore first delta
      if(img->frame_num)  error("frame_num != 0 in idr pix", -1020);
    }
    else if (img->frame_num<img->PreviousFrameNum)
    {             //not first pix of IDRGOP
      img->FrameNumOffset += img->MaxFrameNum;
    }

    // 2nd
    if(img->num_ref_frames_in_pic_order_cnt_cycle) 
      img->AbsFrameNum = img->FrameNumOffset+img->frame_num;
    else 
      img->AbsFrameNum=0;
    if(img->disposable_flag && img->AbsFrameNum>0)
      img->AbsFrameNum--;

    // 3rd
    img->ExpectedDeltaPerPicOrderCntCycle=0;
    if(active_sps->num_ref_frames_in_pic_order_cnt_cycle)
    for(i=0;i<(int) active_sps->num_ref_frames_in_pic_order_cnt_cycle;i++)
      img->ExpectedDeltaPerPicOrderCntCycle += active_sps->offset_for_ref_frame[i];

    if(img->AbsFrameNum)
    {
      img->PicOrderCntCycleCnt = (img->AbsFrameNum-1)/img->num_ref_frames_in_pic_order_cnt_cycle;
      img->FrameNumInPicOrderCntCycle = (img->AbsFrameNum-1)%img->num_ref_frames_in_pic_order_cnt_cycle;
      img->ExpectedPicOrderCnt = img->PicOrderCntCycleCnt*img->ExpectedDeltaPerPicOrderCntCycle;
      for(i=0;i<=(int)img->FrameNumInPicOrderCntCycle;i++)
        img->ExpectedPicOrderCnt += img->offset_for_ref_frame[i];
    }
    else 
      img->ExpectedPicOrderCnt=0;

    if(img->disposable_flag)
      img->ExpectedPicOrderCnt += img->offset_for_non_ref_pic;

    // TIAN DONG: The following processing may need to be updated. POC200301
    // and some codes can move to post_poc().
    // post processing: for management of POC values of pictures in buffer
    if(img->field_pic_flag==0)
    {           //frame pix
      img->toppoc = img->ExpectedPicOrderCnt + img->delta_pic_order_cnt[0];
      img->bottompoc = img->toppoc + img->offset_for_top_to_bottom_field + img->delta_pic_order_cnt[1];
      img->ThisPOC = img->framepoc = (img->toppoc < img->bottompoc)? img->toppoc : img->bottompoc; // POC200301
      if(img->PreviousPOC!=img->ThisPOC)
      {         //new frame detected
        if(img->disposable_flag)
          push_poc(img->toppoc,img->bottompoc,NONREFFRAME);
        else
          push_poc(img->toppoc,img->bottompoc,REFFRAME);
      }
    }
    else if (img->bottom_field_flag==0)
    {  //top field 
      img->ThisPOC = img->toppoc = img->ExpectedPicOrderCnt + img->delta_pic_order_cnt[0];
      img->bottompoc = 0;
      if(img->PreviousPOC!=img->ThisPOC  &&  img->FirstFieldType==img->bottom_field_flag)
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
    {  //bottom field
      img->toppoc = 0;
      img->ThisPOC = img->bottompoc = img->ExpectedPicOrderCnt + img->offset_for_top_to_bottom_field + img->delta_pic_order_cnt[0];
      if(img->PreviousPOC!=img->ThisPOC  &&  img->FirstFieldType==img->bottom_field_flag)
      {           //new frame detected
        if(img->disposable_flag)
          push_poc(0,img->bottompoc,NONREFFRAME);
        else 
          push_poc(0,img->bottompoc,REFFRAME);
      }
      else 
        bottomrefpoc[0] = img->bottompoc;          //2nd field of same frame
    }
  
    // 4th (last) update "Previous" stuff for next slice
//    post_poc( img );
    break;


  case 2: // POC MODE 2
    if(img->idr_flag) // IDR picture
    {
      img->FrameNumOffset=0;     //  first pix of IDRGOP, 
      img->ThisPOC = img->framepoc = img->toppoc = img->bottompoc = 0;
      if(img->frame_num)  error("frame_num != 0 in idr pix", -1020);
    }
    else
    {
      if (img->frame_num<img->PreviousFrameNum)
        img->FrameNumOffset += img->MaxFrameNum;
      img->AbsFrameNum = img->FrameNumOffset+img->frame_num;
      if(img->disposable_flag)
        img->ThisPOC = (2*img->AbsFrameNum - 1);
      else
        img->ThisPOC = (2*img->AbsFrameNum);

      if (img->field_pic_flag==0)
        img->toppoc = img->bottompoc = img->framepoc = img->ThisPOC;
      else if (img->bottom_field_flag==0)
         img->toppoc = img->framepoc = img->ThisPOC;
      else img->bottompoc = img->framepoc = img->ThisPOC;
    }

    break;


  default:
    //error must occurs
    assert( 1==0 );
    break;
  }

  //temp stuff to track tr
  if(img->field_pic_flag)
  {
    img->tr_fld = img->ThisPOC;
    if(img->bottom_field_flag)
    {
      currSlice->picture_id = img->tr = img->bottompoc;
    }
    else
    {   //top field
      currSlice->picture_id = img->tr = img->toppoc;
    }
  }
  else
  {           //frame pix  -  use toppoc/2
    img->tr_frm = img->ThisPOC/2;
    currSlice->picture_id = img->tr = img->toppoc/2;
  }
            //moved from above for stuff that still uses img->tr
            //soon to be obsolete
  if(!img->current_slice_nr)
  { 
    if((img->type != B_SLICE) || !img->disposable_flag) 
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
}

void post_poc(struct img_par *img)
{
  switch ( img->pic_order_cnt_type )
  {
  case 0: // POC MODE 0
    if( !img->disposable_flag )
    {
      img->PrevPicOrderCntLsb = img->pic_order_cnt_lsb;
      img->PicOrderCntMsb = img->CurrPicOrderCntMsb;
    }
    break;


  case 1: // POC MODE 1
    img->PreviousFrameNum=img->frame_num;
    img->Previousfield_pic_flag=img->field_pic_flag;
    img->Previousbottom_field_flag=img->bottom_field_flag;
    img->Previousnal_reference_idc=img->idr_flag;               
    img->Previousdelta_pic_order_cnt[0]=img->delta_pic_order_cnt[0];
    img->Previousdelta_pic_order_cnt[1]=img->delta_pic_order_cnt[1];
    img->PreviousPOC=img->ThisPOC;
    break;


  case 2: // POC MODE 2
    if (!img->disposable_flag)
      img->PreviousFrameNum=img->frame_num;
    break;


  default:
    //error must occurs
    assert( 1==0 );
    break;
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
    printf ("log2_max_frame_num_minus4             %d\n", active_sps->log2_max_frame_num_minus4);         // POC200301
    printf ("log2_max_pic_order_cnt_lsb_minus4     %d\n", active_sps->log2_max_pic_order_cnt_lsb_minus4);
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

// this function is likely to be updated according to how POC values are managed in JM.
int poc_distance( int refa, int refb)
{
  return toprefpoc[refb + 1] - toprefpoc[refa + 1];
}

/*!
 ************************************************************************
 * \brief
 *    return the poc of img as per (8-1) JVT-F100d2
 *  POC200301
 ************************************************************************
 */
int picture_order(struct img_par *img)
{
  if (img->field_pic_flag==0) // is a frame
    return img->framepoc;
  else if (img->bottom_field_flag==0) // top field
    return img->toppoc;
  else // bottom field
    return img->bottompoc;
}
