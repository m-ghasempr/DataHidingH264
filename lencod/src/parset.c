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
 **************************************************************************************
 * \file
 *    parset.c
 * \brief
 *    Picture and Sequence Parameter set generation and handling
 *  \date 25 November 2002
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details) 
 *      - Stephan Wenger        <stewe@cs.tu-berlin.de>
 *
 **************************************************************************************
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <malloc.h>
#include <string.h>
 
#include "contributors.h"
#include "parsetcommon.h"
#include "nalu.h"
#include "parset.h"
#include "fmo.h"
#include "global.h"
#include "vlc.h"

// Local helpers
static int IdentifyProfile();
static int IdentifyLevel();
static int IdentifyNumRefFrames();
static int GenerateVUISequenceParameters();


/*! 
 *************************************************************************************
 * \brief
 *    int GenerateSeq_parameter_set_NALU ();
 *
 * \param 
 *    None.  Uses the global variables through FillParameterSetStructures()
 *
 * \return
 *    A NALU containing the Sequence ParameterSet
 *
 *************************************************************************************
 */
 
NALU_t *GenerateSeq_parameter_set_NALU ()
{
  seq_parameter_set_rbsp_t *sps = NULL; 
  pic_parameter_set_rbsp_t *pps = NULL;
  NALU_t *n = AllocNALU(64000);
  int RBSPlen = 0;
  int NALUlen;
  byte rbsp[MAXRBSPSIZE];

  sps = AllocSPS();
  pps = AllocPPS();

  FillParameterSetStructures (sps, pps);
  RBSPlen = GenerateSeq_parameter_set_rbsp (sps, rbsp);
  NALUlen = RBSPtoNALU (rbsp, n, RBSPlen, NALU_TYPE_SPS, NALU_PRIORITY_HIGHEST, 0, 1);
  n->startcodeprefix_len = 4;
  
  FreeSPS (sps);
  FreePPS (pps);

  return n;
}
  

/*! 
 *************************************************************************************
 * \brief
 *    NALU_t *GeneratePic_parameter_set_NALU ();
 *
 * \param 
 *    None.  Uses the global variables through FillParameterSetStructures()
 *
 * \return
 *    A NALU containing the Picture Parameter Set
 *
 *************************************************************************************
 */
 
NALU_t *GeneratePic_parameter_set_NALU()
{
  seq_parameter_set_rbsp_t *sps = NULL; 
  pic_parameter_set_rbsp_t *pps = NULL;
  NALU_t *n = AllocNALU(64000);
  int RBSPlen = 0;
  int NALUlen;
  byte rbsp[MAXRBSPSIZE];

  sps = AllocSPS();
  pps = AllocPPS();

  FillParameterSetStructures (sps, pps);
  RBSPlen = GeneratePic_parameter_set_rbsp (pps, rbsp);
  NALUlen = RBSPtoNALU (rbsp, n, RBSPlen, NALU_TYPE_PPS, NALU_PRIORITY_HIGHEST, 0, 1);
  n->startcodeprefix_len = 4;
  FreeSPS (sps);
  FreePPS (pps);

  return n;
}


/*!
 ************************************************************************
 * \brief
 *    FillParameterSetStructures: extracts info from global variables and
 *    generates a picture and sequence parameter set structure
 *
 * \par Input:
 *    seq_parameter_set_rbsp_t *sps,  Sequence parameter set to be filled
 *    pic_parameter_set_rbsp_t *pps   Picture parameter set to be filled
 * \par
 *    Function reads all kinds of values from several global variables,
 *    including input-> and image-> and fills in the sps and pps.  Many
 *    values are current hard-coded to defaults, especially most of the
 *    VUI stuff.  Currently, the sps and pps structures are fixed length
 *    with the exception of the fully flexible FMO map (mode 6).  This
 *    mode is not supported.  Hence, the function does not need to
 *    allocate memory for the FMOmap, the pointer slice_group_id is
 *    always NULL.  If one wants to implement FMO mode 6, one would need
 *    to malloc the memory for the map here, and the caller would need
 *    to free it after use.
 *
 * \par 
 *    Limitations
 *    Currently, the encoder does not support multiple parameter sets,
 *    primarily because the config file does not support it.  Hence the
 *    pic_parameter_set_id and the seq_parameter_set_id are always zero.
 *    If one day multiple parameter sets are implemented, it would
 *    make sense to break this function into two, one for the picture and
 *    one for the sequence.
 *    Currently, FMO is not supported
 *    The following pps and sps elements seem not to be used in the encoder 
 *    or decoder and, hence, a guessed default value is conveyed:
 *
 *    pps->num_ref_idx_l0_active_minus1 = img->num_ref_pic_active_fwd_minus1;
 *    pps->num_ref_idx_l1_active_minus1 = img->num_ref_pic_active_bwd_minus1;
 *    pps->chroma_qp_index_offset = 0;
 *    sps->required_frame_num_update_behaviour_flag = FALSE;
 *    sps->direct_temporal_constrained_flag = FALSE;
 *
 * \par
 *    Regarding the QP
 *    The previous software versions coded the absolute QP only in the 
 *    slice header.  This is kept, and the offset in the PPS is coded 
 *    even if we could save bits by intelligently using this field.
 *
 ************************************************************************
 */

void FillParameterSetStructures (seq_parameter_set_rbsp_t *sps, 
                                 pic_parameter_set_rbsp_t *pps)
{
  unsigned i;
  // *************************************************************************
  // Sequence Parameter Set
  // *************************************************************************
  assert (sps != NULL);
  assert (pps != NULL);
  // Profile and Level should be calculated using the info from the config
  // file.  Calculation is hidden in IndetifyProfile() and IdentifyLevel()
  sps->profile_idc = IdentifyProfile();
  sps->level_idc = IdentifyLevel();

	sps->more_than_one_slice_group_allowed_flag = (input->FmoType!=0);
  // no out of order slice support
	sps->arbitrary_slice_order_allowed_flag = 0;
	sps->redundant_slices_allowed_flag = (input->redundant_slice_flag);

  // Parameter Set ID hardcoded to zero
  sps->seq_parameter_set_id = 0;

  //! POC stuff:
  //! The following values are hard-coded in init_poc().  Apparently,
  //! the poc implementation covers only a subset of the poc functionality.
  //! Here, the same subset is implemented.  Changes in the POC stuff have
  //! also to be reflected here
  sps->log2_max_frame_num_minus4 = LOG2_MAX_FRAME_NUM_MINUS4;
  sps->log2_max_pic_order_cnt_lsb_minus4 = LOG2_MAX_PIC_ORDER_CNT_LSB_MINUS4;   // POC200301
  sps->pic_order_cnt_type = img->pic_order_cnt_type;
  sps->num_ref_frames_in_pic_order_cnt_cycle = img->num_ref_frames_in_pic_order_cnt_cycle;
  sps->delta_pic_order_always_zero_flag = img->delta_pic_order_always_zero_flag;
  sps->offset_for_non_ref_pic = img->offset_for_non_ref_pic;
  sps->offset_for_top_to_bottom_field = img->offset_for_top_to_bottom_field;
  // This is the only one used, because num_ref_frames_in_pic_order_cnt_cycle 
  // is hard coded to 1.  
  sps->offset_for_ref_frame[0] = img->offset_for_ref_frame[0];
  // End of POC stuff

  // Number of Reference Frames
  sps->num_ref_frames = IdentifyNumRefFrames();

  //required_frame_num_update_behaviour_flag hardcoded to zero
  sps->required_frame_num_update_behaviour_flag = FALSE;    // double check

  // Picture size, finally a simple one :-)
  sps->frame_width_in_mbs_minus1 = (input->img_width/16) -1;
  sps->frame_height_in_mbs_minus1 = (input->img_height/16) -1;

  // a couple of flags, simple
  sps->frame_mbs_only_flag = (input->InterlaceCodingOption == FRAME_CODING);
  sps->mb_adaptive_frame_field_flag = (input->InterlaceCodingOption == MB_CODING);
  sps->direct_8x8_inference_flag = FALSE;

  // Sequence VUI not implemented, signalled as not present
  sps->vui_parameters_present_flag = FALSE;

  // *************************************************************************
  // Picture Parameter Set 
  // *************************************************************************

  pps->seq_parameter_set_id = 0;
  pps->pic_parameter_set_id = 0;
  pps->entropy_coding_mode = (input->symbol_mode==UVLC?0:1);

  // JVT-Fxxx (by Stephan Wenger, make this flag unconditional
  pps->pic_order_present_flag = img->pic_order_present_flag;


  // Begin FMO stuff
  pps->num_slice_groups_minus1 = input->num_slice_groups_minus1;

  if (pps->num_slice_groups_minus1 > 0)
    switch (input->mb_allocation_map_type)
    {
    case 0:
      pps->mb_slice_group_map_type = 6;       // This implementation always uses the fully flexible map
      printf ("Param.c: FMO type 0 not yet signalled, using 6 instead\n");
      goto FMOTYPE6;
      // This code should work as soon as the re-generation of the MBAmap in
      // the decoder does its job
      // pps->mb_slice_group_map_type = 0;    // according to the draft
      // for (i=0; i< pps->num_slice_groups_minus1+1; i++)
      //   pps->run_length[i] = img->width/16;    // Line interleaving
      break;
    case 1:
      pps->mb_slice_group_map_type = 6;       // This implementation always uses the fully flexible map
      printf ("Param.c: FMO type 1 not yet signalled, using 6 instead\n");
      goto FMOTYPE6;
      // This code should work as soon as the re-generation of the MBAmap in
      // the decoder does its job
      // pps->mb_slice_group_map_type = 1;    // according to the draft
      break;
    case 2:
      assert(pps->num_slice_groups_minus1 == 1);      // restriction of the config file
      for(i=0; i<pps->num_slice_groups_minus1; i++)
      {
        pps->top_left_mb[i] = input->top_left_mb;
        pps->bottom_right_mb[i] = input->bottom_right_mb;
      }
     break;
    case 3:
    case 4:
    case 5:
      pps->slice_group_change_direction_flag = input->slice_group_change_direction;
      pps->slice_group_change_rate_minus1 = input->slice_group_change_rate_minus1;
      break;
    case 6:
FMOTYPE6:
      pps->mb_slice_group_map_type = 6;       // This implementation always uses the fully flexible map
      pps->slice_group_id_cnt_minus1 = (img->height/16) * (img->width/16) - 1;
      for (i=0;i<=pps->slice_group_id_cnt_minus1; i++)
        pps->slice_group_id[i++] = FmoMB2SliceGroup (i);
      break;
    default:
      printf ("Param.c: FMO type invalid, default\n");
      assert (0==1);
    }
// End FMO stuff

  pps->num_ref_idx_l0_active_minus1 = img->num_ref_pic_active_fwd_minus1;   // check this
  pps->num_ref_idx_l1_active_minus1 = img->num_ref_pic_active_bwd_minus1;   // check this
  
  pps->weighted_pred_flag = input->WeightedPrediction;
  pps->weighted_bipred_idc = input->WeightedBiprediction;

  pps->pic_init_qp_minus26 = 0;         // hard coded to zero, QP lives in the slice header
  pps->pic_init_qs_minus26 = 0;

  pps->chroma_qp_index_offset = 0;      // double check: is this chroma fidelity thing already implemented???

  pps->deblocking_filter_parameters_present_flag = input->LFSendParameters;
  pps->constrained_intra_pred_flag = input->UseConstrainedIntraPred;
  
  pps->redundant_pic_cnt_present_flag = 0;

  // the picture vui consists currently of the cropping rectangle, which cannot
  // used by the current decoder and hence is never sent.
  pps->frame_cropping_flag = FALSE;
};



/*! 
 *************************************************************************************
 * \brief
 *    int GenerateSeq_parameter_set_rbsp (seq_parameter_set_rbsp_t *sps, char *rbsp);
 *
 * \param 
 *    sps:  sequence parameter structure
 *    rbsp:  buffer to be filled with the rbsp, size should be at least MAXIMUMPARSETRBSPSIZE
 *
 * \return
 *    size of the RBSP in bytes
 *
 * \note
 *    Sequence Parameter VUI function is called, but the function implements
 *    an exit (-1)
 *************************************************************************************
 */
 
int GenerateSeq_parameter_set_rbsp (seq_parameter_set_rbsp_t *sps, char *rbsp)
{
  DataPartition *partition;
  int len = 0, LenInBytes;
  unsigned i;

  assert (rbsp != NULL);
  // In order to use the entropy coding functions from golomb.c we need 
  // to allocate a partition structure.  It will be freed later in this
  // function
  if ((partition=calloc(1,sizeof(DataPartition)))==NULL) no_mem_exit("SeqParameterSet:partition");
  if ((partition->bitstream=calloc(1, sizeof(Bitstream)))==NULL) no_mem_exit("SeqParameterSet:bitstream");
  // .. and use the rbsp provided (or allocated above) for the data
  partition->bitstream->streamBuffer = rbsp;
  partition->bitstream->bits_to_go = 8;

  len+=u_v  (8, "SPS: profile_idc",                             sps->profile_idc,                               partition);
  len+=u_v  (8, "SPS: level_idc",                               sps->level_idc,                                 partition);
  len+=u_1  ("SPS: more_than_one_slice_group_allowed_flag",     sps->more_than_one_slice_group_allowed_flag,    partition);
  len+=u_1  ("SPS: arbitrary_slice_order_allowed_flag",         sps->arbitrary_slice_order_allowed_flag,        partition);
  len+=u_1  ("SPS: redundant_slices_allowed_flag",              sps->redundant_slices_allowed_flag,             partition);
  len+=ue_v ("SPS: seq_parameter_set_id",                    sps->seq_parameter_set_id,                      partition);
  len+=ue_v ("SPS: log2_max_frame_num_minus4",               sps->log2_max_frame_num_minus4,                 partition);
  len+=ue_v ("SPS: pic_order_cnt_type",                      sps->pic_order_cnt_type,                        partition);
  // POC200301
  if (sps->pic_order_cnt_type == 0)
    len+=ue_v ("SPS: log2_max_pic_order_cnt_lsb_minus4",     sps->log2_max_pic_order_cnt_lsb_minus4,         partition);
  else if (sps->pic_order_cnt_type == 1)
  {
    len+=u_1  ("SPS: delta_pic_order_always_zero_flag",        sps->delta_pic_order_always_zero_flag,          partition);
    len+=se_v ("SPS: offset_for_non_ref_pic",                  sps->offset_for_non_ref_pic,                    partition);
    len+=se_v ("SPS: offset_for_top_to_bottom_field",          sps->offset_for_top_to_bottom_field,            partition);
    len+=ue_v ("SPS: num_ref_frames_in_pic_order_cnt_cycle",   sps->num_ref_frames_in_pic_order_cnt_cycle,     partition);
    for (i=0; i<sps->num_ref_frames_in_pic_order_cnt_cycle; i++)
      len+=se_v ("SPS: offset_for_ref_frame",                  sps->offset_for_ref_frame[i],                      partition);
  }
  len+=ue_v ("SPS: num_ref_frames",                          sps->num_ref_frames,                            partition);
  len+=u_1  ("SPS: required_frame_num_update_behaviour_flag",sps->required_frame_num_update_behaviour_flag,  partition);
  len+=ue_v ("SPS: frame_width_in_mbs_minus1",               sps->frame_width_in_mbs_minus1,                 partition);
  len+=ue_v ("SPS: frame_height_in_mbs_minus1",              sps->frame_height_in_mbs_minus1,                partition);
  len+=u_1  ("SPS: frame_mbs_only_flag",                     sps->frame_mbs_only_flag,                       partition);
  if (!sps->frame_mbs_only_flag)
  {
    len+=u_1  ("SPS: mb_adaptive_frame_field_flag",            sps->mb_adaptive_frame_field_flag,              partition);
  }
  len+=u_1  ("SPS: direct_8x8_inference_flag",               sps->direct_8x8_inference_flag,                 partition);
  len+=u_1  ("SPS: vui_parameters_present_flag",             sps->vui_parameters_present_flag,               partition);
  if (sps->vui_parameters_present_flag)
    len+=GenerateVUISequenceParameters();    // currently a dummy, asserting

  SODBtoRBSP(partition->bitstream);     // copies the last couple of bits into the byte buffer
  
  LenInBytes=partition->bitstream->byte_pos;

  free (partition->bitstream);
  free (partition);
  
  return LenInBytes;
}


/*! 
 *************************************************************************************
 * \brief
 *    int GeneratePic_parameter_set_rbsp (pic_parameter_set_rbsp_t *sps, char *rbsp);
 *
 * \param 
 *    pps: sequence parameter structure
 *    rbsp:  buffer to be filled with the rbsp, size should be at least MAXIMUMPARSETRBSPSIZE
 *
 * \return
 *    size of the RBSP in bytes, negative in case of an error
 *
 * \note
 *    Picture Parameter VUI function is called, but the function implements
 *    an exit (-1)
 *************************************************************************************
 */
 
int GeneratePic_parameter_set_rbsp (pic_parameter_set_rbsp_t *pps, char *rbsp)
{
  DataPartition *partition;
  int len = 0, LenInBytes;
  unsigned i;
  unsigned NumberBitsPerSliceGroupId;

  assert (rbsp != NULL);

  // In order to use the entropy coding functions from golomb.c we need 
  // to allocate a partition structure.  It will be freed later in this
  // function
  if ((partition=calloc(1,sizeof(DataPartition)))==NULL) no_mem_exit("PicParameterSet:partition");
  if ((partition->bitstream=calloc(1, sizeof(Bitstream)))==NULL) no_mem_exit("PicParameterSet:bitstream");
  // .. and use the rbsp provided (or allocated above) for the data
  partition->bitstream->streamBuffer = rbsp;
  partition->bitstream->bits_to_go = 8;

  len+=ue_v ("PPS: pic_parameter_set_id",                    pps->pic_parameter_set_id,                      partition);
  len+=ue_v ("PPS: seq_parameter_set_id",                    pps->seq_parameter_set_id,                      partition);
  len+=u_1  ("PPS: entropy_coding_mode",                     pps->entropy_coding_mode,                       partition);
  len+=u_1  ("PPS: pic_order_present_flag",                  pps->pic_order_present_flag,                    partition);
  len+=ue_v ("PPS: num_slice_groups_minus1",                 pps->num_slice_groups_minus1,                   partition);

  // FMO stuff
  if(pps->num_slice_groups_minus1 > 0 )
  {
    len+=ue_v ("PPS: mb_slice_group_map_type",                 pps->mb_slice_group_map_type,                   partition);
    if (pps->mb_slice_group_map_type == 0)
      for (i=0; i<pps->num_slice_groups_minus1; i++)
        len+=ue_v ("PPS: run_length[i]",                           pps->run_length[i],                             partition);
    else if (pps->mb_slice_group_map_type==2)
      for (i=0; i<pps->num_slice_groups_minus1; i++)
      {
        len+=ue_v ("PPS: top_left_mb[i]",                          pps->top_left_mb[i],                           partition);
        len+=ue_v ("PPS: bottom_right_mb[i]",                      pps->bottom_right_mb[i],                       partition);
      }
    else if (pps->mb_slice_group_map_type == 3 ||
             pps->mb_slice_group_map_type == 4 ||
             pps->mb_slice_group_map_type == 5) 
    {
      len+=u_1  ("PPS: slice_group_change_direction_flag",         pps->slice_group_change_direction_flag,         partition);
      len+=ue_v ("PPS: slice_group_change_rate_minus1",            pps->slice_group_change_rate_minus1,            partition);
    } 
    else if (pps->mb_slice_group_map_type == 6)
    {
      if (pps->num_slice_groups_minus1>=4)
        NumberBitsPerSliceGroupId=3;
      else if (pps->num_slice_groups_minus1>=2)
        NumberBitsPerSliceGroupId=2;
      else if (pps->num_slice_groups_minus1>=1)
        NumberBitsPerSliceGroupId=1;
      else
        NumberBitsPerSliceGroupId=0;
        
      len+=ue_v ("PPS: slice_group_id_cnt_minus1",                  pps->slice_group_id_cnt_minus1,                 partition);
      for(i=0; i<=pps->slice_group_id_cnt_minus1; i++)
        len+= u_v  (NumberBitsPerSliceGroupId, "PPS: >slice_group_id[i]",                            pps->slice_group_id[i],                         partition);
    }
  }
  // End of FMO stuff

  len+=ue_v ("PPS: num_ref_idx_l0_active_minus1",             pps->num_ref_idx_l0_active_minus1,              partition);
  len+=ue_v ("PPS: num_ref_idx_l1_active_minus1",             pps->num_ref_idx_l1_active_minus1,              partition);
  len+=u_1  ("PPS: weighted_pred_flag",                       pps->weighted_pred_flag,                        partition);
  len+=u_v  (2, "PPS: weighted_bipred_idc",                   pps->weighted_bipred_idc,                       partition);
  len+=se_v ("PPS: pic_init_qp_minus26",                      pps->pic_init_qp_minus26,                       partition);
  len+=se_v ("PPS: pic_init_qs_minus26",                      pps->pic_init_qs_minus26,                       partition);
  len+=se_v ("PPS: chroma_qp_index_offset",                   pps->chroma_qp_index_offset,                    partition);
  len+=u_1  ("PPS: deblocking_filter_parameters_present_flag",pps->deblocking_filter_parameters_present_flag, partition);
  len+=u_1  ("PPS: constrained_intra_pred_flag",              pps->constrained_intra_pred_flag,               partition);
  len+=u_1  ("PPS: redundant_pic_cnt_present_flag",           pps->redundant_pic_cnt_present_flag,            partition);
  len+=u_1  ("PPS: frame_cropping_flag",                      pps->frame_cropping_flag,                       partition);
  if (pps->frame_cropping_flag)
  {
    len+=ue_v ("PPS: frame_cropping_rect_left_offset",          pps->frame_cropping_rect_left_offset,           partition);
    len+=ue_v ("PPS: frame_cropping_rect_right_offset",         pps->frame_cropping_rect_right_offset,          partition);
    len+=ue_v ("PPS: frame_cropping_rect_top_offset",           pps->frame_cropping_rect_top_offset,            partition);
    len+=ue_v ("PPS: frame_cropping_rect_bottom_offset",        pps->frame_cropping_rect_bottom_offset,         partition);
  }

  SODBtoRBSP(partition->bitstream);     // copies the last couple of bits into the byte buffer
  
  LenInBytes=partition->bitstream->byte_pos;

  // Get rid of the helper structures
  free (partition->bitstream);
  free (partition);

  return LenInBytes;
}



/*! 
 *************************************************************************************
 * \brief
 *    Returns the Profile
 *
 * \param none
 *
 * \return
 *    Profile according to Annex A
 *
 * \note
 *    Function is currently a dummy.  Should "calculate" the profile from those
 *    config file parameters.  E.g.
 *
 *    Profile = Baseline;
 *    if (CABAC Used || Interlace used) Profile=Main;
 *    if (!Cabac Used) && (Bframes | SPframes) Profile = Streaming;
 *
 *************************************************************************************
 */
int IdentifyProfile()
{
  return PROFILE_IDC;       // Baseline
};

/*! 
 *************************************************************************************
 * \brief
 *    Returns the Level
 *
 * \param none
 *
 * \return
 *    Level according to Annex A
 *
 * \note
 *    This function is currently a dummy, but should calculate the level out of 
 *    the config file parameters (primarily the picture size)
 *************************************************************************************
 */
int IdentifyLevel()
{
  return LEVEL_IDC;
};


/*! 
 *************************************************************************************
 * \brief
 *    Returns the number of reference frame buffers
 *
 * \param none
 *
 * \return
 *    Number of reference frame buffers used
 *
 * \note
 *    This function currently maps to input->no_multpred.  With all this interlace
 *    stuff this may or may not be correct.  If you determine a problem with the
 *    memory management for Interlace, then this could be one possible problem.
 *    However, so far no problem have been determined by my limited testing of
 *    a stupid 1950's technology :-)  StW, 11/27/02
 *************************************************************************************
 */

int IdentifyNumRefFrames()
{
  return input->no_multpred;
}


static int GenerateVUISequenceParameters()
{
  printf ("Sequence Parameter VUI not yet implemented, this should never happen, exit\n");
  exit (-1);
}

