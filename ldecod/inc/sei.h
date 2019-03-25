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
 * \file sei.h
 *
 * \brief
 *    Prototypes for sei.h
 *************************************************************************************
 */

#ifndef SEI_H
#define SEI_H

//! definition of SEI payload type
typedef enum {
  SEI_BUFFERING_PERIOD = 0,
  SEI_PIC_TIMING,
  SEI_PAN_SCAN_RECT,
  SEI_FILLER_PAYLOAD,
  SEI_USER_DATA_REGISTERED_ITU_T_T35,
  SEI_USER_DATA_UNREGISTERED,
  SEI_RANDOM_ACCESS_POINT,
  SEI_DEC_REF_PIC_MARKING_REPETITION,
  SEI_SPARE_PIC,
  SEI_SCENE_INFO,
  SEI_SUB_SEQ_INFO,
  SEI_SUB_SEQ_LAYER_CHARACTERISTICS,
  SEI_SUB_SEQ_CHARACTERISTICS,
  SEI_FULL_FRAME_FREEZE,
  SEI_FULL_FRAME_FREEZE_RELEASE,
  SEI_FULL_FRAME_SNAPSHOT,
  SEI_PROGRESSIVE_REFINEMENT_SEGMENT_START,
  SEI_PROGRESSIVE_REFINEMENT_SEGMENT_END,
  SEI_MOTION_CONSTRAINED_SLICE_GROUP_SET,

  SEI_MAX_ELEMENTS  //!< number of maximum syntax elements
} SEI_type;

#define MAX_FN 256

void InterpretSEIMessage(byte* msg, int size, ImageParameters *img);
void interpret_spare_pic( byte* payload, int size, ImageParameters *img );
void interpret_subsequence_info( byte* payload, int size, ImageParameters *img );
void interpret_subsequence_layer_characteristics_info( byte* payload, int size, ImageParameters *img );
void interpret_subsequence_characteristics_info( byte* payload, int size, ImageParameters *img );
void interpret_scene_information( byte* payload, int size, ImageParameters *img ); // JVT-D099
void interpret_user_data_registered_itu_t_t35_info( byte* payload, int size, ImageParameters *img );
void interpret_user_data_unregistered_info( byte* payload, int size, ImageParameters *img );
void interpret_pan_scan_rect_info( byte* payload, int size, ImageParameters *img );
void interpret_random_access_info( byte* payload, int size, ImageParameters *img );
void interpret_filler_payload_info( byte* payload, int size, ImageParameters *img );
void interpret_dec_ref_pic_marking_repetition_info( byte* payload, int size, ImageParameters *img );
void interpret_full_frame_freeze_info( byte* payload, int size, ImageParameters *img );
void interpret_full_frame_freeze_release_info( byte* payload, int size, ImageParameters *img );
void interpret_full_frame_snapshot_info( byte* payload, int size, ImageParameters *img );
void interpret_progressive_refinement_start_info( byte* payload, int size, ImageParameters *img );
void interpret_progressive_refinement_end_info( byte* payload, int size, ImageParameters *img );
void interpret_motion_constrained_slice_group_set_info( byte* payload, int size, ImageParameters *img );
void interpret_reserved_info( byte* payload, int size, ImageParameters *img );
void interpret_buffering_period_info( byte* payload, int size, ImageParameters *img );
void interpret_picture_timing_info( byte* payload, int size, ImageParameters *img );

#endif
