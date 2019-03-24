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
  SEI_ZERO,        //!< 0 is undefined, useless
  SEI_TEMPORAL_REF,
  SEI_CLOCK_TIMESTAMP,
  SEI_PANSCAN_RECT,
  SEI_BUFFERING_PERIOD,
  SEI_HRD_PICTURE,
  SEI_FILLER_PAYLOAD,
  SEI_USER_DATA_REGISTERED_ITU_T_T35,
  SEI_USER_DATA_UNREGISTERED,
  SEI_RANDOM_ACCESS_POINT,
  SEI_REF_PIC_BUFFER_MANAGEMENT_REPETITION,
  SEI_SPARE_PICTURE,
  SEI_SCENE_INFORMATION,
  SEI_SUBSEQ_INFORMATION,
  SEI_SUBSEQ_LAYER_CHARACTERISTICS,
  SEI_SUBSEQ_CHARACTERISTICS,
  SEI_MAX_ELEMENTS  //!< number of maximum syntax elements
} SEI_type;

#define MAX_FN 256

#define AGGREGATION_PACKET_TYPE 4
#define SEI_PACKET_TYPE 5  // Tian Dong: See VCEG-N72, it need updates

#define NORMAL_SEI 0
#define AGGREGATION_SEI 1

void InterpretSEIMessage(byte* msg, int size, struct img_par *img);
void interpret_spare_picture( byte* payload, int size, struct img_par *img );
void interpret_subsequence_info( byte* payload, int size, struct img_par *img );
void interpret_subsequence_layer_info( byte* payload, int size, struct img_par *img );
void interpret_subsequence_characteristics_info( byte* payload, int size, struct img_par *img );

#endif