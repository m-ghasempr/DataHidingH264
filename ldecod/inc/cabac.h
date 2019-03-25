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
 ***************************************************************************
 * \file
 *    cabac.h
 *
 * \brief
 *    Headerfile for entropy coding routines
 *
 * \author
 *    Detlev Marpe                                                         \n
 *    Copyright (C) 2000 HEINRICH HERTZ INSTITUTE All Rights Reserved.
 *
 * \date
 *    21. Oct 2000 (Changes by Tobias Oelbaum 28.08.2001)
 ***************************************************************************
 */

#ifndef _CABAC_H_
#define _CABAC_H_

#include "global.h"

MotionInfoContexts* create_contexts_MotionInfo(void);
TextureInfoContexts* create_contexts_TextureInfo(void);
void init_contexts_MotionInfo(struct img_par *img, MotionInfoContexts *enco_ctx);
void init_contexts_TextureInfo(struct img_par *img, TextureInfoContexts *enco_ctx);
void delete_contexts_MotionInfo(MotionInfoContexts *enco_ctx);
void delete_contexts_TextureInfo(TextureInfoContexts *enco_ctx);


void readMB_typeInfo_CABAC(SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);
void readB8_typeInfo_CABAC(SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);
void readIntraPredMode_CABAC(SyntaxElement *se, struct inp_par *inp,struct img_par *img, DecodingEnvironmentPtr dep_dp);
void readRefFrame_CABAC(SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);
void readMVD_CABAC(SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);
void readCBP_CABAC(SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);
void readRunLevel_CABAC(SyntaxElement *se, struct inp_par *inp, struct img_par *img,  DecodingEnvironmentPtr dep_dp);
void readDquant_CABAC(SyntaxElement *se,struct inp_par *inp,struct img_par *img,DecodingEnvironmentPtr dep_dp);
void readCIPredMode_CABAC(SyntaxElement *se,struct inp_par *inp,struct img_par *img,DecodingEnvironmentPtr dep_dp);
void readMB_skip_flagInfo_CABAC( SyntaxElement *se, struct inp_par *inp, struct img_par *img, DecodingEnvironmentPtr dep_dp);
void readFieldModeInfo_CABAC(SyntaxElement *se,struct inp_par *inp,struct img_par *img,DecodingEnvironmentPtr dep_dp); 

int  readSyntaxElement_CABAC(SyntaxElement *se, struct img_par *img, struct inp_par *inp, DataPartition *this_dataPart);

int  check_next_mb_and_get_field_mode_CABAC(SyntaxElement *se,struct img_par *img,struct inp_par *inp,DataPartition  *act_dp);
void CheckAvailabilityOfNeighborsCABAC();


#endif  // _CABAC_H_

