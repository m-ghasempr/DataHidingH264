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
 ************************************************************************
 * \file image.h
 *
 * \author
 *  Inge Lille-Langøy               <inge.lille-langoy@telenor.com>
 *  Copyright (C) 1999  Telenor Satellite Services, Norway
 ************************************************************************
 */
#ifndef _IMAGE_H_
#define _IMAGE_H_

#include "mbuffer.h"

extern StorablePicture *enc_picture;
extern StorablePicture *enc_frame_picture;
extern StorablePicture *enc_top_picture;
extern StorablePicture *enc_bottom_picture;

int encode_one_frame ();
Boolean dummy_slice_too_big(int bits_slice);
void copy_rdopt_data (int field_type);    //!< For MB level field/frame coding tools

#endif

