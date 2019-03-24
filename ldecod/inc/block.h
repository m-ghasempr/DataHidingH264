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
 * \file block.h
 *
 * \author
 *  Inge Lille-Langøy               <inge.lille-langoy@telenor.com>    \n
 *  Telenor Satellite Services                                         \n
 *  P.O.Box 6914 St.Olavs plass                                        \n
 *  N-0130 Oslo, Norway
 *
 ************************************************************************
 */

#ifndef _BLOCK_H_
#define _BLOCK_H_

#include "global.h"

#define DQ_BITS         6
#define DQ_ROUND        (1<<(DQ_BITS-1))

const int JQQ = 1048576;
const int JQQ2 = 524288;
const int JQQ3= 349525;
const int JQQ4= 174762;

int const MAP[4][4]=
{
  {0,2,4,5},
  {1,3,5,5},
  {2,4,5,5},
  {3,5,5,5},
};

extern const int JQ1[];
extern const int JQ[32];
extern const byte FILTER_STR[32][4];//!< defined in image.h
extern const byte QP_SCALE_CR[40] ;
extern const int dequant_coef[6][4][4];

#endif

