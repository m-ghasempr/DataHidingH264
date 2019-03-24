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
 * \file erc_globals.h
 *
 * \brief
 *      global header file for error concealment module
 *
 * \author
 *      - Viktor Varsa                     <viktor.varsa@nokia.com>
 *      - Ye-Kui Wang                   <wangy@cs.tut.fi>
 ************************************************************************
 */

#ifndef _ERC_GLOBALS_H_
#define _ERC_GLOBALS_H_


#include <memory.h>

/* "block" means an 8x8 pixel area */

/* Region modes */
#define REGMODE_INTER_COPY       0  /* Copy region */
#define REGMODE_INTER_PRED       1  /* Inter region with motion vectors */
#define REGMODE_INTRA            2  /* Intra region */
#define REGMODE_SPLITTED         3  /* Any region mode higher than this indicates that the region 
                                       is splitted which means 8x8 block */
#define REGMODE_INTER_COPY_8x8   4
#define REGMODE_INTER_PRED_8x8   5
#define REGMODE_INTRA_8x8        6

/* YUV pixel domain image arrays for a video frame */
typedef struct
{
  byte *yptr;
  byte *uptr;
  byte *vptr;
} frame;

/* region structure stores information about a region that is needed for concealment */
typedef struct 
{
  byte regionMode;  /* region mode as above */
  int xMin;         /* X coordinate of the pixel position of the top-left corner of the region */
  int yMin;         /* Y coordinate of the pixel position of the top-left corner of the region */
  int32 mv[3];      /* motion vectors in 1/4 or 1/8 pixel units: mvx = mv[0], mvy = mv[1], 
                              and ref_frame = mv[2] */
} objectBuffer_t;

#endif

