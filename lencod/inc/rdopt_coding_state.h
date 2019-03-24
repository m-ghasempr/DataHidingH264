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
 *    rdopt_coding_state.h
 *
 * \author
 *    Heiko Schwarz
 *
 * \date
 *    17. April 2001
 *
 * \brief
 *    Headerfile for storing/restoring coding state
 *    (for rd-optimized mode decision)
 **************************************************************************
 */

#ifndef _RD_OPT_CS_H_
#define _RD_OPT_CS_H_

#include "global.h"



typedef struct {

  // important variables of data partition array
  int                   no_part;
  EncodingEnvironment  *encenv;
  Bitstream            *bitstream;

  // contexts for binary arithmetic coding
  int                   symbol_mode;
  MotionInfoContexts   *mot_ctx;
  TextureInfoContexts  *tex_ctx;

  // syntax element number and bitcounters
  int                   currSEnr;
  int                   bitcounter[MAX_BITCOUNTER_MB];

  // elements of current macroblock
  int                   mvd[2][BLOCK_MULTIPLE][BLOCK_MULTIPLE][2];
  int                   coeffs_count[BLOCK_MULTIPLE][BLOCK_MULTIPLE];


} CSobj;
typedef CSobj* CSptr;


void  delete_coding_state  (CSptr);  //!< delete structure
CSptr create_coding_state  ();       //!< create structure

void  store_coding_state   (CSptr);  //!< store parameters
void  reset_coding_state   (CSptr);  //!< restore parameters


//a few vars exported for routines in abt.c
extern int mpr8x8[16][16];
extern int ****cofAC,****cofAC8x8;
extern int ***cofDC;
extern int **cofAC4x4,****cofAC4x4intern;
extern int cbp, cbp8x8, cnt_nonz_8x8;
extern CSptr cs_mb,cs_b8,cs_cm,cs_imb,cs_ib8,cs_ib4,cs_pc;


#endif

