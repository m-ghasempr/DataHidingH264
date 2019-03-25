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
 * \file mb_access.h
 *
 * \brief
 *    Functions for macroblock neighborhoods
 *
 * \author
 *     Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Karsten Sühring          <suehring@hhi.de>
 *************************************************************************************
 */

#ifndef _MB_ACCESS_H_
#define _MB_ACCESS_H_

void CheckAvailabilityOfNeighbors();

void getNeighbour(int curr_mb_nr, int xN, int yN, int luma, PixelPos *pix);
void getLuma4x4Neighbour (int curr_mb_nr, int block_x, int block_y, int rel_x, int rel_y, PixelPos *pix);
void getChroma4x4Neighbour (int curr_mb_nr, int block_x, int block_y, int rel_x, int rel_y, PixelPos *pix);

int  mb_is_available(int mbAddr, int currMbAddr);
void get_mb_pos (int mb_addr, int *x, int*y);
void get_mb_block_pos (int mb_addr, int *x, int*y);



#endif
