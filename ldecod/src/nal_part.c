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
 * \file  nal_part.c
 *
 * \brief
 *    Network Adaptation layer for partition file
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Tobias Oelbaum <oelbaum@hhi.de, oelbaum@drehvial.de>
 ************************************************************************
 */

#include <string.h>

#include "contributors.h"
#include "global.h"
#include "elements.h"

int assignSE2partition[][SE_MAX_ELEMENTS] =
{
  // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19  // elementnumber (no not uncomment)
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },   //!< all elements in one partition no data partitioning
  {  0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 0, 0 }    //!< three partitions per slice
};

int PartitionMode;

/*!
 ************************************************************************
 * \brief
 *    Resets the entries in the bitstream struct
 ************************************************************************
 */
void free_Partition(Bitstream *currStream)
{
  byte *buf = currStream->streamBuffer;

  currStream->bitstream_length = 0;
  currStream->frame_bitoffset = 0;
  currStream->ei_flag =0;
  memset (buf, 0x00, MAX_CODED_FRAME_SIZE);
}
