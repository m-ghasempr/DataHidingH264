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

#include "contributors.h"

#include <assert.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <string.h>

#include "global.h"
#include "elements.h"

extern void tracebits(const char *trace_str,  int len,  int info,int value1, int value2) ;


/*!
 ************************************************************************
 * \brief
 *    Input file containing Partition structures
 * \return
 *    TRUE if a startcode follows (meaning that all symbols are used).  \n
 *    FALSE otherwise.
 * \para
 *     Looks like this is dead code.  StW, 7.7.02
 ************************************************************************/
/*
int slice_startcode_follows(struct img_par *img, struct inp_par *inp)
{
  Slice *currSlice = img->currentSlice;
  int dp_Nr = assignSE2partition[currSlice->dp_mode][SE_MBTYPE];
  DataPartition *dP = &(currSlice->partArr[dp_Nr]);
  Bitstream   *currStream = dP->bitstream;
  byte *buf = currStream->streamBuffer;
  int frame_bitoffset = currStream->frame_bitoffset;
  int info;
assert (!currStream->ei_flag);
printf ("slice_startcode_follows returns %d\n", (-1 == GetVLCSymbol (buf, frame_bitoffset, &info, currStream->bitstream_length)));

  if (currStream->ei_flag)
    return (img->current_mb_nr == currSlice->last_mb_nr);
  else
  {
    if (-1 == GetVLCSymbol (buf, frame_bitoffset, &info, currStream->bitstream_length))
      return TRUE;
    else
      return FALSE;
  }
}
*/
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
