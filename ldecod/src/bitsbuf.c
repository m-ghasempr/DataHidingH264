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
 * \file bitsbuf.c
 *
 * \brief
 *    Bit Stream Buffer Management
 *
 *  This module is introduced to facilitate the handling of the bit stream
 *  reading.  It assumes a start code taht is a UVLC code word with len 31
 *  and info 0 for Picture Start code, and info = 1 for Slice start code.
 *  The module should be used for all entropy coding schemes, including UVLC
 *  and CABAC.
 *
 * \author
 *    Stephan Wenger   stewe@cs,tu-berlin.de
 *************************************************************************************
 */

#include <stdlib.h>
#include <assert.h>

#include "global.h"
#include "bitsbuf.h"

FILE *bits;                //!< the bit stream file


/*!
 ************************************************************************
 * \brief
 *    returns if new start code is found at byte aligned position buf.
 *    new-startcode is of form N 0x00 bytes, followed by a 0x01 byte.
 *
 *  \return
 *     1 if start-code is found or                      \n
 *     0, indicating that there is no start code
 *
 *  \param Buf
 *         pointer to byte-stream
 *  \param zeros_in_startcode
 *         indicates number of 0x00 bytes in start-code.
 ************************************************************************
 */
static int FindStartCode (byte *Buf, int zeros_in_startcode)
{
  int info;
  int i;

  info = 1;
  for (i = 0; i < zeros_in_startcode; i++)
    if(Buf[i] != 0)
      info = 0;

  if(Buf[i] != 1)
    info = 0;
  return info;
}

/*!
 ************************************************************************
 * \brief
 *    Returns the number of bytes in copied into the buffer
 *    The buffer includes the own start code, but excludes the
 *    next slice's start code.
 *    Important: Works with the new start-code. See below.
 * \return
 *     0 if there is nothing any more to read (EOF)                         \n
 *    -1 in case of any error
 *
 *
 * \note
 *   GetOneSliceIntoSourceBitBuffer () expects start codes as follows:       \n
 *   Start code: either 2 or 3 0x00 bytes followed by one 0x01 byte.                           \n
 *  \note Side-effect: Returns length of start-code in bytes. \n
 *
 * \note
 *   getOneSliceIntoSourceBitBuffer expects start codes
 *   at byte aligned positions in the file
 *
 ************************************************************************
 */
int GetOneSliceIntoSourceBitBuffer(struct img_par *img, struct inp_par *inp, byte *Buf, int *startcodeprefix_len)
{
  int info2, info3, pos = 0;
  int StartCodeFound, rewind;
  Slice *currSlice = img->currentSlice;

  *startcodeprefix_len=3;

  // read the first 32 bits (which MUST contain a start code, or eof)
  info2 = 0;
  info3 = 0;
  
  if (3 != fread (Buf, 1, 3, bits))
  {
    return 0;
  }

  info2 = FindStartCode (Buf, 2);
  if(info2 != 1) {
    if(1 != fread(Buf+3, 1, 1, bits))
      return 0;
    info3 = FindStartCode (Buf, 3);
  }

  if (info2 != 1 && info3 != 1)
  {
    printf ("GetOneSliceIntoSourceBitBuffer: no Start Code at the begin of the slice, return -1\n");
    return -1;
  }

  if( info2 == 1) {
    *startcodeprefix_len = 3;
    pos = 3;
  }
  else if(info3 ==1 ) {
    pos = 4;
    *startcodeprefix_len = 4;
  }
  else
    printf( " Panic: Error \n");

  StartCodeFound = 0;
  info2 = 0;
  info3 = 0;

  while (!StartCodeFound)
  {
    if (feof (bits))
    // WYK: Oct. 8, 2001. for detection of frame and slice loss.
    {
      currSlice->next_header = EOS;
      return pos-1; // modified to "pos-1" instead of "pos" 
    }
    Buf[pos++] = fgetc (bits);
    info3 = FindStartCode(&Buf[pos-4], 3);
    if(info3 != 1)
      info2 = FindStartCode(&Buf[pos-3], 2);
    StartCodeFound = (info2 == 1 || info3 == 1);
  }

 
  // Here, we have found another start code (and read length of startcode bytes more than we should
  // have.  Hence, go back in the file
  rewind = 0;
  if(info3 == 1)
    rewind = -4;
  else if (info2 == 1)
    rewind = -3;
  else
    printf(" Panic: Error in next start code search \n");

  if (0 != fseek (bits, rewind, SEEK_CUR))
  {
    snprintf (errortext, ET_SIZE, "GetOneSliceIntoSourceBitBuffer: Cannot fseek %d in the bit stream file", rewind);
    error(errortext, 600);
  }

  return (pos+rewind);
}




/*!
 ************************************************************************
 * \brief
 *    Attempts to open the bit stream file named fn
 * \return
 *    0 on success,                                                    \n
 *    -1 in case of error
 ************************************************************************
 */
int OpenBitstreamFile (char *fn)
{
  if (NULL == (bits=fopen(fn, "rb")))
    return -1;
  return 0;
}


/*!
 ************************************************************************
 * \brief
 *    Closes the bit stream file
 ************************************************************************
 */
void CloseBitstreamFile()
{
  fclose (bits);
}

