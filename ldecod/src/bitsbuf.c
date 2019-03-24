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

static int SourceBitBufferPtr;    //!< the position for the next byte-write into the bit buffer
FILE *bits;                //!< the bit stream file


/*!
 ************************************************************************
 * \brief
 *    return the byte position of the bitstream file
 ************************************************************************
 */
int getBitsPos()
{
  return ftell(bits);
}
/*!
 ************************************************************************
 * \brief
 *    set the position of the bitsream file
 ************************************************************************
 */
int setBitsPos(int offset)
{
printf ("in setBitsPos (%d)\n", offset);
  return fseek(bits, offset, SEEK_SET);
}


/*!
 ************************************************************************
 * \brief
 *    Initializes the Source Bit Buffer
 ************************************************************************
 */
void InitializeSourceBitBuffer()
{
  SourceBitBufferPtr = 0;
}


/*!
 ************************************************************************
 * \brief
 *    returns the type of the start code at byte aligned position buf.
 *
 *  \return
 *     the info field of a start code (len == 31) or                      \n
 *     -1, indicating that there is no start code here
 *
 *  This function could be optimized considerably by checking for zero bytes first,
 *  but the current implementation gives us more freedom to define what a
 *  start code actually is
 *  Note that this function could be easily extended to search for non
 *  byte aligned start codes, by simply checking all 8 bit positions
 *  (and not only the zero position
 ************************************************************************
 */
static int TypeOfStartCode (byte *Buf)
{
  int info;

#ifdef _EXP_GOLOMB
  // A bit of optimization first: an EXP Golomb start code starts always with 2 zero bytes and a '1' byte
  if ((Buf[0] != 0) || (Buf[1] != 1) || (Buf[2] != 0))
    return -1;
#else
  // A bit of optimization first: a start code starts always with 3 zero bytes
  if ((Buf[0] != 0) || (Buf[1] != 0) || (Buf[2] != 0))
    return -1;
#endif
  if (31 != GetVLCSymbol (Buf, 0, &info, MAX_CODED_FRAME_SIZE))
  {
    return -1;
  }
  if (info != 0 && info != 1)   // the only two start codes currently allowed
    return -1;
  return info;
}

/*!
 ************************************************************************
 * \brief
 *    Returns the number of bytes in copied into the buffer
 *    The buffer includes the own start code, but excludes the
 *    next slice's start code
 * \return
 *     0 if there is nothing any more to read (EOF)                         \n
 *    -1 in case of any error
 *
 * \note
 *   GetOneSliceIntoSourceBitBuffer() expects start codes as follows:       \n
 *   Slice start code: UVLC, len == 31, info == 1                           \n
 *   Picture Start code: UVLC, len == 31, info == 0                         \n
 * \note
 *   getOneSliceIntoSourceBitBuffer expects Slice and Picture start codes
 *   at byte aligned positions in the file
 *
 ************************************************************************
 */
int GetOneSliceIntoSourceBitBuffer(struct img_par *img, struct inp_par *inp, byte *Buf)
{
  int info, pos;
  int StartCodeFound;
  Slice *currSlice = img->currentSlice;

  InitializeSourceBitBuffer(); // WYK: Useless, can be erased. Why use this?
  // read the first 32 bits (which MUST contain a start code, or eof)
  if (4 != fread (Buf, 1, 4, bits))
  {
    return 0;
  }
  info = TypeOfStartCode (Buf);
  if (info < 0)
  {
    printf ("GetOneSliceIntoSourceBitBuffer: no Start Code at the begin of the slice, return -1\n");
    return -1;
  }
  if (info != 0 && info != 1)
  {
    printf ("GetOneSliceIntoSourceBitBuffer: found start code with invalid info %d, return -1\n", info);
    return -1;
  }
  pos = 4;
  StartCodeFound = 0;

  while (!StartCodeFound)
  {
    if (feof (bits))
    // WYK: Oct. 8, 2001. for detection of frame and slice loss.
    {
      currSlice->next_header = EOS;
      return pos-1; // modified to "pos-1" instead of "pos" 
    }
    Buf[pos++] = fgetc (bits);

    info = TypeOfStartCode(&Buf[pos-4]);
    StartCodeFound = (info == 0 || info == 1);
  }

  // Here, we have found another start code (and read four bytes more than we should
  // have.  Hence, go back in the file

  if (0 != fseek (bits, -4, SEEK_CUR))
  {
    snprintf (errortext, ET_SIZE, "GetOneSliceIntoSourceBitBuffer: Cannot fseek -4 in the bit stream file");
    error(errortext, 600);
  }

  return (pos-4);
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

/*!
 ************************************************************************
 * \brief
 *    read bytes from input
 * \return
 *    Number of bytes read from partition
 ************************************************************************
 */
int GetOnePartitionIntoSourceBitBuffer(int PartitionSize, byte *Buf)
{
  int pos;
  InitializeSourceBitBuffer();
  for (pos=0; pos<PartitionSize; pos++)
  {
    if (feof (bits))
      return pos;
    Buf[pos] = fgetc (bits);
  }
  return pos;
}


