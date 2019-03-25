
/*!
 *************************************************************************************
 * \file annexb.c
 *
 * \brief
 *    Annex B Byte Stream format
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Stephan Wenger                  <stewe@cs.tu-berlin.de>
 *************************************************************************************
 */

#include "global.h"
#include "annexb.h"
#include "memalloc.h"

static int  BitStreamFile = -1;                //!< the bit stream file
static const unsigned int IOBUFFERSIZE = 65536;
static byte *iobuffer;
static byte *iobufferread;
static unsigned int bytesinbuffer;
static int is_eof;

static int IsFirstByteStreamNALU = 1;


/*!
************************************************************************
* \brief
*    fill IO buffer
************************************************************************
*/
static inline int getChunk()
{
  unsigned int readbytes = read (BitStreamFile, iobuffer, IOBUFFERSIZE);
  if (0==readbytes)
  {
    is_eof = TRUE;
    return 0;
  }

  bytesinbuffer = readbytes;
  iobufferread  = iobuffer;
  return readbytes;
}


/*!
************************************************************************
* \brief
*    returns a byte from IO buffer
************************************************************************
*/
static inline byte getfbyte()
{
  if (0==bytesinbuffer)
  {
    if (0 == getChunk())
      return 0;
  }
  bytesinbuffer--;
  return (*iobufferread++);
}


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
 *     pointer to byte-stream
 *  \param zeros_in_startcode
 *     indicates number of 0x00 bytes in start-code.
 ************************************************************************
 */
static inline int FindStartCode (unsigned char *Buf, int zeros_in_startcode)
{
  int i;

  for (i = 0; i < zeros_in_startcode; i++)
  {
    if(*(Buf++) != 0)
    {
      return 0;
    }
  }

  if(*Buf != 1)
    return 0;

  return 1;
}


/*!
 ************************************************************************
 * \brief
 *    Returns the size of the NALU (bits between start codes in case of
 *    Annex B.  nalu->buf and nalu->len are filled.  Other field in
 *    nalu-> remain uninitialized (will be taken care of by NALUtoRBSP.
 *
 * \return
 *     0 if there is nothing any more to read (EOF)
 *    -1 in case of any error
 *
 *  \note Side-effect: Returns length of start-code in bytes.
 *
 * \note
 *   GetAnnexbNALU expects start codes at byte aligned positions in the file
 *
 ************************************************************************
 */

static int nextstartcodebytes = 0;

int GetAnnexbNALU (NALU_t *nalu)
{
  int i;
  int info2 = 0, info3 = 0, pos = 0;
  int StartCodeFound = 0;
  static byte *Buf, *pBuf; 
  int LeadingZero8BitsCount = 0;

  if ((Buf = (byte*) calloc (nalu->max_size , sizeof(char))) == NULL) 
    no_mem_exit("GetAnnexbNALU: Buf");

  pBuf = Buf;

  if (nextstartcodebytes != 0)
  {
    for (i=0; i<nextstartcodebytes-1; i++)
    {
      (*pBuf++) = 0;
      pos++;
    }
    (*pBuf++) = 1;
    pos++;
  }
  else
  {
    while(!is_eof)
    {
      pos++;
      if ((*(pBuf++)= getfbyte())!= 0)
        break;
    }
  }
  if(is_eof == TRUE)
  {
    free (Buf);
    if(pos==0)
    {
      return 0;
    }
    else
    {
      printf( "GetAnnexbNALU can't read start code\n");
      return -1;
    }
  }  

  if(*(pBuf - 1) != 1 || pos < 3)
  {
    printf ("GetAnnexbNALU: no Start Code at the beginning of the NALU, return -1\n");
    free (Buf);
    return -1;
  }

  if (pos == 3)
  {
    nalu->startcodeprefix_len = 3;
  }
  else
  {
    LeadingZero8BitsCount = pos - 4;
    nalu->startcodeprefix_len = 4;
  }

  //the 1st byte stream NAL unit can has leading_zero_8bits, but subsequent ones are not
  //allowed to contain it since these zeros(if any) are considered trailing_zero_8bits
  //of the previous byte stream NAL unit.
  if(!IsFirstByteStreamNALU && LeadingZero8BitsCount > 0)
  {
    printf ("GetAnnexbNALU: The leading_zero_8bits syntax can only be present in the first byte stream NAL unit, return -1\n");
    free (Buf);
    return -1;
  }

  LeadingZero8BitsCount = pos;
  IsFirstByteStreamNALU = 0;

  while (!StartCodeFound)
  {
    if (is_eof == TRUE)
    {
      pBuf -= 2;
      while(*(pBuf--)==0)
        pos--;

      nalu->len = (pos - 1) - LeadingZero8BitsCount;
      memcpy (nalu->buf, Buf + LeadingZero8BitsCount, nalu->len);
      nalu->forbidden_bit     = (*(nalu->buf) >> 7) & 1;
      nalu->nal_reference_idc = (*(nalu->buf) >> 5) & 3;
      nalu->nal_unit_type     = (*(nalu->buf)) & 0x1f;
      nextstartcodebytes = 0;

      // printf ("GetAnnexbNALU, eof case: pos %d nalu->len %d, nalu->reference_idc %d, nal_unit_type %d \n", pos, nalu->len, nalu->nal_reference_idc, nalu->nal_unit_type);

#if TRACE
      fprintf (p_trace, "\n\nLast NALU in File\n\n");
      fprintf (p_trace, "Annex B NALU w/ %s startcode, len %d, forbidden_bit %d, nal_reference_idc %d, nal_unit_type %d\n\n",
        nalu->startcodeprefix_len == 4?"long":"short", nalu->len, nalu->forbidden_bit, nalu->nal_reference_idc, nalu->nal_unit_type);
      fflush (p_trace);
#endif

      free (Buf);
      return (pos - 1);
    }

    pos++;
    *(pBuf ++)  = getfbyte();    
    info3 = FindStartCode(pBuf - 4, 3);
    if(info3 != 1)
    {
      info2 = FindStartCode(pBuf - 3, 2);
      StartCodeFound = (info2 == 1);
    }
    else
      StartCodeFound = 1;
  }

  // Here, we have found another start code (and read length of startcode bytes more than we should
  // have.  Hence, go back in the file
  if(info3 == 1)  //if the detected start code is 00 00 01, trailing_zero_8bits is sure not to be present
  {
    pBuf -= 5;
    while(*(pBuf--) == 0)
      pos--;
    nextstartcodebytes = 4;
  }
  else if (info2 == 1)
    nextstartcodebytes = 3;
  else
  {
    printf(" Panic: Error in next start code search \n");
    free (Buf);
    return -1;
  }

  pos -= nextstartcodebytes;

  // Here the leading zeros(if any), Start code, the complete NALU, trailing zeros(if any)
  // and the next start code is in the Buf.
  // The size of Buf is pos - rewind, pos are the number of bytes excluding the next
  // start code, and (pos) - LeadingZero8BitsCount
  // is the size of the NALU.

  nalu->len = pos - LeadingZero8BitsCount;
  memcpy (nalu->buf, Buf + LeadingZero8BitsCount, nalu->len);
  nalu->forbidden_bit     = (*(nalu->buf) >> 7) & 1;
  nalu->nal_reference_idc = (*(nalu->buf) >> 5) & 3;
  nalu->nal_unit_type     = (*(nalu->buf)) & 0x1f;
  nalu->lost_packets = 0;

  //printf ("GetAnnexbNALU, regular case: pos %d nalu->len %d, nalu->reference_idc %d, nal_unit_type %d \n", pos, nalu->len, nalu->nal_reference_idc, nalu->nal_unit_type);
#if TRACE
  fprintf (p_trace, "\n\nAnnex B NALU w/ %s startcode, len %d, forbidden_bit %d, nal_reference_idc %d, nal_unit_type %d\n\n",
    nalu->startcodeprefix_len == 4?"long":"short", nalu->len, nalu->forbidden_bit, nalu->nal_reference_idc, nalu->nal_unit_type);
  fflush (p_trace);
#endif

  free (Buf);
  return (pos);
}


/*!
 ************************************************************************
 * \brief
 *    Opens the bit stream file named fn
 * \return
 *    none
 ************************************************************************
 */
void OpenAnnexBFile (char *fn)
{
  if (NULL != iobuffer)
  {
    error ("OpenAnnexBFile: tried to open Annex B file twice",500);
  }
  if ((BitStreamFile = open(fn, OPENFLAGS_READ)) == -1)
  {
    snprintf (errortext, ET_SIZE, "Cannot open Annex B ByteStream file '%s'", fn);
    error(errortext,500);
  }
  iobuffer = malloc (IOBUFFERSIZE * sizeof (byte));
  if (NULL == iobuffer)
  {
    error ("OpenAnnexBFile: cannot allocate IO buffer",500);
  }
  is_eof = FALSE;
  getChunk();
}


/*!
 ************************************************************************
 * \brief
 *    Closes the bit stream file
 ************************************************************************
 */
void CloseAnnexBFile()
{
  if (BitStreamFile != -1)
  {
    close(BitStreamFile);
    BitStreamFile = - 1;
  }
  free (iobuffer);
  iobuffer = NULL;
}

