
/*!
 ************************************************************************
 * \file  nalu.c
 *
 * \brief
 *    Decoder NALU support functions
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Stephan Wenger   <stewe@cs.tu-berlin.de>
 ************************************************************************
 */

#include "global.h"
#include "nalu.h"
#include "annexb.h"
#include "rtp.h"

BitsFile bitsfile;

static int LastAccessUnitExists  = 0;
static int NALUCount = 0;


/*!
*************************************************************************************
* \brief
*    Initialize bitstream reading structure
*
* \param
*    filemode: 
*
*************************************************************************************
*/

void initBitsFile (int filemode)
{
  switch (filemode)
  {
  case PAR_OF_ANNEXB:
    bitsfile.OpenBitsFile     = OpenAnnexBFile;
    bitsfile.CloseBitsFile    = CloseAnnexBFile;
    bitsfile.GetNALU          = GetAnnexbNALU;
    break;
  case PAR_OF_RTP:
    bitsfile.OpenBitsFile     = OpenRTPFile;
    bitsfile.CloseBitsFile    = CloseRTPFile;
    bitsfile.GetNALU          = GetRTPNALU;
    break;
  default:
    error ("initBitsFile: Unknown bitstream file mode", 255);
    break;
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Converts a NALU to an RBSP
 *
 * \param
 *    nalu: nalu structure to be filled
 *
 * \return
 *    length of the RBSP in bytes
 *************************************************************************************
 */

static int NALUtoRBSP (NALU_t *nalu)
{
  assert (nalu != NULL);

  nalu->len = EBSPtoRBSP (nalu->buf, nalu->len, 1) ;

  return nalu->len ;
}

/*!
************************************************************************
* \brief
*    Read the next NAL unit (with error handling)
************************************************************************
*/
int read_next_nalu(NALU_t *nalu)
{
  int ret;

  ret = bitsfile.GetNALU(nalu);

  if (ret < 0)
  {
    snprintf (errortext, ET_SIZE, "Error while getting the NALU in file format %s, exit\n", params->FileFormat==PAR_OF_ANNEXB?"Annex B":"RTP");
    error (errortext, 601);
  }
  if (ret == 0)
  {
    FreeNALU(nalu);
    return 0;
  }

  //In some cases, zero_byte shall be present. If current NALU is a VCL NALU, we can't tell
  //whether it is the first VCL NALU at this point, so only non-VCL NAL unit is checked here.
  CheckZeroByteNonVCL(nalu);

  ret = NALUtoRBSP(nalu);

  if (ret < 0)
    error ("Invalid startcode emulation prevention found.", 602);


  // Got a NALU
  if (nalu->forbidden_bit)
  {
    error ("Found NALU with forbidden_bit set, bit error?", 603);
  }

  return nalu->len;
}

void CheckZeroByteNonVCL(NALU_t *nalu)
{
  int CheckZeroByte=0;

  //This function deals only with non-VCL NAL units
  if(nalu->nal_unit_type>=1&&nalu->nal_unit_type<=5)
    return;

  //for SPS and PPS, zero_byte shall exist
  if(nalu->nal_unit_type==NALU_TYPE_SPS || nalu->nal_unit_type==NALU_TYPE_PPS)
    CheckZeroByte=1;
  //check the possibility of the current NALU to be the start of a new access unit, according to 7.4.1.2.3
  if(nalu->nal_unit_type==NALU_TYPE_AUD  || nalu->nal_unit_type==NALU_TYPE_SPS ||
    nalu->nal_unit_type==NALU_TYPE_PPS || nalu->nal_unit_type==NALU_TYPE_SEI ||
    (nalu->nal_unit_type>=13 && nalu->nal_unit_type<=18))
  {
    if(LastAccessUnitExists)
    {
      LastAccessUnitExists=0;    //deliver the last access unit to decoder
      NALUCount=0;
    }
  }
  NALUCount++;
  //for the first NAL unit in an access unit, zero_byte shall exists
  if(NALUCount==1)
    CheckZeroByte=1;
  if(CheckZeroByte && nalu->startcodeprefix_len==3)
  {
    printf("Warning: zero_byte shall exist\n");
    //because it is not a very serious problem, we do not exit here
  }
}

void CheckZeroByteVCL(NALU_t *nalu)
{
  int CheckZeroByte=0;

  //This function deals only with VCL NAL units
  if(!(nalu->nal_unit_type>=1&&nalu->nal_unit_type<=5))
    return;

  if(LastAccessUnitExists)
  {
    NALUCount=0;
  }
  NALUCount++;
  //the first VCL NAL unit that is the first NAL unit after last VCL NAL unit indicates
  //the start of a new access unit and hence the first NAL unit of the new access unit.           (sounds like a tongue twister :-)
  if(NALUCount == 1)
    CheckZeroByte = 1;
  LastAccessUnitExists = 1;
  if(CheckZeroByte && nalu->startcodeprefix_len==3)
  {
    printf("warning: zero_byte shall exist\n");
    //because it is not a very serious problem, we do not exit here
  }
}
