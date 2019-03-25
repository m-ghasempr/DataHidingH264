
/*!
 ************************************************************************
 * \file  nalu.c
 *
 * \brief
 *    Common NALU support functions
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Stephan Wenger   <stewe@cs.tu-berlin.de>
 ************************************************************************
 */

#include "global.h"
#include "nalu.h"
#include "nal.h"

/*!
 *************************************************************************************
 * \brief
 *    Converts an RBSP to a NALU
 *
 * \param rbsp
 *    byte buffer with the rbsp
 * \param nalu
 *    nalu structure to be filled
 * \param rbsp_size
 *    size of the rbsp in bytes
 * \param nal_unit_type
 *    as in JVT doc
 * \param nal_reference_idc
 *    as in JVT doc
 * \param min_num_bytes
 *    some incomprehensible CABAC stuff
 * \param UseAnnexbLongStartcode
 *    when 1 and when using AnnexB bytestreams, then use a long startcode prefix
 *
 * \return
 *    length of the NALU in bytes
 *************************************************************************************
 */

int RBSPtoNALU (unsigned char *rbsp, NALU_t *nalu, int rbsp_size, int nal_unit_type, int nal_reference_idc,
                int min_num_bytes, int UseAnnexbLongStartcode)
{
  int len;

  assert (nalu != NULL);
  assert (nal_reference_idc <=3 && nal_reference_idc >=0);
  assert (nal_unit_type > 0 && nal_unit_type <= 10);
  assert (rbsp_size < MAXRBSPSIZE);
  
  nalu->startcodeprefix_len = UseAnnexbLongStartcode ? 4 : 3;
  nalu->forbidden_bit       = 0;  
  nalu->nal_reference_idc   = (NalRefIdc) nal_reference_idc;
  nalu->nal_unit_type       = (NaluType) nal_unit_type;    

  len = RBSPtoEBSP (nalu->buf, rbsp, rbsp_size);
  nalu->len = len;

  return len;
}

/*!
 ************************************************************************
 *  \brief
 *     write AUD NALU
 ************************************************************************
 */

int Write_AUD_NALU( void )
{  
  int     RBSPlen = 0;
  int     NALUlen, len;
  byte    rbsp[MAXRBSPSIZE];
  NALU_t *nalu = AllocNALU( 64000 );

  switch( img->type )
  {
  case I_SLICE:
    img->primary_pic_type = 0;
    break;
  case P_SLICE:
    img->primary_pic_type = 1;
    break;
  case B_SLICE:
    img->primary_pic_type = 2;
    break;
  }
  RBSPlen = 1;
  rbsp[0] = img->primary_pic_type << 5;
  rbsp[0] |= (1 << 4);

  // write RBSP into NALU
  NALUlen = RBSPtoNALU( rbsp, nalu, RBSPlen, NALU_TYPE_AUD, NALU_PRIORITY_HIGHEST, 0, 1 );
  // write NALU into bitstream
  len     = WriteNALU( nalu );

  FreeNALU( nalu );

  return len;
}
