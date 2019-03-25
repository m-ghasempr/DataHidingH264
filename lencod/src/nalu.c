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

#include <stdio.h>
#include <assert.h>
#include <malloc.h>
#include "memory.h"

#include "global.h"
#include "nalu.h"

/*! 
 *************************************************************************************
 * \brief
 *    Converts an RBSP to a NALU
 *
 * \param 
 *    rbsp: byte buffer with the rbsp
 *    nalu: nalu structure to be filled
 *    rbsp_size: size of the rbsp in bytes
 *    nal_unit_type, nal_reference_idc: as in JVT doc
 *    min_num_bytes: some incomprehensible CABAC stuff
 *    UseAnnexbLongStartcode: when 1 and when using AnnexB bytestreams, then use a long startcode prefix
 *
 * \return
 *    length of the NALU in bytes
 *************************************************************************************
 */

int RBSPtoNALU (char *rbsp, NALU_t *nalu, int rbsp_size, int nal_unit_type, int nal_reference_idc, 
                int min_num_bytes, int UseAnnexbLongStartcode)
{
  int len;

  assert (nalu != NULL);
  assert (nal_reference_idc <=3 && nal_reference_idc >=0);
  assert (nal_unit_type > 0 && nal_unit_type <= 10);
  assert (rbsp_size < MAXRBSPSIZE);

  nalu->forbidden_bit = 0;
  nalu->nal_reference_idc = nal_reference_idc;
  nalu->nal_unit_type = nal_unit_type;
  nalu->startcodeprefix_len = UseAnnexbLongStartcode?4:3;
  nalu->buf[0] =
    nalu->forbidden_bit << 7      |
    nalu->nal_reference_idc << 5  |
    nalu->nal_unit_type;

  memcpy (&nalu->buf[1], rbsp, rbsp_size);
// printf ("First Byte %x\n", nalu->buf[0]);
// printf ("RBSPtoNALU: Before: NALU len %d\t RBSP %x %x %x %x\n", rbsp_size, (unsigned) nalu->buf[1], (unsigned) nalu->buf[2], (unsigned) nalu->buf[3], (unsigned) nalu->buf[4]);

  len = 1 + RBSPtoEBSP (&nalu->buf[1], 0, rbsp_size, min_num_bytes);

// printf ("RBSPtoNALU: After : NALU len %d\t EBSP %x %x %x %x\n", rbsp_size, (unsigned) nalu->buf[1], (unsigned) nalu->buf[2], (unsigned) nalu->buf[3], (unsigned) nalu->buf[4]);
// printf ("len %d\n\n", len);
  nalu->len = len;

  return len;
}


