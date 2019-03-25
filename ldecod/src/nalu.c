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
 *    Converts a NALU to an RBSP
 *
 * \param 
 *    rbsp: byte buffer with the rbsp
 *    nalu: nalu structure to be filled
 *
 * \return
 *    length of the RBSP in bytes
 *************************************************************************************
 */

int NALUtoRBSP (NALU_t *nalu, char *rbsp)

{
  int len;

  assert (nalu != NULL);
  assert (rbsp != NULL);

  nalu->forbidden_bit = (nalu->buf[0]<<7) & 1;
  nalu->nal_reference_idc = (nalu->buf[0]<<5) & 3;
  nalu->nal_unit_type = (nalu->buf[0]) & 0x1f;


  memcpy (rbsp, &nalu->buf[1], nalu->len-1);

// printf ("First Byte %x\n", nalu->buf[0]);
// printf ("RBSPtoNALU: Before: NALU len %d\t RBSP %x %x %x %x\n", nalu->len, (unsigned) nalu->buf[1], (unsigned) nalu->buf[2], (unsigned) nalu->buf[3], (unsigned) nalu->buf[4]);

  len = EBSPtoRBSP (rbsp, nalu->len-1, 0);

// printf ("RBSPtoNALU: After : NALU len %d\t EBSP %x %x %x %x\n", nalu->len, (unsigned) nalu->buf[1], (unsigned) nalu->buf[2], (unsigned) nalu->buf[3], (unsigned) nalu->buf[4]);
// printf ("len %d\n\n", len);

  return len;
}

