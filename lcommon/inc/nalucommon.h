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
 **************************************************************************************
 * \file
 *    nalucommon.h.h
 * \brief
 *    NALU handling common to encoder and decoder
 *  \date 25 November 2002
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details) 
 *      - Stephan Wenger        <stewe@cs.tu-berlin.de>
 ***************************************************************************************
 */


#ifndef _NALUCOMMON_H_
#define _NALUCOMMON_H_

#define MAXNALUSIZE 64000     /*! just some number */
#define MAXRBSPSIZE 64000

#define NALU_TYPE_SLICE   1
#define NALU_TYPE_DPA     2
#define NALU_TYPE_DPB     3
#define NALU_TYPE_DPC     4
#define NALU_TYPE_IDR     5
#define NALU_TYPE_SEI     6
#define NALU_TYPE_SPS     7
#define NALU_TYPE_PPS     8
#define NALU_TYPE_PD      9
#define NALU_TYPE_FILL    10

#define NALU_PRIORITY_HIGHEST     3
#define NALU_PRIORITY_HIGH        2
#define NALU_PRIRITY_LOW          1
#define NALU_PRIORITY_DISPOSABLE  0


typedef struct 
{
  int startcodeprefix_len;      //! 4 for parameter sets and first slice in picture, 3 for everything else (suggested)
  unsigned len;                 //! Length of the NAL unit (Excluding the start code, which does not belong to the NALU)
  int nal_unit_type;            //! NALU_TYPE_xxxx
  int nal_reference_idc;        //! NALU_PRIORITY_xxxx
  int forbidden_bit;            //! should be always FALSE
  char buf[MAXNALUSIZE];        //! conjtains the first byte followed by the EBSP
} NALU_t;


NALU_t *AllocNALU();
void FreeNALU(NALU_t *n);

#endif