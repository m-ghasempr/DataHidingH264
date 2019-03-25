
/*!
 **************************************************************************************
 * \file
 *    nalu.h
 * \brief
 *    Common NALU support functions
 *
 * \date 25 November 2002
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Stephan Wenger        <stewe@cs.tu-berlin.de>
 ***************************************************************************************
 */


#ifndef _NALU_H_
#define _NALU_H_

#include "nalucommon.h"

typedef struct sBitsFile
{
  void (*OpenBitsFile)    (ImageParameters *p_Img, char *filename);
  void (*CloseBitsFile)   (ImageParameters *p_Img);
  int  (*GetNALU)         (ImageParameters *p_Img, NALU_t *nalu);
} BitsFile;

extern void initBitsFile (ImageParameters *p_Img, int filemode);
extern void CheckZeroByteNonVCL(ImageParameters *p_Img, NALU_t *nalu);
extern void CheckZeroByteVCL   (ImageParameters *p_Img, NALU_t *nalu);

extern int read_next_nalu(ImageParameters *p_Img, NALU_t *nalu);

#endif
