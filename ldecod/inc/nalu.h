
/*!
 **************************************************************************************
 * \file
 *    parset.h
 * \brief
 *    Picture and Sequence Parameter Sets, encoder operations
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
  void (*OpenBitsFile)    (char *filename);
  void (*CloseBitsFile)   (void);
  int  (*GetNALU)         (NALU_t *nalu);
} BitsFile;

extern BitsFile bitsfile;

extern void initBitsFile (int filemode);

extern void CheckZeroByteNonVCL(NALU_t *nalu);
extern void CheckZeroByteVCL(NALU_t *nalu);

extern int read_next_nalu(NALU_t *nalu);

#endif
