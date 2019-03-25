
/*!
 ***************************************************************************
 * \file
 *    q_offsets.h
 *
 * \brief
 *    Headerfile for q_offsets array
 *
 * \date
 *    18. Nov 2004
 ***************************************************************************
 */

#ifndef _Q_OFFSETS_H_
#define _Q_OFFSETS_H_

static const int OffsetBits = 11;

extern void Init_QOffsetMatrix      (ImageParameters *p_Img, InputParameters *p_Inp);
extern void CalculateOffset4x4Param (ImageParameters *p_Img, InputParameters *p_Inp);
extern void CalculateOffset8x8Param (ImageParameters *p_Img, InputParameters *p_Inp);
extern void free_QOffsets           (QuantParameters *p_Quant, InputParameters *p_Inp);

#endif
