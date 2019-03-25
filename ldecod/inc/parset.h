
/*!
 **************************************************************************************
 * \file
 *    parset.h
 * \brief
 *    Picture and Sequence Parameter Sets, decoder operations
 * 
 * \date 25 November 2002
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Stephan Wenger        <stewe@cs.tu-berlin.de>
 ***************************************************************************************
 */


#ifndef _PARSET_H_
#define _PARSET_H_

#include "parsetcommon.h"
#include "nalucommon.h"

static const byte ZZ_SCAN[16]  =
{  0,  1,  4,  8,  5,  2,  3,  6,  9, 12, 13, 10,  7, 11, 14, 15
};

static const byte ZZ_SCAN8[64] =
{  0,  1,  8, 16,  9,  2,  3, 10, 17, 24, 32, 25, 18, 11,  4,  5,
   12, 19, 26, 33, 40, 48, 41, 34, 27, 20, 13,  6,  7, 14, 21, 28,
   35, 42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23, 30, 37, 44, 51,
   58, 59, 52, 45, 38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63
};

extern void Scaling_List(int *scalingList, int sizeOfScalingList, Boolean *UseDefaultScalingMatrix, Bitstream *s);

extern void InitVUI(seq_parameter_set_rbsp_t *sps);
extern int  ReadVUI(DataPartition *p, seq_parameter_set_rbsp_t *sps);
extern int  ReadHRDParameters(DataPartition *p, hrd_parameters_t *hrd);

extern void PPSConsistencyCheck (pic_parameter_set_rbsp_t *pps);
extern void SPSConsistencyCheck (seq_parameter_set_rbsp_t *sps);

extern void MakePPSavailable (ImageParameters *p_Img, int id, pic_parameter_set_rbsp_t *pps);
extern void MakeSPSavailable (ImageParameters *p_Img, int id, seq_parameter_set_rbsp_t *sps);

extern void ProcessSPS (ImageParameters *p_Img, NALU_t *nalu);
extern void ProcessPPS (ImageParameters *p_Img, NALU_t *nalu);

extern void UseParameterSet (Slice *currSlice, int PicParsetId);

extern void CleanUpPPS(ImageParameters *p_Img);

extern void activate_sps (ImageParameters *p_Img, seq_parameter_set_rbsp_t *sps);
extern void activate_pps (ImageParameters *p_Img, pic_parameter_set_rbsp_t *pps);

#endif
