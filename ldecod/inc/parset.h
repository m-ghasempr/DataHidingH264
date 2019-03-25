
/*!
 **************************************************************************************
 * \file
 *    parset.h
 * \brief
 *    Picture and Sequence Parameter Sets, decoder operations
 *    This code reflects JVT version xxx
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

void PPSConsistencyCheck (pic_parameter_set_rbsp_t *pps);
void SPSConsistencyCheck (seq_parameter_set_rbsp_t *sps);

void MakePPSavailable (int id, pic_parameter_set_rbsp_t *pps);
void MakeSPSavailable (int id, seq_parameter_set_rbsp_t *sps);

void ProcessSPS (NALU_t *nalu);
void ProcessPPS (NALU_t *nalu);

void UseParameterSet (int PicParsetId);

void activate_sps (seq_parameter_set_rbsp_t *sps);
void activate_pps (pic_parameter_set_rbsp_t *pps);

#endif
