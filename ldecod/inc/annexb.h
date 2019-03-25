
/*!
 *************************************************************************************
 * \file annexb.h
 *
 * \brief
 *    Annex B byte stream buffer handling.
 *
 *************************************************************************************
 */

#ifndef _ANNEXB_H_
#define _ANNEXB_H_

#include "nalucommon.h"

int  GetAnnexbNALU (NALU_t *nalu);
void OpenAnnexBFile (char *fn);
void CloseAnnexBFile(void);

#endif

