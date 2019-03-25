/*!
 **************************************************************************************
 * \file
 *    annexb.h
 * \brief
 *    Byte stream operations support
 *
 *  \date 7 December 2002
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Stephan Wenger        <stewe@cs.tu-berlin.de>
 ***************************************************************************************
 */

#ifndef _ANNEXB_H_
#define _ANNEXB_H_

#include "nalucommon.h"

extern int WriteAnnexbNALU (VideoParameters *p_Vid, NALU_t *n);
extern void OpenAnnexbFile (VideoParameters *p_Vid, char *fn);
extern void CloseAnnexbFile(VideoParameters *p_Vid);


#endif //_ANNEXB_H_
