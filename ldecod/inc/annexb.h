
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

typedef struct annex_b_struct 
{
  int  BitStreamFile;                //!< the bit stream file
  byte *iobuffer;
  byte *iobufferread;
  unsigned int bytesinbuffer;
  int is_eof;

  int IsFirstByteStreamNALU;
  int nextstartcodebytes;
  byte *Buf;  
} ANNEXB_t;

extern int  GetAnnexbNALU  (ImageParameters *p_Img, NALU_t *nalu);
extern void OpenAnnexBFile (ImageParameters *p_Img, char *fn);
extern void CloseAnnexBFile(ImageParameters *p_Img);
extern void init_annex_b(ANNEXB_t *annex_b);
extern void malloc_annex_b(ImageParameters *p_Img);
extern void free_annex_b(ImageParameters *p_Img);

#endif

