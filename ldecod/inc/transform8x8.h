/*!
 ***************************************************************************
 *
 * \file transform8x8.h
 *
 * \brief
 *    prototypes of 8x8 transform functions
 *
 * \date
 *    9. October 2003
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Yuri Vatis
 **************************************************************************/

#ifndef _TRANSFORM8X8_H_
#define _TRANSFORM8X8_H_

int    intrapred8x8(Macroblock *currMB, ColorPlane pl, struct img_par *img, int b8);
void   itrans8x8(ColorPlane pl, struct img_par *img, int ioff, int joff);

#endif
