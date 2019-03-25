
/*!
 ***************************************************************************
 *
 * \file intrarefresh.h
 *
 * \brief
 *    Pseudo-Raqndom Intra macroblock refresh support
 *
 * \date
 *    16 June 2002
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 **************************************************************************/

#ifndef _INTRAREFRESH_H_
#define _INTRAREFRESH_H_

extern void RandomIntraInit(ImageParameters *p_Img, int xsize, int ysize, int refresh);
extern void RandomIntraUninit(ImageParameters *p_Img);
extern int  RandomIntra (ImageParameters *p_Img, int mb);   //! returns 1 for MBs that need forced Intra
extern void RandomIntraNewPicture (ImageParameters *p_Img);  //! to be called once per picture


#endif //_INTRAREFRESH_H_
