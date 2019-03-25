
/*!
 ***************************************************************************
 *
 * \file fmo.h
 *
 * \brief
 *    Support for Flexilble Macroblock Ordering (FMO)
 *
 * \date
 *    19 June, 2002
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 **************************************************************************/

#ifndef _FMO_H_
#define _FMO_H_


extern int fmo_init (ImageParameters *p_Img);
extern int FmoFinit (ImageParameters *p_Img);

extern int FmoGetNumberOfSliceGroup(ImageParameters *p_Img);
extern int FmoGetLastMBOfPicture   (ImageParameters *p_Img);
extern int FmoGetLastMBInSliceGroup(ImageParameters *p_Img, int SliceGroup);
extern int FmoGetSliceGroupId      (ImageParameters *p_Img, int mb);
extern int FmoGetNextMBNr          (ImageParameters *p_Img, int CurrentMbNr);

#endif
