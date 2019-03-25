
/*!
 ***************************************************************************
 *
 * \file fmo.h
 *
 * \brief
 *    Support for Flexible Macroblock Ordering
 *
 * \date
 *    16 June 2002
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 **************************************************************************/

#ifndef _FMO_H_
#define _FMO_H_

extern int  FmoInit                       (ImageParameters * p_Img, pic_parameter_set_rbsp_t * pps, seq_parameter_set_rbsp_t * sps);
extern void FmoUninit                     (ImageParameters *p_Img);
extern int  FmoFinit                      (seq_parameter_set_rbsp_t * sps);
extern int  FmoMB2SliceGroup              (ImageParameters *p_Img, int mb);
extern int  FmoGetFirstMBOfSliceGroup     (ImageParameters *p_Img, int SliceGroupID);
extern int  FmoGetFirstMacroblockInSlice  (ImageParameters *p_Img, int SliceGroup);
extern int  FmoGetNextMBNr                (ImageParameters *p_Img, int CurrentMbNr);
extern int  FmoGetPreviousMBNr            (ImageParameters *p_Img, int CurrentMbNr);
extern int  FmoGetLastCodedMBOfSliceGroup (ImageParameters *p_Img, int SliceGroupID);
extern int  FmoStartPicture               (ImageParameters *p_Img);
extern int  FmoEndPicture                 (void);
extern int  FmoSliceGroupCompletelyCoded  (ImageParameters *p_Img, int SliceGroupID);
extern void FmoSetLastMacroblockInSlice   (ImageParameters *p_Img, int mb);


#endif
