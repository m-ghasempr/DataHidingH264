
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

#define MAXIMUM_FMO_INITIALIZATION_INTEGERS 1000

int FmoInit (int xs, int ys, int NumSliceGroups, int FmoMode, int* MapData);
void FmoUninit ();
int FmoFinit ();
int FmoMB2SliceGroup (int mb);
int FmoGetFirstMBOfSliceGroup (int SliceGroupID);
int FmoGetFirstMacroblockInSlice (int SliceGroup);
int FmoGetNextMBNr (int CurrentMbNr);
int FmoGetLastCodedMBOfSliceGroup (int SliceGroupID);
int FmoStartPicture ();
int FmoEndPicture();
int FmoSliceGroupCompletelyCoded(int SliceGroupID);
void FmoSetLastMacroblockInSlice (int mb);

// JVT-D097
int FmoInitEvolvingMBAmap (int FmoMode, int XSize, int YSize, int *MBAmap);
int FmoUpdateEvolvingMBAmap (int FmoMode, int XSize, int YSize, int *MBAmap);

extern int *MBAmap; 
extern int fmo_evlv_NewPeriod;
extern int slice_group_change_cycle;
// End JVT-D097

#endif
