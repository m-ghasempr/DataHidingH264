/*
***********************************************************************
* COPYRIGHT AND WARRANTY INFORMATION
*
* Copyright 2001, International Telecommunications Union, Geneva
*
* DISCLAIMER OF WARRANTY
*
* These software programs are available to the user without any
* license fee or royalty on an "as is" basis. The ITU disclaims
* any and all warranties, whether express, implied, or
* statutory, including any implied warranties of merchantability
* or of fitness for a particular purpose.  In no event shall the
* contributor or the ITU be liable for any incidental, punitive, or
* consequential damages of any kind whatsoever arising from the
* use of these programs.
*
* This disclaimer of warranty extends to the user of these programs
* and user's customers, employees, agents, transferees, successors,
* and assigns.
*
* The ITU does not represent or warrant that the programs furnished
* hereunder are free of infringement of any third-party patents.
* Commercial implementations of ITU-T Recommendations, including
* shareware, may be subject to royalty fees to patent holders.
* Information regarding the ITU-T patent policy is available from
* the ITU Web site at http://www.itu.int.
*
* THIS IS NOT A GRANT OF PATENT RIGHTS - SEE THE ITU-T PATENT POLICY.
************************************************************************
*/

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
