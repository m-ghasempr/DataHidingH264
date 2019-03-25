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
 *****************************************************************************
 *
 * \file fmo.c
 *
 * \brief
 *    Support for Flexible Macroblock Ordering: MBAmap handling
 *
 * \date
 *    16 June, 2002
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <memory.h>
#include <malloc.h>

#include "fmo.h"
#include "global.h"
#include "defines.h"
#include "header.h"

#define MAXSLICEGROUPIDS 8

int *MBAmap = NULL;   
static int PictureXSize, PictureYSize, PictureSizeInMBs;

// JVT-D097
int fmo_evlv_NewPeriod;
int slice_group_change_cycle;

static int fmo_evlv_x;
static int fmo_evlv_y;
static int fmo_evlv_left;
static int fmo_evlv_right;
static int fmo_evlv_top;
static int fmo_evlv_bottom;
static int fmo_evlv_directx;
static int fmo_evlv_directy;

static int FmoGenerateType3MBAmap (int NumSliceGroups, int XSize, int YSize, int *MBAmap);
static int FmoBoxoutCounterClockwise (int XSize, int YSize, int *MBAmap);
static int FmoBoxoutClockwise (int XSize, int YSize, int *MBAmap);
static int FmoRasterScan(int XSize, int YSize, int *MBAmap);
static int FmoInverseRasterScan(int XSize, int YSize, int *MBAmap);
static int FmoWipeRight(int XSize, int YSize, int *MBAmap);
static int FmoWipeLeft(int XSize, int YSize, int *MBAmap);

// End JVT-D097

static int FmoGenerateDefaultMap (int XSize, int YSize, int *MBAmap);

static int FmoGenerateType0MBAmap (int NumSliceGroupIDs, int *run_length,int XSize, int YSize, int *MBAmap);
static int FmoGenerateType1MBAmap (int NumSliceGroups, int XSize, int YSize, int *MBAmap);
int FmoGenerateType2MBAmap (int NumSliceGroups, int *MapData, int xs, int ys, int *MBAmap);

static int FirstMBInSlice[MAXSLICEGROUPIDS];

/*!
 *****************************************************************************
 *  How does a MBAmap look like?
 *
 *  An MBAmap is a one-diemnsional array of ints.  Each int 
 *  represents an MB in scan order.  A zero or positive value represents
 *  a slice group ID.  Negative values are reserved for future extensions.
 *  The numbering range for the SliceGroupIDs is 0..7 as per JVT-C167.
 *
 *  This module contains a static variable MBAmap.  This is the MBAmap of the
 *  picture currently coded.  It can be accessed only through the access
 *  functions.
 *****************************************************************************
*/



/*!
 ************************************************************************
 * \brief
 *    FmoInit: Initializes the FMO Module (called once per sequence)
 *
 * \par Input:
 *    currParSet: the ParameterSet to be used to initialize the FMO 
 *                module
 *
 * \par
 *    if we ever decide to use some hierarchy for Parameter Sets, this
 *    would be (at least) the picture level.
 ************************************************************************
 */

int FmoInit (int xs, int ys, int NumSliceGroups, int FmoMode, int* MapData)

{
  int i;

  for (i=0;i<MAXSLICEGROUPIDS; i++)
    FirstMBInSlice[i] = -1;
  PictureXSize = xs;
  PictureYSize = ys;
  PictureSizeInMBs = PictureXSize*PictureYSize;
  if (MBAmap != NULL)
    free (MBAmap);
  if ((MBAmap = malloc (sizeof (int) * PictureSizeInMBs)) == NULL)
  {
    printf ("Cannot allocate %d bytes for the MBAmap, exit\n", PictureSizeInMBs * sizeof (int));
    exit (-1);
  }

  // a couple of asserts to check out things:
  assert (PictureXSize == input->img_width/16);
  assert (PictureYSize == input->img_height/16);

  assert (NumSliceGroups >= 0);    // must have at least one slice group
  assert (NumSliceGroups < 8);     // cannot have more than 8 slice groups

  if (NumSliceGroups == 0)         // One Slice Group == No FMO, scan order slices
  {
//    printf ("No FMO\n\n");
    FmoGenerateDefaultMap (PictureXSize, PictureYSize, MBAmap);
  }
  else
  {
    switch (FmoMode)
    {
    case 0:
      FmoGenerateType0MBAmap (NumSliceGroups, MapData, xs, ys, MBAmap);
      break;
    case 1:
      FmoGenerateType1MBAmap (NumSliceGroups, xs, ys, MBAmap);
      break;
    case 2:
      FmoGenerateType2MBAmap (NumSliceGroups, MapData, xs, ys, MBAmap);
      break;

    // JVT-D095, JVT-D097
    case 3:
      FmoGenerateType3MBAmap (NumSliceGroups, xs, ys, MBAmap);
      break;
    case 4:
    case 5:
    case 6:
      assert(NumSliceGroups == 1);
      FmoInitEvolvingMBAmap (FmoMode, PictureXSize, PictureYSize, MBAmap); // need to be updated before coding of each picture
      break;
    // End JVT-D095, JVT-D097

    default:
      printf ("Illegal FmoMode %d , exit \n", FmoMode);
      exit (-1);
    }
  }

/*
{
int x, y;
printf ("FMOInit: Using MBAmap as follows\n");
for (y=0;y<PictureYSize; y++) {
for (x=0; x<PictureXSize;x++) printf ("%d ", MBAmap [y*xs+x]);
printf ("\n"); }
printf ("\n");
}
*/

  return 0;
}

void FmoUninit()
{
  free(MBAmap);
}

/*!
 ************************************************************************
 * \brief
 *    FmoStartPicture: initializes FMO at the begin of each new picture
 *
 * \par Input:
 *    None
 ************************************************************************
 */

int FmoStartPicture ()
{
  int i;

  assert (MBAmap != NULL);
  assert (PictureXSize == input->img_width/16);
  assert (PictureYSize == input->img_height/16);

  for (i=0; i<MAXSLICEGROUPIDS; i++)
    FirstMBInSlice[i] = FmoGetFirstMBOfSliceGroup (i);
  return 0;
}



/*!
 ************************************************************************
 * \brief
 *    FmoEndPicture: Ends the Scattered Slices Module (called once
 *    per picture)
 *
 * \par Input:
 *    None
 ************************************************************************
 */

int FmoEndPicture ()
{
  // Do nothing
  return 0;
}


/*!
 ************************************************************************
 * \brief
 *    FmoMB2Slice: Returns SliceID for a given MB
 *
 * \par Input:
 *    Macroblock Nr (in scan order)
 ************************************************************************
 */


int FmoMB2SliceGroup (int mb)
{
  assert (mb < PictureSizeInMBs);
  assert (MBAmap != NULL);
  return MBAmap[mb];
}

/*!
 ************************************************************************
 * \brief
 *    FmoGetNextMBBr: Returns the MB-Nr (in scan order) of the next
 *    MB in the (FMO) Slice, -1 if the SliceGroup is finished
 *
 * \par Input:
 *    None
 ************************************************************************
 */


int FmoGetNextMBNr (int CurrentMbNr)
{
  int SliceGroupID=FmoMB2SliceGroup (CurrentMbNr);
  
  while (++CurrentMbNr<PictureSizeInMBs && MBAmap[CurrentMbNr] != SliceGroupID)
    ;
  if (mref==mref_fld) //KS: dirty hack - field coding
  {
    if (CurrentMbNr >= PictureSizeInMBs / 2)
      return -1;    // No further MB in this slice group 
    else
      return CurrentMbNr;
  }
  else
  {
    if (CurrentMbNr >= PictureSizeInMBs)
      return -1;    // No further MB in this slice group 
    else
      return CurrentMbNr;
  }
}

/*!
 ************************************************************************
 * \brief
 *    FmoGetFirstMBOfSliceGroup: Returns the MB-Nr (in scan order) of the 
 *    next first MB of the Slice group, -1 if no such MB exists
 *
 * \par Input:
 *    SliceGroupID: Id of SliceGroup
 ************************************************************************
 */

int FmoGetFirstMBOfSliceGroup (int SliceGroupID)
{
  int i = 0;
  while ((i<PictureSizeInMBs) && (FmoMB2SliceGroup (i) != SliceGroupID))
    i++;
  if (i < PictureSizeInMBs)
    return i;
  else
    return -1;
}


/*!
 ************************************************************************
 * \brief
 *    FmoGetLastCodedMBOfSlice: Returns the MB-Nr (in scan order) of 
 *    the last MB of the slice group
 *
 * \par Input:
 *    SliceGroupID
 * \par Return
 *    MB Nr in case of success (is always >= 0)
 *    -1 if the SliceGroup doesn't exist
 ************************************************************************
 */

int FmoGetLastCodedMBOfSliceGroup (int SliceGroupID)
{
  int i;
  int LastMB = -1;

  for (i=0; i<PictureSizeInMBs; i++)
    if (FmoMB2SliceGroup (i) == SliceGroupID)
      LastMB = i;
  return LastMB;
}


void FmoSetLastMacroblockInSlice (int mb)
{
  // called by terminate_slice(), writes the last processed MB into the
  // FirstMBInSlice[MAXSLICEGROUPIDS] array.  FmoGetFirstMacroblockInSlice()
  // uses this info to identify the first uncoded MB in each slice group
  
  int currSliceGroup = FmoMB2SliceGroup (mb);
  assert (mb >= 0);
  mb = FmoGetNextMBNr (mb);   // The next (still uncoded) MB, or -1 if SG is finished
  FirstMBInSlice[currSliceGroup] = mb;
}

int FmoGetFirstMacroblockInSlice (int SliceGroup)
{
  return FirstMBInSlice[SliceGroup];
  // returns the first uncoded MB in each slice group, -1 if there is no
  // more to do in this slice group
}


int FmoSliceGroupCompletelyCoded(int SliceGroupID)
{
  if (FmoGetFirstMacroblockInSlice (SliceGroupID) < 0)  // slice group completelty coded or not present
    return TRUE;
  else
    return FALSE;
}

/*
 ************************************************************************
 ************************************************************************
 ************************************************************************
 * 
 * MBAmap generation support
 *
 * The current WD (JVT-C167) defines three types of MBAmap generation
 * mechanisms: slice interleaving (type 0), dispersed MB allocation
 * (type 1) and fully flexible, explicit MBAmap transmission (type 2).
 *
 * This section contains code for the generation of all the MBAmaps 
 * types.
 *
 * The config file (in case of the encoder) or the binary transmission of
 * the ParameterSets (in case of the decoder) contains information about 
 * the type, the run/length for type 0, and an explicit MBAmap for type 2.  
 * And, of course, also the number of slice groups.  This information is 
 * sufficient to build the MBAmap for the internal operation of the 
 * encoder and decoder.
 * 
 * I'm unsure whether this could should reside in the common section.
 *
 ************************************************************************
 ************************************************************************
 ************************************************************************
*/



/*      
 ************************************************************************
 * \brief
 *    FmoGenerateDefaultMBAmap: generates a default MBAmap (scan order)
 *
 * \par Input:
 *    XSize, Ysize:       Size of MBAmap (== picture) in macroblocks
 *    MBAMap:             MBAmap to be filled
 *
 * \par Return
 *    always 0 (success)
 *
 ************************************************************************
 */


int FmoGenerateDefaultMap (int XSize, int YSize, int *MBAmap)
{
  assert (MBAmap != NULL);
  memset (MBAmap, 0, sizeof (int) * XSize * YSize);
  return 0;
}

  
    
/*      
 ************************************************************************
 * \brief
 *    FmoGenerateType0MBAmap: generates a type 0 MBAmap (Slice 
 *    interleaveing)
 *
 * \par Input:
 *    NumSliceGroupIDs:   the number of slice groups to be used
 *    run_length:         an array of run_len structures, see JVT-C167
 *    XSize, Ysize:       Size of MBAmap (== picture) in macroblocks
 *    MBAMap:             MBAmap to be filled
 *
 * \par Return
 *    always 0 (success)
 *
 * \par Side effects
 *    
 *
 * Syntax of FMO config run-lenght file:
 *
 * a maximum MAXIMUM_FMO_INITIALIZATION_INTEGERS lines that can be 
 * fscanned() by the string "(%d,%d)\n".  The first int is the run,
 * the second the length.  The following two line example yields slice
 * interlaving for a QCIF picture:
 * (0,11)
 * (0,11)
 ************************************************************************
 */
 
int FmoGenerateType0MBAmap (int NumSliceGroupIDs, 
                            int *run_length,
                            int XSize, int YSize,
                            int *MBAmap)
{
  // in the encoder:
  // open file
  // read trun-length
  // interpret them
  // side effect: generate the bit stream (or at least the integer) for 
  // later transmission in the Parset NALU
/*
  FILE *f;
  int run [MAXIMUM_FMO_INITIALIZATION_INTEGERS];
  int len [MAXIMUM_FMO_INITIALIZATION_INTEGERS];
  int ValidRunLength = 0, i;
*/
  printf ("Type 0 MBAmap generation not yet implemented\n");
  exit (-1);
/*
  if ((f = fopen (input->FmoConfigFileName, "r")) == NULL)
  {
    printf ("Cannot open FMO config file %d for reading, FMO type %d\n",
      input->FmoConfigFileName, input->FmoType);
    exit (-1);
  }
  while (2 == fscanf (f, "(%d,%d)\n", run[ValidRunLength], len[ValidRunLength]))
    ValidRunLength++;
  fclose (f);

{
printf ("FMO config file of type 0 contains %d run-len info as follows\n", ValidRunLength);
for (i=0; i<ValidRunLength;i++)
printf ("(%d,%d)\n", run[i], len[i]);
printf ("\n");
}
 
  // Now all the info is present to generate MBAmap (and the bit stream for the
  // Parameter Set


*/
  return 0;
}


/*      
 ************************************************************************
 * \brief
 *    FmoGenerateType1MBAmap: generates a type 1 MBAmap
 *
 * \par
 *    Type 1 MBAmaps are scattered slices.  See JVT-C167 section 8.2.3 
 *    for the algorith description and formula
 *
 * \par Input:
 *    NumSliceGroupIDs:   the number of slice groups to be used
 *    XSize, Ysize:       Size of MBAmap (== picture) in macroblocks
 *    MBAMap:             MBAmap to be filled
 * \par Return
 *    always 0 (success)
 ************************************************************************
 */

int FmoGenerateType1MBAmap (int NumSliceGroups, 
                            int XSize, int YSize,
                            int *MBAmap)
{
  int x;                      // Macroblock Number
  int n = XSize;              // Number of columns
  int p = NumSliceGroups+1;     // Number of Slice Groups

  for (x=0; x<XSize*YSize; x++)
  {
    MBAmap[x] = ((x/XSize)&1)?      // even or odd?
      ((x%n)+1+p/2)%p:              // odd
      ((x%n)+1) % p;                // even
  }

/*
{
int i;
printf ("FMO MBAmap of type 1 as follows\n");
for (i=0; i<XSize*YSize;i++){
if (i%XSize==0) printf ("\n");
printf ("%d ", MBAmap[i]);}
printf ("\n");
}
*/
  return 0;
}


/*      
 ************************************************************************
 * \brief
 *    FmoGenerateType2MBAmap: generates a type 2 MBAmap
 *
 * \par
 *    Type 2 MBAmaps are fully flexible.  See JVT-C167 
 * \par Input:
 *    XSize, Ysize:       Size of MBAmap (== picture) in macroblocks
 *    MBAMap:             MBAmap to be filled
 * \par Return
 *    always 0 (success)
 ************************************************************************
 */

int FmoGenerateType2MBAmap (int NumSliceGroups, int *MapData, int xs, int ys, int *MBAmap)
{
  printf ("Type 2 MBAmap generation not yet implemented\n");
  exit (-1);

/*
  FILE *f;

  if ((f = fopen (input->FmoConfigFileName, "r")) == NULL)
  {
    printf ("Cannot open FMO config file %d for reading, FMO type %d\n",
      input->FmoConfigFileName, input->FmoType);
    exit (-1);
  }
*/
  return 0;
}

// JVT-D095
int FmoGenerateType3MBAmap (int NumSliceGroups, int XSize, int YSize, int *MBAmap)
{
  int x, y, xx;
  int n = XSize;              // Number of columns
  int p = NumSliceGroups+1;     // Number of Slice Groups

  int rx0, rx1, ry0, ry1;   // coordinates of the rectangule

  assert (NumSliceGroups == 1);

  rx0 = input->top_left_mb%n;
  ry0 = input->top_left_mb/n;
  rx1 = input->bottom_right_mb%n;
  ry1 = input->bottom_right_mb/n;

  for (y=0; y<YSize; y++)
  for (x=0; x<XSize; x++)
  {
    xx = y*XSize+x;
    if(x >= rx0 && x <= rx1 && y >= ry0 && y<= ry1) // within the rectangular slice group
      MBAmap[xx] = 0;
    else
      MBAmap[xx] = 1;
  }
  return 0;
}
// End JVT-D095

// JVT-D097
int FmoInitEvolvingMBAmap (int FmoMode, int XSize, int YSize, int *MBAmap)
{
  int i;

  assert (MBAmap != NULL);
  for(i=0; i< XSize*YSize; i++) 
    MBAmap[i] = 1;

  fmo_evlv_NewPeriod = 0;
  slice_group_change_cycle = -1;

  switch(FmoMode)
  {
  case 4:
    if(input->slice_group_change_direction == 0)
    {
      fmo_evlv_x = XSize / 2;
      fmo_evlv_y = YSize / 2;
      fmo_evlv_directx = -1;
      fmo_evlv_directy = 0;
    }
    else 
    {
      fmo_evlv_x = ( XSize - 1 ) / 2;
      fmo_evlv_y = ( YSize - 1) / 2;
      fmo_evlv_directx = 0;
      fmo_evlv_directy = 1;
    }
    fmo_evlv_left = fmo_evlv_x;
    fmo_evlv_right = fmo_evlv_x;
    fmo_evlv_top = fmo_evlv_y;
    fmo_evlv_bottom = fmo_evlv_y;
    break;
  case 5:
    if(input->slice_group_change_direction == 0)
    {
      fmo_evlv_x = 0;
      fmo_evlv_y = 0;
    }
    else
    {
      fmo_evlv_x = XSize-1;
      fmo_evlv_y = YSize-1;
    }
    break;
  case 6:
    if(input->slice_group_change_direction == 0)
    {
      fmo_evlv_x = 0;
      fmo_evlv_y = 0;
    }
    else
    {
      fmo_evlv_x = XSize-1;
      fmo_evlv_y = YSize-1;
    }
    break;
  }

  return 0;
}

int FmoUpdateEvolvingMBAmap (int FmoMode, int XSize, int YSize, int *MBAmap)
{
  switch(FmoMode)
  {
  case 4:
    if(input->slice_group_change_direction == 0)
      FmoBoxoutClockwise (XSize, YSize, MBAmap);
    else 
      FmoBoxoutCounterClockwise (XSize, YSize, MBAmap);
    break;
  case 5:
    if(input->slice_group_change_direction == 0)
      FmoRasterScan (XSize, YSize, MBAmap);
    else
      FmoInverseRasterScan (XSize, YSize, MBAmap);
    break;
  case 6:
    if(input->slice_group_change_direction == 0)
      FmoWipeRight (XSize, YSize, MBAmap);
    else
      FmoWipeLeft (XSize, YSize, MBAmap);
    break;
  }

/*
{
int xx, yy;
for (yy=0;yy<YSize; yy++) {
for (xx=0; xx<XSize;xx++) printf ("%d ", MBAmap [yy*XSize+xx]);
printf ("\n"); }
printf ("\n");
}
*/

  return 0;
}

int FmoWipeLeft(int XSize, int YSize, int *MBAmap)
{
  int i;
  int x = fmo_evlv_x;
  int y = fmo_evlv_y;

  slice_group_change_cycle++;

  for(i=0; i<input->slice_group_change_rate_minus1+1; i++)
  {
    // update the MBAmap unit of the MB (x,y)
    MBAmap[y*XSize+x] = 0;

    // go to the next MB
    if(y > 0) y--;
    else if(x > 0)
    {
      y = YSize-1;
      x--;
    }
    else 
    {
      fmo_evlv_NewPeriod = 1;
      break;
    }
  }
  fmo_evlv_x = x;
  fmo_evlv_y = y;

  return 0;
}

int FmoWipeRight(int XSize, int YSize, int *MBAmap)
{
  int i;
  int x = fmo_evlv_x;
  int y = fmo_evlv_y;

  slice_group_change_cycle++;

  for(i=0; i<input->slice_group_change_rate_minus1+1; i++)
  {
    // update the MBAmap unit of the MB (x,y)
    MBAmap[y*XSize+x] = 0;

    // go to the next MB
    if(y < YSize-1) y++;
    else if(x < XSize-1)
    {
      y = 0;
      x++;
    }
    else 
    {
      fmo_evlv_NewPeriod = 1;
      break;
    }
  }
  fmo_evlv_x = x;
  fmo_evlv_y = y;

  return 0;
}

int FmoInverseRasterScan(int XSize, int YSize, int *MBAmap)
{
  int i;
  int nextMBnum = fmo_evlv_y * XSize + fmo_evlv_x;

  slice_group_change_cycle++;

  for(i=0; i<input->slice_group_change_rate_minus1+1; i++)
  {
    // update the next MBAmap unit
    MBAmap[nextMBnum] = 0;

    // go to the next MB
    nextMBnum--;
    // check whether passed already the last MB in the evolving period
    if( nextMBnum < 0 ) 
    {
      fmo_evlv_NewPeriod = 1;
      break;
    }
  }
  fmo_evlv_x = nextMBnum%XSize;
  fmo_evlv_y = nextMBnum/XSize;

  return 0;
}

int FmoRasterScan(int XSize, int YSize, int *MBAmap)
{
  int i;
  int nextMBnum = fmo_evlv_y * XSize + fmo_evlv_x;

  slice_group_change_cycle++;

  for(i=0; i<input->slice_group_change_rate_minus1+1; i++)
  {
    // update the next MBAmap unit
    MBAmap[nextMBnum] = 0;

    // go to the next MB
    nextMBnum++;
    // check whether passed already the last MB in the evolving period
    if( nextMBnum >= XSize*YSize ) 
    {
      fmo_evlv_NewPeriod = 1;
      break;
    }
  }
  fmo_evlv_x = nextMBnum%XSize;
  fmo_evlv_y = nextMBnum/XSize;

  return 0;
}

int FmoBoxoutCounterClockwise (int XSize, int YSize, int *MBAmap)
{
  int i;
  int W = XSize, H = YSize;
  
  int x = fmo_evlv_x;
  int y = fmo_evlv_y;
  int left = fmo_evlv_left;
  int right = fmo_evlv_right;
  int top = fmo_evlv_top;
  int bottom = fmo_evlv_bottom;
  int directx = fmo_evlv_directx;
  int directy = fmo_evlv_directy;

  slice_group_change_cycle++;

  for(i=0; i<input->slice_group_change_rate_minus1+1; i++)
  {
    // update the MBAmap unit of the MB (x,y)
    MBAmap[y*XSize+x] = 0;

    // go to the next mb (x, y)
    if ( directx == -1 && directy == 0 )
    {
      if (x > left) x--;
      else if (x == 0)
      {
        y = bottom + 1;
        bottom++;
        directx = 1;
        directy = 0;
      }
      else if (x == left)
      {
        x--;
        left--;
        directx = 0;
        directy = 1;
      }
    }
    else if ( directx == 1 && directy == 0 )
    {
      if (x < right) x++;
      else if (x == W - 1)
      {
        y = top - 1;
        top--;
        directx = -1;
        directy = 0;
      }
      else if (x == right)
      {
        x++;
        right++;
        directx = 0;
        directy = -1;
      }
    }
    else if ( directx == 0 && directy == -1 )
    {
      if ( y > top) y--;
      else if (y == 0)
      {
        x = left - 1;
        left--;
        directx = 0;
        directy = 1;
      }
      else if (y == top)
      {
        y--;
        top--;
        directx = -1;
        directy = 0;
      }
    }
    else if ( directx == 0 && directy == 1 )
    {
      if (y < bottom) y++;
      else if (y == H - 1)
      {
        x = right+1;
        right++;
        directx = 0;
        directy = -1;
      }
      else if (y == bottom)
      {
        y++;
        bottom++;
        directx = 1;
        directy = 0;
      }
    }

    // check whether passed already the last MB in the evolving period
    if( !(left >= 0 && right < W && top >= 0 && bottom < H) ) 
    {
      fmo_evlv_NewPeriod = 1;
      break;
    }
  }

  fmo_evlv_x = x;
  fmo_evlv_y = y;
  fmo_evlv_left = left;
  fmo_evlv_right = right;
  fmo_evlv_top = top;
  fmo_evlv_bottom = bottom;
  fmo_evlv_directx = directx;
  fmo_evlv_directy = directy;

  return 0;
}

int FmoBoxoutClockwise (int XSize, int YSize, int *MBAmap)
{
  int i;
  int W = XSize, H = YSize;
  
  int x = fmo_evlv_x;
  int y = fmo_evlv_y;
  int left = fmo_evlv_left;
  int right = fmo_evlv_right;
  int top = fmo_evlv_top;
  int bottom = fmo_evlv_bottom;
  int directx = fmo_evlv_directx;
  int directy = fmo_evlv_directy;

  slice_group_change_cycle++;

  for(i=0; i<input->slice_group_change_rate_minus1+1; i++)
  {
    // update the MBAmap unit of the MB (x,y)
    MBAmap[y*XSize+x] = 0;

    // go to the next mb (x, y)
    if ( directx == -1 && directy == 0 )
    {
      if (x > left) x--;
      else if (x == 0)
      {
        y = top - 1;
        top--;
        directx = 1;
        directy = 0;
      }
      else if (x == left)
      {
        x--;
        left--;
        directx = 0;
        directy = -1;
      }
    }
    else if ( directx == 1 && directy == 0 )
    {
      if (x < right) x++;
      else if (x == W - 1)
      {
        y = bottom + 1;
        bottom++;
        directx = -1;
        directy = 0;
      }
      else if (x == right)
      {
        x++;
        right++;
        directx = 0;
        directy = 1;
      }
    }
    else if ( directx == 0 && directy == -1 )
    {
      if ( y > top) y--;
      else if (y == 0)
      {
        x = right + 1;
        right++;
        directx = 0;
        directy = 1;
      }
      else if (y == top)
      {
        y--;
        top--;
        directx = 1;
        directy = 0;
      }
    }
    else if ( directx == 0 && directy == 1 )
    {
      if (y < bottom) y++;
      else if (y == H - 1)
      {
        x = left - 1;
        left--;
        directx = 0;
        directy = -1;
      }
      else if (y == bottom)
      {
        y++;
        bottom++;
        directx = -1;
        directy = 0;
      }
    }

    // check whether passed already the last MB in the evolving period
    if( !(left >= 0 && right < W && top >= 0 && bottom < H) ) 
    {
      fmo_evlv_NewPeriod = 1;
      break;
    }
  }

  fmo_evlv_x = x;
  fmo_evlv_y = y;
  fmo_evlv_left = left;
  fmo_evlv_right = right;
  fmo_evlv_top = top;
  fmo_evlv_bottom = bottom;
  fmo_evlv_directx = directx;
  fmo_evlv_directy = directy;

  return 0;
}
// End JVT-D097
