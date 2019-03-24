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

static int *MBAmap = NULL;   
static int PictureXSize, PictureYSize, PictureSizeInMBs;

static int FmoGenerateDefaultMap (int XSize, int YSize, int *MBAmap);

static int FmoGenerateType0MBAmap (int NumSliceGroupIDs, int *run_length,int XSize, int YSize, int *MBAmap);
static int FmoGenerateType1MBAmap (int NumSliceGroups, 
                                   int XSize, int YSize,
                                   int *MBAmap);
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
    //printf ("No FMO\n\n");
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



