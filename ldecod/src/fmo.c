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
 * \file scatter.c
 *
 * \brief
 *    Support for Flexible Macroblock Ordering (FMO)
 *
 * \date
 *    19 June, 2002
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


#include "global.h"
#include "elements.h"
#include "defines.h"
#include "header.h"
#include "fmo.h"

static int *MBAmap;   

static int PictureXSize, PictureYSize, PictureSizeInMBs;
static int NumberOfSliceGroups;		// the number of slice groups -1 (0 == scan order, 7 == maximum)



/*!
 ************************************************************************
 * \brief
 *    FmoInit: Initializes (called always when a Parameter Set with a
 *    new MBAmap is used).  Current simplification: we have only one
 *    ParameterSet, hence it is sufficient to run this function whenever
 *    a new PUP is received.  this needs to be fixed when using multiple
 *    Parameter Sets (including a consistsency check that the MBAmap
 *    stays identical within a picture)
 *
 * \par Input:
 *    xs, ys: size of picture (in MBs, 11, 9 for QCIF)
 *    pattern: scatter pattern, or NULL for default pattern
 *
 ************************************************************************
 */

int FmoInit (int xs, int ys, int NewMBAmap[], int SizeOfNewMBAmap)
{
  int i;

  NumberOfSliceGroups = 0;
  PictureXSize = xs;
  PictureYSize = ys;
  PictureSizeInMBs = xs*ys;

// printf ("In FMOInit\n");

  if (MBAmap != NULL)
    free (MBAmap);
  if ((MBAmap = malloc (PictureSizeInMBs * sizeof (int))) == NULL)
  {
    printf ("cannot allocated %d bytes for MBAmap, exit\n", PictureSizeInMBs * sizeof (int));
    exit (-1);
  }

  if (NewMBAmap == NULL || SizeOfNewMBAmap == 0)    // no new MBAmap ->default MBAmap
  {
    memset (MBAmap, 0, PictureSizeInMBs * sizeof (int));
    return 0;
  }

  for (i=0; i<PictureSizeInMBs; i++)
  {
    MBAmap[i] = NewMBAmap[i%SizeOfNewMBAmap];
	if (MBAmap[i] > NumberOfSliceGroups)
      NumberOfSliceGroups = MBAmap[i];
  }

/*
printf ("FmoInit: Using MBAmap as follows\n");
for (y=0;y<ys; y++) {
for (x=0; x<xs;x++) printf ("%d ", MBAmap [y*xs+x]);
printf ("\n"); }
printf ("\n");
*/
  return 0;
}


/*!
 ************************************************************************
 * \brief
 *    FmoFinit: Terminates FMO
 *
 ************************************************************************
 */

int FmoFinit()
{
  if (MBAmap != NULL)
    free (MBAmap);
  return 0;
}



/*!
 ************************************************************************
 * \brief
 *    FmoStartPicture: 
 *
 * \par Input:
 *    None
 ************************************************************************
 */

void FmoStartPicture()
{
  //nothing so far
}

/*!
 ************************************************************************
 * \brief
 *    FmoEndPicture: 
 *
 * \par Input:
 *    None
 ************************************************************************
 */

void FmoEndPicture()
{
  //nothing so far
}


/*!
 ************************************************************************
 * \brief
 *    FmoGetNumberOfSliceGroup() 
 *
 * \par Input:
 *    None
 ************************************************************************
 */
int FmoGetNumberOfSliceGroup()
{
  return NumberOfSliceGroups;
}


/*!
 ************************************************************************
 * \brief
 *    FmoGetLastMBOfPicture() 
 *    returns the macroblock number of the last MB in a picture.  This
 *    mb happens to be the last macroblock of the picture if there is only
 *    one slice group
 *
 * \par Input:
 *    None
 ************************************************************************
 */

int FmoGetLastMBOfPicture(int structure)
{
  return FmoGetLastMBInSliceGroup (FmoGetNumberOfSliceGroup(), structure);
}


/*!
 ************************************************************************
 * \brief
 *    FmoGetLastMBInSliceGroup: Returns MB number of last MB in SG
 *
 * \par Input:
 *    SliceGroupID (0 to 7)
 ************************************************************************
 */

int FmoGetLastMBInSliceGroup (int SliceGroup, int structure)
{
  int i;
  
  // KS: dirty hack for interlace
  int PictureSize = (structure?PictureSizeInMBs/2:PictureSizeInMBs);

  for (i=PictureSize-1; i>=0; i--)
    if (FmoMB2Slice (i) == SliceGroup)
      return i;
  return -1;

};


/*!
 ************************************************************************
 * \brief
 *    FmoMB2Slice: Returns SliceGroupID for a given MB
 *
 * \par Input:
 *    Macroblock Nr (in scan order)
 ************************************************************************
 */

int FmoMB2Slice (int mb)
{
  assert (mb < PictureSizeInMBs);
  assert (MBAmap != NULL);
  return MBAmap[mb];
}

/*!
 ************************************************************************
 * \brief
 *    FmoGetNextMBBr: Returns the MB-Nr (in scan order) of the next
 *    MB in the (scattered) Slice, -1 if the slice is finished
 *
 * \par Input:
 *    CurrentMbNumber
 ************************************************************************
 */


int FmoGetNextMBNr (int CurrentMbNr, int structure)
{
  int SliceGroup = FmoMB2Slice (CurrentMbNr);
  // KS: dirty hack for interlace
  int PictureSize = (structure?PictureSizeInMBs/2:PictureSizeInMBs);
  
  while (++CurrentMbNr<PictureSize && MBAmap [CurrentMbNr] != SliceGroup)
    ;
  if (CurrentMbNr >= PictureSize)
    return -1;    // No further MB in this slice (could be end of picture)
  else
    return CurrentMbNr;
}

#if 0

/*!
 ************************************************************************
 * \brief
 *    ScatterGetFirstMBOfSlice: Returns the MB-Nr (in scan order) of the 
 *    next first MB of the scattered Slice, -1 if no such slice is exists
 *
 * \par Input:
 *    SliceID: Id of Slice
 ************************************************************************
 */

int ScatterGetFirstMBOfSlice (int SliceID, int structure)
{
  int i = 0;
  // KS: dirty hack for interlace
  int PictureSize = (structure?PictureSizeInMBs/2:PictureSizeInMBs);

  while ((i<PictureSize) && (ScatterMB2Slice (i) != SliceID))
    i++;
  if (i < PictureSize)
    return i;
  else
    return -1;
}


/*!
 ************************************************************************
 * \brief
 *    ScatterFindFirstUnassignedMB: Returns the MB-Nr (in scan order) of the 
 *    first MB that is not yet assigned to a slice-1 if no such MB exists
 *
 * \par Input:
 *    None
 ************************************************************************
 */

int ScatterFindFirstUnassignedMB(int structure)
{
  return ScatterGetFirstMBOfSlice (-1, structure);
}


/*!
 ************************************************************************
 * \brief
 *    ScatterGetLastCodedMBOfSlice: Returns the MB-Nr (in scan order) of 
 *    the last MB of the slice
 *
 * \par Input:
 *    SliceID
 * \par Return
 *    MB Nr in case of success (is always >= 0)
 *    -1 if the slice doesn't exist
 *    -2 if last MB is unknown (not yet used)
 ************************************************************************
 */

int ScatterGetLastCodedMBOfSlice (int SliceID, int structure)
{
  int i;
  int LastMB = -1;
  // KS: dirty hack for interlace
  int PictureSize = (structure?PictureSizeInMBs/2:PictureSizeInMBs);

  for (i=0; i<PictureSize; i++)
    if (ScatterMB2Slice (i) == SliceID)
      LastMB = i;
  return LastMB;
}




/************************************************************************
/************************************************************************
 *
 *  Pseudo-Random Intra MB walk-around refresh
 *
/************************************************************************
/************************************************************************
/************************************************************************/

static int *RefreshPattern;
static int *IntraMBs;
static int WalkAround = 0;
static int NumberOfMBs = 0;
static int NumberIntraPerPicture;
 
/*!
 ************************************************************************
 * \brief
 *    RandomIntraInit: Initializes Random Intra module.  Should be called
 *    only after initialization (or changes) of the picture size or the
 *    random intra refresh value.  In version jm1.7 it is impossible to
 *    change those values on-the-fly, hence RandomIntraInit should be
 *    called immediately after the parsing of the config file
 *
 * \par Input:
 *    xsize, ysize: size of the picture (in MBs)
 *    refresh     : refresh rate in MBs per picture
 ************************************************************************
 */

void RandomIntraInit(int xsize, int ysize, int refresh)
{
  int i, pos;

  srand (1);      // A fixed random initializer to make things reproducable
  NumberOfMBs = xsize * ysize;
  NumberIntraPerPicture = refresh;

  RefreshPattern = malloc (sizeof (int) * NumberOfMBs);
  assert (RefreshPattern != NULL);

  IntraMBs = malloc (sizeof (int) * refresh);
  assert (IntraMBs != NULL);

  for (i= 0; i<NumberOfMBs; i++)
    RefreshPattern[i] = -1;

   for (i=0; i<NumberOfMBs; i++)
   {
     do
     {
       pos = rand() % NumberOfMBs;
     } while (RefreshPattern [pos] != -1);
     RefreshPattern [pos] = i;
   }
/*
for (i=0; i<NumberOfMBs; i++) printf ("%d\t", RefreshPattern[i]);
getchar();
*/
}

/*!
 ************************************************************************
 * \brief
 *    RandomIntra: Code an MB as Intra?
 *
 * \par Input
 *    MacroblockNumberInScanOrder
 * \par Output
 *    1 if an MB should be forced to Intra, according the the 
 *      RefreshPattern
 *    0 otherwise
 *
 ************************************************************************
 */

int RandomIntra (int mb)
{
  int i;

  for (i=0; i<NumberIntraPerPicture; i++)
    if (IntraMBs[i] == mb)
      return 1;
  return 0;
}


/*!
 ************************************************************************
 * \brief
 *    RandomIntraNewPicture: Selects new MB for forced Intra
 *
 * \par
 *    This function should be called exactly once per picture, and 
 *    requires a finished initialization 
 *
 ************************************************************************
 */

void RandomIntraNewPicture ()
{
  int i, j;

  WalkAround += NumberIntraPerPicture;
  for (j=0,i=WalkAround; j<NumberIntraPerPicture; j++, i++)
    IntraMBs[j] = RefreshPattern [i%NumberOfMBs];
}
#endif
