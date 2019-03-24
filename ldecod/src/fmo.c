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

int *MBAmap;   

static int PictureXSize, PictureYSize, PictureSizeInMBs;
static int NumberOfSliceGroups;    // the number of slice groups -1 (0 == scan order, 7 == maximum)

// JVT-D095, JVT-D097
static int FmoGenerateType3MBAmap (struct img_par *img, struct inp_par *inp, int NumSliceGroups, int XSize, int YSize, int *MBAmap);
static int FmoBoxoutCounterClockwise (struct img_par *img, struct inp_par *inp, int XSize, int YSize, int *MBAmap);
static int FmoBoxoutClockwise (struct img_par *img, struct inp_par *inp, int XSize, int YSize, int *MBAmap);
static int FmoRasterScan(struct img_par *img, struct inp_par *inp, int XSize, int YSize, int *MBAmap);
static int FmoInverseRasterScan(struct img_par *img, struct inp_par *inp, int XSize, int YSize, int *MBAmap);
static int FmoWipeRight(struct img_par *img, struct inp_par *inp, int XSize, int YSize, int *MBAmap);
static int FmoWipeLeft(struct img_par *img, struct inp_par *inp, int XSize, int YSize, int *MBAmap);
// End JVT-D095, JVT-D097


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

int FmoInit (struct img_par *img, struct inp_par *inp, int xs, int ys, int NewMBAmap[], int SizeOfNewMBAmap)
{
  int i, FmoMode;

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

  // JVT-D095, JVT-D097
  NumberOfSliceGroups = img->num_slice_groups_minus1;
  FmoMode = img->mb_allocation_map_type;
  if (NumberOfSliceGroups && FmoMode >= 3)
  {
    switch (FmoMode)
    {
    case 3:
      FmoGenerateType3MBAmap (img, inp, NumberOfSliceGroups, xs, ys, MBAmap);
      break;
    case 4:
    case 5:
    case 6:
      assert(NumberOfSliceGroups == 1);
      //FmoInitEvolvingMBAmap (img, inp, FmoMode, PictureXSize, PictureYSize, MBAmap); // need to be updated before coding of each picture
      break;
    default:
      printf ("Illegal FmoMode %d , exit \n", FmoMode);
      exit (-1);
    }
  }
  // End JVT-D095, JVT-D097

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
  int PictureSize = (structure?(structure == 3 ? PictureSizeInMBs:PictureSizeInMBs/2):PictureSizeInMBs);

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
  int PictureSize = PictureSize = (structure?(structure == 3 ? PictureSizeInMBs:PictureSizeInMBs/2):PictureSizeInMBs);
  
  while (++CurrentMbNr<PictureSize && MBAmap [CurrentMbNr] != SliceGroup)
    ;
  if (CurrentMbNr >= PictureSize)
    return -1;    // No further MB in this slice (could be end of picture)
  else
    return CurrentMbNr;
}


// JVT-D095
int FmoGenerateType3MBAmap (struct img_par *img, struct inp_par *inp, int NumSliceGroups, int XSize, int YSize, int *MBAmap)
{
  int x, y, xx;
  int n = XSize;              // Number of columns
  int p = NumSliceGroups+1;     // Number of Slice Groups

  int rx0, rx1, ry0, ry1;   // coordinates of the rectangule

  assert (NumSliceGroups == 1);

  rx0 = img->top_left_mb%n;
  ry0 = img->top_left_mb/n;
  rx1 = img->bottom_right_mb%n;
  ry1 = img->bottom_right_mb/n;

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
int FmoUpdateEvolvingMBAmap (struct img_par *img, struct inp_par *inp, int *MBAmap)
{
  int FmoMode = img->mb_allocation_map_type;
  int XSize = img->width/16;
  int YSize = img->height/16;
  int i;

  for (i=0; i<YSize*XSize; i++)
      MBAmap[i] = 1;

  switch(FmoMode)
  {
  case 4:
    if(img->slice_group_change_direction == 0)
      FmoBoxoutClockwise (img, inp, XSize, YSize, MBAmap);
    else 
      FmoBoxoutCounterClockwise (img, inp, XSize, YSize, MBAmap);
    break;
  case 5:
    if(img->slice_group_change_direction == 0)
      FmoRasterScan (img, inp, XSize, YSize, MBAmap);
    else
      FmoInverseRasterScan (img, inp, XSize, YSize, MBAmap);
    break;
  case 6:
    if(img->slice_group_change_direction == 0)
      FmoWipeRight (img, inp, XSize, YSize, MBAmap);
    else
      FmoWipeLeft (img, inp, XSize, YSize, MBAmap);
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

int FmoWipeLeft(struct img_par *img, struct inp_par *inp, int XSize, int YSize, int *MBAmap)
{
  int i, x, y, n;

  x = XSize-1; 
  y = YSize-1;

  n = (img->slice_group_change_cycle+1) * (img->slice_group_change_rate_minus1+1);

  for(i=0; i<n; i++)
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
      break;
  }

  return 0;
}

int FmoWipeRight(struct img_par *img, struct inp_par *inp, int XSize, int YSize, int *MBAmap)
{
  int i, x, y, n;

  x = 0; 
  y = 0;

  n = (img->slice_group_change_cycle+1) * (img->slice_group_change_rate_minus1+1);

  for(i=0; i<n; i++)
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
      break;
  }

  return 0;
}

int FmoInverseRasterScan(struct img_par *img, struct inp_par *inp, int XSize, int YSize, int *MBAmap)
{
  int i, x, y, n, nextMBnum;

  x = XSize -1; 
  y = YSize -1;

  n = (img->slice_group_change_cycle+1) * (img->slice_group_change_rate_minus1+1);
  nextMBnum = y * XSize + x;

  for(i=0; i<n; i++)
  {
    // update the next MBAmap unit
    MBAmap[nextMBnum] = 0;

    // go to the next MB
    nextMBnum--;
    // check whether passed already the last MB in the evolving period
    if( nextMBnum < 0 ) 
      break;
  }

  return 0;
}

int FmoRasterScan(struct img_par *img, struct inp_par *inp, int XSize, int YSize, int *MBAmap)
{
  int i, x, y, n, nextMBnum;

  x = 0; 
  y = 0;

  n = (img->slice_group_change_cycle+1) * (img->slice_group_change_rate_minus1+1);
  nextMBnum = y * XSize + x;

  for(i=0; i<n; i++)
  {
    // update the next MBAmap unit
    MBAmap[nextMBnum] = 0;

    // go to the next MB
    nextMBnum++;
    // check whether passed already the last MB in the evolving period
    if( nextMBnum >= XSize*YSize ) 
      break;
  }

  return 0;
}

int FmoBoxoutCounterClockwise (struct img_par *img, struct inp_par *inp, int XSize, int YSize, int *MBAmap)
{
  int i, n;
  int W = XSize, H = YSize;
  
  int x = ( XSize - 1 ) / 2;
  int y = ( YSize - 1) / 2;
  int directx = 0;
  int directy = 1;
  int left = x;
  int right = x;
  int top = y;
  int bottom = y;

  n = (img->slice_group_change_cycle+1) * (img->slice_group_change_rate_minus1+1);

  for(i=0; i<n; i++)
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
      break;
  }

  return 0;
}

int FmoBoxoutClockwise (struct img_par *img, struct inp_par *inp, int XSize, int YSize, int *MBAmap)
{
  int i, n;
  int W = XSize, H = YSize;
  
  int x = XSize / 2;
  int y = YSize / 2;
  int directx = -1;
  int directy = 0;
  int left = x;
  int right = x;
  int top = y;
  int bottom = y;

  n = (img->slice_group_change_cycle+1) * (img->slice_group_change_rate_minus1+1);

  for(i=0; i<n; i++)
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
      break;
  }

  return 0;
}
// End JVT-D097


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
  int PictureSize = PictureSize = (structure?(structure == 3 ? PictureSizeInMBs:PictureSizeInMBs/2):PictureSizeInMBs);;

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
  int PictureSize = PictureSize = (structure?(structure == 3 ? PictureSizeInMBs:PictureSizeInMBs/2):PictureSizeInMBs);;

  for (i=0; i<PictureSize; i++)
    if (ScatterMB2Slice (i) == SliceID)
      LastMB = i;
  return LastMB;
}




/************************************************************************
 ************************************************************************
 *
 *  Pseudo-Random Intra MB walk-around refresh
 *
 ************************************************************************
 ************************************************************************
 ************************************************************************/

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
