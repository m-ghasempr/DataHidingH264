/*!
 *****************************************************************************
 *
 * \file fmo.c
 *
 * \brief
 *    Support for Flexible Macroblock Ordering for different Slice Group Modes: MBAmap handling
 *
 * \date
 *    16 June, 2002  Modified April 25, 2004
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 *    Dong Wang        Dong.Wang@bristol.ac.uk
 *
 *****************************************************************************/

/*!
 ****************************************************************************
 *   Notes by Dong Wang (April 25 2004)
 *
 *  Source codes are modified to support 7 slice group types (fmo modes).
 *  The functions for generating map are very similar to that in decoder, but have
 *  a little difference.
 *
 *  The MB map is calculated at the beginning of coding of each picture (frame or field).
 *
 *  'slice_group_change_cycle' in structure 'ImageParameters' is the syntax in the slice
 *  header. It's set to be 1 before the initialization of FMO in function code_a_picture().
 *  It can be changed every time if needed.
 *
 ****************************************************************************
 */

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

//#define PRINT_FMO_MAPS  1


#include "global.h"

#include "fmo.h"


static void FmoGenerateType0MapUnitMap (ImageParameters * p_Img, pic_parameter_set_rbsp_t * pps);
static void FmoGenerateType1MapUnitMap (ImageParameters * p_Img, pic_parameter_set_rbsp_t * pps);
static void FmoGenerateType2MapUnitMap (ImageParameters * p_Img, pic_parameter_set_rbsp_t * pps);
static void FmoGenerateType3MapUnitMap (ImageParameters * p_Img, pic_parameter_set_rbsp_t * pps);
static void FmoGenerateType4MapUnitMap (ImageParameters * p_Img, pic_parameter_set_rbsp_t * pps);
static void FmoGenerateType5MapUnitMap (ImageParameters * p_Img, pic_parameter_set_rbsp_t * pps);
static void FmoGenerateType6MapUnitMap (ImageParameters * p_Img, pic_parameter_set_rbsp_t * pps);


static int FmoGenerateMapUnitToSliceGroupMap (ImageParameters * p_Img, pic_parameter_set_rbsp_t * pps);
static int FmoGenerateMBAmap                 (ImageParameters * p_Img, seq_parameter_set_rbsp_t* sps);


/*!
 ************************************************************************
 * \brief
 *    Generates p_Img->MapUnitToSliceGroupMap 
 *
 * \param p_Img
 *    Image Parameter to be used for map generation
 * \param pps
 *    Picture Parameter set to be used for map generation
 *
 ************************************************************************
 */
static int FmoGenerateMapUnitToSliceGroupMap (ImageParameters * p_Img, pic_parameter_set_rbsp_t * pps)
{
  p_Img->PicSizeInMapUnits = p_Img->PicHeightInMapUnits * p_Img->PicWidthInMbs;


  if (pps->slice_group_map_type == 6)
  {
    if ((pps->pic_size_in_map_units_minus1+1) != p_Img->PicSizeInMapUnits)
    {
      error ("wrong pps->pic_size_in_map_units_minus1 for used SPS and FMO type 6", 500);
    }
  }

  // allocate memory for p_Img->MapUnitToSliceGroupMap 
  if (p_Img->MapUnitToSliceGroupMap)
    free (p_Img->MapUnitToSliceGroupMap);

  if ((p_Img->MapUnitToSliceGroupMap = malloc ((p_Img->PicSizeInMapUnits) * sizeof (byte))) == NULL)
  {
    printf ("cannot allocated %d bytes for p_Img->MapUnitToSliceGroupMap , exit\n", (int) ( p_Img->PicSizeInMapUnits * sizeof (byte)));
    exit (-1);
  }

  if (pps->num_slice_groups_minus1 == 0)    // only one slice group
  {
    memset (p_Img->MapUnitToSliceGroupMap , 0,  p_Img->PicSizeInMapUnits * sizeof (byte));
    return 0;
  }

  switch (pps->slice_group_map_type)
  {
  case 0:
    FmoGenerateType0MapUnitMap (p_Img, pps);
    break;
  case 1:
    FmoGenerateType1MapUnitMap (p_Img, pps);
    break;
  case 2:
    FmoGenerateType2MapUnitMap (p_Img, pps);
    break;
  case 3:
    FmoGenerateType3MapUnitMap (p_Img, pps);
    break;
  case 4:
    FmoGenerateType4MapUnitMap (p_Img, pps);
    break;
  case 5:
    FmoGenerateType5MapUnitMap (p_Img, pps);
    break;
  case 6:
    FmoGenerateType6MapUnitMap (p_Img, pps);
    break;
  default:
    printf ("Illegal slice_group_map_type %d , exit \n", pps->slice_group_map_type);
    exit (-1);
  }
  return 0;
}


/*!
 ************************************************************************
 * \brief
 *    Generates MBAmap from p_Img->MapUnitToSliceGroupMap 
 *
 * \param p_Img
 *    Image Parameter to be used for map generation
 * \param sps
 *    Sequence Parameter set to be used for map generation
 *
 ************************************************************************
 */
static int FmoGenerateMBAmap (ImageParameters * p_Img, seq_parameter_set_rbsp_t* sps)
{
  unsigned i;

  // allocate memory for p_Img->MBAmap
  if (p_Img->MBAmap)
    free (p_Img->MBAmap);


  if ((p_Img->MBAmap = malloc ((p_Img->PicSizeInMbs) * sizeof (byte))) == NULL)
  {
    printf ("cannot allocated %d bytes for p_Img->MBAmap, exit\n", (int) ((p_Img->PicSizeInMbs) * sizeof (byte)));
    exit (-1);
  }

  if ((sps->frame_mbs_only_flag) || p_Img->field_picture)
  {
    for (i=0; i<p_Img->PicSizeInMbs; i++)
    {
      p_Img->MBAmap[i] = p_Img->MapUnitToSliceGroupMap [i];
    }
  }
  else
    if (sps->mb_adaptive_frame_field_flag  &&  (! p_Img->field_picture))
    {
      for (i=0; i<p_Img->PicSizeInMbs; i++)
      {
        p_Img->MBAmap[i] = p_Img->MapUnitToSliceGroupMap [i >> 1];
      }
    }
    else
    {
      for (i=0; i<p_Img->PicSizeInMbs; i++)
      {
        p_Img->MBAmap[i] = p_Img->MapUnitToSliceGroupMap [(i/(2*p_Img->PicWidthInMbs))*p_Img->PicWidthInMbs+(i%p_Img->PicWidthInMbs)];
      }
    }
    return 0;
}


/*!
 ************************************************************************
 * \brief
 *    FMO initialization: Generates p_Img->MapUnitToSliceGroupMap  and p_Img->MBAmap.
 *
 * \param p_Img
 *    Image Parameter to be used for map generation
 * \param pps
 *    Picture Parameter set to be used for map generation
 * \param sps
 *    Sequence Parameter set to be used for map generation
 ************************************************************************
 */
int FmoInit(ImageParameters * p_Img, pic_parameter_set_rbsp_t * pps, seq_parameter_set_rbsp_t * sps)
{

#ifdef PRINT_FMO_MAPS
  unsigned i,j;
  int bottom;
#endif

  int k;
  for (k = 0; k < MAXSLICEGROUPIDS; k++)
    p_Img->FirstMBInSlice[k] = -1;

  FmoGenerateMapUnitToSliceGroupMap(p_Img, pps);
  FmoGenerateMBAmap(p_Img, sps);

#ifdef PRINT_FMO_MAPS
  printf("\n");
  printf("FMO Map (Units):\n");

  for (j=0; j<p_Img->PicHeightInMapUnits; j++)
  {
    for (i=0; i<p_Img->PicWidthInMbs; i++)
    {
      printf("%d ",p_Img->MapUnitToSliceGroupMap [i+j*p_Img->PicWidthInMbs]);
    }
    printf("\n");
  }
  printf("\n");

  if(sps->mb_adaptive_frame_field_flag==0)
  {
    printf("FMO Map (Mb):\n");
    for (j=0; j<(p_Img->PicSizeInMbs/p_Img->PicWidthInMbs); j++)
    {
      for (i=0; i<p_Img->PicWidthInMbs; i++)
      {
        printf("%d ",p_Img->MBAmap[i+j*p_Img->PicWidthInMbs]);
      }
      printf("\n");
    }
    printf("\n");
  }
  else
  {
    printf("FMO Map (Mb in scan order for MBAFF):\n");
    for (j=0; j<(p_Img->PicSizeInMbs/p_Img->PicWidthInMbs); j++)
    {
      for (i=0; i<p_Img->PicWidthInMbs; i++)
      {
        bottom=(j & 0x01);
        printf("%d ",p_Img->MBAmap[(j-bottom)*p_Img->PicWidthInMbs+i*2+bottom]);
      }
      printf("\n");

    }
    printf("\n");

  }

#endif

  return 0;
}


/*!
 ************************************************************************
 * \brief
 *    Free memory if allocated by FMO functions
 ************************************************************************
 */
void FmoUninit(ImageParameters *p_Img)
{
  if (p_Img->MBAmap)
  {
    free (p_Img->MBAmap);
    p_Img->MBAmap = NULL;
  }
  if (p_Img->MapUnitToSliceGroupMap )
  {
    free (p_Img->MapUnitToSliceGroupMap );
    p_Img->MapUnitToSliceGroupMap  = NULL;
  }

}


/*!
 ************************************************************************
 * \brief
 *    Generate interleaved slice group map type MapUnit map (type 0)
 *
 * \param p_Img
 *    Image Parameter to be used for map generation
 * \param pps
 *    Picture Parameter set to be used for map generation
 ************************************************************************
 */
static void FmoGenerateType0MapUnitMap (ImageParameters * p_Img, pic_parameter_set_rbsp_t * pps )
{
  unsigned iGroup, j;
  unsigned i = 0;
  do
  {
    for( iGroup = 0;
    (iGroup <= pps->num_slice_groups_minus1) && (i < p_Img->PicSizeInMapUnits);
    i += pps->run_length_minus1[iGroup++] + 1)
    {
      for( j = 0; j <= pps->run_length_minus1[ iGroup ] && i + j < p_Img->PicSizeInMapUnits; j++ )
        p_Img->MapUnitToSliceGroupMap [i+j] = (byte)  iGroup;
    }
  }
  while( i < p_Img->PicSizeInMapUnits );
}


/*!
 ************************************************************************
 * \brief
 *    Generate dispersed slice group map type MapUnit map (type 1)
 *
 * \param p_Img
 *    Image Parameter to be used for map generation
 * \param pps
 *    Picture Parameter set to be used for map generation
 ************************************************************************
 */
static void FmoGenerateType1MapUnitMap (ImageParameters * p_Img, pic_parameter_set_rbsp_t * pps )
{
  unsigned i;
  for( i = 0; i < p_Img->PicSizeInMapUnits; i++ )
  {
    p_Img->MapUnitToSliceGroupMap [i] = (byte) (((i%p_Img->PicWidthInMbs)+(((i/p_Img->PicWidthInMbs)*(pps->num_slice_groups_minus1+1))>>1))
      %(pps->num_slice_groups_minus1+1));
  }
}

/*!
 ************************************************************************
 * \brief
 *    Generate foreground with left-over slice group map type MapUnit map (type 2)
 *
 * \param p_Img
 *    Image Parameter to be used for map generation
 * \param pps
 *    Picture Parameter set to be used for map generation
 ************************************************************************
 */
static void FmoGenerateType2MapUnitMap (ImageParameters * p_Img, pic_parameter_set_rbsp_t * pps )
{
  int iGroup;
  unsigned i, x, y;
  unsigned yTopLeft, xTopLeft, yBottomRight, xBottomRight;

  for( i = 0; i < p_Img->PicSizeInMapUnits; i++ )
    p_Img->MapUnitToSliceGroupMap [ i ] = (byte) pps->num_slice_groups_minus1;

  for( iGroup = pps->num_slice_groups_minus1 - 1 ; iGroup >= 0; iGroup-- )
  {
    yTopLeft = pps->top_left[ iGroup ] / p_Img->PicWidthInMbs;
    xTopLeft = pps->top_left[ iGroup ] % p_Img->PicWidthInMbs;
    yBottomRight = pps->bottom_right[ iGroup ] / p_Img->PicWidthInMbs;
    xBottomRight = pps->bottom_right[ iGroup ] % p_Img->PicWidthInMbs;
    for( y = yTopLeft; y <= yBottomRight; y++ )
      for( x = xTopLeft; x <= xBottomRight; x++ )
        p_Img->MapUnitToSliceGroupMap [ y * p_Img->PicWidthInMbs + x ] = (byte) iGroup;
  }
}


/*!
 ************************************************************************
 * \brief
 *    Generate box-out slice group map type MapUnit map (type 3)
 *
 * \param p_Img
 *    Image Parameter to be used for map generation
 * \param pps
 *    Picture Parameter set to be used for map generation
 ************************************************************************
 */
static void FmoGenerateType3MapUnitMap (ImageParameters * p_Img, pic_parameter_set_rbsp_t * pps )
{
  unsigned i, k;
  int leftBound, topBound, rightBound, bottomBound;
  int x, y, xDir, yDir;
  int mapUnitVacant;

  unsigned mapUnitsInSliceGroup0 = imin((pps->slice_group_change_rate_minus1 + 1) * p_Img->slice_group_change_cycle, p_Img->PicSizeInMapUnits);

  for( i = 0; i < p_Img->PicSizeInMapUnits; i++ )
    p_Img->MapUnitToSliceGroupMap [ i ] = 2;

  x = ( p_Img->PicWidthInMbs - pps->slice_group_change_direction_flag ) / 2;
  y = ( p_Img->PicHeightInMapUnits - pps->slice_group_change_direction_flag ) / 2;

  leftBound   = x;
  topBound    = y;
  rightBound  = x;
  bottomBound = y;

  xDir =  pps->slice_group_change_direction_flag - 1;
  yDir =  pps->slice_group_change_direction_flag;

  for( k = 0; k < p_Img->PicSizeInMapUnits; k += mapUnitVacant )
  {
    mapUnitVacant = ( p_Img->MapUnitToSliceGroupMap [ y * p_Img->PicWidthInMbs + x ]  ==  2 );
    if( mapUnitVacant )
      p_Img->MapUnitToSliceGroupMap [ y * p_Img->PicWidthInMbs + x ] = (byte) ( k >= mapUnitsInSliceGroup0 );

    if( xDir  ==  -1  &&  x  ==  leftBound )
    {
      leftBound = imax( leftBound - 1, 0 );
      x = leftBound;
      xDir = 0;
      yDir = 2 * pps->slice_group_change_direction_flag - 1;
    }
    else
      if( xDir  ==  1  &&  x  ==  rightBound )
      {
        rightBound = imin( rightBound + 1, (int)p_Img->PicWidthInMbs - 1 );
        x = rightBound;
        xDir = 0;
        yDir = 1 - 2 * pps->slice_group_change_direction_flag;
      }
      else
        if( yDir  ==  -1  &&  y  ==  topBound )
        {
          topBound = imax( topBound - 1, 0 );
          y = topBound;
          xDir = 1 - 2 * pps->slice_group_change_direction_flag;
          yDir = 0;
        }
        else
          if( yDir  ==  1  &&  y  ==  bottomBound )
          {
            bottomBound = imin( bottomBound + 1, (int)p_Img->PicHeightInMapUnits - 1 );
            y = bottomBound;
            xDir = 2 * pps->slice_group_change_direction_flag - 1;
            yDir = 0;
          }
          else
          {
            x = x + xDir;
            y = y + yDir;
          }
  }

}

/*!
 ************************************************************************
 * \brief
 *    Generate raster scan slice group map type MapUnit map (type 4)
 *
 * \param p_Img
 *    Image Parameter to be used for map generation
 * \param pps
 *    Picture Parameter set to be used for map generation
 ************************************************************************
 */
static void FmoGenerateType4MapUnitMap (ImageParameters * p_Img, pic_parameter_set_rbsp_t * pps )
{

  unsigned mapUnitsInSliceGroup0 = imin((pps->slice_group_change_rate_minus1 + 1) * p_Img->slice_group_change_cycle, p_Img->PicSizeInMapUnits);
  unsigned sizeOfUpperLeftGroup = pps->slice_group_change_direction_flag ? ( p_Img->PicSizeInMapUnits - mapUnitsInSliceGroup0 ) : mapUnitsInSliceGroup0;

  unsigned i;

  for( i = 0; i < p_Img->PicSizeInMapUnits; i++ )
    if( i < sizeOfUpperLeftGroup )
      p_Img->MapUnitToSliceGroupMap [ i ] = (byte) pps->slice_group_change_direction_flag;
    else
      p_Img->MapUnitToSliceGroupMap [ i ] = (byte) (1 - pps->slice_group_change_direction_flag);

}

/*!
 ************************************************************************
 * \brief
 *    Generate wipe slice group map type MapUnit map (type 5)
 *
 * \param p_Img
 *    Image Parameter to be used for map generation
 * \param pps
 *    Picture Parameter set to be used for map generation
 ************************************************************************
*/
static void FmoGenerateType5MapUnitMap (ImageParameters * p_Img, pic_parameter_set_rbsp_t * pps )
{

  unsigned mapUnitsInSliceGroup0 = imin((pps->slice_group_change_rate_minus1 + 1) * p_Img->slice_group_change_cycle, p_Img->PicSizeInMapUnits);
  unsigned sizeOfUpperLeftGroup = pps->slice_group_change_direction_flag ? ( p_Img->PicSizeInMapUnits - mapUnitsInSliceGroup0 ) : mapUnitsInSliceGroup0;

  unsigned i,j, k = 0;

  for( j = 0; j < p_Img->PicWidthInMbs; j++ )
    for( i = 0; i < p_Img->PicHeightInMapUnits; i++ )
      if( k++ < sizeOfUpperLeftGroup )
        p_Img->MapUnitToSliceGroupMap [ i * p_Img->PicWidthInMbs + j ] = (byte) pps->slice_group_change_direction_flag;
      else
        p_Img->MapUnitToSliceGroupMap [ i * p_Img->PicWidthInMbs + j ] = (byte) (1 - pps->slice_group_change_direction_flag);

}

/*!
 ************************************************************************
 * \brief
 *    Generate explicit slice group map type MapUnit map (type 6)
 *
 * \param p_Img
 *    Image Parameter to be used for map generation
 * \param pps
 *    Picture Parameter set to be used for map generation
 ************************************************************************
 */
static void FmoGenerateType6MapUnitMap (ImageParameters * p_Img, pic_parameter_set_rbsp_t * pps )
{
  unsigned i;
  for (i=0; i<p_Img->PicSizeInMapUnits; i++)
  {
    p_Img->MapUnitToSliceGroupMap [i] = (byte) pps->slice_group_id[i];
  }
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
int FmoStartPicture (ImageParameters *p_Img)
{
  int i;

  assert (p_Img->MBAmap != NULL);

  for (i=0; i<MAXSLICEGROUPIDS; i++)
    p_Img->FirstMBInSlice[i] = FmoGetFirstMBOfSliceGroup (p_Img, i);
  return 0;
}



/*!
 ************************************************************************
 * \brief
 *    FmoEndPicture: Ends the Scattered Slices Module (called once
 *    per picture).
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
int FmoMB2SliceGroup (ImageParameters *p_Img, int mb)
{
  assert (mb < (int)p_Img->PicSizeInMbs);
  assert (p_Img->MBAmap != NULL);
  return p_Img->MBAmap[mb];
}

/*!
 ************************************************************************
 * \brief
 *    FmoGetNextMBBr: Returns the MB-Nr (in scan order) of the next
 *    MB in the (FMO) Slice, -1 if the SliceGroup is finished
 *
 * \par Input:
 *    CurrentMbNr
 ************************************************************************
 */
int FmoGetNextMBNr (ImageParameters *p_Img, int CurrentMbNr)
{

  int  SliceGroupID = FmoMB2SliceGroup (p_Img, CurrentMbNr);

  while (++CurrentMbNr<(int)p_Img->PicSizeInMbs &&  p_Img->MBAmap[CurrentMbNr] != SliceGroupID)
    ;

  if (CurrentMbNr >= (int)p_Img->PicSizeInMbs)
    return -1;    // No further MB in this slice (could be end of picture)
  else
    return CurrentMbNr;
}


/*!
 ************************************************************************
 * \brief
 *    FmoGetNextMBBr: Returns the MB-Nr (in scan order) of the next
 *    MB in the (FMO) Slice, -1 if the SliceGroup is finished
 *
 * \par Input:
 *    CurrentMbNr
 ************************************************************************
 */
int FmoGetPreviousMBNr (ImageParameters *p_Img, int CurrentMbNr)
{

  int  SliceGroupID = FmoMB2SliceGroup (p_Img, CurrentMbNr);
  CurrentMbNr--;
  while (CurrentMbNr>=0 &&  p_Img->MBAmap[CurrentMbNr] != SliceGroupID)
    CurrentMbNr--;

  if (CurrentMbNr < 0)
    return -1;    // No previous MB in this slice
  else
    return CurrentMbNr;
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
int FmoGetFirstMBOfSliceGroup (ImageParameters *p_Img, int SliceGroupID)
{
  int i = 0;
  while ((i<(int)p_Img->PicSizeInMbs) && (FmoMB2SliceGroup (p_Img, i) != SliceGroupID))
    i++;

  if (i < (int)p_Img->PicSizeInMbs)
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
int FmoGetLastCodedMBOfSliceGroup (ImageParameters *p_Img, int SliceGroupID)
{
  int i;
  int LastMB = -1;

  for (i=0; i<(int)p_Img->PicSizeInMbs; i++)
    if (FmoMB2SliceGroup (p_Img, i) == SliceGroupID)
      LastMB = i;
  return LastMB;
}


void FmoSetLastMacroblockInSlice ( ImageParameters *p_Img, int mb)
{
  // called by terminate_slice(), writes the last processed MB into the
  // FirstMBInSlice[MAXSLICEGROUPIDS] array.  FmoGetFirstMacroblockInSlice()
  // uses this info to identify the first uncoded MB in each slice group

  int currSliceGroup = FmoMB2SliceGroup (p_Img, mb);
  assert (mb >= 0);
  mb = FmoGetNextMBNr (p_Img, mb);   // The next (still uncoded) MB, or -1 if SG is finished
  p_Img->FirstMBInSlice[currSliceGroup] = mb;
}

int FmoGetFirstMacroblockInSlice ( ImageParameters *p_Img, int SliceGroup)
{
  return p_Img->FirstMBInSlice[SliceGroup];
  // returns the first uncoded MB in each slice group, -1 if there is no
  // more to do in this slice group
}


int FmoSliceGroupCompletelyCoded( ImageParameters *p_Img, int SliceGroupID)
{
  if (FmoGetFirstMacroblockInSlice (p_Img, SliceGroupID) < 0)  // slice group completelty coded or not present
    return TRUE;
  else
    return FALSE;
}



