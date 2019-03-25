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
 ************************************************************************
 *  \file
 *     decodeiff.c
 *  \brief
 *     Impelentation of the interim file format 
 *     as defined in Appendix III of WD-1 (JVT-A003).
 *
 *     Instructions for a new release of the Joint Model:
 *     1. Remember to update WORKING_DRAFT_MAJOR_NO and WORKING_DRAFT_MINOR_NO
 *        in defines.h appropriately
 *     2. Remember to modify the encoder as well
 *
 *     Instructions for modifications to the interim file format:
 *     1. Update INTERIM_FILE_MAJOR_NO and INTERIM_FILE_MINOR_NO
 *        in defines.h appropriately. It is up to the implementor
 *        to decide, which modification deserves an update in the major
 *        version number.
 *
 *     The box header only permits compact size (32 bits mode)
 *     data partition is not supported.
 *     The source code supports multiple segments per file, but
 *     the feature has not been tested.
 *      
 *  \author
 *      - Dong Tian                              <tian@cs.tut.fi>
 *      - Miska M. Hannuksela                    <miska.hannuksela@nokia.com>
 ************************************************************************
 */

/*!
 ************************************************************************
 *  \flowchart
 * main()
 * {
 *   init
 *   IFFSequenceHeader(); // Read the first Boxes, and set the initial parameter set 
 *                        // according to the input file
 *   if ( inp->of_mode == PAR_OF_IFF )
 *     while ( parse_one_box() != -1 );
 *   else
 *     while (decode_one_frame() != EOS);
 * }
 * 
 * parse_one_box()
 * {
 *   call different funcitons according to the box type
 *   switch( box type )
 *   {
 *     ......
 *     case segment:
 *       parse_one_segment
 *       {
 *         while ( decode_one_frame() );
 *       }
 *     break;
 *     ......
 *   }
 * }
 ************************************************************************
 */

#include <stdio.h>
#include <assert.h>
#include <memory.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>

#include "global.h"
#include "decodeiff.h"
#include "mbuffer.h"
#include "vlc.h"

FileTypeBox box_ft;
FileHeaderBox box_fh;
ContentInfoBox box_ci;
AlternateTrackInfoBox box_ati;
ParameterSetBox box_ps;
SegmentBox box_s;
AlternateTrackHeaderBox box_ath;
PictureInformationBox box_pi;
PayloadInfo currPayloadInfo;
PictureInfo currPictureInfo, oldPictureInfo;
LayerBox box_layr[MAX_LAYER_NUMBER];
SubSequenceBox box_sseq[MAX_LAYER_NUMBER];
AlternateTrackMediaBox box_atm;
SwitchPictureBox box_sp;

int CurrentParameterSet=-1;  // Tian Dong: must not be set to 0 here.
static int BeginOfPictureOrSlice; // SOP or SOS
int IFFEndOfFile = FALSE;
int ThisAlternateTrack = 0;   // indicate which track in the file will be decoded.

int isBigEndian=0;

/*!
 ************************************************************************
 * \brief
 *      checks if the System is big- or little-endian
 * \return
 *      0, little-endian
 *      1, big-endian
 ************************************************************************
 */
int testEndian()
{
  short s;
  byte *p;

  p=(byte*)&s;

  s=1;

  return (*p==0);
}

/*!
 ************************************************************************
 * \brief
 *      Test the file if it is in Interim File Format by reading FileTypeBox
 * \return
 *      0, if it is in Interim File Format
 *      -1, otherwise
  * \param fp
 *      input file pointer
************************************************************************
 */
int testFileTypeBox(FILE* fp)
{
  int left;
  int size, x;

  assert( fp != NULL );

  box_ft.compatibleBrands = NULL;

  x = readfile( &box_ft.type.size, 4, 1, fp );
  if ( -1 == x ) return -1;
  x = readfile( &box_ft.type.type, 4, 1, fp );
  if ( -1 == x ) return -1;
  size = box_ft.type.size;
  if ( box_ft.type.type != BOX_FTYP ) return -1;

  if ( box_ft.type.size > 0 )
  {
    x = readfile( box_ft.majorBrand, 4, 1, fp );
    if ( -1 == x ) return -1;
    x = readfile( &box_ft.jmMajorVersion, 2, 1, fp );
    if ( -1 == x ) return -1;
    x = readfile( &box_ft.jmMinorVersion, 2, 1, fp );
    if ( -1 == x ) return -1;

    if ( box_ft.jmMajorVersion != WORKING_DRAFT_MAJOR_NO ||
      box_ft.jmMinorVersion != WORKING_DRAFT_MINOR_NO)    
    {
      printf( "Error: The working draft version is not supported.\n" );
      return -1;
    }

    left = (int)size - SIZEOF_BOXTYPE - 8;
    box_ft.compatibleBrands = calloc( 1, left );
    box_ft.numCompatibleBrands = left / 4;
    x = fread( box_ft.compatibleBrands, 1, left, fp );
    if ( left != x ) return -1;
  }
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      Free the memory allocated for File Type Box
 * \return
 *      None
 ************************************************************************
 */
void freeFileTypeBox()
{
  if ( box_ft.compatibleBrands != NULL )
    free( box_ft.compatibleBrands );
}

/*!
 ************************************************************************
 * \brief
 *      Read the File Header Box
 * \return
 *      0, if success
 *      -1, otherwise
 * \param fp
 *      input file pointer
 * \param size
 *      size of the box
 ************************************************************************
 */
int rdFileHeaderBox( FILE* fp, size_t size )
{
  byte cd;
  int x;
  assert( fp != NULL );

  x = readfile(&box_fh.majorVersion, 1, 1, fp);
  if ( -1 == x ) return -1;
  x = readfile(&box_fh.minorVersion, 1, 1, fp);
  if ( -1 == x ) return -1;

  if ( 
    box_fh.majorVersion != INTERIM_FILE_MAJOR_NO ||
    box_fh.minorVersion != INTERIM_FILE_MINOR_NO )
  {
    printf( "Error: The interim file format version is not supported.\n" );
    return -1;
  }

  x = readfile(&box_fh.timescale, 4, 1, fp);
  if ( -1 == x ) return -1;
  x = readfile(&box_fh.numUnitsInTick, 4, 1, fp);
  if ( -1 == x ) return -1;
  x = readfile(&box_fh.duration, 8, 1, fp);
  if ( -1 == x ) return -1;
  x = readfile(&box_fh.pixAspectRatioX, 2, 1, fp);
  if ( -1 == x ) return -1;
  x = readfile(&box_fh.pixAspectRatioY, 2, 1, fp);
  if ( -1 == x ) return -1;
  x = readfile(&box_fh.maxPicId, 2, 1, fp);
  if ( -1 == x ) return -1;
  x = readfile(&box_fh.numAlternateTracks, 1, 1, fp);
  if ( -1 == x ) return -1;
  x = readfile(&cd, 1, 1, fp);
  if ( -1 == x ) return -1;
  box_fh.numBytesInPayloadCountMinusOne = (cd >> 6);
  box_fh.numBytesInPictureOffsetMinusTwo = ((cd & 0x3F) >> 4);
  box_fh.numBytesInPictureDisplayTimeMinusOne = ((cd & 0x0f) >> 2);
  box_fh.numBytesInPictureCountMinusOne = (cd & 0x03);
  x = readfile(&cd, 1, 1, fp);
  if ( -1 == x ) return -1;
  box_fh.numBytesInPayloadSizeMinusOne = (cd >> 6);

  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      Free the memory allocated for File Header Box
 *      Do nothing
 * \return
 *      None
 ************************************************************************
 */
void freeFileHeaderBox()
{
}

/*!
 ************************************************************************
 * \brief
 *      Tian Dong, Feb 10, 2002:
 *      Just move the file pointer forward to skip the content info box
 * \return
 *      0, if success
 *      -1, otherwise
 * \param fp
 *      input file pointer
 * \param size
 *      size of the box
 ************************************************************************
 */
int rdContentInfoBox( FILE* fp, size_t size )
{
  // do nothing
  if ( 0 != fseek( fp, size - SIZEOF_BOXTYPE, SEEK_CUR ) ) return -1;
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      Free the memory allocated for Content Info Box
 *      Do nothing
 * \return
 *      None
 ************************************************************************
 */
void freeContentInfoBox()
{
}

/*!
 ************************************************************************
 * \brief
 *      Read the Alterate Track Info Box
 *      Tian Dong, Feb 10, 2002:
 *      Only one track in the input file
 * \return
 *      0, if success
 *      -1, otherwise
 * \param fp
 *      input file pointer
 * \param size
 *      size of the box
 ************************************************************************
 */
int rdAlternateTrackInfoBox( FILE* fp, unsigned long size )
{
  assert( fp != NULL );
  assert( box_fh.numAlternateTracks == 1 );
  box_ati.info = NULL;
  box_ati.info = calloc( box_fh.numAlternateTracks, sizeof(AlternateTrackInfo) );
  if ( box_ati.info == NULL ) return -1;
  if ( -1 == readfile( &box_ati.info[0].displayWindowWidth, 2, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ati.info[0].displayWindowHeight, 2, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ati.info[0].maxSDUSize, 2, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ati.info[0].avgSDUSize, 2, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ati.info[0].avgBitRate, 4, 1, fp ) ) return -1;

  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      Free the memory allocated for Alternate Track Info Box
 * \return
 *      None
 ************************************************************************
 */
void freeAlternateTrackInfoBox()
{
  if ( box_ati.info != NULL ) 
    free(box_ati.info);
}

/*!
 ************************************************************************
 * \brief
 *      Read the ParameterSet Box
 *      Tian Dong, Feb 10, 2002:
 *      Only one Parameter Set Box in the input file
 * \return
 *      0, if success
 *      -1, otherwise
 * \param fp
 *      input file pointer
 * \param size
 *      size of the box
 ************************************************************************
 */
int rdParameterSetBox( FILE* fp, unsigned long size )
{
  unsigned char cd;
  assert( fp );

  if ( -1 == readfile( &box_ps.parameterSetID, 2, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ps.profile, 1, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ps.level, 1, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ps.version, 1, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ps.pictureWidthInMBs, 2, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ps.pictureHeightInMBs, 2, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ps.displayRectangleOffsetTop, 2, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ps.displayRectangleOffsetLeft, 2, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ps.displayRectangleOffsetBottom, 2, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ps.displayRectangleOffsetRight, 2, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ps.displayMode, 1, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ps.displayRectangleOffsetFromWindowTop, 2, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ps.displayRectangleOffsetFromWindowLeftBorder, 2, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ps.loopFilterParametersFlag, 1, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ps.entropyCoding, 1, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ps.partitioningType, 1, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ps.intraPredictionType, 1, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ps.bufCycle, 1, 1, fp ) ) return -1;

  if ( -1 == readfile( &cd, 1, 1, fp ) ) return -1;
  if ( cd == 0x80 ) box_ps.requiredPictureNumberUpdateBehavior= TRUE;
  else box_ps.requiredPictureNumberUpdateBehavior = FALSE;

  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      From the current position, find the parameter set which id is given.
 *      The file pointer must point to the beginning of a box.
 * \return
 *      0, if success
 *      -1, otherwise
 * \param fp
 *      input file pointer
 * \param pos
 *      used to store the parameter set box's position
 * \param pid
 *      give the id of the wanted parameter set
 ************************************************************************
 */
int findParameterSetBox( FILE* fp, long *pos, int pid )
{
  int id, size, type;

  assert( fp );
  while( !feof( fp ) )
  {
    if ( readfile( &size, 4, 1, fp ) == -1 ) return -1;
    if ( readfile( &type, 4, 1, fp ) == -1 ) return -1;
    if ( type == BOX_PRMS )
    {
      if ( readfile( &id, 4, 1, fp ) == -1 ) return -1;
      if ( id == pid ) 
      {
        *pos = ftell( fp ) - 4;
        return 0;   // succeed
      }
      if ( 0 != fseek( fp, size-12, SEEK_CUR ) ) return -1;
    }
    else
    {
      if ( 0 != fseek( fp, size-8, SEEK_CUR ) ) return -1;
    }
  }
  return -1;
}

/*!
 ************************************************************************
 * \brief
 *      Do nothing
 * \return
 *      None
 ************************************************************************
 */
void freeParameterSetBox()
{
}

/*!
 ************************************************************************
 * \brief
 *      parse one segment
 *      Tian Dong, Feb 10, 2002:
 *      Only one Parameter Set Box in the input file
 * \return
 *      0, if success
 *      -1, otherwise
 * \param img
 *      image pointer
 * \param inp
 *      input parameter pointer
 * \param snr
 *      snr pointer
 * \param fp
 *      input file pointer
 * \param size
 *      size of the box
 ************************************************************************
 */
int parse_one_segment( struct img_par* img, struct inp_par* inp, struct snr_par *snr, FILE* fp, size_t size )
{
  unsigned long limited, storedpos;

  limited = ftell( fp );
  limited += (size - SIZEOF_BOXTYPE);

  if ( -1 == readfile( &box_s.fileSize, 8, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_s.startTick, 8, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_s.segmentDuration, 8, 1, fp ) ) return -1;
  // storedpos will save the file pointer to point to the current track header
  storedpos = ftell( fp );
  if ( -1 == find_track_meta( fp, limited, &box_ath.storedpos, ThisAlternateTrack ) ) return -1;
  if ( 0 != fseek( fp, storedpos, SEEK_SET ) ) return -1;
  if ( -1 == find_track_media( fp, limited, &box_atm.storedpos, ThisAlternateTrack ) ) return -1;

  if ( -1 == rdOneTrack( img, inp, snr, fp ) ) return -1;

  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      Free all the resources allocated for the segment
 * \return
 *      None
 ************************************************************************
 */
void freeSegmentBox()
{
  freeAlternateTrackHeaderBox();
  freeAlternateTrackMediaBox();
}

/*!
 ************************************************************************
 * \brief
 *      read and decode all the frames in one track of a segment
 * \return
 *      0, if success
 *      -1, otherwise
 * \param img
 *      image pointer
 * \param inp
 *      input parameter pointer
 * \param snr
 *      snr pointer
 * \param fp
 *      input file pointer
 ************************************************************************
 */
int rdOneTrack( struct img_par* img, struct inp_par* inp, struct snr_par *snr, FILE *fp )
{
  unsigned long numPictures, pos;
  unsigned long size, type, boxsize, readsize;
  int ret = EOS + 100;
  int layerno = 0;

  numPictures=0;

  if ( 0 != fseek( fp, box_ath.storedpos, SEEK_SET ) ) return -1; // ath
  if ( -1 == readfile( &boxsize, 4, 1, fp ) ) return -1;
  if ( boxsize == 0 ) return -1;
  if ( -1 == readfile( &type, 4, 1, fp ) ) return -1;
  assert( type == BOX_ATRH );

  if ( -1 == readfile( &box_ath.numLayers, 1, 1, fp ) ) return -1;

  pos = ftell( fp );
  boxsize = boxsize - SIZEOF_BOXTYPE - 1;
  readsize = 0;
  while ( readsize < boxsize )
  {
    if ( 0 != fseek( fp, pos, SEEK_SET ) ) return -1;
    if ( -1 == readfile( &size, 4, 1, fp ) ) return -1;
    if ( size == 0 ) return -1;
    if ( -1 == readfile( &type, 4, 1, fp ) ) return -1;
    
    switch ( type )
    {
    case BOX_PICI:
      if ( -1 == readfile_s( &numPictures, sizeof(unsigned long), box_fh.numBytesInPictureCountMinusOne+1, 1, fp ) ) return -1;
        box_ath.storedpos = ftell( fp );
      while ( numPictures-- > 0 && ret != EOS )
        ret = decode_one_frame( img, inp, snr );
      break;
    case BOX_LAYR:
      rdLayerBox(layerno++, size, fp);
      break;
    default:
      printf("Unknow boxes found in the track!\n");
      exit(0);
      break;
    }
    readsize += size;
    pos += size;
  }
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      try to find the ATRH box from the input file
 * \return
 *      0, if success
 *      -1, otherwise
 * \param fp
 *      input file pointer
 * \param limited
 *      how far will go to find it
 * \param storedpos
 *      if found, save the position 
 * \param tracknr
 *      which track to find
 ************************************************************************
 */
int find_track_meta( FILE *fp, unsigned long limited, unsigned long* storedpos, int tracknr )
{
  int nr = -1;
  unsigned long size, type, pos;

  *storedpos = 0;

  assert( fp != NULL );
  assert( tracknr < box_fh.numAlternateTracks );

  pos = 0;
  while ( tracknr > nr && pos < limited && pos != (unsigned long)-1 )
  {
    pos = (unsigned)ftell(fp);
    if ( -1 == readfile( &size, 4, 1, fp ) ) return -1;     // in 32 bits mode
    if ( -1 == readfile( &type, 4, 1, fp ) ) return -1;
    if ( type == BOX_ATRH ) nr++;
    if ( tracknr > nr )
    { if ( 0 != fseek( fp, size - SIZEOF_BOXTYPE, SEEK_CUR ) ) return -1; }
  }

  if ( tracknr != nr ) return -1;     // have not found the track needed

  *storedpos = pos;           // save the file pointer to the track header
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      try to find the ATRM box from the input file
 * \return
 *      0, if success
 *      -1, otherwise
 * \param fp
 *      input file pointer
 * \param limited
 *      how far will go to find it
 * \param storedpos
 *      if found, save the position 
 * \param tracknr
 *      which track to find
 ************************************************************************
 */
int find_track_media( FILE *fp, unsigned long limited, unsigned long* storedpos, int tracknr )
{
  int nr = -1;
  unsigned long size, type, pos;

  *storedpos = 0;

  assert( fp != NULL );
  assert( tracknr < box_fh.numAlternateTracks );

  pos = 0;
  while ( tracknr > nr && pos < limited && pos != (unsigned long)-1 )
  {
    if ( -1 == readfile( &size, 4, 1, fp ) ) return -1;     // in 32 bits mode
    if ( -1 == readfile( &type, 4, 1, fp ) ) return -1;
    pos = (unsigned)ftell(fp);
    if ( type == BOX_ATRM ) nr++;
    if ( tracknr > nr )
    { if ( 0 != fseek( fp, size - SIZEOF_BOXTYPE, SEEK_CUR ) ) return -1; }
  }

  if ( tracknr != nr ) return -1;     // have not found the track needed

  *storedpos = ftell( fp );           // save the file pointer to the track header
  box_atm.currPictureOffset = *storedpos;   // maintain the media data pointer
  box_atm.currPayloadOffset = *storedpos;
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      Do nothing
 * \return
 *      None
 ************************************************************************
 */
void freeAlternateTrackHeaderBox()
{
}

/*!
 ************************************************************************
 * \brief
 *      Read the picture info before it is decoded
 * \return
 *      0, if success
 *      -1, otherwise
 * \param fp
 *      input file pointer
 ************************************************************************
 */
int rdPictureInfo( FILE* fp )
{
  byte cd;
  size_t num = 0;
  
  assert( fp != NULL );

  memcpy( &oldPictureInfo, &currPictureInfo, sizeof(PictureInfo) );

  fseek( fp, box_ath.storedpos, SEEK_SET );

  num += readfile( &cd, 1, 1, fp );
  if ( (cd&0x80) != 0 ) currPictureInfo.intraPictureFlag = TRUE;
  else currPictureInfo.intraPictureFlag = FALSE;
  if ( (cd&0x40) != 0 ) currPictureInfo.syncPictureFlag = TRUE;
  else currPictureInfo.syncPictureFlag = FALSE;

  currPictureInfo.pictureOffset=0;
  currPictureInfo.pictureDisplayTime=0;
  currPictureInfo.numPayloads=0;

  num += readfile_s( &currPictureInfo.pictureOffset, sizeof(currPictureInfo.pictureOffset), box_fh.numBytesInPictureOffsetMinusTwo + 2, 1, fp );
  num += readfile_s( &currPictureInfo.pictureDisplayTime, sizeof(currPictureInfo.pictureDisplayTime), box_fh.numBytesInPictureDisplayTimeMinusOne+1, 1, fp );

  if (box_ath.numLayers != 0)
  {
    if ( -1 == readfile( &currPictureInfo.layerNumber, 1, 1, fp ) ) return -1;
    if ( -1 == readfile( &currPictureInfo.subSequenceIdentifier, 2, 1, fp ) ) return -1;
    if ( currPictureInfo.syncPictureFlag )
    {
      if ( -1 == readfile( &currPictureInfo.originLayerNumber, 1, 1, fp ) ) return -1;
      if ( -1 == readfile( &currPictureInfo.originSubSequenceIdentifier, 2, 1, fp ) ) return -1;
    }
  }

  num += readfile_s( &currPictureInfo.numPayloads, sizeof(currPictureInfo.numPayloads), box_fh.numBytesInPayloadCountMinusOne+1, 1, fp );

  if ( num != (unsigned)(1+box_fh.numBytesInPictureOffsetMinusTwo + 2 +
    box_fh.numBytesInPictureDisplayTimeMinusOne+1 +
    box_fh.numBytesInPayloadCountMinusOne+1) ) 
    return -1;

  currPayloadInfo.storedpos = ftell( fp );  // update the pointer, it will be increased when read one payload info
  currPayloadInfo.payloadnr = 0;
  BeginOfPictureOrSlice = SOP;
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      read the payloadinfo data and set the decoder parameter
 * \return
 *      0, if success
 *      -1, otherwise
 * \param img
 *      image pointer
 * \param inp
 *      input parameter pointer
 * \param pp
 *      read the data to pp
 * \param fp
 *      input file pointer
 ************************************************************************
 */
int rdPayloadInfo( struct img_par *img, struct inp_par* inp, PayloadInfo* pp, FILE *fp )
{
  byte cd;

  assert( fp != NULL );

  fseek( fp, pp->storedpos, SEEK_SET );

  pp->payloadSize=0;
  if ( -1 == readfile_s( &pp->payloadSize, sizeof (pp->payloadSize), box_fh.numBytesInPayloadSizeMinusOne+1, 1, fp ) ) return -1;
  if ( -1 == readfile( &pp->headerSize, 1, 1, fp ) ) return -1;
  if ( -1 == readfile( &cd, 1, 1, fp ) ) return -1;
  pp->payloadType = (cd >> 4);
  pp->errorIndication = ((cd&0x0F) >> 3);
  pp->reserved = (cd & 0x07);

  switch ( pp->payloadType )
  {
  case 0:
  case PAYLOAD_TYPE_IDERP:
    pp->buffer.bitstream_length = pp->headerSize - (box_fh.numBytesInPayloadSizeMinusOne+1 + 2);
    pp->buffer.streamBuffer = alloca( pp->buffer.bitstream_length );
    assert( pp->buffer.streamBuffer != NULL );
    if ( fread( pp->buffer.streamBuffer, 1, pp->buffer.bitstream_length, fp ) != (unsigned)pp->buffer.bitstream_length )
      return -1;
    pp->buffer.frame_bitoffset = 0;
    decomposeSliceHeader( img, inp, pp );
    break;
  default:
    // not supported now.
    assert( 0 == 1 );
    break;
  }

  pp->payloadnr++;
  if ( currPictureInfo.numPayloads == pp->payloadnr ) // the last payloadinfo is read, set the position of next pic
    box_ath.storedpos = ftell( fp );
  else // next payloadinfo 
    pp->storedpos = ftell( fp );

  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      set the decoder parameter according to a payloadinfo
 * \return
 *      None
 * \param img
 *      image pointer
 * \param inp
 *      input parameter pointer
 * \param pp
 *      read the data to pp
 ************************************************************************
 */
void decomposeSliceHeader( struct img_par *img, struct inp_par* inp, PayloadInfo* pp )
{
  SyntaxElement sym;
  Bitstream* buf;
  Slice *currSlice = img->currentSlice;
  int bitptr = 0;
  int tmp1;
  RMPNIbuffer_t *tmp_rmpni,*tmp_rmpni2;
  MMCObuffer_t *tmp_mmco,*tmp_mmco2;
  int done;
  static int last_imgtr_frm=0,modulo_ctr_frm=0,last_imgtr_fld=0,modulo_ctr_fld=0;
  static int last_imgtr_frm_b=0,modulo_ctr_frm_b=0,last_imgtr_fld_b=0,modulo_ctr_fld_b=0;

  buf = &(pp->buffer);

  sym.type = SE_HEADER;       // This will be true for all symbols generated here
  sym.mapping = linfo;        // Mapping rule: Simple code number to len/info
  
  sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
  sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
  pp->parameterSet = sym.value1;
  IFFUseParameterSet( pp->parameterSet, img, inp );
  buf->frame_bitoffset += sym.len;
  bitptr += sym.len;

  sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
  sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
  pp->structure = sym.value1;
  currSlice->structure = img->structure = sym.value1;
  buf->frame_bitoffset += sym.len;
  bitptr += sym.len;

  sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
  sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
  pp->pictureID = sym.value1;
  currSlice->picture_id = img->tr = sym.value1;
  buf->frame_bitoffset += sym.len;
  bitptr += sym.len;

  sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
  sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
  pp->sliceType = sym.value1;
  switch (sym.value1)
  {
  //! Potential BUG: do we need to distinguish between INTER_IMG_MULT and INTER_IMG?
  //!    similar with B_IMG_! and B_IMG_MULT
  //! also: need to define Slice types for SP images
  //! see VCEG-N72r1 for the Slice types, which are mapped here to img->type
  case 0:
    img->type = INTER_IMG;
    break;
  case 1:
    img->type = B_IMG;
    break;
  case 2:
    img->type = SP_IMG;
    break;
  case 3:
    img->type = currSlice->picture_type = INTRA_IMG;
    break;
  case 4:
    img->type = currSlice->picture_type = SI_IMG;
  default:
    printf ("Panic: unknown Slice type %d, conceal by loosing slice\n", pp->sliceType);
    currSlice->ei_flag = 1;
  }  
  buf->frame_bitoffset += sym.len;
  bitptr += sym.len;

  sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &(sym.inf), buf->bitstream_length );
  sym.mapping(sym.len, sym.inf, &(img->disposable_flag), &(sym.value2));
  buf->frame_bitoffset += sym.len;
  bitptr+=sym.len;

  if(img->type==INTER_IMG)
  {
    sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &(sym.inf), buf->bitstream_length );
    sym.mapping(sym.len, sym.inf, &(img->num_ref_pic_active_fwd), &(sym.value2));
    buf->frame_bitoffset += sym.len;
    bitptr+=sym.len;
    img->num_ref_pic_active_fwd++;
  }
  else if(img->type==B_IMG)
  {
 
    sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &(sym.inf), buf->bitstream_length );
    sym.mapping(sym.len, sym.inf, &(img->num_ref_pic_active_fwd), &(sym.value2));
    buf->frame_bitoffset += sym.len;
    bitptr+=sym.len;
    img->num_ref_pic_active_fwd++;

    sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &(sym.inf), buf->bitstream_length );
    sym.mapping(sym.len, sym.inf, &(img->num_ref_pic_active_bwd), &(sym.value2));
    buf->frame_bitoffset += sym.len;
    bitptr+=sym.len;
    img->num_ref_pic_active_bwd++;
  }

  if (img->type <= INTRA_IMG || img->type >= SP_IMG || !img->disposable_flag) 
  {
    if (img->structure == FRAME)
    {     
      if(img->tr <last_imgtr_frm) 
        modulo_ctr_frm++;
      
      last_imgtr_frm = img->tr;
      img->tr_frm = img->tr + (256*modulo_ctr_frm);
    }
    else
    {
      if(img->tr <last_imgtr_fld) 
        modulo_ctr_fld++;
      
      last_imgtr_fld = img->tr;
      img->tr_fld = img->tr + (256*modulo_ctr_fld);
    }
  }
  else
  {
    if (img->structure == FRAME)
    {     
      if(img->tr <last_imgtr_frm_b) 
        modulo_ctr_frm_b++;
      
      last_imgtr_frm_b = img->tr;
      img->tr_frm = img->tr + (256*modulo_ctr_frm_b);
    }
    else
    {
      if(img->tr <last_imgtr_fld_b) 
        modulo_ctr_fld_b++;
      
      last_imgtr_fld_b = img->tr;
      img->tr_fld = img->tr + (256*modulo_ctr_fld_b);
    }
  }
  if(img->type != B_IMG) {
    img->pstruct_next_P = img->structure;
    if(img->structure == TOP_FIELD)
    {
      img->imgtr_last_P = img->imgtr_next_P;
      img->imgtr_next_P = img->tr_fld;
    }
    else if(img->structure == FRAME)
    {
      img->imgtr_last_P = img->imgtr_next_P;
      img->imgtr_next_P = 2*img->tr_frm;
    }
  }
  
  if(img->type==B_IMG)
  {
    if(img->disposable_flag==0) 
    {
      if(img->structure == TOP_FIELD)
      {
        img->imgtr_last_P = img->imgtr_next_P;
        img->imgtr_next_P = img->tr_fld;
      }
      else if(img->structure == FRAME)
      {
        img->imgtr_last_P = img->imgtr_next_P;
        img->imgtr_next_P = 2*img->tr_frm;
      }
    }
  }

  img->mb_frame_field_flag = 0;

  if (img->structure == 3)
    img->mb_frame_field_flag = 1;

  sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
  sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
  pp->firstMBInSliceX = sym.value1;
  buf->frame_bitoffset += sym.len;
  bitptr += sym.len;

  sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
  sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
  pp->firstMBInSliceY = sym.value1;
  buf->frame_bitoffset += sym.len;
  bitptr += sym.len;

  currSlice->start_mb_nr = pp->firstMBInSliceY * box_ps.pictureWidthInMBs + pp->firstMBInSliceX;
  
  if ( img->type == B_IMG)
  {
    sym.len = GetfixedSymbol( buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length,1 );
    img->direct_type= pp->directType = sym.inf ;
    buf->frame_bitoffset += sym.len;
    bitptr += sym.len;
  }
  


  sym.mapping = linfo_dquant;
  sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
  sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
  pp->initialQP = sym.value1 + (MAX_QP - MIN_QP + 1)/2;
  currSlice->qp = img->qp = pp->initialQP;
  buf->frame_bitoffset += sym.len;
  bitptr += sym.len;
  if ( img->type == SP_IMG || img->type == SI_IMG)
  {
    if ( img->type != SI_IMG )
    {
      sym.len = GetfixedSymbol( buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length,1 );
      img->sp_switch=sym.inf ;
      buf->frame_bitoffset += sym.len;
      bitptr += sym.len;
    }
    sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
    sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
    img->qpsp = sym.value1 + (MAX_QP - MIN_QP + 1)/2;
    buf->frame_bitoffset += sym.len;
    bitptr += sym.len;
  }

  if (inp->LFParametersFlag)
  {
    sym.len = GetfixedSymbol( buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length,1 );
    currSlice->LFDisableIdc = sym.inf;
    buf->frame_bitoffset += sym.len;
    bitptr+=1;

    if (!currSlice->LFDisableIdc)
    {
      sym.mapping = linfo_dquant;
      sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
      sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
      currSlice->LFAlphaC0Offset = sym.value1 << 1;
      buf->frame_bitoffset += sym.len;
      bitptr+=sym.len;

      sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
      sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
      currSlice->LFBetaOffset = sym.value1 << 1;
      buf->frame_bitoffset += sym.len;
      bitptr+=sym.len;
    }
  }

  sym.mapping = linfo;        // Mapping rule: Simple code number to len/info
  // Multi-Picture Buffering Syntax, Feb 27, 2002
  // RPSF: Reference Picture Selection Flags
  sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &(sym.inf), buf->bitstream_length );
  sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
  buf->frame_bitoffset += sym.len;
  bitptr+=sym.len;

  // PN: Picture Number
  sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &(sym.inf), buf->bitstream_length );
  sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
  buf->frame_bitoffset += sym.len;
  bitptr+=sym.len;
  img->pn=sym.value1;

  // Tian Dong: payloadType 8 is IDERP packet, in which PN is reset,
  // so a PN gap is expected. Be sure not try to fill such a gap.
  if ( pp->payloadType != PAYLOAD_TYPE_IDERP ) fill_PN_gap(img);
  // Tian Dong: The following condition will ensure the buffer reset will
  // not be done before the successive B frames are decoded.
  if ( pp->payloadType == PAYLOAD_TYPE_IDERP && fb != NULL && frm != NULL )
  {
    reset_buffers();
    frm->picbuf_short[0]->used=0;
    frm->picbuf_short[0]->picID=-1;
    frm->picbuf_short[0]->lt_picID=-1;
    frm->short_used=0;
  }

  // RPSL: Reference picture selection layer
  sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &(sym.inf), buf->bitstream_length );
  sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
  buf->frame_bitoffset += sym.len;
  bitptr+=sym.len;

  if (sym.value1)
  {
    // read Reference Picture Selection Layer
    // free old buffer content
    while (img->currentSlice->rmpni_buffer)
    { 
      tmp_rmpni=img->currentSlice->rmpni_buffer;
 
      img->currentSlice->rmpni_buffer=tmp_rmpni->Next;
      free (tmp_rmpni);
    } 
    done=0;
    // if P or B frame RMPNI

    if ((img->type==INTER_IMG)||(img->type==B_IMG)||(img->type==SP_IMG))
    {
      do
      {
    
        sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &(sym.inf), buf->bitstream_length );
        sym.mapping(sym.len, sym.inf, &tmp1, &(sym.value2));
        buf->frame_bitoffset += sym.len;
        bitptr+=sym.len;


        // check for illegal values
        if ((tmp1<0)||(tmp1>3))
          error ("Invalid RMPNI operation specified",400);

        if (tmp1!=3)
        {
          tmp_rmpni=(RMPNIbuffer_t*)calloc (1,sizeof (RMPNIbuffer_t));
          tmp_rmpni->Next=NULL;
          tmp_rmpni->RMPNI=tmp1;

          // get the additional parameter
          sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &(sym.inf), buf->bitstream_length );
          sym.mapping(sym.len, sym.inf, &(tmp_rmpni->Data), &(sym.value2));
          buf->frame_bitoffset += sym.len;
          bitptr+=sym.len;

          // add RMPNI to list
          if (img->currentSlice->rmpni_buffer==NULL) 
          {
            img->currentSlice->rmpni_buffer=tmp_rmpni;
          }
          else
          {
            tmp_rmpni2=img->currentSlice->rmpni_buffer;
            while (tmp_rmpni2->Next!=NULL) 
              tmp_rmpni2=tmp_rmpni2->Next;
            tmp_rmpni2->Next=tmp_rmpni;
          }
        } else
        {
          // free temporary memory (no need to save end loop operation)
          done=1;
        }
      } while (!done);
    }
  }

  // RBPT 
  sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &(sym.inf), buf->bitstream_length );
  sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
  buf->frame_bitoffset += sym.len;
  bitptr+=sym.len;

  // free old buffer content
  while (img->mmco_buffer)
  { 
    tmp_mmco=img->mmco_buffer;

    img->mmco_buffer=tmp_mmco->Next;
    free (tmp_mmco);
  } 
  
  // read Memory Management Control Operation
  if (sym.value1)
  {
    // read Memory Management Control Operation 
    do
    {

      tmp_mmco=(MMCObuffer_t*)calloc (1,sizeof (MMCObuffer_t));
      tmp_mmco->Next=NULL;
    
      sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &(sym.inf), buf->bitstream_length );
      sym.mapping(sym.len, sym.inf, &(tmp_mmco->MMCO), &(sym.value2));
      buf->frame_bitoffset += sym.len;
      bitptr+=sym.len;

      switch (tmp_mmco->MMCO)
      {
      case 0:
      case 5:
        break;
      case 1:
        sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &(sym.inf), buf->bitstream_length );
        sym.mapping(sym.len, sym.inf, &(tmp_mmco->DPN), &(sym.value2));
        buf->frame_bitoffset += sym.len;
        bitptr+=sym.len;
        break;
      case 2:
        sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &(sym.inf), buf->bitstream_length );
        sym.mapping(sym.len, sym.inf, &(tmp_mmco->LPIN), &(sym.value2));
        buf->frame_bitoffset += sym.len;
        bitptr+=sym.len;
        break;
      case 3:
        sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &(sym.inf), buf->bitstream_length );
        sym.mapping(sym.len, sym.inf, &(tmp_mmco->DPN), &(sym.value2));
        buf->frame_bitoffset += sym.len;
        bitptr+=sym.len;
        sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &(sym.inf), buf->bitstream_length );
        sym.mapping(sym.len, sym.inf, &(tmp_mmco->LPIN), &(sym.value2));
        buf->frame_bitoffset += sym.len;
        bitptr+=sym.len;
        break;
      case 4:
        sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &(sym.inf), buf->bitstream_length );
        sym.mapping(sym.len, sym.inf, &(tmp_mmco->MLIP1), &(sym.value2));
        buf->frame_bitoffset += sym.len;
        bitptr+=sym.len;
        break;
      default:
        error ("Invalid MMCO operation specified",400);
        break;
      }

      // add MMCO to list
      if (img->mmco_buffer==NULL) 
      {
        img->mmco_buffer=tmp_mmco;
      }
      else
      {
        tmp_mmco2=img->mmco_buffer;
        while (tmp_mmco2->Next!=NULL) tmp_mmco2=tmp_mmco2->Next;
        tmp_mmco2->Next=tmp_mmco;
      }
      
    }while (tmp_mmco->MMCO!=0);
  }
}

/*!
 ************************************************************************
 * \brief
 *      read the payload data
 * \return
 *      0, if success
 *      -1, otherwise
 * \param img
 *      image pointer
 * \param inp
 *      input parameter pointer
 * \param pp
 *      corresponding payloadinfo
 * \param fp
 *      input file pointer
 ************************************************************************
 */
int rdOnePayload( struct img_par *img, struct inp_par* inp, PayloadInfo *pp, FILE* fp )
{
  int *read_len;        
  Slice* currSlice = img->currentSlice;
  Bitstream* currStream = currSlice->partArr[0].bitstream;
  byte* buf = currStream->streamBuffer;
  DecodingEnvironmentPtr dep = &((currSlice->partArr[0]).de_cabac);

  currStream->frame_bitoffset = 0;
  currStream->code_len = 0;
  fseek( fp, box_atm.currPayloadOffset, SEEK_SET );
  memset( buf, 0, MAX_CODED_FRAME_SIZE );
  if ( pp->payloadSize == fread( buf, 1, (unsigned int)pp->payloadSize, fp ) )
  {
    currStream->bitstream_length = (int)pp->payloadSize;
    box_atm.currPayloadOffset += (long)pp->payloadSize;

    currStream->bitstream_length = (int) pp->payloadSize = EBSPtoRBSP( buf, (int) pp->payloadSize, 0);
    currStream->bitstream_length = (int) pp->payloadSize = RBSPtoSODB( buf, (int) pp->payloadSize);

    currSlice->ei_flag = 0;
    currSlice->dp_mode = PAR_DP_1;
    currSlice->max_part_nr=1;
    read_len = &(currSlice->partArr[0].bitstream->read_len);
    *read_len = (int)pp->payloadSize;

    if(active_pps->entropy_coding_mode  == CABAC)
    {
      dep = &((currSlice->partArr[0]).de_cabac);
      arideco_start_decoding(dep, buf, 0, read_len, img->type);
    }

    // At this point the slice is ready for decoding. 
    currSlice->next_header = IFFGetFollowingSliceHeader(img, pp); // no use for the info in nextp, nextsh yet. 
    return 0;
  }
  IFFEndOfFile = TRUE;
  return -1;
}

/*!
 ************************************************************************
 * \brief
 *      To check if the last payload in the current picture is read
 * \return
 *      SOP, if TRUE
 *      SOS, if FALSE, means there are more slices or payloads
 * \param img
 *      image pointer
 * \param pp
 *      corresponding payloadinfo
 ************************************************************************
 */
int IFFGetFollowingSliceHeader( struct img_par *img, PayloadInfo* pp )
{
  if ( currPictureInfo.numPayloads == pp->payloadnr ) // the last payloadinfo is read, set the position of next pic
    return SOP;
  else // next payloadinfo 
    return SOS;
}

/*!
 ************************************************************************
 * \brief
 *      Do nothing
 * \return
 *      None
 ************************************************************************
 */
void freeAlternateTrackMediaBox()
{
}

/*!
 ************************************************************************
 * \brief
 *      use parameter set n to set the decoder parameters
 * \return
 *      0, if success
 *      -1, otherwise
 * \param n
 *      indicate which parameter set will be used
 * \param img
 *      image pointer
 * \param inp
 *      input parameter pointer
 ************************************************************************
 */
int IFFUseParameterSet( int n, struct img_par* img, struct inp_par* inp )
{
  if ( n == CurrentParameterSet )
    return 0;

  CurrentParameterSet = n;
  img->width = box_ps.pictureWidthInMBs * MB_BLOCK_SIZE;
  img->width_cr = box_ps.pictureWidthInMBs * MB_BLOCK_SIZE / 2;

  img->height = box_ps.pictureHeightInMBs * MB_BLOCK_SIZE;
  img->height_cr = box_ps.pictureHeightInMBs * MB_BLOCK_SIZE / 2;

  if ( box_ps.entropyCoding == 0 ) active_pps->entropy_coding_mode  = UVLC;
  else active_pps->entropy_coding_mode  = CABAC;

//!  inp->partition_mode = box_ps.partitioningType;
//!  inp->UseConstrainedIntraPred = box_ps.intraPredictionType;
  inp->buf_cycle = box_ps.bufCycle;
  inp->LFParametersFlag = box_ps.loopFilterParametersFlag;
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      read one slice from the input file
 * \return
 *      EOS, if met with the end of sequence
 *      SOP, if this is the start of a picture
 *      SOS, if this is the start of a slice
 * \param img
 *      image pointer
 * \param inp
 *      input parameter pointer
 ************************************************************************
 */
int readSliceIFF( struct img_par* img, struct inp_par* inp )
{
  extern FILE *bits;

  if ( IFFEndOfFile ) 
    return EOS;

  if ( -1 == rdPayloadInfo(img, inp, &currPayloadInfo, bits) ) return EOS;
  if ( -1 == rdOnePayload(img, inp, &currPayloadInfo, bits ) ) return EOS;

  if ( BeginOfPictureOrSlice == SOP )
  {
    BeginOfPictureOrSlice = SOS;
    return SOP;
  }

  return SOS;
}

/*!
 ************************************************************************
 * \brief
 *      read the data from file, bytes are in big Endian order.
 * \return
 *      how many bytes being read, if success
 *      -1, otherwise
 * \param outf
 *      output file pointer
 ************************************************************************
 */
int readfile( void* buf, size_t size, size_t count, FILE* fp ) 
{
  byte* p;
  int num = 0;
    
  if (isBigEndian)
    p = (byte*)buf;
  else
    p = (byte*)buf+size-1;

  assert( fp != NULL );
  assert( buf != NULL );
  assert( count == 1 );

  while ( size > 0 )
  {
    if (isBigEndian)
    {
      if ( fread( p++, 1, 1, fp ) != 1 ) return -1;
    }
    else
    {
      if ( fread( p--, 1, 1, fp ) != 1 ) return -1;
    }
    size--;
    num++;
  }
  return num;
}

/*!
 ************************************************************************
 * \brief
 *      read the data from file, bytes are in big Endian order.
 *      to be used if buffers size differs from number of written bytes
 * \return
 *      how many bytes being read, if success
 *      -1, otherwise
 * \param outf
 *      output file pointer
 ************************************************************************
 */
int readfile_s( void* buf, size_t bufsize, size_t size, size_t count, FILE* fp ) 
{
  byte* p;
  int num = 0;
    
  if (isBigEndian)
    p = (byte*)buf+(bufsize-size);
  else
    p = (byte*)buf+size-1;

  assert( fp != NULL );
  assert( buf != NULL );
  assert( count == 1 );

  while ( size > 0 )
  {
    if (isBigEndian)
    {
      if ( fread( p++, 1, 1, fp ) != 1 ) return -1;
    }
    else
    {
      if ( fread( p--, 1, 1, fp ) != 1 ) return -1;
    }
    size--;
    num++;
  }
  return num;
}

/*!
 ************************************************************************
 * \brief
 *      read one box, and then parse it.
 * \return
 *      0, if success
 *      -1, otherwise
 * \param img
 *      image pointer
 * \param inp
 *      input parameter pointer
 * \param snr
 *      snr pointer
 * \param bits
 *      input file pointer
 ************************************************************************
 */
int parse_one_box( struct img_par* img, struct inp_par* inp, struct snr_par* snr, FILE* fp )
{
  unsigned int type, size, storedpos;
  assert( fp );

  storedpos = ftell( fp );

  size = -1;
  if ( -1 == readfile( &size, 4, 1, fp ) ) return -1;
  if ( size == 0 ) return -1;
  if ( -1 == readfile( &type, 4, 1, fp ) ) return -1;

  switch ( type )
  {
  case BOX_JVTH: //<! 
    if ( -1 == rdFileHeaderBox( fp, size ) ) return -1;
    break;
  case BOX_CINF: //<! 
    if ( -1 == rdContentInfoBox( fp, size ) ) return -1;
    break;
  case BOX_ATIN: //<! 
    if ( -1 == rdAlternateTrackInfoBox( fp, size ) ) return -1;
    break;
  case BOX_PRMS: //<! 
    if ( -1 == rdParameterSetBox( fp, size ) ) return -1;
    break;
  case BOX_SEGM: //<! 
    if ( -1 == parse_one_segment( img, inp, snr, fp, size ) ) return -1;
    break;
  case BOX_ATRH: //<! 
  case BOX_SWPC: //<! 
  case BOX_ATRM:  //<! 
  default:
    printf( "Unexpected boxs %d at %ld\n", type, ftell( fp ));
    return -1;
    break;
  }

  // point to the next box
  if ( 0 != fseek( fp, storedpos+size, SEEK_SET ) ) return -1;
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      Read the first Boxes, and set the initial parameter set according to the input file
 * \return
 *      0, if success
 *      -1, otherwise
 * \param img
 *      image pointer
 * \param inp
 *      input parameter pointer
 * \param bits
 *      input file pointer
 ************************************************************************
 */
int IFFSequenceHeader( struct img_par *img, struct inp_par *inp, FILE *bits )
{
  unsigned long size, type;
  long storedpos, pos;

  assert( bits != NULL );

  memset(&currPictureInfo, 0, sizeof(PictureInfo));
  memset(&oldPictureInfo, 0, sizeof(PictureInfo));

  // read FileTypeBox
  if ( -1 == testFileTypeBox( bits ) ) return -1;
  size = -1;
  if ( -1 == readfile( &size, 4, 1, bits ) ) return -1;
  if ( size == 0 ) return -1;
  if ( -1 == readfile( &type, 4, 1, bits ) ) return -1;
  if ( type != BOX_JVTH ) return -1;
  if ( -1 == rdFileHeaderBox( bits, size ) ) return -1;

  // now, try to find the parameter set box, which parametersetID is equal to CurrentParameterSet
  storedpos = ftell( bits );
  if ( 0 != findParameterSetBox( bits, &pos, 0 ) ) return -1;
  fseek( bits, pos, SEEK_SET );
  if ( -1 == rdParameterSetBox( bits, 0 ) ) return -1;
  IFFUseParameterSet( 0, img, inp );    // use parameter set 0

  fseek( bits, storedpos, SEEK_SET );
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      finished, free all resources allocated for interim file format
 * \return
 *      None
 * \param img
 *      image pointer
 * \param inp
 *      input parameter pointer
 * \param bits
 *      input file pointer
 ************************************************************************
 */
void terminateInterimFile()
{ 
  freeSegmentBox();
  freeParameterSetBox();
  freeAlternateTrackInfoBox();
  freeContentInfoBox();
  freeFileHeaderBox();
  freeFileTypeBox();
}

int rdLayerBox( int no, unsigned long boxsize, FILE *bits )
{
  unsigned long size, type, readsize;
  unsigned INT32 avgBitRate, avgFrameRate;
  time_t ltime;
  FILE *fp = fopen("enhanced_layers.txt", "at");
  if ( fp == NULL )
  {
    printf("Cannot open enhanced_layers.txt for reporting the information in Layer Box.\n");
    return -1;
  }

  if ( no == 0 )
  {
    time( &ltime );
    fprintf( fp, "\n\n*****************************************\n");
    fprintf( fp, "Reporting the information in Layer Boxes\nGenerated by the JVT decoder\nDate & Time: %s", ctime( &ltime ));
    fprintf( fp, "*****************************************\n");
  }
  fprintf(fp, "\nLayer Box %d:\n", no );

  if (-1 == readfile( &avgBitRate, 4, 1, bits ) ) return -1;
  if (-1 == readfile( &avgFrameRate, 4, 1, bits ) ) return -1;

  fprintf(fp, "Average bit-rate is %ld, average frame rate is %ld.\n", avgBitRate, avgFrameRate );
  fprintf(fp, "The sub-sequences in this layer:\n");

  readsize = 0;
  boxsize = boxsize - SIZEOF_BOXTYPE - 8;
  while (readsize < boxsize)
  {
    if ( -1 == readfile( &size, 4, 1, bits ) ) return -1;
    if ( size == 0 ) return -1;
    if ( -1 == readfile( &type, 4, 1, bits ) ) return -1;
    if ( type != BOX_SSEQ ) return -1;
    if ( -1 == rdSubSequence(fp, bits) ) return -1;
    readsize += size;
  }
  fclose( fp );
  return 0;
}

int rdSubSequence(FILE* fp, FILE* bits)
{
  int i;
  unsigned INT16 subSequenceIdentifier;
  Boolean continuationFromPreviousSegmentFlag=FALSE, continuationToNextSegmentFlag=FALSE;
  Boolean startTickAvailableFlag=FALSE;
  unsigned INT64 ssStartTick, ssDuration;
  unsigned INT32 avgBitRate, avgFrameRate;
  unsigned INT16 numReferencedSubSequences;
  unsigned INT8  layerNumber;
  unsigned INT16 depSubSequenceIdentifier;
  unsigned INT16 cd;

  if (-1 == readfile( &subSequenceIdentifier, 2, 1, bits ) ) return -1;
  fprintf(fp, "\nsubSequenceIdentifier = %d\n", subSequenceIdentifier );
  if (-1 == readfile( &cd, 2, 1, bits ) ) return -1;
  if ( cd&0x8000 ) continuationFromPreviousSegmentFlag = TRUE;
  if ( cd&0x4000 ) continuationToNextSegmentFlag = TRUE;
  if ( cd&0x2000 ) startTickAvailableFlag = TRUE;
  fprintf(fp, "continuationFromPreviousSegmentFlag = %d\n", continuationFromPreviousSegmentFlag );
  fprintf(fp, "continuationToNextSegmentFlag = %d\n", continuationToNextSegmentFlag);
  fprintf(fp, "startTickAvailableFlag = %d\n", startTickAvailableFlag);
  if (-1 == readfile( &ssStartTick, 8, 1, bits ) ) return -1;
  fprintf(fp, "ssStartTick = %lld\n", ssStartTick);
  if (-1 == readfile( &ssDuration, 8, 1, bits ) ) return -1;
  fprintf(fp, "ssDuration = %lld\n", ssDuration);
  if (-1 == readfile( &avgBitRate, 4, 1, bits ) ) return -1;
  fprintf(fp, "avgBitRate = %ld\n", avgBitRate);
  if (-1 == readfile( &avgFrameRate, 4, 1, bits ) ) return -1;
  fprintf(fp, "avgFrameRate = %ld\n", avgFrameRate);
  if (-1 == readfile( &numReferencedSubSequences, 2, 1, bits ) ) return -1;
  fprintf(fp, "numReferencedSubSequences = %d\n", numReferencedSubSequences);

  for (i=0; i<numReferencedSubSequences; i++)
  {
    if (-1 == readfile( &layerNumber, 1, 1, bits ) ) return -1;
    fprintf(fp, "layerNumber = %d\n", layerNumber);
    if (-1 == readfile( &depSubSequenceIdentifier, 2, 1, bits ) ) return -1;
    fprintf(fp, "depSubSequenceIdentifier= %d\n", depSubSequenceIdentifier);
  }
  return 0;
}
