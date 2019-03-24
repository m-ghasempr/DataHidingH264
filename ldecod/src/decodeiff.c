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
#include "global.h"
#include "decodeiff.h"

FileTypeBox box_ft;
FileHeaderBox box_fh;
ContentInfoBox box_ci;
AlternateTrackInfoBox box_ati;
ParameterSetBox box_ps;
SegmentBox box_s;
AlternateTrackHeaderBox box_ath;
PayloadInfo currPayloadInfo;
PictureInfo currPictureInfo;
AlternateTrackMediaBox box_atm;
SwitchPictureBox box_sp;

extern int CurrentParameterSet = -1;  // Tian Dong: must not be set to 0 here.
static int BeginOfPictureOrSlice; // SOP or SOS
int IFFEndOfFile = FALSE;
int ThisAlternateTrack = 0;   // indicate which track in the file will be decoded.

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
  if ( -1 == readfile( &box_ps.entropyCoding, 1, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ps.motionResolution, 1, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ps.partitioningType, 1, 1, fp ) ) return -1;
  if ( -1 == readfile( &box_ps.intraPredictionType, 1, 1, fp ) ) return -1;

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
  unsigned long numPictures;
  int ret = EOS + 100;
//  FILE *f;

//  f = fopen("dummy2.txt", "at");

  if ( 0 != fseek( fp, box_ath.storedpos, SEEK_SET ) ) return -1;
  if ( -1 == readfile( &numPictures, box_fh.numBytesInPictureCountMinusOne+1, 1, fp ) ) return -1;
  currPictureInfo.storedpos = ftell( fp );

  while ( numPictures-- > 0 && ret != EOS )
    ret = decode_one_frame( img, inp, snr );
/*  {
    rdPictureInfo( fp );
    fprintf( f, "picture payload, total %d:\n", currPictureInfo.numPayloads );
    i=9;
    while (i-- > 0)
    {
      ret = readSliceIFF( img, inp );
      fprintf( f, "picture %d payload nr %d 's payloadsize: %d\n", currPayloadInfo.pictureID, currPayloadInfo.payloadnr, currPayloadInfo.payloadSize);
    }
  }
  fclose( f ); */
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
    if ( -1 == readfile( &size, 4, 1, fp ) ) return -1;     // in 32 bits mode
    if ( -1 == readfile( &type, 4, 1, fp ) ) return -1;
    pos = (unsigned)ftell(fp);
    if ( type == BOX_ATRH ) nr++;
    if ( tracknr > nr )
    { if ( 0 != fseek( fp, size - SIZEOF_BOXTYPE, SEEK_CUR ) ) return -1; }
  }

  if ( tracknr != nr ) return -1;     // have not found the track needed

  *storedpos = ftell( fp );           // save the file pointer to the track header
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

  fseek( fp, currPictureInfo.storedpos, SEEK_SET );

  num += readfile( &cd, 1, 1, fp );
  if ( cd == 0x80 ) currPictureInfo.intraPictureFlag = 1;
  else currPictureInfo.intraPictureFlag = 0;

  num += readfile( &currPictureInfo.pictureOffset, box_fh.numBytesInPictureOffsetMinusTwo + 2, 1, fp );
  num += readfile( &currPictureInfo.pictureDisplayTime, box_fh.numBytesInPictureDisplayTimeMinusOne+1, 1, fp );
  num += readfile( &currPictureInfo.numPayloads, box_fh.numBytesInPayloadCountMinusOne+1, 1, fp );

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
  Slice *currSlice = img->currentSlice;
  byte cd;

  assert( fp != NULL );

  fseek( fp, pp->storedpos, SEEK_SET );

  if ( -1 == readfile( &pp->payloadSize, box_fh.numBytesInPayloadSizeMinusOne+1, 1, fp ) ) return -1;
  if ( -1 == readfile( &pp->headerSize, 1, 1, fp ) ) return -1;
  if ( -1 == readfile( &cd, 1, 1, fp ) ) return -1;
  pp->payloadType = (cd >> 4);
  pp->errorIndication = ((cd&0x0F) >> 3);
  pp->reserved = (cd & 0x07);

  switch ( pp->payloadType )
  {
  case 0:
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
    currPictureInfo.storedpos = ftell( fp );
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
  int dP_nr = 0;
  Slice *currSlice = img->currentSlice;
  Bitstream *currStream = (currSlice->partArr[dP_nr]).bitstream;
  DataPartition *partition = &(currSlice->partArr[dP_nr]);
  int bitptr = 0;

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
    img->type = currSlice->picture_type = (inp->buf_cycle > 1)?INTER_IMG_MULT:INTER_IMG_1;
    break;
  case 1:
    img->type = currSlice->picture_type = (inp->buf_cycle > 1)?B_IMG_MULT:B_IMG_1;
    break;
  case 2:
    img->type = currSlice->picture_type = (inp->buf_cycle > 1)?SP_IMG_MULT:SP_IMG_1;
    break;
  case 3:
    img->type = currSlice->picture_type = INTRA_IMG;
    break;
  default:
    printf ("Panic: unknown Slice type %d, conceal by loosing slice\n", pp->sliceType);
    currSlice->ei_flag = 1;
  }  
//  currSlice->picture_type = img->type = sym.value1;
  buf->frame_bitoffset += sym.len;
  bitptr += sym.len;

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

  sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
  sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
  pp->initialQP = 31 - sym.value1;
  currSlice->qp = img->qp = pp->initialQP;
  buf->frame_bitoffset += sym.len;
  bitptr += sym.len;

  if ( inp->symbol_mode == CABAC )
  {
    sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
    sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
    pp->lastMBnr = sym.value1;
    currSlice->last_mb_nr = sym.value1;
    buf->frame_bitoffset += sym.len;
    bitptr += sym.len;
  }
  if ( img->type == SP_IMG_1 || img->type == SP_IMG_MULT )
  {
    sym.len = GetVLCSymbol( buf->streamBuffer, buf->frame_bitoffset, &sym.inf, buf->bitstream_length );
    sym.mapping(sym.len, sym.inf, &(sym.value1), &(sym.value2));
    img->qpsp = 31 - sym.value1;
    buf->frame_bitoffset += sym.len;
    bitptr += sym.len;
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

    currSlice->ei_flag = 0;
    currSlice->dp_mode = PAR_DP_1;
    currSlice->max_part_nr=1;
    read_len = &(currSlice->partArr[0].bitstream->read_len);
    *read_len = (int)pp->payloadSize;

    if(inp->symbol_mode == CABAC)
    {
      dep = &((currSlice->partArr[0]).de_cabac);
      arideco_start_decoding(dep, buf, 0, read_len);
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
  Slice* currSlice = img->currentSlice;
  Bitstream* currStream = currSlice->partArr[0].bitstream;
  byte* buf = currStream->streamBuffer;
  
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
  Slice* currSlice = img->currentSlice;

  if ( n == CurrentParameterSet )
    return 0;

  CurrentParameterSet = n;
  img->width = box_ps.pictureWidthInMBs * MB_BLOCK_SIZE;
  img->width_cr = box_ps.pictureWidthInMBs * MB_BLOCK_SIZE / 2;

  img->height = box_ps.pictureHeightInMBs * MB_BLOCK_SIZE;
  img->height_cr = box_ps.pictureHeightInMBs * MB_BLOCK_SIZE / 2;

  if ( box_ps.entropyCoding == 0 ) inp->symbol_mode = UVLC;
  else inp->symbol_mode = CABAC;

  switch ( box_ps.motionResolution )
  {
  case 2:
    img->mv_res = 0;
    break;
  case 3:
    img->mv_res = 1;
    break;
  case 0:
    img->mv_res = 2;
    break;
  case 1:
    img->mv_res = 3;
    break;
  default:
    assert( 1==0 );
  }

  inp->partition_mode = box_ps.partitioningType;

  inp->UseConstrainedIntraPred = box_ps.intraPredictionType;

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
  byte* p = (byte*)buf+size-1;
  int num = 0;

  assert( fp != NULL );
  assert( buf != NULL );
  assert( count == 1 );

  while ( size > 0 )
  {
    if ( fread( p--, 1, 1, fp ) != 1 ) return -1;
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
    printf( "Unexpected boxs %d at %d\n", type, ftell( fp ));
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