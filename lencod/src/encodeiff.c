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
 *     encodeiff.c
 *  \brief
 *     Implementation of the interim file format, as defined
 *     in Appendix III of WD-1 (JVT-A003r1).
 *
 *     Instructions for a new release of the Joint Model:
 *     1. Remember to update WORKING_DRAFT_MAJOR_NO and WORKING_DRAFT_MINOR_NO
 *        in defines.h appropriately
 *     2. Remember to modify the decoder as well
 *
 *     Instructions for modifications to the interim file format:
 *     1. Update INTERIM_FILE_MAJOR_NO and INTERIM_FILE_MINOR_NO
 *        in defines.h appropriately. It is up to the implementor
 *        to decide, which modification deserves an update in the major
 *        version number.
 *     2. If you modify one of the box types that are of constant size,
 *        such as the Parameter Set Box, update the hard-coded size of 
 *        the box appropriately.
 *
 *     Feb 10, 2002:
 *     Until now, the box header only permits the compact size (32 bits mode),
 *     and the Data Partition mode is not supported. One segment per file
 *     is created.
 *      
 *  \author
 *      - Dong Tian                             <tian@cs.tut.fi>
 *      - Miska M. Hannuksela                   <miska.hannuksela@nokia.com>
 *
 ************************************************************************
 */

 /*!
 ************************************************************************
 *  Tian Dong:
 *  \flowchart of how the output module works:
 *  main()
 *  {
 *    start_sequence()
 *    {
 *      SequenceHeader(out)
 *      {
 *        initInterimFile()
 *      }
 *    }
 *    initSegmentBox()
 *    {
 *      ......
 *      initAlternateTrackHeaderBox();
 *      initSwitchPictureBox();
 *      initAlternateTrackMediaBox();
 *    }
 *    loop over pictures in the same segment
 *    {
 *      encode_one_frame()
 *      {
 *        init_frame()
 *        update header of AlternateTrackHeaderBox 
 *        initPictureInfo 
 *        while (end_of_frame == FALSE) // loop over slices
 *        {
 *          encode_one_slice(&sym)
 *          {
 *            start_slice()
 *            {
 *              newPayloadInfo
 *              addOnePayloadInfo
 *            }
 *            loop over MB encoding // write the MBData to mem buffer
 *            terminate_slice()  file handle
 *            {
 *              updating the payloadInfo
 *              writing media data to the tmp file
 *            }
 *          }
 *        }
 *        wrPictureInfo // write pictureInfo & all payloadInfos to a temporary file
 *        {
 *          write the pictureInfo header to file
 *          loop to write every payloadInfo to file
 *          {
 *            wrPayloadInfo
 *            update info about SDU in AlternateTrackInfo Box
 *          }
 *        }
 *        freePictureInfo
 *      }
 *    }
 *    updateAlternateTrackHeaderBox();
 *    updateAlternateTrackMediaBox();
 *    updateSegmentBox();
 *    terminate_sequence()
 *    {
 *      TerminateInterimFile( outf )
 *      {
 *        wrFileTypeBox( outf );
 *        wrFileHeaderBox( outf );
 *        wrContentInfoBox( outf );
 *        wrAlternateTrackInfoBox( outf );
 *        wrParameterSetBox( outf );
 *        wrSegmentBox( outf );
 *        freeAll
 *      }
 *    }
 *  }
 ************************************************************************
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <malloc.h>
#include "global.h"
#include "mbuffer.h"
#include "encodeiff.h"

FileTypeBox box_ft;
FileHeaderBox box_fh;
ContentInfoBox box_ci;
AlternateTrackInfoBox box_ati;
ParameterSetBox box_ps;
SegmentBox box_s;
AlternateTrackHeaderBox box_ath;
PayloadInfo* pCurrPayloadInfo;
PictureInfo currPictureInfo;
AlternateTrackMediaBox box_atm;
SwitchPictureBox box_sp;

/*!
 ************************************************************************
 * \brief
 *      Initiate the File Type Box. 
 * \return
 *      0, if success
 *      -1, otherwise
 ************************************************************************
 */
int initFileTypeBox()
{
  box_ft.type.size = SIZEOF_BOXTYPE + 12;
  box_ft.type.type = BOX_FTYP;

  memcpy( box_ft.majorBrand, "jvt ", 4 );
  box_ft.jmMajorVersion = WORKING_DRAFT_MAJOR_NO;
  box_ft.jmMinorVersion = WORKING_DRAFT_MINOR_NO;
  
  box_ft.numCompatibleBrands = 1;
  box_ft.compatibleBrands = calloc( box_ft.numCompatibleBrands*4, 1 );
  if ( box_ft.compatibleBrands == NULL ) return -1;
  memcpy( box_ft.compatibleBrands, "jvt ", 4 );

  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      The function to write the File Type Box to output file
 * \return
 *      how many bytes been written, if success
 *      -1, if failed
 * \param fp
 *      output file pointer
 ************************************************************************
 */
size_t wrFileTypeBox(FILE* fp)
{
  size_t num = 0, ret;
  assert( fp != NULL );

  // to write the structure to file
  num += writefile( &box_ft.type.size, 4, 1, fp );
  num += writefile( &box_ft.type.type, 4, 1, fp );
  num += writefile( box_ft.majorBrand, 4, 1, fp );
  num += writefile( &box_ft.jmMajorVersion, 2, 1, fp );
  num += writefile( &box_ft.jmMinorVersion, 2, 1, fp );
  ret = fwrite( box_ft.compatibleBrands, 4, box_ft.numCompatibleBrands, fp );
  num += (ret * 4);

  if ( num == box_ft.type.size ) return num;
  return -1;
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
  free( box_ft.compatibleBrands );
}

/*!
 ************************************************************************
 * \brief
 *      Initiate the File Header Box
 * \return
 *      0, if success
 *      -1, otherwise
 ************************************************************************
 */
int initFileHeaderBox( )
{
  // to write the file header
  box_fh.type.size = SIZEOF_BOXTYPE + 27;
  box_fh.type.type = BOX_JVTH;

  box_fh.majorVersion = INTERIM_FILE_MAJOR_NO;
  box_fh.minorVersion = INTERIM_FILE_MINOR_NO;
  box_fh.timescale = 30000;
  box_fh.numUnitsInTick = 1001;
  box_fh.duration = 0; // to be updated at the end of coding
  box_fh.pixAspectRatioX = 1;
  box_fh.pixAspectRatioY = 1;
  box_fh.maxPicId = input->PicIdModulus - 1;
  box_fh.numAlternateTracks = 1;
  box_fh.numBytesInPayloadCountMinusOne = 1;
  box_fh.numBytesInPictureOffsetMinusTwo = 1;
  box_fh.numBytesInPictureDisplayTimeMinusOne = 1;
  box_fh.numBytesInPictureCountMinusOne = 1;
  box_fh.numBytesInPayloadSizeMinusOne = 1;
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      The function to write the File Header Box to output file
 * \return
 *      how many bytes been written, if success
 *      -1, if failed
 * \param fp
 *      output file pointer
 ************************************************************************
 */
size_t wrFileHeaderBox( FILE* fp )
{
  byte cd;
  size_t num = 0;
  assert( fp != NULL );

  box_fh.duration = box_ati.info[0].last_frame + 1;     // update

  num += writefile( &box_fh.type.size, 4, 1, fp );
  num += writefile( &box_fh.type.type, 4, 1, fp );
  num += writefile( &box_fh.majorVersion, 1, 1, fp );
  num += writefile( &box_fh.minorVersion, 1, 1, fp );
  num += writefile( &box_fh.timescale, 4, 1, fp );
  num += writefile( &box_fh.numUnitsInTick, 4, 1, fp );
  num += writefile( &box_fh.duration, 8, 1, fp );
  num += writefile( &box_fh.pixAspectRatioX, 2, 1, fp );
  num += writefile( &box_fh.pixAspectRatioY, 2, 1, fp );
  num += writefile( &box_fh.maxPicId, 2, 1, fp );
  num += writefile( &box_fh.numAlternateTracks, 1, 1, fp );

  cd = 0;
  cd = (box_fh.numBytesInPayloadCountMinusOne << 6) |
    (box_fh.numBytesInPictureOffsetMinusTwo << 4) |
    (box_fh.numBytesInPictureDisplayTimeMinusOne << 2) |
    (box_fh.numBytesInPictureCountMinusOne << 0);
  num += writefile( &cd, 1, 1, fp );

  cd = 0 +
    (box_fh.numBytesInPayloadSizeMinusOne << 6);
  num += writefile( &cd, 1, 1, fp );

  if ( num == box_fh.type.size ) return num;
  return -1;
}

/*!
 ************************************************************************
 * \brief
 *      Free the memory allocated for File Type Box
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
 *      Initiate the Content Info Box
 *      Tian Dong, Feb 10, 2002:
 *      Do nothing, The Box is skipped in current implementation.
 * \return
 *      0, if success
 *      -1, otherwise
 ************************************************************************
 */
int initContentInfoBox()
{
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      The function to write the Content Info Box to output file
 * \return
 *      how many bytes been written, if success
 *      -1, if failed
 * \param fp
 *      output file pointer
 ************************************************************************
 */
size_t wrContentInfoBox( FILE* fp )
{
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
void freeContentInfoBox()
{
}

// Functions on AlternateTrackInfoBox
/*!
 ************************************************************************
 * \brief
 *      Initiate the Alternate Track Info Box
 *      Tian Dong, Feb 10, 2002:
 *      Only one track media contained in the output file
 * \return
 *      0, if success
 *      -1, otherwise
 ************************************************************************
 */
int initAlternateTrackInfoBox()
{
  assert(box_fh.numAlternateTracks == 1);   // Only one track media contained in the output file

  box_ati.type.size = SIZEOF_BOXTYPE + 12 * box_fh.numAlternateTracks;
  box_ati.type.type = BOX_ATIN;

  box_ati.info = calloc( sizeof(AlternateTrackInfo), box_fh.numAlternateTracks );
  if ( box_ati.info == NULL ) return -1;

  box_ati.info[0].displayWindowWidth = input->img_width;
  box_ati.info[0].displayWindowHeight = input->img_height;
  box_ati.info[0].maxSDUSize = 0; // to be updated
  box_ati.info[0].avgSDUSize = 0; // to be updated
  box_ati.info[0].avgBitRate = 0; // to be updated

  box_ati.info[0].numSDU = 0;
  box_ati.info[0].sumSDUSize = 0;
  box_ati.info[0].last_frame = 0;
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      The function to write the Alternate Track Info Box to output file
 * \return
 *      how many bytes been written, if success
 *      -1, if failed
 * \param fp
 *      output file pointer
 ************************************************************************
 */
size_t wrAlternateTrackInfoBox( FILE* fp )
{
  float t;
  size_t num = 0;

  assert( fp != NULL );

  // update avgSDUSize
  if ( box_ati.info[0].numSDU != 0 )
    box_ati.info[0].avgSDUSize = (unsigned INT16)(box_ati.info[0].sumSDUSize / box_ati.info[0].numSDU);
  t = (float)(box_ati.info[0].last_frame + 1) * (float)box_fh.numUnitsInTick / (float)box_fh.timescale;
  box_ati.info[0].avgBitRate = (INT32) (box_ati.info[0].sumSDUSize / t );

  // write them to the file
  num += writefile( &box_ati.type.size, 4, 1, fp );
  num += writefile( &box_ati.type.type, 4, 1, fp );

  num += writefile( &box_ati.info[0].displayWindowWidth, 2, 1, fp );
  num += writefile( &box_ati.info[0].displayWindowHeight, 2, 1, fp );
  num += writefile( &box_ati.info[0].maxSDUSize, 2, 1, fp );
  num += writefile( &box_ati.info[0].avgSDUSize, 2, 1, fp );
  num += writefile( &box_ati.info[0].avgBitRate, 4, 1, fp );

  if ( num == box_ati.type.size ) return num;
  return -1;
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
  free(box_ati.info);
}

/*!
 ************************************************************************
 * \brief
 *      Initiate the Parameter Set Box
 *      Tian Dong, Feb 10, 2002:
 *      Only one parameter set, whose ID is set to 0, contained in the output file
 * \return
 *      0, if success
 *      -1, otherwise
 ************************************************************************
 */
int initParameterSetBox()
{
  box_ps.type.size = SIZEOF_BOXTYPE + 27;     // 26 => 27, add bufCycle
  box_ps.type.type = BOX_PRMS;

  box_ps.parameterSetID = 0;
  box_ps.profile = 0;
  box_ps.level = 0;
  box_ps.version = 0;
  box_ps.pictureWidthInMBs = input->img_width / MB_BLOCK_SIZE;  
  box_ps.pictureHeightInMBs = input->img_height / MB_BLOCK_SIZE;
  box_ps.displayRectangleOffsetTop = 0;
  box_ps.displayRectangleOffsetLeft = 0;
  box_ps.displayRectangleOffsetBottom = 0;
  box_ps.displayRectangleOffsetRight = 0;
  box_ps.displayMode = 1;
  box_ps.displayRectangleOffsetFromWindowTop = 0;
  box_ps.displayRectangleOffsetFromWindowLeftBorder = 0;
  if ( input->symbol_mode == UVLC ) box_ps.entropyCoding = 0;
  else box_ps.entropyCoding = 1;
  switch ( input->mv_res )  // the value defined in VCEG-O58 is different from that in .cfg
  {
  case 2:
    box_ps.motionResolution = 0;  // full-pixel
    break;
  case 3:
    box_ps.motionResolution = 1;  // half-pixel
    break;
  case 0:
    box_ps.motionResolution = 2;  // 1/4-pixel
    break;
  case 1:
    box_ps.motionResolution = 3;  // 1/8-pixel
    break;
  default:
    break;
  }
  box_ps.partitioningType = input->partition_mode;
  box_ps.intraPredictionType = input->UseConstrainedIntraPred;
  box_ps.bufCycle = input->no_multpred;
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      The function to write the Parameter Set Box to output file
 * \return
 *      how many bytes been written, if success
 *      -1, if failed
 * \param fp
 *      output file pointer
 ************************************************************************
 */
size_t wrParameterSetBox( FILE* fp )
{
  size_t num = 0;

  assert( fp != NULL );
  num += writefile( &box_ps.type.size, 4, 1, fp );
  num += writefile( &box_ps.type.type, 4, 1, fp );

  num += writefile( &box_ps.parameterSetID, 2, 1, fp );
  num += writefile( &box_ps.profile, 1, 1, fp );
  num += writefile( &box_ps.level, 1, 1, fp );
  num += writefile( &box_ps.version, 1, 1, fp );
  num += writefile( &box_ps.pictureWidthInMBs, 2, 1, fp );
  num += writefile( &box_ps.pictureHeightInMBs, 2, 1, fp );
  num += writefile( &box_ps.displayRectangleOffsetTop, 2, 1, fp );
  num += writefile( &box_ps.displayRectangleOffsetLeft, 2, 1, fp );
  num += writefile( &box_ps.displayRectangleOffsetBottom, 2, 1, fp );
  num += writefile( &box_ps.displayRectangleOffsetRight, 2, 1, fp );
  num += writefile( &box_ps.displayMode, 1, 1, fp );
  num += writefile( &box_ps.displayRectangleOffsetFromWindowTop, 2, 1, fp );
  num += writefile( &box_ps.displayRectangleOffsetFromWindowLeftBorder, 2, 1, fp );
  num += writefile( &box_ps.entropyCoding, 1, 1, fp );
  num += writefile( &box_ps.motionResolution, 1, 1, fp );
  num += writefile( &box_ps.partitioningType, 1, 1, fp );
  num += writefile( &box_ps.intraPredictionType, 1, 1, fp );
  num += writefile( &box_ps.bufCycle, 1, 1, fp );

  if ( num == box_ps.type.size ) return num;
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

// Functions on SegmentBox
/*!
 ************************************************************************
 * \brief
 *      Initiate the Segment Box and all sub-Boxes in it.
 *      Tian Dong, Feb 10, 2002:
 *      Only one Segment Box in the output file
 * \return
 *      0, if success
 *      -1, otherwise
 ************************************************************************
 */
int initSegmentBox()
{
  box_s.type.size = 0;  // to be updated
  box_s.type.type = BOX_SEGM;
  
  box_s.fileSize = 0; // to be updated
  box_s.startTick = 0;
  box_s.firstFrameNr = box_s.lastFrameNr;
  box_s.lastFrameNr = 0;  // serve as the duration of the segment.
  box_s.segmentDuration = 0; // to be updated

  // since we now only deal with the ONE track case, we assert:
  assert( box_fh.numAlternateTracks == 1 );

  if ( -1 == initAlternateTrackHeaderBox() ) return -1;
  if ( -1 == initSwitchPictureBox() ) return -1;
  if ( -1 == initAlternateTrackMediaBox() ) return -1;

  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      Update the data in Segment Box
 * \return
 *      None
 ************************************************************************
 */
void updateSegmentBox()
{
  box_s.type.size = SIZEOF_BOXTYPE + 24 + box_ath.type.size + box_atm.type.size;
  box_s.fileSize = box_s.type.size;
  box_s.segmentDuration = box_s.lastFrameNr - box_s.firstFrameNr + 1;
}

/*!
 ************************************************************************
 * \brief
 *      The function to write the Segment Box to output file
 * \return
 *      how many bytes been written, if success
 *      -1, if failed
 * \param fp
 *      output file pointer
 ************************************************************************
 */
size_t wrSegmentBox( FILE *fp )
{
  size_t num = 0;
  assert( fp );

  // since we now only deal with the ONE track case, we assert:
  assert( box_fh.numAlternateTracks == 1 );

  num += writefile( &box_s.type.size, 4, 1, fp );
  num += writefile( &box_s.type.type, 4, 1, fp );
  num += writefile( &box_s.fileSize, 8, 1, fp );
  num += writefile( &box_s.startTick, 8, 1, fp );
  num += writefile( &box_s.segmentDuration, 8, 1, fp );

  num += mergeAlternateTrackHeaderBox( fp );
  num += mergeAlternateTrackMediaBox( fp );

  if ( num == box_s.type.size ) return num;
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      Free all resource allocated for this segment
 * \return
 *      None
 ************************************************************************
 */
void freeSegmentBox()
{
  freeAlternateTrackHeaderBox();
  freeAlternateTrackMediaBox();
}

// Functions on AlternateTrackHeaderBox
/*!
 ************************************************************************
 * \brief
 *      Initiate the Alternate Track Header Box & some data in picture info.
 *      Tian Dong, Feb 10, 2002:
 *      Only one Alternate Track Header Box in the output file.
 * \return
 *      0, if success
 *      -1, otherwise
 ************************************************************************
 */
int initAlternateTrackHeaderBox()
{
  box_ath.type.size = 0;  // to be updated
  box_ath.type.type = BOX_ATRH;
  box_ath.numPictures = 0;  // set to 0

  box_ath.fpMeta = tmpfile();
  if ( box_ath.fpMeta == NULL ) return -1;

  currPictureInfo.lastFrameNr = 0;    // do this for init of picture info
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      Update the data in this Alternate Track Header
 * \return
 *      0, if success
 *      -1, if failed
 ************************************************************************
 */
int updateAlternateTrackHeaderBox()
{
  int pictureDataSize;
  
  assert( box_ath.fpMeta != NULL );
  // update the head data
  fseek( box_ath.fpMeta, 0, SEEK_END );
  pictureDataSize = ftell( box_ath.fpMeta );
  box_ath.type.size = SIZEOF_BOXTYPE + box_fh.numBytesInPictureCountMinusOne+1 + pictureDataSize;
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      Merge the Alternate Track Header Data to the output file
 * \return
 *      how many bytes been appended, if success
 *      -1, if failed
 * \param fp
 *      output file pointer
 ************************************************************************
 */
size_t mergeAlternateTrackHeaderBox( FILE* fp )
{
  FILE* sourcef;
  FILE* destf;
  unsigned char c;
  size_t num = 0;

  sourcef = box_ath.fpMeta;
  destf = fp;

  assert( sourcef != NULL );
  assert( destf != NULL );

  // write the head of Alternate Track Header Box
  num += writefile( &box_ath.type.size, 4, 1, fp );
  num += writefile( &box_ath.type.type, 4, 1, fp );
//  writefile( &box_ath.type.largesize, 8, 1, fp );

  num += writefile( &box_ath.numPictures, box_fh.numBytesInPictureCountMinusOne+1, 1, destf );
  
  // append the data in box_ath.fpMeta into fp:
  fseek( sourcef, 0L, SEEK_SET );
  
  c = fgetc(sourcef);
  while ( !feof( sourcef ) )
  {
    fputc( c, destf );
    num++;
    c = fgetc(sourcef);
  }

  if ( num == box_ath.type.size ) return num;
  return -1;
}

/*!
 ************************************************************************
 * \brief
 *      Free all resource allocated for this Alternate Track Header Box
 * \return
 *      None
 ************************************************************************
 */
void freeAlternateTrackHeaderBox()
{
  fclose( box_ath.fpMeta );
}

// Functions on PictureInfo
/*!
 ************************************************************************
 * \brief
 *      Initiate the picture info, before one frame is encoded
 * \return
 *      0, if success
 *      -1, otherwise
 ************************************************************************
 */
int initPictureInfo()
{
  static int prev_frame_no = 0;

  if ( img->type == INTRA_IMG )
    currPictureInfo.intraPictureFlag = 0x80; // if the picture is INTRA
  else
    currPictureInfo.intraPictureFlag = 0x00; // if the picture is not INTRA
  currPictureInfo.pictureOffset = currPictureInfo.currPictureSize;
  currPictureInfo.currPictureSize = 0;  // reset

  // to set the picture display time (relative to the previous coded frame)
  currPictureInfo.pictureDisplayTime = frame_no - prev_frame_no;
  prev_frame_no = frame_no;

  if ( img->type != B_IMG )
    currPictureInfo.lastFrameNr = frame_no;

  currPictureInfo.numPayloads = 0;  // numPayloads must be set to zero here
  currPictureInfo.payloadData = NULL;
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      The function to write the meta data of the current picture
 * \return
 *      how many bytes been written, if success
 *      -1, if failed
 * \param fp
 *      output file pointer
 ************************************************************************
 */
size_t wrPictureInfo( FILE* fp )
{
  int i;
  INT64 size;
  PayloadInfo* pp;
  size_t num = 0;
//  FILE* f;

//  f = fopen( "dummy1.txt", "at" );
  assert( fp != NULL );

  num += writefile( &currPictureInfo.intraPictureFlag, 1, 1, fp );
  num += writefile( &currPictureInfo.pictureOffset, box_fh.numBytesInPictureOffsetMinusTwo + 2, 1, fp );
  num += writefile( &currPictureInfo.pictureDisplayTime, box_fh.numBytesInPictureDisplayTimeMinusOne + 1, 1, fp );
  num += writefile( &currPictureInfo.numPayloads, box_fh.numBytesInPayloadCountMinusOne + 1, 1, fp );

  if ( num != (unsigned)(1+box_fh.numBytesInPictureOffsetMinusTwo + 2+box_fh.numBytesInPictureDisplayTimeMinusOne + 1+box_fh.numBytesInPayloadCountMinusOne + 1) ) return -1;
//  fprintf( f, "picture %d 's payload, total %d:\n", frame_no, currPictureInfo.numPayloads );

  pp = currPictureInfo.payloadData;
  for ( i = 0; i < currPictureInfo.numPayloads; i++ )
  { 
    assert( pp != NULL );

//    fprintf( f, "payloadsize of %d is: %d\n", i, pp->payloadSize );

    // update and write the payload to file
    if ( -1 == wrPayloadInfo( pp, fp ) ) return -1;

    // then update the parameter in Alternate Track Info Box
    size = pp->payloadSize + pp->headerSize;
    box_ati.info[0].numSDU++;
    box_ati.info[0].sumSDUSize += (long double)size;
    if ( box_ati.info[0].maxSDUSize < size ) box_ati.info[0].maxSDUSize = (unsigned INT16)size;

    pp = pp->next;
  }

//  fclose ( f );
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      Free all resource allocated for this Picture Info
 * \return
 *      None
 ************************************************************************
 */
void freePictureInfo()
{
  int i;
  PayloadInfo* next;

  for ( i = 0; i < currPictureInfo.numPayloads; i++ )
  {
    assert( currPictureInfo.payloadData != NULL );
    next = currPictureInfo.payloadData->next;
    free( currPictureInfo.payloadData );
    currPictureInfo.payloadData = next;
  }
}

// Functions on payloadInfo

/*!
 ************************************************************************
 * \brief
 *      To new a payload info to save the meta data of a slice
 *      Tian Dong, Feb 10, 2002:
 *      No data partition supported.
 * \return
 *      a pointer to the allocated mem, if success
 *      NULL, if failed
 ************************************************************************
 */
PayloadInfo* newPayloadInfo()
{
  PayloadInfo* pli;

  pli = calloc(1, sizeof(PayloadInfo));
  if ( pli == NULL ) return NULL;

  pli->payloadSize = 0;
  pli->headerSize = box_fh.numBytesInPayloadSizeMinusOne+1 + 1 + 1;
  pli->payloadType = 0;  // single slice
  pli->errorIndication = 0;  // no known error
  pli->reserved = 0;
  pli->parameterSet = 0;
  pli->pictureID = 0;
  pli->sliceID = 0;
  pli->sliceType = 0;
  pli->firstMBInSliceX = 0;
  pli->firstMBInSliceY = 0;
  pli->next = NULL;
  pli->pn = img->pn;
  pli->type = img->type;
  pli->max_lindex = img->max_lindex;
  pli->lindex = img->lindex;

  // set values
  if ( input->partition_mode == PAR_DP_1 )
    pli->payloadType = 0;
  else
    assert(0==1);
  pli->parameterSet = box_ps.parameterSetID;
  pli->pictureID = img->currentSlice->picture_id;
  pli->sliceID = img->current_slice_nr;
  
  // copied from rtp.c, and add pli->sliceType2
  switch (img->type)
  {
    case INTER_IMG:
      if (img->types==SP_IMG)
      {
        pli->sliceType = 2;
        pli->sliceType2 = SP_IMG;
      }
      else
      {
        pli->sliceType = 0;
        pli->sliceType2 = INTER_IMG;
      }
      break;
    case INTRA_IMG:
      pli->sliceType = 3;
      pli->sliceType2 = INTRA_IMG;
      break;
    case B_IMG:
      pli->sliceType = 1;
      pli->sliceType2 = B_IMG;
      break;
    case SP_IMG:
      pli->sliceType = 2;
      pli->sliceType2 = SP_IMG;
      break;
    default:
      printf ("Panic: unknown picture type %d, exit\n", img->type);
      assert (0==1);
  }

  pli->firstMBInSliceX = img->mb_x;
  pli->firstMBInSliceY = img->mb_y;
  pli->initialQP = img->qp;
  pli->qpsp = img->qpsp;

  // begin of ERPS
#ifdef _CHECK_MULTI_BUFFER_1_

  /* RPSL: Reference Picture Selection Layer */
  if(img->type!=INTRA_IMG)
  {
    // let's mix some reference frames
    if ((img->pn==5)&&(img->type==INTER_IMG))
    {
      RMPNIbuffer_t *r;
      r = (RMPNIbuffer_t*)calloc (1,sizeof(RMPNIbuffer_t));
      r->RMPNI=0;
      r->Data=2;
      r->Next=NULL;
      img->currentSlice->rmpni_buffer=r;
    }
  }
  reorder_mref(img);
  
#endif

#ifdef _CHECK_MULTI_BUFFER_2_

  // some code to check operation
  if ((img->pn==3) && (img->type==INTER_IMG))
  {
    // check in this frame as long term picture
//    if (img->max_lindex==0)
    {
      // set long term buffer size = 2
      img->max_lindex=2;
    }

    // assign local long term
    init_long_term_buffer(2,img);
    init_mref(img);
    init_Refbuf(img);

    if (img->current_slice_nr==0)
    {
      assign_long_term_id(3,img->lindex,img);
      img->lindex=(img->lindex+1)%img->max_lindex;
    }
  }

#endif 
  // end of ERPS

  return pli;
}

/*!
 ************************************************************************
 * \brief
 *      Add a payload info to the Payloadinfo List, which head is stored in PictureInfo
 *      Tian Dong, Feb 10, 2002:
 *      No data partition supported.
 * \return
 *      1, if success
 *      0, if failed
 * \param pi
 *      The PayloadInfo will be added to the list
 ************************************************************************
 */
int addOnePayloadInfo(PayloadInfo* pi)
{
  PayloadInfo* p;
  PayloadInfo* last;

  last = p = currPictureInfo.payloadData;

  assert( pi != NULL );

  while ( p )
  {
    last = p;
    p = p->next;
  }

  if (last == NULL)
    currPictureInfo.payloadData = pi;  // this is the first payloadInfo
  else
    last->next = pi;  // add the payloadInfo to the end of the payloadInfo list

  currPictureInfo.numPayloads++;
  
  return 1;
}

/*!
 ************************************************************************
 * \brief
 *      The function to write one payloadinfo to file
 * \return
 *      how many bytes been written, if success
 *      -1, if failed
 * \param pp
 *      PayloadInfo pointer
 * \param fp
 *      output file pointer
 ************************************************************************
 */
size_t wrPayloadInfo( PayloadInfo* pp, FILE *fp )
{
  byte cd; 
  Bitstream* bitstream; // used as a MEM buffer of slice header
  SyntaxElement sym;
  size_t num = 0, bytes_written;

  assert( pp != NULL );
  assert( fp != NULL );

  // Initialize the bitsteam:
  bitstream = alloca(sizeof(Bitstream));
  assert( bitstream != NULL );
  bitstream->streamBuffer = alloca(BUFSIZE_FOR_PAYLOADINFO);
  memset( bitstream->streamBuffer, 0, BUFSIZE_FOR_PAYLOADINFO);
  
  bitstream->bits_to_go  = 8;
  bitstream->byte_pos    = 0;
  bitstream->byte_buf    = 0;
  
  // First write the element to the MEM bitstream buffer
  sym.type = SE_HEADER;       // This will be true for all symbols generated here
  sym.mapping = n_linfo2;       // Mapping rule: Simple code number to len/info
  if ( pp->payloadType == 0 )
  {
    // write the parameter set
    sym.value1 = pp->parameterSet;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);
    // write slice header;
    sym.value1 = pp->pictureID;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);
    sym.value1 = pp->sliceType;
//    select_picture_type (&sym);
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);
    sym.value1 = pp->firstMBInSliceX;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);
    sym.value1 = pp->firstMBInSliceY;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);
    sym.value1 = MAX_QP - pp->initialQP;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);
    if ( input->symbol_mode == CABAC )
    {
      sym.value1 = pp->lastMBnr;
      writeSyntaxElement2Buf_UVLC(&sym, bitstream);
    }
    if ( pp->sliceType2 == SP_IMG )
    {
      sym.value1 = MAX_QP - pp->qpsp;
      writeSyntaxElement2Buf_UVLC(&sym, bitstream);
    }

    iff_writeERPS(&sym, pp, bitstream);       // Tian: to support ERPS (Annex U), Feb 27, 2002
  }
  else if ( pp->payloadType == 1 )
  {
    // write the parameter set
    sym.value1 = pp->parameterSet;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);
    // write slice header;
    sym.value1 = pp->pictureID;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);
    sym.value1 = pp->sliceType;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);
    sym.value1 = pp->firstMBInSliceX;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);
    sym.value1 = pp->firstMBInSliceY;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);
    sym.value1 = MAX_QP - pp->initialQP;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);
    // write the slice id
    sym.value1 = pp->sliceID;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);

    iff_writeERPS(&sym, pp, bitstream);       // Tian: to support ERPS (Annex U), Feb 27, 2002
  }
  else if ( pp->payloadType == 2 || pp->payloadType == 3 )
  {
    // write pictureID
    sym.value1 = pp->pictureID;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);
    // write sliceID
    sym.value1 = pp->sliceID;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);

    iff_writeERPS(&sym, pp, bitstream);       // Tian: to support ERPS (Annex U), Feb 27, 2002
  }
  else if ( pp->payloadType == 5 )
  {
    // no additional codewords
  }

  // finishing the MEM buffer
  if (bitstream->bits_to_go < 8)
  { // trailing bits to process
    bitstream->byte_buf <<= bitstream->bits_to_go;
    bitstream->streamBuffer[bitstream->byte_pos++]=bitstream->byte_buf;
    bitstream->bits_to_go = 8;
  }

  bytes_written = bitstream->byte_pos;

  // update & write headerSize
  num += writefile( &pp->payloadSize, box_fh.numBytesInPayloadSizeMinusOne+1, 1, fp );
  pp->headerSize += bytes_written;
  num += writefile( &pp->headerSize, 1, 1, fp );
  cd = (pp->payloadType << 4) | (pp->errorIndication << 3) | (pp->reserved);
  num += writefile( &cd, 1, 1, fp );
  if ( num != (unsigned)(box_fh.numBytesInPayloadSizeMinusOne+1+2) ) return -1;

  // Then write the bitstream to FILE
  if ( bytes_written != fwrite (bitstream->streamBuffer, 1, bytes_written, fp) ) return -1;
  return num+bytes_written;
}

/*!
 ************************************************************************
 * \brief
 *      writes the ERPS syntax elements to a bitstream
 *      imitate write_ERPS(), and change the function calls:
 *    from:
 *      len += writeSyntaxElement_UVLC (sym, partition);
 *    to:
 *      len += writeSyntaxElement2Buf_UVLC(sym, bitstream);
 * \return
 *      how many bytes been written, if success
 *      -1, if failed
 * \param bitstream
 *      destination: where the code to be written
 ************************************************************************
 */
size_t iff_writeERPS(SyntaxElement *sym, PayloadInfo* pp, Bitstream* bitstream)
{
  size_t len=0;

  // RPSF: Reference Picture Selection Flags
  sym->value1 = 0;
  len += writeSyntaxElement2Buf_UVLC(sym, bitstream);

  // PN: Picture Number
  sym->value1 = pp->pn;
  len += writeSyntaxElement2Buf_UVLC(sym, bitstream);

#ifdef _CHECK_MULTI_BUFFER_1_

  // RPSL: Reference Picture Selection Layer
  sym->value1 = 1;
  len += writeSyntaxElement2Buf_UVLC(sym, bitstream);

  if(pp->type!=INTRA_IMG)
  {
    // let's mix some reference frames
    if ((pp->pn==5)&&(pp->type==INTER_IMG))
    {
      // negative ADPN follows
      // RMPNI
      sym->value1 = 0;
      len += writeSyntaxElement2Buf_UVLC(sym, bitstream);
      // ADPN
      sym->value1 = 2;
      len += writeSyntaxElement2Buf_UVLC(sym, bitstream);
    }
    // RMPNI
    sym->value1 = 3;
    len += writeSyntaxElement2Buf_UVLC(sym, bitstream);
  }
  
#else
  // RPSL: Reference Picture Selection Layer
  sym->value1 = 0;
  len += writeSyntaxElement2Buf_UVLC(sym, bitstream);

#endif

#ifdef _CHECK_MULTI_BUFFER_2_

  // Reference Picture Bufering Type
  sym->value1 = 1;
  len += writeSyntaxElement2Buf_UVLC(sym, bitstream);

  // some code to check operation
  if ((pp->pn==3) && (pp->type==INTER_IMG))
  {
    // set long term buffer size = 2
    // MMCO Specify Max Long Term Index
    // command
    sym->value1 = 4;
    len += writeSyntaxElement2Buf_UVLC(sym, bitstream);
    // size = 2+1 (MLP1)
    sym->value1 = 2+1;
    len += writeSyntaxElement2Buf_UVLC(sym, bitstream);

    // assign a long term index to actual frame
    // MMCO Assign Long Term Index to a Picture
    // command
    sym->value1 = 3;
    len += writeSyntaxElement2Buf_UVLC(sym, bitstream);
    // DPN=0 for actual frame 
    sym->value1 = 0;
    len += writeSyntaxElement2Buf_UVLC(sym, bitstream);
    //long term ID
    sym->value1 = pp->lindex;
    len += writeSyntaxElement2Buf_UVLC(sym, bitstream);
  } 
  if ((pp->pn==4) && (pp->type==INTER_IMG))
  {
     if (pp->max_lindex>0)
     {
      // delete long term picture again
      // MMCO Mark a Long-Term Picture as Unused
      // command
      sym->value1 = 2;
      len += writeSyntaxElement2Buf_UVLC(sym, bitstream);
      // MMCO LPIN
      // command
      sym->value1 = (pp->max_lindex+pp->lindex-1)%pp->max_lindex;
      len += writeSyntaxElement2Buf_UVLC(sym, bitstream);
    }
  } 

  // end MMCO loop
  // end loop
  sym->value1 = 0;
  len += writeSyntaxElement2Buf_UVLC(sym, bitstream);
#else
    // RPBT: Reference Picture Bufering Type
    sym->value1 = 0;
    len += writeSyntaxElement2Buf_UVLC(sym, bitstream);
#endif 

  return len;
}

// Functions on SwitchPictureBox
/*!
 ************************************************************************
 * \brief
 *      Initiate the Switch Picture Box. Do nothing
 *      Tian Dong, Feb 10, 2002
 *      The switch picture box is skipped in the output file
 * \return
 *      0, if success
 *      -1, otherwise
 ************************************************************************
 */
int initSwitchPictureBox()
{
  return 0;
}

// Functions on AlternateMediaBox
/*!
 ************************************************************************
 * \brief
 *      Initiate the Alternate Track Media Box & some data in picture info.
 *      Tian Dong, Feb 10, 2002:
 *      Only one Alternate Track Media Box in the output file.
 * \return
 *      0, if success
 *      -1, otherwise
 ************************************************************************
 */
int initAlternateTrackMediaBox()
{
  box_atm.type.size = 0;    // to be updated
  box_atm.type.type = BOX_ATRM;

  box_atm.fpMedia = tmpfile();
  assert( box_atm.fpMedia != NULL );

  // to maintain a picture pointer to point to the beginning 
  // of the latest picture, relative to the beginning of the ATMC.
  currPictureInfo.currPictureSize = SIZEOF_BOXTYPE;  // the first time to set its value. in 64 bits mode
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      Update the data in this Alternate Track Media
 * \return
 *      0, if success
 *      -1, if failed
 ************************************************************************
 */
int updateAlternateTrackMediaBox()
{
  int mediaDataSize;

  assert( box_atm.fpMedia != NULL );
  fseek( box_atm.fpMedia, 0, SEEK_END );
  mediaDataSize = ftell( box_atm.fpMedia );
  box_atm.type.size = SIZEOF_BOXTYPE + mediaDataSize;
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      Merge the Alternate Track Media Data to the output file
 * \return
 *      how many bytes been appended, if success
 *      -1, if failed
 * \param fp
 *      output file pointer
 ************************************************************************
 */
int mergeAlternateTrackMediaBox( FILE* fp )
{
  FILE* sourcef;
  FILE* destf;
  unsigned char c;
  size_t num = 0;

  sourcef = box_atm.fpMedia;
  destf = fp;

  assert( sourcef != NULL );
  assert( destf != NULL );
  
  // write the head
  num += writefile( &box_atm.type.size, 4, 1, fp );
  num += writefile( &box_atm.type.type, 4, 1, fp );
  if ( num != 8 ) return -1;

  // append the data in box_ath.fpMeta into fp:
  fseek( sourcef, 0L, SEEK_SET );
  c = fgetc(sourcef);
  while ( !feof( sourcef ) )
  {
    fputc( c, destf );
    num++;
    c = fgetc(sourcef);
  }
  return num;
}

/*!
 ************************************************************************
 * \brief
 *      Free all resource allocated for this Alternate Track Media Box
 * \return
 *      None
 ************************************************************************
 */
void freeAlternateTrackMediaBox()
{
  fclose( box_atm.fpMedia );
}

/*!
 ************************************************************************
 * \brief
 *      The init function for Interim File Format
 * \return
 *      0, if success
 *      -1, otherwise
 ************************************************************************
 */
int initInterimFile()
{
  box_s.lastFrameNr = 0;
  if ( -1 == initFileTypeBox() ) return -1;
  if ( -1 == initFileHeaderBox() ) return -1;
  if ( -1 == initContentInfoBox() ) return -1;
  if ( -1 == initAlternateTrackInfoBox() ) return -1;
  if ( -1 == initParameterSetBox() ) return -1;
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      The close function for Interim File Format
 * \return
 *      how many bytes being written into the file, if success
 *      -1, otherwise
 * \param outf
 *      output file pointer
 ************************************************************************
 */
size_t terminateInterimFile(FILE* outf)
{
  size_t num = 0, len;

  assert( outf != NULL );

  if ( (len = wrFileTypeBox( outf )) == -1 ) return -1;
  num += len;
  if ( (len = wrFileHeaderBox( outf )) == -1 ) return -1;
  num += len;
  if ( (len = wrContentInfoBox( outf )) == -1 ) return -1;
  num += len;
  if ( (len = wrAlternateTrackInfoBox( outf )) == -1 ) return -1;
  num += len;
  if ( (len = wrParameterSetBox( outf )) == -1 ) return -1;
  num += len;
  if ( (len = wrSegmentBox( outf )) == -1 ) return -1;
  num += len;

  freeSegmentBox();
  freeParameterSetBox();
  freeAlternateTrackInfoBox();
  freeContentInfoBox();
  freeFileHeaderBox();
  freeFileTypeBox();

  return num;
}

/*!
 ************************************************************************
 * \brief
 *      write the data to file, bytes are in big Endian order.
 * \return
 *      how many bytes being written into the file, if success
 *      -1, otherwise
 * \param outf
 *      output file pointer
 ************************************************************************
 */
size_t writefile( void* buf, size_t size, size_t count, FILE* fp )
{
  byte* p = (byte*)buf+size-1;
  int num = 0;

  assert( fp != NULL );
  assert( buf != NULL );
  assert( count == 1 );

  while ( size > 0 )
  {
    fwrite( p--, 1, 1, fp );
    size--;
    num++;
  }
  return num;
}
