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
 *  Tian Dong (Last Updated in June 15, 2002)
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
 *      initAlternateTrackHeaderBox()
 *      {
 *         initPictureInformationBox()
 *         initLayerBox()
 *      }
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
PictureInformationBox box_pi;
PayloadInfo* pCurrPayloadInfo;
PictureInfo currPictureInfo;
LayerBox box_layr[MAX_LAYER_NUMBER];
SubSequenceBox box_sseq[MAX_LAYER_NUMBER];
AlternateTrackMediaBox box_atm;
SwitchPictureBox box_sp;

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
  box_ps.type.size = SIZEOF_BOXTYPE + 28;     // 26 => 27 => 28, add bufCycle, temporal scalability
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
  if (input->NumFramesInELSubSeq!=0) box_ps.requiredPictureNumberUpdateBehavior=1;
  else box_ps.requiredPictureNumberUpdateBehavior=0;
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
  unsigned char cd;

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
  if (box_ps.requiredPictureNumberUpdateBehavior==1) cd=0x80;
  else cd=0;
  num += writefile( &cd, 1, 1, fp );

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
 *      Tian Dong, May 30, 2002:
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
  
  if ( box_ps.requiredPictureNumberUpdateBehavior != 0 ) box_ath.numLayers = 2;
  else box_ath.numLayers = 0;

  assert((box_ath.numLayers <= MAX_LAYER_NUMBER));

  if ( initPictureInformationBox() == -1 ) return -1;
  if ( initLayerBox() == -1 ) return -1;
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      Update the data in this Alternate Track Header Box
 *      Only one PicureInformationBox in ATH, May 30, 2002
 * \return
 *      0, if success
 *      -1, if failed
 ************************************************************************
 */
int updateAlternateTrackHeaderBox()
{
  int i;
  assert(box_ath.numLayers <= MAX_LAYER_NUMBER );

  updatePictureInformationBox();
  if ( input->NumFramesInELSubSeq != 0 )
    updateLayerBox();

  // update the head data
  box_ath.type.size = SIZEOF_BOXTYPE + 1 + box_pi.type.size;
  for (i=0; i<box_ath.numLayers; i++)
    box_ath.type.size += box_layr[i].type.size;
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
  FILE* destf;
  size_t num = 0, ret;

  destf = fp;
  assert( destf != NULL );

  // write the head of Picture Information
  num += writefile( &box_ath.type.size, 4, 1, destf );
  num += writefile( &box_ath.type.type, 4, 1, destf );
//  writefile( &box_ath.type.largesize, 8, 1, destf );

  num += writefile( &box_ath.numLayers, 1, 1, destf );
  
  ret = mergePictureInformationBox( destf );
  if ( ret == -1 ) return -1;
  num += ret;

  if ( input->NumFramesInELSubSeq != 0 ) 
  {
    ret = mergeLayerBox( destf );
    if ( ret == -1 ) return -1;
    num += ret;
  }

  if ( num != box_ath.type.size ) return -1;
  return num;
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
  freePictureInformationBox();
  freeLayerBox();
}

/*!
 ************************************************************************
 * \brief
 *      Initiate the Picture Information Box & some data in picture info.
 *      Tian Dong, Feb 10, 2002:
 *      Only one Alternate Track Header Box in the output file.
 * \return
 *      0, if success
 *      -1, otherwise
 ************************************************************************
 */
int initPictureInformationBox()
{
  box_pi.type.size = 0;  // to be updated
  box_pi.type.type = BOX_PICI;
  box_pi.numPictures = 0;  // set to 0

  box_pi.fpMeta = tmpfile();
  if ( box_pi.fpMeta == NULL ) return -1;

  currPictureInfo.lastFrameNr = 0;    // do this for init of picture info
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      Update the data in this Picture Information Box
 * \return
 *      0, if success
 *      -1, if failed
 ************************************************************************
 */
int updatePictureInformationBox()
{
  int pictureDataSize;
  
  assert( box_pi.fpMeta != NULL );
  // update the head data
  fseek( box_pi.fpMeta, 0, SEEK_END );
  pictureDataSize = ftell( box_pi.fpMeta );
  box_pi.type.size = SIZEOF_BOXTYPE + box_fh.numBytesInPictureCountMinusOne+1 + pictureDataSize;
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      Merge the PictureInformation Data to the output file
 * \return
 *      how many bytes been appended, if success
 *      -1, if failed
 * \param fp
 *      output file pointer
 ************************************************************************
 */
size_t mergePictureInformationBox( FILE* fp )
{
  FILE* sourcef;
  FILE* destf;
  unsigned char c;
  size_t num = 0;

  sourcef = box_pi.fpMeta;
  destf = fp;

  assert( sourcef != NULL );
  assert( destf != NULL );

  // write the head of Picture Information
  num += writefile( &box_pi.type.size, 4, 1, fp );
  num += writefile( &box_pi.type.type, 4, 1, fp );
//  writefile( &box_pi.type.largesize, 8, 1, fp );

  num += writefile_s( &box_pi.numPictures, sizeof(box_pi.numPictures), box_fh.numBytesInPictureCountMinusOne+1, 1, destf );
  
  // append the data in box_ath.fpMeta into fp:
  fseek( sourcef, 0L, SEEK_SET );
  
  c = fgetc(sourcef);
  while ( !feof( sourcef ) )
  {
    fputc( c, destf );
    num++;
    c = fgetc(sourcef);
  }

  if ( num == box_pi.type.size ) return num;
  return -1;
}

/*!
 ************************************************************************
 * \brief
 *      Free all resource allocated for this PictureInformation Box
 * \return
 *      None
 ************************************************************************
 */
void freePictureInformationBox()
{
  fclose( box_pi.fpMeta );
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

  if ( img->type == INTRA_IMG ) currPictureInfo.intraPictureFlag = TRUE;
  else currPictureInfo.intraPictureFlag = FALSE;
  if ( img->type == SP_IMG ) currPictureInfo.syncPictureFlag = TRUE; // ???
  else currPictureInfo.syncPictureFlag = FALSE; // ???

  currPictureInfo.pictureOffset = currPictureInfo.currPictureSize;
  currPictureInfo.currPictureSize = 0;  // reset

  // to set the picture display time (relative to the previous coded frame)
  currPictureInfo.pictureDisplayTime = frame_no - prev_frame_no;
  prev_frame_no = frame_no;

  if ( img->type != B_IMG )
    currPictureInfo.lastFrameNr = frame_no;

  currPictureInfo.numPayloads = 0;  // numPayloads must be set to zero here
  currPictureInfo.payloadData = NULL;

  // the follow values may be updated.
  currPictureInfo.layerNumber = 0;
  currPictureInfo.subSequenceIdentifier = 0;
  currPictureInfo.originLayerNumber = 0;
  currPictureInfo.originSubSequenceIdentifier = 0;

  if ( input->NumFramesInELSubSeq != 0 )
  {
    if (IMG_NUMBER%(input->NumFramesInELSubSeq+1)==0 && img->type!=B_IMG)
      currPictureInfo.layerNumber = 0;
    else
    {
      currPictureInfo.layerNumber = 1;
    }
    currPictureInfo.subSequenceIdentifier = box_sseq[currPictureInfo.layerNumber].subSequenceIdentifier;
    if (img->type!=B_IMG)
    {
      currPictureInfo.refFromLayerNumber = currPictureInfo.layerNumber;
      currPictureInfo.refFromSubSequenceIdentifier = currPictureInfo.subSequenceIdentifier;
    }
  }

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
  unsigned char cd;
  size_t num = 0;

  assert( fp != NULL );

  if ( currPictureInfo.intraPictureFlag == TRUE ) cd = 0x80;
  else cd = 0x00;
  if ( currPictureInfo.syncPictureFlag == TRUE ) cd |= 0x40;
  num += writefile( &cd, 1, 1, fp );

  num += writefile_s( &currPictureInfo.pictureOffset, sizeof(currPictureInfo.pictureOffset), box_fh.numBytesInPictureOffsetMinusTwo + 2, 1, fp );
  num += writefile_s( &currPictureInfo.pictureDisplayTime, sizeof(currPictureInfo.pictureDisplayTime), box_fh.numBytesInPictureDisplayTimeMinusOne + 1, 1, fp );

  if ( box_ath.numLayers )
  {
    if ( -1 == writefile( &currPictureInfo.layerNumber, 1, 1, fp ) ) return -1;
    if ( -1 == writefile( &currPictureInfo.subSequenceIdentifier, 2, 1, fp ) ) return -1;
    if ( currPictureInfo.syncPictureFlag )
    {
      if ( -1 == writefile( &currPictureInfo.originLayerNumber, 1, 1, fp ) ) return -1;
      if ( -1 == writefile( &currPictureInfo.originSubSequenceIdentifier, 2, 1, fp ) ) return -1;
    }
  }

  num += writefile_s( &currPictureInfo.numPayloads, sizeof(currPictureInfo.numPayloads), box_fh.numBytesInPayloadCountMinusOne + 1, 1, fp );

  // check if the operations to write file are successful.
  if ( num != (unsigned)(1+box_fh.numBytesInPictureOffsetMinusTwo + 2+box_fh.numBytesInPictureDisplayTimeMinusOne + 1+box_fh.numBytesInPayloadCountMinusOne + 1) ) return -1;

  pp = currPictureInfo.payloadData;
  for ( i = 0; i < currPictureInfo.numPayloads; i++ )
  { 
    assert( pp != NULL );

    // update and write the payload to file
    if ( -1 == wrPayloadInfo( pp, fp ) ) return -1;

    // then update the parameter in Alternate Track Info Box
    size = pp->payloadSize + pp->headerSize;
    box_ati.info[0].numSDU++;
    box_ati.info[0].sumSDUSize += (long double)size;
    if ( box_ati.info[0].maxSDUSize < size ) box_ati.info[0].maxSDUSize = (unsigned INT16)size;

    pp = pp->next;
  }

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
  assert ( input->partition_mode == PAR_DP_1 );

  // Tian Dong: JVT-C083. June 15, 2002
  // Calculating payload type, the conditions may be updated when more
  // coding options are supported later.
  if ( FirstFrameIn2ndIGOP == img->number )
    pli->payloadType = PAYLOAD_TYPE_IDERP;

  pli->parameterSet = box_ps.parameterSetID;
  pli->pictureID = img->currentSlice->picture_id;
  pli->sliceID = img->current_slice_nr;
  pli->pstruct = img->pstruct;

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
#ifdef _ABT_FLAG_IN_SLICE_HEADER_
  pli->abtMode   = input->abt;
#endif
  pli->initialQP = img->qp;
  pli->qpsp = img->qpsp;

  // Tian: begin of ERPS
  pli->numRMPNI = 0;
  if ( img->type != INTRA_IMG )
  {
    // do the possible remapping 
    remap_ref_short_term(pli);
  }

  // begin of ERPS
#ifdef _CHECK_MULTI_BUFFER_1_

  // Tian Dong: It is suggested to delete the lines 
  // for test purpose. June 7, 2002
  printf("_CHECK_MULTI_BUFFER_1_ cannot be defined!\n");
  exit(0);

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

  // Tian Dong: It is suggested to delete such check coding lines
  // for test purpose. June 7, 2002
  assert(0);

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
  int start_mb_nr;

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
  if ( pp->payloadType == 0 || pp->payloadType == PAYLOAD_TYPE_IDERP )
  {
    // write the parameter set
    sym.value1 = pp->parameterSet;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);
    // write picture structure
    sym.value1 = pp->pstruct;
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
    start_mb_nr = (img->width/16)*pp->firstMBInSliceY+pp->firstMBInSliceX;
#ifdef _ABT_FLAG_IN_SLICE_HEADER_
    sym.value1 = pp->abtMode;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);
#endif
    sym.mapping = dquant_linfo;
    sym.value1 = pp->initialQP - (MAX_QP - MIN_QP +1)/2;
    writeSyntaxElement2Buf_UVLC (&sym, bitstream);
 
    if ( pp->sliceType2 ==SP_IMG)
    {
      sym.value1 = pp->qpsp - (MAX_QP - MIN_QP +1)/2;
      writeSyntaxElement2Buf_UVLC (&sym, bitstream);
    }
    sym.mapping = n_linfo2;
    iff_writeERPS(&sym, pp, bitstream); // Tian: to support ERPS (Annex U), Feb 27, 2002

    // Tian Dong: June 10, 2002
    // Update: a differential value is conveyed rather than an absolute value. 
  }
  else if ( pp->payloadType == 1 )
  {
    // write the parameter set
    sym.value1 = pp->parameterSet;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);
    // write picture structure
    sym.value1 = pp->pstruct;
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
#ifdef _ABT_FLAG_IN_SLICE_HEADER_
    sym.value1 = pp->abtMode;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);
#endif
    sym.mapping = dquant_linfo;
    sym.value1 = pp->initialQP - (MAX_QP - MIN_QP + 1)/2;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);
    sym.mapping = n_linfo2;
    // write the slice id
    sym.value1 = pp->sliceID;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);

    iff_writeERPS(&sym, pp, bitstream);       // Tian: to support ERPS (Annex U), Feb 27, 2002
  }
  else if ( pp->payloadType == 2 || pp->payloadType == 3 )
  {
    // write picture structure
    sym.value1 = pp->pstruct;
    writeSyntaxElement2Buf_UVLC(&sym, bitstream);
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
  num += writefile_s( &pp->payloadSize, sizeof(pp->payloadSize), box_fh.numBytesInPayloadSizeMinusOne+1, 1, fp );
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
  int i; 

  // RPSF: Reference Picture Selection Flags
  sym->value1 = 0;
  len += writeSyntaxElement2Buf_UVLC(sym, bitstream);

  // PN: Picture Number
  sym->value1 = pp->pn;
  len += writeSyntaxElement2Buf_UVLC(sym, bitstream);

  if (pp->numRMPNI)
  {
    sym->value1 = 1;
    len += writeSyntaxElement2Buf_UVLC(sym, bitstream);

    // now write the data:
    for (i=0; i<pp->numRMPNI; i++)
    {
      assert( pp->rmpni_RMPNI[i] >= 0 && pp->rmpni_RMPNI[i] <= 3 );
      sym->value1 = pp->rmpni_RMPNI[i];
//      printf("write RMPNI %d\n", pp->rmpni_RMPNI[i]);
      len += writeSyntaxElement2Buf_UVLC(sym, bitstream);
      if ( pp->rmpni_RMPNI[i] != 3 )
      {
        sym->value1 = pp->rmpni_Data[i];
        len += writeSyntaxElement2Buf_UVLC(sym, bitstream);
      }
    }
  }
  else
  {
    sym->value1 = 0;
    len += writeSyntaxElement2Buf_UVLC(sym, bitstream);
  }

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
//  sym->value1 = 0;
//  len += writeSyntaxElement2Buf_UVLC(sym, bitstream);

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

/*!
 ************************************************************************
 * \brief
 *      Initiate the Layer Box
 *      Tian Dong, May 30, 2002:
 * \return
 *      0, if success
 *      -1, otherwise
 ************************************************************************
 */
int initLayerBox()
{
  int i;
  if ( input->NumFramesInELSubSeq == 0 ) return 0;
  if ( box_ath.numLayers > MAX_LAYER_NUMBER ) return -1;
  for (i=0; i<box_ath.numLayers; i++)
  {
    box_layr[i].type.size = 0;  // to be updated
    box_layr[i].type.type = BOX_LAYR;
    box_layr[i].avgBitRate = 0;
    box_layr[i].avgFrameRate = 0;

    box_layr[i].fp = tmpfile();
    if ( box_layr[i].fp == NULL ) return -1;
  }
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      Update the average bit rate and frame rate in this Layer
 * \return
 *      0, if success
 *      -1, if failed
 ************************************************************************
 */
int updateLayerBox()
{
  int i;
  int pictureDataSize;

  if ( input->NumFramesInELSubSeq == 0 ) return 0;
  for (i=0; i<box_ath.numLayers; i++)
  {
    assert( box_layr[i].fp != NULL );
    // update the head data
    fseek( box_layr[i].fp, 0, SEEK_END );
    pictureDataSize = ftell( box_layr[i].fp );
    box_layr[i].type.size = SIZEOF_BOXTYPE + 8 + pictureDataSize;

    // replace the following two lines to calculate the average values
    box_layr[i].avgBitRate = 0;
    box_layr[i].avgFrameRate = 0;
  }
  return 0;
}

size_t mergeLayerBox( FILE* fp )
{
  size_t num = 0, num2 = 0;
  int layr;
  FILE* sourcef, *destf;
  unsigned char c;

  if ( input->NumFramesInELSubSeq == 0 ) return 0;
  for (layr=0; layr<box_ath.numLayers; layr++)
  {
    num = 0;
    num += writefile( &box_layr[layr].type.size, 4, 1, fp );
    num += writefile( &box_layr[layr].type.type, 4, 1, fp );

    num += writefile( &box_layr[layr].avgBitRate, 4, 1, fp );
    num += writefile( &box_layr[layr].avgFrameRate, 4, 1, fp );

    // append the data in box_layr[layr].fp into fp:
    sourcef = box_layr[layr].fp;
    destf = fp;

    fseek( sourcef, 0L, SEEK_SET );
  
    c = fgetc(sourcef);
    while ( !feof( sourcef ) )
    {
      fputc( c, destf );
      num++;
      c = fgetc(sourcef);
    }

    if ( num != box_layr[layr].type.size ) return -1;
    num2 += num;
  }
  return num2;
}

void freeLayerBox()
{
  int layr;
  if ( input->NumFramesInELSubSeq == 0 ) return;
  for (layr=0; layr<box_ath.numLayers; layr++)
  {
    fclose( box_layr[layr].fp );
  }
}

/*!
 ************************************************************************
 * \brief
 *      Initiate the Sub Sequence Box
 *      Tian Dong, May 30, 2002:
 * \return
 *      0, if success
 *      -1, otherwise
 ************************************************************************
 */
int initSubSequenceBox(int layr)
{
  static unsigned INT16 id=0;

  if ( input->NumFramesInELSubSeq == 0 ) return 0;
  if (layr<0 || layr>box_ath.numLayers) return -1;

  box_sseq[layr].type.size = 0;  // to be updated
  box_sseq[layr].type.type = BOX_SSEQ;

  // set the data in current sub-sequence here...
  box_sseq[layr].subSequenceIdentifier = id++;
  box_sseq[layr].continuationFromPreviousSegmentFlag = FALSE;
  box_sseq[layr].continuationToNextSegmentFlag = FALSE;
  box_sseq[layr].startTickAvailableFlag = 0;
  box_sseq[layr].ssStartTick = 0;
  box_sseq[layr].ssDuration = 0;
  box_sseq[layr].avgBitRate = 0;
  box_sseq[layr].avgFrameRate = 0;
  box_sseq[layr].numReferenceSubSequences = 0; // to be updated. 

  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      Update the Sub Sequence Box
 *      Tian Dong, May 30, 2002:
 * \return
 *      number of bytes being written
 *      -1, if failed
 ************************************************************************
 */
int updateSubSequenceBox( int layr )
{
  if ( input->NumFramesInELSubSeq == 0 ) return 0;
  assert( layr >= 0 && layr < box_ath.numLayers );

  box_sseq[layr].type.size = SIZEOF_BOXTYPE + 30 + 3*box_sseq[layr].numReferenceSubSequences;
  box_sseq[layr].ssDuration = 0;
  box_sseq[layr].avgBitRate = 0;
  box_sseq[layr].avgFrameRate = 0;
  return 0;
}

/*!
 ************************************************************************
 * \brief
 *      Write the Sub Sequence Box
 *      Tian Dong, May 30, 2002:
 * \input
 *      layr: the layor number
 *      fp: the destination file, normally it is box_layr[layr].fp
 * \return
 *      number of bytes being written
 *      -1, if failed
 ************************************************************************
 */
size_t wrSubSequenceBox( int layr )
{
  size_t num = 0, num2;
  unsigned INT16 cd;
  unsigned int i;
  FILE *fp;

  if ( input->NumFramesInELSubSeq == 0 ) return 0;
  fp = box_layr[layr].fp;
  assert( fp );
  if ( box_sseq[layr].numReferenceSubSequences < 0 || box_sseq[layr].numReferenceSubSequences > MAX_DEPENDENT_SUBSEQ )
    return -1;

  // since we now only deal with the ONE track case, we assert:
  assert( box_fh.numAlternateTracks == 1 );
  assert( layr >= 0 && layr < box_ath.numLayers );

  num += writefile( &box_sseq[layr].type.size, 4, 1, fp );
  num += writefile( &box_sseq[layr].type.type, 4, 1, fp );
  num += writefile( &box_sseq[layr].subSequenceIdentifier, 2, 1, fp );

  if ( box_sseq[layr].continuationFromPreviousSegmentFlag ) cd = 0x8000;
  else cd = 0;
  if ( box_sseq[layr].continuationToNextSegmentFlag) cd |= 0x4000;
  if ( box_sseq[layr].startTickAvailableFlag) cd |= 0x2000;
  num += writefile( &cd, 2, 1, fp );

  num += writefile( &box_sseq[layr].ssStartTick, 8, 1, fp );
  num += writefile( &box_sseq[layr].ssDuration, 8, 1, fp );
  num += writefile( &box_sseq[layr].avgBitRate, 4, 1, fp );
  num += writefile( &box_sseq[layr].avgFrameRate, 4, 1, fp );
  num += writefile( &box_sseq[layr].numReferenceSubSequences, 2, 1, fp );
  
  if ( num != SIZEOF_BOXTYPE+30 ) return -1;

  for (i=0; i<box_sseq[layr].numReferenceSubSequences; i++)
  {
    num2 = writefile( &box_sseq[layr].dependencyData[i].layerNumber, 1, 1, fp );
    num2 += writefile( &box_sseq[layr].dependencyData[i].subSequenceIdentifier, 2, 1, fp );
    num += num2;
    if ( num2 != 3 ) return -1;
  }

  return num;
}

/*!
 ************************************************************************
 * \brief
 *      Free the resources allocated for the Sub Sequence Box
 *      Tian Dong, May 30, 2002:
 ************************************************************************
 */
void freeSubSequenceBox( int layr )
{
  if ( input->NumFramesInELSubSeq == 0 ) return;
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
      fwrite( p++, 1, 1, fp );
    else
      fwrite( p--, 1, 1, fp );
    size--;
    num++;
  }
  return num;
}

/*!
 ************************************************************************
 * \brief
 *      write the data to file, bytes are in big Endian order.
 *      to be used if buffers size differs from number of written bytes
 * \return
 *      how many bytes being written into the file, if success
 *      -1, otherwise
 * \param outf
 *      output file pointer
 ************************************************************************
 */
size_t writefile_s( void* buf, size_t bufsize, size_t size, size_t count, FILE* fp )
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
      fwrite( p++, 1, 1, fp );
    else
      fwrite( p--, 1, 1, fp );
    size--;
    num++;
  }
  return num;
}

/*!
 ************************************************************************
 * \date:
 *      June 15, 2002
 * \brief
 *      If the first frame in a new sub-sequence is to be encoded,
 *        this function will do some initialization for the new sub-seq.
 ************************************************************************
 */
void begin_sub_sequence()
{
  if ( input->of_mode != PAR_OF_IFF || input->NumFramesInELSubSeq == 0 ) return;

  // begin to encode the base layer subseq?
  if ( IMG_NUMBER == 0 )
  {
//    printf("begin to encode the base layer subseq\n");
    initSubSequenceBox(0);  // init the sub-sequence in the base layer
  }
  // begin to encode the enhanced layer subseq?
  if ( IMG_NUMBER % (input->NumFramesInELSubSeq+1) == 1 )
  {
//    printf("begin to encode the enhanced layer subseq\n");
    initSubSequenceBox(1);  // init the sub-sequence in the enhanced layer
    add_dependent_subseq(1);
  }
}

/*!
 ************************************************************************
 * \date:
 *      June 15, 2002
 * \brief
 *      If the last frame of current sub-sequence has been encoded,
 *        this function will do some finalization for the sub-seq.
 ************************************************************************
 */
void end_sub_sequence()
{
  if ( input->of_mode != PAR_OF_IFF || input->NumFramesInELSubSeq == 0 ) return;

  // end of the base layer:
  if ( img->number == input->no_frames-1 )
  {
//    printf("end of encoding the base layer subseq\n");
    updateSubSequenceBox(0);
    if ( -1 == wrSubSequenceBox(0) )
    {
      printf("Serious warning: error occurs when trying to write sub sequence box 0.\n");
    }
  }
  // end of the enhanced layer:
  if ( ((IMG_NUMBER%(input->NumFramesInELSubSeq+1)==0) && (input->successive_Bframe != 0) && (IMG_NUMBER>0)) || // there are B frames
    ((IMG_NUMBER%(input->NumFramesInELSubSeq+1)==input->NumFramesInELSubSeq) && (input->successive_Bframe==0))   // there are no B frames
    )
  {
//    printf("end of encoding the enhanced layer subseq\n");
    add_dependent_subseq(1);
    updateSubSequenceBox(1);
    if ( -1 == wrSubSequenceBox(1) )
    {
      printf("Serious warning: error occurs when trying to write sub sequence box 1.\n");
    }
  }
}


/*!
 ************************************************************************
 * \date:
 *      June 15, 2002
 * \brief
 *    Hide some frames in the short term frame buffer
 *      The number of frames can be used for the forward prediction will
 *      be set in fb->num_short_used.
 ************************************************************************
 */
void remap_ref_short_term(PayloadInfo* pp)
{
  int i;
  RMPNIbuffer_t *r, *t;
  int currLayerNo, currSubSeqNo;
  int pnp, pnq;
  int delta, mdelta;
  Boolean need_rmpni = FALSE;

  currLayerNo = currPictureInfo.layerNumber;
  currSubSeqNo = currPictureInfo.subSequenceIdentifier;
  pnp = img->pn;

  // check if the remapping is need or not?
  for (i=0; i<fb->short_used; i++)
  {
    if ( fb->picbuf_short[i]->layer_no > currLayerNo ||
      (fb->picbuf_short[i]->layer_no == currLayerNo && fb->picbuf_short[i]->sub_seq_no != currSubSeqNo)
      )
    {
      need_rmpni = TRUE;
      break;
    }
  }

  if ( !need_rmpni ) return;

  fb->num_short_used = 0;

  // Tian Dong: we use the reverse order for re-mapping to let the latest frame
  // to have the relative index 0. See our implementation accompanying document, section 1.2.5
  t = img->currentSlice->rmpni_buffer;
  for (i=fb->short_used-1; i>=0; i--)
  {
    if ( fb->picbuf_short[i]->layer_no < currLayerNo ||
      (fb->picbuf_short[i]->layer_no == currLayerNo && fb->picbuf_short[i]->sub_seq_no == currSubSeqNo)
      )
    {
      r = (RMPNIbuffer_t*)calloc (1,sizeof(RMPNIbuffer_t));
      r->Next=NULL;
      fb->num_short_used++;

      // caculate the abs_diff_pic_numbers from picID. TBD...
      pnq = fb->picbuf_short[i]->picID;
      delta = pnq - pnp;
      if ( delta < 0 )
      {
        if ( delta < -fb->short_size/2 -1 )
          mdelta = delta+fb->short_size;
        else
          mdelta = delta;
      }
      else
      {
        if ( delta > fb->short_size/2 )
          mdelta = delta - fb->short_size;
        else
          mdelta = delta;
      }

      r->Data = abs( mdelta );
      if ( mdelta < 0 ) r->RMPNI = 0;
      else  r->RMPNI=1;

      pp->rmpni_Data[pp->numRMPNI] = r->Data;
      pp->rmpni_RMPNI[pp->numRMPNI] = r->RMPNI;
      pp->numRMPNI++;
      assert( pp->numRMPNI < 6 );

      // add a new node to the list
      if ( img->currentSlice->rmpni_buffer == NULL )
        img->currentSlice->rmpni_buffer = r;
      else
        t->Next = r;
      t = r;

      pnp=pnq;
    }
  }
  // the end loop code number is 3
  if ( img->currentSlice->rmpni_buffer != NULL )
  {
    r = (RMPNIbuffer_t*)calloc (1,sizeof(RMPNIbuffer_t));
    r->RMPNI=3;
    r->Data=0;
    r->Next=NULL;
    t->Next = r;

    pp->rmpni_Data[pp->numRMPNI] = r->Data;
    pp->rmpni_RMPNI[pp->numRMPNI] = r->RMPNI;
    pp->numRMPNI++;
    assert( pp->numRMPNI < 6 );
  }

  reorder_mref(img);

  for (i=fb->num_short_used;i<fb->short_used;i++)
  {
    mref[i]=NULL;
    mcef[i]=NULL;
    Refbuf11[i]=NULL;
  }

}

/*!
 ************************************************************************
 * \date:
 *      June 15, 2002
 * \brief
 *      To add the dependent sub-sequence to current sub-sequence.
 *      (stored in dependencyData[*])
 *      Note: It only check the possible dependent frames for the frame
 *      to-be-encoded, thus to let all the dependent sub-sequence be 
 *      collected, this function should be called twice: in the 
 *      begin_sub_sequence() and end_sub_sequence().
 * \param layr
 *    the layer number of current sub-sequence: box_sseq[this_layr].
 ************************************************************************
 */
void add_dependent_subseq(int layr)
{
  int i;
  for (i=0; i<frm->short_used; i++)
  {
    if (!in_dependency_set(layr, frm->picbuf_short[i]->sub_seq_no, frm->picbuf_short[i]->layer_no) &&
      can_predict_from(layr, frm->picbuf_short[i]->sub_seq_no, frm->picbuf_short[i]->layer_no) 
      )
    {
      assert(box_sseq[layr].numReferenceSubSequences+1 < MAX_DEPENDENT_SUBSEQ);
      box_sseq[layr].dependencyData[box_sseq[layr].numReferenceSubSequences].subSequenceIdentifier = 
        frm->picbuf_short[i]->sub_seq_no;
      box_sseq[layr].dependencyData[box_sseq[layr].numReferenceSubSequences].layerNumber = 
        frm->picbuf_short[i]->layer_no;
      box_sseq[layr].numReferenceSubSequences++;
    }
  }
}

/*!
 ************************************************************************
 * \date:
 *      June 15, 2002
 * \brief
 *      To judge if any frame from sub-sequence 1 is in the dependency 
 *      data of sub-sequence 2.
 * \return
 *    FALSE: if not.
 *    TRUE: if yes.
 * \param this_layr
 *    the layer number of sub-sequence 2, box_sseq[this_layr] is sub-seq 2
 * \param sub_seq_no
 *    the id of sub-sequence 1
 * \param layer_no
 *    the layer no of sub-sequence 1
 ************************************************************************
 */
Boolean in_dependency_set(int this_layr, int sub_seq_no, int layer_no)
{
  unsigned int i;
  for (i=0; i<box_sseq[this_layr].numReferenceSubSequences; i++)
  {
    if ( box_sseq[this_layr].dependencyData[i].subSequenceIdentifier == sub_seq_no &&
      box_sseq[this_layr].dependencyData[i].layerNumber == layer_no 
      )
      return TRUE;
  }
  return FALSE;
}

/*!
 ************************************************************************
 * \date:
 *      June 15, 2002
 * \brief
 *      To judge if sub-sequence 1 can predict from sub-sequence 2
 * \return
 *    FALSE: if it can not or
 *           sub-sequnece 1 and sub-sequence 2 are the same one;
 *    TRUE: otherwise.
 * \param this_layr
 *    the layer number of sub-sequence 1, box_sseq[this_layr] is sub-seq 1
 * \param sub_seq_no
 *    the id of sub-sequence 2
 * \param layer_no
 *    the layer no of sub-sequence 2
 ************************************************************************
 */
Boolean can_predict_from(int this_layr, int sub_seq_no, int layer_no)
{
  if (this_layr>layer_no) return TRUE;
  if (this_layr<layer_no) return FALSE;
//  if (box_sseq[this_layr].subSequenceIdentifier == sub_seq_no) return TRUE;
  return FALSE;
}
