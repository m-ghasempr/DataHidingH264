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
 *     encodeiff.h
 *  \brief
 *     definitions for H.26L interim file format, as defined in VCEG-O58
 *  \author
 *      - Tian, Dong                             <tian@cs.tut.fi>
 *
 * ************************************************************************
 */

#ifndef ENCODEIFF_H
#define ENCODEIFF_H

#define SIZEOF_BOXTYPE 8  // 8: compact size, 32 bits mode, 16: extended size, 64 bits mode
#define BUFSIZE_FOR_PAYLOADINFO 2048

typedef byte __int2;

//! Box Types
typedef enum
{
    BOX_FTYP, //<! 
    BOX_JVTH, //<! 
    BOX_CINF, //<! 
    BOX_ATIN, //<! 
    BOX_PRMS, //<! 
    BOX_SEGM, //<! 
    BOX_ATRH, //<! 
    BOX_SWPC, //<! 
    BOX_ATRM  //<! 
} TYPE_OF_BOX;

typedef struct
{
  unsigned __int32 size;
  unsigned __int32 type;
//  unsigned __int64 largesize;
} BoxType;

// 1
typedef struct
{
  BoxType type;

  unsigned char    majorBrand[4];
  unsigned __int16 jmMajorVersion;
  unsigned __int16 jmMinorVersion;
  unsigned int     numCompatibleBrands;
  unsigned char*   compatibleBrands;
} FileTypeBox;

// 2
typedef struct
{
  BoxType type;

  unsigned __int8 majorVersion;
  unsigned __int8 minorVersion;
  unsigned __int32 timescale;
  unsigned __int32 numUnitsInTick;
  unsigned __int64 duration;
  unsigned __int16 pixAspectRatioX;
  unsigned __int16 pixAspectRatioY;
  unsigned __int16 maxPicId;
  unsigned __int8 numAlternateTracks;
  __int2 numBytesInPayloadCountMinusOne;
  __int2 numBytesInPictureOffsetMinusTwo;
  __int2 numBytesInPictureDisplayTimeMinusOne;
  __int2 numBytesInPictureCountMinusOne;
  __int2 numBytesInPayloadSizeMinusOne;
} FileHeaderBox;

// 3
typedef struct
{
  BoxType type;

  unsigned __int64 creationTime;
  unsigned __int64 modificationTime;

  unsigned __int8 titleNumBytes;
  unsigned char* title;

  unsigned __int8 authorNumBytes;
  unsigned char* author;

  unsigned __int8 copyrightNumBytes;
  unsigned char* copyright;

  unsigned __int16 descriptionNumBytes;
  unsigned char* description;

  unsigned __int16 URINumBytes;
  unsigned char* URI;
} ContentInfoBox;

// 4
typedef struct
{
  unsigned __int16 displayWindowWidth;
  unsigned __int16 displayWindowHeight;
  unsigned __int16 maxSDUSize;
  unsigned __int16 avgSDUSize;
  unsigned __int32 avgBitRate;
  long double      sumSDUSize;
  long double      numSDU;
  int              last_frame;
} AlternateTrackInfo;

typedef struct
{
  BoxType type;
  AlternateTrackInfo *info;
} AlternateTrackInfoBox;

// 5
typedef struct
{
  BoxType type;

  unsigned __int16 parameterSetID;
  unsigned __int8 profile;
  unsigned __int8 level;
  unsigned __int8 version;
  unsigned __int16 pictureWidthInMBs;
  unsigned __int16 pictureHeightInMBs;
  unsigned __int16 displayRectangleOffsetTop;
  unsigned __int16 displayRectangleOffsetLeft;
  unsigned __int16 displayRectangleOffsetBottom;
  unsigned __int16 displayRectangleOffsetRight;
  unsigned __int8 displayMode;
  unsigned __int16 displayRectangleOffsetFromWindowTop;
  unsigned __int16 displayRectangleOffsetFromWindowLeftBorder;
  unsigned __int8 entropyCoding;
  unsigned __int8 motionResolution;
  unsigned __int8 partitioningType;
  unsigned __int8 intraPredictionType;
} ParameterSetBox;

// 6
typedef struct
{
  BoxType type;

  unsigned __int64 fileSize;
  unsigned __int64 startTick;
  unsigned __int64 segmentDuration;
  unsigned __int64 firstFrameNr;
  unsigned __int64 lastFrameNr;
} SegmentBox;

// 7

typedef struct sPayloadInfo
{
  unsigned __int64 payloadSize;
  unsigned __int8  headerSize;
  unsigned __int8  payloadType;
  unsigned __int8  errorIndication;
  unsigned __int8  reserved;

  unsigned __int16 parameterSet;

  unsigned __int8  pictureID;
  unsigned __int8  sliceID;

  unsigned __int8  sliceType;
  unsigned __int8  firstMBInSliceX;
  unsigned __int8  firstMBInSliceY;
  signed __int8  initialQP;

  int sliceType2;   // save according to original value
  int qpsp;

  int lastMBnr;
  
  struct sPayloadInfo* next;

} PayloadInfo;

typedef struct
{
  Boolean intraPictureFlag;
  __int64 pictureOffset;
  __int64 currPictureSize;
  __int64 pictureDisplayTime;
  unsigned __int64 numPayloads;

  __int64 lastFrameNr;
  
  PayloadInfo* payloadData;

} PictureInfo;

typedef struct
{
  BoxType type;
  unsigned __int64 numPictures;

  FILE* fpMeta;   // save the PictureInfo to tmp files
} AlternateTrackHeaderBox;

typedef struct
{
  BoxType type;
  // more attributes:

} SwitchPictureBox;

typedef struct
{
  BoxType type;
  FILE* fpMedia;
} AlternateTrackMediaBox;

extern FileTypeBox box_ft;
extern FileHeaderBox box_fh;
extern ContentInfoBox box_ci;
extern AlternateTrackInfoBox box_ati;
extern ParameterSetBox box_ps;
extern SegmentBox box_s;
extern AlternateTrackHeaderBox box_ath;
extern PayloadInfo* pCurrPayloadInfo;
extern PictureInfo currPictureInfo;
extern AlternateTrackMediaBox box_atm;
extern SwitchPictureBox box_sp;

// functions

// Functions on FileTypeBox
int initFileTypeBox();
size_t wrFileTypeBox(FILE* fp);
void freeFileTypeBox();

// Functions on FileHeaderBox
int initFileHeaderBox();
size_t wrFileHeaderBox( FILE* fp );
void freeFileHeaderBox();

// Functions on ContentInfoBox
int initContentInfoBox();
size_t wrContentInfoBox( FILE* fp );
void freeContentInfoBox();

// Functions on AlternateTrackInfoBox
int initAlternateTrackInfoBox();
size_t wrAlternateTrackInfoBox( FILE* fp );
void freeAlternateTrackInfoBox();

// Functions on ParameterSetBox
int initParameterSetBox();
size_t wrParameterSetBox( FILE* fp );
void freeParameterSetBox();

// Functions on SegmentBox
int initSegmentBox();
void updateSegmentBox();
size_t wrSegmentBox( FILE *fp );
void freeSegmentBox();

// Functions on AlternateTrackHeaderBox
int initAlternateTrackHeaderBox();
int updateAlternateTrackHeaderBox();
size_t mergeAlternateTrackHeaderBox( FILE* fp );
void freeAlternateTrackHeaderBox();

// Functions on PictureInfo
int initPictureInfo();
size_t wrPictureInfo( FILE* fp );
void freePictureInfo();


// Functions on payloadInfo
PayloadInfo* newPayloadInfo();
int addOnePayloadInfo(PayloadInfo* pi);
size_t wrPayloadInfo( PayloadInfo* pp, FILE *fp );

// Functions on SwitchPictureBox
int initSwitchPictureBox();

// Functions on AlternateMediaBox
int initAlternateTrackMediaBox();
int updateAlternateTrackMediaBox();
int mergeAlternateTrackMediaBox( FILE* fp );
void freeAlternateTrackMediaBox();

// Other Functions 
int initInterimFile();
size_t terminateInterimFile(FILE* outf);
size_t writefile( void* buf, size_t size, size_t count, FILE* fp );

#endif
