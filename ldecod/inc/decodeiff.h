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
 *     decodeiff.h
 *  \brief
 *     definitions for H.26L interim file format, as defined in VCEG-O58
 *  \author
 *      - Dong Tian                             <tian@cs.tut.fi>
 ************************************************************************
 */

#ifndef DECODEIFF_H
#define DECODEIFF_H

#ifdef _WIN32
#define INT2 __int8
#define INT8 __int8
#define INT16 __int16
#define INT32 __int32
#define INT64 __int64
#else
#define INT2 char
#define INT8 char
#define INT16 short
#define INT32 long
#define INT64 long long int   // this may be not 64 bit on some compilers
#endif

#define SIZEOF_BOXTYPE 8  // 8: 32 bits mode, 16: 64 bits mode

#define MAX_LAYER_NUMBER 2
#define MAX_DEPENDENT_SUBSEQ 5

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
    BOX_PICI, //<! 
    BOX_LAYR, //<! 
    BOX_SSEQ, //<! 
    BOX_SWPC, //<! 
    BOX_ATRM  //<! 
} TYPE_OF_BOX;

// 0
typedef struct
{
  unsigned INT32 size;
  unsigned INT32 type;
  unsigned INT64 largesize;
} BoxType;


#define ALTERNATETRACK_MAXNUNBER 10 
#define BUFSIZE_FOR_PAYLOADINFO 2048

// 1
typedef struct
{
  BoxType type;

  unsigned char    majorBrand[4];
  unsigned INT16 jmMajorVersion;
  unsigned INT16 jmMinorVersion;
  unsigned int     numCompatibleBrands;
  unsigned char*   compatibleBrands;
} FileTypeBox;

// 2
typedef struct
{
  BoxType type;

  unsigned INT8 majorVersion;
  unsigned INT8 minorVersion;
  unsigned INT32 timescale;
  unsigned INT32 numUnitsInTick;
  unsigned INT64 duration;
  unsigned INT16 pixAspectRatioX;
  unsigned INT16 pixAspectRatioY;
  unsigned INT16 maxPicId;  
  unsigned INT8 numAlternateTracks;
  INT2 numBytesInPayloadCountMinusOne;
  INT2 numBytesInPictureOffsetMinusTwo;
  INT2 numBytesInPictureDisplayTimeMinusOne;
  INT2 numBytesInPictureCountMinusOne;
  INT2 numBytesInPayloadSizeMinusOne;
} FileHeaderBox;

// 3
typedef struct
{
  BoxType type;

  unsigned INT64 creationTime;
  unsigned INT64 modificationTime;

  unsigned INT8 titleNumBytes;
  unsigned char* title;

  unsigned INT8 authorNumBytes;
  unsigned char* author;

  unsigned INT8 copyrightNumBytes;
  unsigned char* copyright;

  unsigned INT16 descriptionNumBytes;
  unsigned char* description;

  unsigned INT16 URINumBytes;
  unsigned char* URI;
} ContentInfoBox;

// 4
typedef struct
{
  unsigned INT16 displayWindowWidth;
  unsigned INT16 displayWindowHeight;
  unsigned INT16 maxSDUSize;
  unsigned INT16 avgSDUSize;
  unsigned INT32 avgBitRate;
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

  unsigned INT16 parameterSetID;
  unsigned INT8 profile;
  unsigned INT8 level;
  unsigned INT8 version;
  unsigned INT16 pictureWidthInMBs;
  unsigned INT16 pictureHeightInMBs;
  unsigned INT16 displayRectangleOffsetTop;
  unsigned INT16 displayRectangleOffsetLeft;
  unsigned INT16 displayRectangleOffsetBottom;
  unsigned INT16 displayRectangleOffsetRight;
  unsigned INT8 displayMode;
  unsigned INT16 displayRectangleOffsetFromWindowTop;
  unsigned INT16 displayRectangleOffsetFromWindowLeftBorder;
  unsigned INT8 entropyCoding;
  unsigned INT8 motionResolution;
  unsigned INT8 partitioningType;
  unsigned INT8 intraPredictionType;
  unsigned INT8 bufCycle;
  Boolean requiredPictureNumberUpdateBehavior;
} ParameterSetBox;

// 6
typedef struct
{
  BoxType type;

  unsigned INT64 fileSize;
  unsigned INT64 startTick;
  unsigned INT64 segmentDuration;
  unsigned INT64 firstFrameNr;
  unsigned INT64 lastFrameNr;
} SegmentBox;

// 7
typedef struct
{
  BoxType type;
  unsigned INT8 numLayers;
  long storedpos;
} AlternateTrackHeaderBox;  // 042

// 8

typedef struct sPayloadInfo
{
  unsigned INT64 payloadSize;
  unsigned INT8  headerSize;
  unsigned INT8  payloadType;
  unsigned INT8  errorIndication;
  unsigned INT8  reserved;

  unsigned INT16 parameterSet;

  unsigned INT8  pictureID;
  unsigned INT8  structure;
  unsigned INT8  sliceID;

  unsigned INT8  sliceType;
  unsigned INT8  firstMBInSliceX;
  unsigned INT8  firstMBInSliceY;
  unsigned INT8  directType;
  INT8  initialQP;

  int lastMBnr;

  long    storedpos;
  unsigned int     payloadnr;

  Bitstream        buffer;
} PayloadInfo;

typedef struct
{
  Boolean intraPictureFlag;
  Boolean syncPictureFlag;    // 042
  INT64 pictureOffset;
  INT64 currPictureSize;
  INT64 pictureDisplayTime;

  unsigned INT8 layerNumber;  // 042
  unsigned INT16 subSequenceIdentifier;
  unsigned INT8 originLayerNumber;
  unsigned INT16 originSubSequenceIdentifier;

  unsigned INT64 numPayloads;

  INT64 lastFrameNr;
  long picPos;
} PictureInfo;

typedef struct
{
  BoxType type;
  unsigned INT64 numPictures;
  unsigned long storedpos;
} PictureInformationBox;

// 9
typedef struct
{
  BoxType type;
  unsigned INT32 avgBitRate;
  unsigned INT32 avgFrameRate;
  FILE* fp;
} LayerBox;

// 10
typedef struct
{
  unsigned INT8 layerNumber;
  unsigned INT16 subSequenceIdentifier;
} DependencyInfo;

typedef struct
{
  BoxType type;

  unsigned INT16 subSequenceIdentifier;
  Boolean continuationFromPreviousSegmentFlag;
  Boolean continuationToNextSegmentFlag;
  Boolean startTickAvailableFlag;
  unsigned INT64 ssStartTick;
  unsigned INT64 ssDuration;
  unsigned INT32 avgBitRate;
  unsigned INT32 avgFrameRate;
  unsigned INT32 numReferenceSubSequences;
  DependencyInfo dependencyData[MAX_DEPENDENT_SUBSEQ];
} SubSequenceBox;

// 11
typedef struct
{
  BoxType type;
  // more attributes:

} SwitchPictureBox;

// 12
typedef struct
{
  BoxType type;

  long  currPictureOffset;        // media data pointer
  long  currPayloadOffset;
  unsigned long storedpos;
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

extern int isBigEndian;
// functions

int testEndian();

// Functions on FileTypeBox
int testFileTypeBox(FILE* fp);
void freeFileTypeBox();

// Functions on FileHeaderBox
int rdFileHeaderBox( FILE* fp, size_t size );
void freeFileHeaderBox();

// Functions on ContentInfoBox
int rdContentBox( FILE* fp, size_t size );
void freeContentInfoBox();

// Functions on AlternateTrackInfoBox
int rdAlternateTrackInfoBox( FILE* fp, unsigned long size );
void freeAlternateTrackInfoBox();

// Functions on ParameterSetBox
int findParameterSetBox( FILE* fp, long *pos, int pid );
void freeParameterSetBox();

// Functions on SegmentBox
int parse_one_segment( struct img_par* img, struct inp_par* inp, struct snr_par *snr, FILE* fp, size_t size );
void freeSegmentBox();

// Functions on AlternateTrackHeaderBox
void freeAlternateTrackHeaderBox();

// Functions on PictureInfo
int rdPictureInfo( FILE* fp );

// Functions on payloadInfo
void decomposeSliceHeader( struct img_par *img, struct inp_par* inp, PayloadInfo* pp );

// Functions on layer box
int rdLayerBox( int no, unsigned long boxsize, FILE *bits );

// Functions on sub sequence box
int rdSubSequence(FILE* fp, FILE* bits);

// Functions on SwitchPictureBox

// Functions on AlternateMediaBox
void freeAlternateTrackMediaBox();

// Other Functions 
int readSliceIFF( struct img_par* img, struct inp_par* inp );
int IFFSequenceHeader( struct img_par *img, struct inp_par *inp, FILE *bits );
int IFFUseParameterSet( int n, struct img_par* img, struct inp_par* inp );
int initInterimFile();
void terminateInterimFile();

int readfile( void* buf, size_t size, size_t count, FILE* fp );
int readfile_s( void* buf, size_t bufsize, size_t size, size_t count, FILE* fp );
int find_track_meta( FILE *fp, unsigned long limited, unsigned long* storedpos, int tracknr );
int find_track_media( FILE *fp, unsigned long limited, unsigned long* storedpos, int tracknr );
int rdOneTrack( struct img_par* img, struct inp_par* inp, struct snr_par *snr, FILE *fp );
int rdPayloadInfo( struct img_par *img, struct inp_par* inp, PayloadInfo* pp, FILE *fp );
int parse_one_box( struct img_par* img, struct inp_par* inp, struct snr_par *snr, FILE *fp );
int rdOnePayload( struct img_par *img, struct inp_par* inp, PayloadInfo *pp, FILE* fp );
int IFFGetFollowingSliceHeader( struct img_par *img, PayloadInfo* pp );

#endif
