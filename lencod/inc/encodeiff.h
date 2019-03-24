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
 *      - Dong Tian                             <tian@cs.tut.fi>
 *
 * ************************************************************************
 */

#ifndef ENCODEIFF_H
#define ENCODEIFF_H

#ifdef WIN32
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
#define INT64 long long int  // This may not be 64 bit on some compilers
#endif

#define SIZEOF_BOXTYPE 8  // 8: compact size, 32 bits mode, 16: extended size, 64 bits mode
#define BUFSIZE_FOR_PAYLOADINFO 2048

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

typedef struct
{
  unsigned INT32 size;
  unsigned INT32 type;
//  unsigned INT64 largesize;
} BoxType;

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
  unsigned INT8 loopFilterParametersFlag;
  unsigned INT8 entropyCoding;
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
  unsigned INT8  pstruct;
  unsigned INT8  sliceID;

  unsigned INT8  sliceType;
  unsigned INT8  firstMBInSliceX;
  unsigned INT8  firstMBInSliceY;
  unsigned INT8  directType;
#ifdef _ABT_FLAG_IN_SLICE_HEADER_
  unsigned INT8  abtMode;
#endif
  unsigned INT8  disposable_flag;
  unsigned INT8  explicit_B_prediction;
  unsigned INT8  num_ref_pic_active_fwd_minus1;
  unsigned INT8  num_ref_pic_active_bwd_minus1;

  INT8  initialQP;

  int sliceType2;   // save according to original value
  unsigned INT8  spSwitchFlag;
  int qpsp;

  int pn;           // save img->pn
  int type;         // save img->type
  int max_lindex;   // save img->max_lindex
  int lindex;       // save img->lindex

  struct sPayloadInfo* next;

  int numRMPNI;
  int rmpni_Data[6];
  int rmpni_RMPNI[6];

  unsigned INT8   filter_parameters_flag;
  unsigned INT8   lf_disable;
  INT8            lf_alpha_c0_offset_div2;
  INT8            lf_beta_offset_div2;

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

  unsigned INT8 refFromLayerNumber;  // 042
  unsigned INT16 refFromSubSequenceIdentifier;

  unsigned INT64 numPayloads;

  INT64 lastFrameNr;
  
  PayloadInfo* payloadData;

} PictureInfo;

typedef struct
{
  BoxType type;
  unsigned INT64 numPictures;

  FILE* fpMeta;   // save the PictureInfo to tmp files
} PictureInformationBox;   // 042

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
extern PictureInformationBox box_pi;
extern PayloadInfo* pCurrPayloadInfo;
extern PictureInfo currPictureInfo;
extern AlternateTrackMediaBox box_atm;
extern SwitchPictureBox box_sp;

extern int isBigEndian;

// functions

int testEndian();

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

// Functions on PictureInformationBox
int initPictureInformationBox();
int updatePictureInformationBox();
size_t mergePictureInformationBox( FILE* fp );
void freePictureInformationBox();

// Functions on PictureInfo
int initPictureInfo();
size_t wrPictureInfo( FILE* fp );
void freePictureInfo();


// Functions on payloadInfo
PayloadInfo* newPayloadInfo();
int addOnePayloadInfo(PayloadInfo* pi);
size_t wrPayloadInfo( PayloadInfo* pp, FILE *fp );
size_t iff_writeERPS(SyntaxElement *sym, PayloadInfo* pp, Bitstream* bitstream);

// Functions on LayerBox
int initLayerBox();
int updateLayerBox();
size_t mergeLayerBox( FILE* fp );
void freeLayerBox();

// Functions on SubSequenceBox
int initSubSequenceBox( int layr );
int updateSubSequenceBox( int layr );
size_t wrSubSequenceBox( int layr );
void freeSubSequenceBox( int layr );
void begin_sub_sequence();
void end_sub_sequence();

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
size_t writefile_s( void* buf, size_t bufsize, size_t size, size_t count, FILE* fp );
void remap_ref_short_term(PayloadInfo* pp);
void add_dependent_subseq(int layr);
Boolean in_dependency_set(int this_layr, int sub_seq_no, int layer_no);
Boolean can_predict_from(int this_layr, int sub_seq_no, int layer_no);

#endif

