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
 *     global.h
 *  \brief
 *     global definitions for for H.26L encoder.
 *  \author
 *     Copyright (C) 1999  Telenor Satellite Services,Norway
 *                         Ericsson Radio Systems, Sweden
 *
 *     Inge Lille-Langøy               <inge.lille-langoy@telenor.com>
 *
 *     Telenor Satellite Services
 *     Keysers gt.13                       tel.:   +47 23 13 86 98
 *     N-0130 Oslo,Norway                  fax.:   +47 22 77 79 80
 *
 *     Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *
 *     Ericsson Radio Systems
 *     KI/ERA/T/VV
 *     164 80 Stockholm, Sweden
 *
 ************************************************************************
 */
#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include <stdio.h>
#include "defines.h"

#ifndef WIN32
  #include "minmax.h"
#else
  #define  snprintf _snprintf
#endif


/***********************************************************************
 * T y p e    d e f i n i t i o n s    f o r    T M L
 ***********************************************************************
 */

typedef unsigned char byte;    //!< byte type definition
#define pel_t byte

//! Data Partitioning Modes
typedef enum
{
  PAR_DP_1,   //!< no data partitioning is supported
  PAR_DP_3,   //!< data partitioning with 3 partitions
} PAR_DP_TYPE;


//! Output File Types
typedef enum
{
  PAR_OF_26L,    //!< Current TML description
  PAR_OF_RTP,    //!< RTP packets in outfile
  PAR_OF_IFF     //!< Interim File Format
} PAR_OF_TYPE;

//! Boolean Type
typedef enum {
  FALSE,
  TRUE
} Boolean;

typedef enum {
  FRAME_CODING,
  ADAPTIVE_CODING
} CodingType;

//! definition of H.26L syntax elements
typedef enum {
  SE_HEADER,
  SE_PTYPE,
  SE_MBTYPE,
  SE_REFFRAME,
  SE_INTRAPREDMODE,
  SE_MVD,
  SE_CBP_INTRA,
  SE_LUM_DC_INTRA,
  SE_CHR_DC_INTRA,
  SE_LUM_AC_INTRA,
  SE_CHR_AC_INTRA,
  SE_CBP_INTER,
  SE_LUM_DC_INTER,
  SE_CHR_DC_INTER,
  SE_LUM_AC_INTER,
  SE_CHR_AC_INTER,
  SE_DELTA_QUANT_INTER,
  SE_DELTA_QUANT_INTRA,
  SE_BFRAME,
  SE_EOS,
  SE_MAX_ELEMENTS  //!< number of maximum syntax elements
} SE_type;         // substituting the definitions in elements.h


typedef enum {
  INTER_MB,
  INTRA_MB_4x4,
  INTRA_MB_16x16
} IntraInterDecision;


typedef enum {
  BITS_HEADER,
  BITS_TOTAL_MB,
  BITS_MB_MODE,
  BITS_INTER_MB,
  BITS_CBP_MB,
  BITS_COEFF_Y_MB,
  BITS_COEFF_UV_MB,
  BITS_DELTA_QUANT_MB,
  MAX_BITCOUNTER_MB
} BitCountType;


typedef enum {
  SINGLE_SCAN,
  DOUBLE_SCAN
} ScanMode;


typedef enum {
  NO_SLICES,
  FIXED_MB,
  FIXED_RATE,
  CALLBACK,
  FMO
} SliceMode;


typedef enum {
  UVLC,
  CABAC
} SymbolMode;


typedef enum {        // ABT
  NO_ABT,
  INTER_ABT,
  INTER_INTRA_ABT
} ABTCodingMode;      // ABT

typedef enum {        // ABT
  B8x8=0,
  B8x4,
  B4x8,
  B4x4
} ABTBlockMode;       // ABT

// mwi: use same enum as in decoder
typedef enum {
  FRAME,
  TOP_FIELD,
  BOTTOM_FIELD
} PictureType;           //!< New enum for field processing
// ~mwi


/***********************************************************************
 * D a t a    t y p e s   f o r  C A B A C
 ***********************************************************************
 */

//! struct to characterize the state of the arithmetic coding engine
typedef struct
{
  unsigned int  Elow, Erange;
  unsigned int  Ebuffer;
  unsigned int  Ebits_to_go;
  unsigned int  Ebits_to_follow;
  byte          *Ecodestrm;
  int           *Ecodestrm_len;
  // storage in case of recode MB
  unsigned int  ElowS, ErangeS;
  unsigned int  EbufferS;
  unsigned int  Ebits_to_goS;
  unsigned int  Ebits_to_followS;
  byte          *EcodestrmS;
  int           *Ecodestrm_lenS;
} EncodingEnvironment;

typedef EncodingEnvironment *EncodingEnvironmentPtr;

//! struct for context management
typedef struct
{
  unsigned int  cum_freq[2];    //!< cumulated frequency counts
  Boolean       in_use;         //!< flag for context in use
  unsigned int  max_cum_freq;   //!< maximum frequency count

} BiContextType;

typedef BiContextType *BiContextTypePtr;


/**********************************************************************
 * C O N T E X T S   F O R   T M L   S Y N T A X   E L E M E N T S
 **********************************************************************
 */


#define NUM_MB_TYPE_CTX  11
#define NUM_B8_TYPE_CTX  9
#define NUM_MV_RES_CTX   10
#define NUM_REF_NO_CTX   6
#define NUM_DELTA_QP_CTX 4

typedef struct
{
  BiContextTypePtr mb_type_contexts[3];
  BiContextTypePtr b8_type_contexts[2];
  BiContextTypePtr mv_res_contexts [2];
  BiContextTypePtr ref_no_contexts [2];
  BiContextTypePtr delta_qp_inter_contexts;
  BiContextTypePtr delta_qp_intra_contexts;
} MotionInfoContexts;


#define NUM_IPR_CTX    2
#define NUM_CBP_CTX    4
#define NUM_TRANS_TYPE 9
#define NUM_LEVEL_CTX  4
#define NUM_RUN_CTX    2
#define NUM_COEFF_COUNT_CTX 6
#define NUM_TRANS_TYPE_ABT      6
#define NUM_COEFF_COUNT_CTX_ABT 8
#define NUM_RUN_CTX_ABT         5


typedef struct
{
  BiContextTypePtr ipr_contexts [9];
  BiContextTypePtr cbp_contexts [2][3];
  BiContextTypePtr level_context[4*NUM_TRANS_TYPE];
  BiContextTypePtr run_context  [2*NUM_TRANS_TYPE];
  BiContextTypePtr coeff_count_context    [NUM_TRANS_TYPE];
  BiContextTypePtr ABT_run_context        [2*NUM_TRANS_TYPE_ABT];
  BiContextTypePtr ABT_coeff_count_context[NUM_TRANS_TYPE_ABT];
} TextureInfoContexts;

//*********************** end of data type definition for CABAC *******************


/***********************************************************************
 * N e w   D a t a    t y p e s   f o r    T M L
 ***********************************************************************
 */

/*! Buffer structure for RMPNI commands */
typedef struct RMPNIbuffer_s
{
  int RMPNI;
  int Data;
  struct RMPNIbuffer_s *Next;
} RMPNIbuffer_t;


//! Syntaxelement
typedef struct syntaxelement
{
  int                 type;           //!< type of syntax element for data part.
  int                 value1;         //!< numerical value of syntax element
  int                 value2;         //!< for blocked symbols, e.g. run/level
  int                 len;            //!< length of code
  int                 inf;            //!< info part of UVLC code
  unsigned int        bitpattern;     //!< UVLC bitpattern
  int                 context;        //!< CABAC context
  int                 k;              //!< CABAC context for coeff_count,uv
  int                 golomb_grad;    //needed if type is a golomb element (ABT)
  int                 golomb_maxlevels; // if this is zero, do not use the golomb coding. (ABT)

#if TRACE
  #define             TRACESTRING_SIZE 100            //!< size of trace string
  char                tracestring[TRACESTRING_SIZE];  //!< trace string
#endif

  //!< for mapping of syntaxElement to UVLC
  void    (*mapping)(int value1, int value2, int* len_ptr, int* info_ptr);
  //!< used for CABAC: refers to actual coding method of each individual syntax element type
  void    (*writing)(struct syntaxelement *, EncodingEnvironmentPtr);

} SyntaxElement;

//! Macroblock
typedef struct macroblock
{
  int                 currSEnr;                   //!< number of current syntax element
  int                 slice_nr;
  int                 delta_qp;
  int                 qp ;
  int                 bitcounter[MAX_BITCOUNTER_MB];
  struct macroblock   *mb_available[3][3];        /*!< pointer to neighboring MBs in a 3x3 window of current MB, which is located at [1][1] \n
                                                       NULL pointer identifies neighboring MBs which are unavailable */

  // some storage of macroblock syntax elements for global access
  int                 mb_type;
  int                 mvd[2][BLOCK_MULTIPLE][BLOCK_MULTIPLE][2];          //!< indices correspond to [forw,backw][block_y][block_x][x,y]
  int                 intra_pred_modes[BLOCK_MULTIPLE*BLOCK_MULTIPLE];
  int                 coeffs_count[BLOCK_MULTIPLE][BLOCK_MULTIPLE];
  int                 cbp ;
  int                 cbp_blk ;    //!< 1 bit set for every 4x4 block with coefs (not implemented for INTRA)
  int                 b8mode[4];
  int                 b8pdir[4];
  int                 useABT[4];         // flag, indicating application of ABT for each 8x8 block in this MB (1=ABT)
  int                 abt_mode[4];       // transform mode used for ABT for each 8x8 block.
  int                 abt_pred_mode[4];  // mode used for ABT block prediction for each 8x8 block.
} Macroblock;



//! Bitstream
typedef struct
{
  int             byte_pos;           //!< current position in bitstream;
  int             bits_to_go;         //!< current bitcounter
  byte            byte_buf;           //!< current buffer for last written byte
  int             stored_byte_pos;    //!< storage for position in bitstream;
  int             stored_bits_to_go;  //!< storage for bitcounter
  int             header_len;
  byte            header_byte_buffer;
  byte            stored_byte_buf;    //!< storage for buffer of last written byte

  byte            byte_buf_skip;      //!< current buffer for last written byte
  int             byte_pos_skip;      //!< storage for position in bitstream;
  int             bits_to_go_skip;    //!< storage for bitcounter

  byte            *streamBuffer;      //!< actual buffer for written bytes
  int             write_flag;

  int             tmp_byte_pos;       //!< temp storage for position in bitstream
  int             last_startcode;     //!< location of last valid startcode

} Bitstream;

//! DataPartition
typedef struct datapartition
{

  Bitstream           *bitstream;
  Bitstream           *bitstream_frm;   //!< frame stream buffer
  Bitstream           *bitstream_fld;   //!< field stream buffer
  EncodingEnvironment ee_cabac;

  int                 (*writeSyntaxElement)(SyntaxElement *, struct datapartition *);
                      /*!< virtual function;
                           actual method depends on chosen data partition and
                           entropy coding method  */
} DataPartition;

//! Slice
typedef struct
{
  int                 picture_id;
  int                 qp;
  int                 picture_type; //!< picture type
  int                 start_mb_nr;
  int                 max_part_nr;  //!< number of different partitions
  int                 num_mb;       //!< number of MBs in the slice
  DataPartition       *partArr;     //!< array of partitions
  MotionInfoContexts  *mot_ctx;     //!< pointer to struct of context models for use in CABAC
  TextureInfoContexts *tex_ctx;     //!< pointer to struct of context models for use in CABAC

  RMPNIbuffer_t        *rmpni_buffer; //!< stores the slice temporary buffer remapping commands

  Boolean             (*slice_too_big)(int bits_slice); //!< for use of callback functions

} Slice;


// global picture format dependend buffers, mem allocation in image.c
byte   **imgY_frm;               //!< Encoded luma images
byte  ***imgUV_frm;              //!< Encoded croma images
byte   **imgY_org_frm;           //!< Reference luma image
byte  ***imgUV_org_frm;          //!< Reference croma image
int   ***tmp_mv_frm;             //!< motion vector buffer
int    **refFrArr_frm;           //!< Array for reference frames of each block

byte   **imgY;               //!< Encoded luma images
byte  ***imgUV;              //!< Encoded croma images
byte   **imgY_org;           //!< Reference luma image
byte  ***imgUV_org;          //!< Reference croma image
//int    **refFrArr;           //!< Array for reference frames of each block
byte   **imgY_pf;            //!< Post filter luma image
byte  ***imgUV_pf;           //!< Post filter croma image
byte  ***mref;               //!< 1/4 pix luma
byte ****mcef;               //!< pix chroma
int    **img4Y_tmp;          //!< for quarter pel interpolation
int   ***tmp_mv;             //!< motion vector buffer
int    **refFrArr;           //!< Array for reference frames of each block

// B pictures
// motion vector : forward, backward, direct
int  ***tmp_fwMV;
int  ***tmp_bwMV;
int  ***dfMV;
int  ***dbMV;
int   **fw_refFrArr;
int   **bw_refFrArr;
byte  **nextP_imgY;
byte ***nextP_imgUV;
pel_t **Refbuf11;            //!< 1/1th pel (full pel) reference frame buffer
int  ***colB8mode;           //!< ABT: Array containing the modes of the collocated 8x8 Blocks.

// global picture format dependend buffers, mem allocation in image.c (field picture)
byte   **imgY_org_top;
byte  ***imgUV_org_top;
byte   **imgY_org_bot;
byte  ***imgUV_org_bot;
byte   **imgY_top;               //!< Encoded luma images
byte  ***imgUV_top;              //!< Encoded croma images
byte   **imgY_bot;               //!< Encoded luma images
byte  ***imgUV_bot;              //!< Encoded croma images
byte   **imgY_com;               //!< Encoded luma images
byte  ***imgUV_com;              //!< Encoded croma images
pel_t **Refbuf11_fld;            //!< 1/1th pel (full pel) reference frame buffer
int    **refFrArr_top;           //!< Array for reference frames of each block
int    **refFrArr_bot;           //!< Array for reference frames of each block

// global picture format dependend buffers, mem allocation in image.c (field picture)
byte  ***mref_fld;               //!< 1/4 pix luma
byte ****mcef_fld;               //!< pix chroma

// global picture format dependend buffers, mem allocation in image.c (frame buffer)
byte  ***mref_frm;               //!< 1/4 pix luma
byte ****mcef_frm;               //!< pix chroma

// B pictures
// motion vector : forward, backward, direct
int   **fw_refFrArr_top;
int   **bw_refFrArr_top;
int   **fw_refFrArr_bot;
int   **bw_refFrArr_bot;
int   ***tmp_mv_top;             //!< motion vector buffer
int   ***tmp_mv_bot;             //!< motion vector buffer

// global picture format dependend buffers, mem allocation in image.c (frame buffer)
byte   **imgY_org_frm;
byte  ***imgUV_org_frm;
byte   **imgY_frm;               //!< Encoded luma images
byte  ***imgUV_frm;              //!< Encoded croma images
pel_t **Refbuf11_frm;            //!< 1/1th pel (full pel) reference frame buffer
int    **refFrArr_frm;           //!< Array for reference frames of each block
int   direct_mode;

// B pictures
// motion vector : forward, backward, direct
int   **fw_refFrArr_frm;
int   **bw_refFrArr_frm;

// Buffers for rd optimization with packet losses, Dim. Kontopodis
/* int  **resY;             //!< Residue of Luminance
byte ***decY;            //!< Decoded values at the simulated decoders
byte ****decref;         //!< Reference frames of the simulated decoders
byte ***decY_best;       //!< Decoded frames for the best mode for all decoders
byte **RefBlock;
byte **status_map;
byte **dec_mb_mode;
byte **dec_mb_ref;*/
byte **pixel_map;   //!< Shows the latest reference frame that is reliable for each pixel
byte **refresh_map; //!< Stores the new values for pixel_map  
int intras;         //!< Counts the intra updates in each frame.

int  Bframe_ctr, frame_no, nextP_tr_fld, nextP_tr_frm;
int  tot_time;

#define ET_SIZE 300      //!< size of error text buffer
char errortext[ET_SIZE]; //!< buffer for error message for exit with error()

//! Info for the "decoders-in-the-encoder" used for rdoptimization with packet losses
typedef struct
{
  int  **resY;             //!< Residue of Luminance
  byte ***decY;            //!< Decoded values at the simulated decoders
  byte ****decref;         //!< Reference frames of the simulated decoders
  byte ***decY_best;       //!< Decoded frames for the best mode for all decoders
  byte **RefBlock;
  byte **status_map;
  byte **dec_mb_mode;
} Decoders;
extern Decoders *decs;

//! SNRParameters
typedef struct
{
  float snr_y;               //!< current Y SNR
  float snr_u;               //!< current U SNR
  float snr_v;               //!< current V SNR
  float snr_y1;              //!< SNR Y(dB) first frame
  float snr_u1;              //!< SNR U(dB) first frame
  float snr_v1;              //!< SNR V(dB) first frame
  float snr_ya;              //!< Average SNR Y(dB) remaining frames
  float snr_ua;              //!< Average SNR U(dB) remaining frames
  float snr_va;              //!< Average SNR V(dB) remaining frames
} SNRParameters;

                             //! all input parameters
typedef struct
{
  int no_frames;                //!< number of frames to be encoded
  int no_fields;                //!< number of fields to be encoded
  int qp0;                      //!< QP of first frame
  int qpN;                      //!< QP of remaining frames
  int jumpd;                    //!< number of frames to skip in input sequence (e.g 2 takes frame 0,3,6,9...)
  int mv_res;                   //!< motion vector resolution: 0: 1/4-pel, 1: 1/8-pel
  int hadamard;                 /*!< 0: 'normal' SAD in 1/3 pixel search.  1: use 4x4 Haphazard transform and '
                                     Sum of absolute transform difference' in 1/3 pixel search                   */
  int search_range;             /*!< search range - integer pel search and 16x16 blocks.  The search window is
                                     generally around the predicted vector. Max vector is 2xmcrange.  For 8x8
                                     and 4x4 block sizes the search range is 1/2 of that for 16x16 blocks.       */
  int abt;                      //!< Use Adaptive Block Transforms (ABT) 0:Off, 1: Inter ABT, 2: Inter & Intra ABT
  int no_multpred;              /*!< 1: prediction from the last frame only. 2: prediction from the last or
                                     second last frame etc.  Maximum 5 frames                                    */
  int img_width;                //!< GH: if CUSTOM image format is chosen, use this size
  int img_height;               //!< GH: width and height must be a multiple of 16 pels
  int yuv_format;               //!< GH: YUV format (0=4:0:0, 1=4:2:0, 2=4:2:2, 3=4:4:4,currently only 4:2:0 is supported)
  int color_depth;              //!< GH: YUV color depth per component in bit/pel (currently only 8 bit/pel is supported)
  int intra_upd;                /*!< For error robustness. 0: no special action. 1: One GOB/frame is intra coded
                                     as regular 'update'. 2: One GOB every 2 frames is intra coded etc.
                                     In connection with this intra update, restrictions is put on motion vectors
                                     to prevent errors to propagate from the past                                */
  int blc_size[8][2];           //!< array for different block sizes
  int slice_mode;               //!< Indicate what algorithm to use for setting slices
  int slice_argument;           //!< Argument to the specified slice algorithm
  int UseConstrainedIntraPred;  //!< 0: Inter MB pixels are allowed for intra prediction 1: Not allowed
  int  infile_header;           //!< If input file has a header set this to the length of the header
  char infile[100];             //!< YUV 4:2:0 input format
  char outfile[100];            //!< H.26L compressed output bitstream
  char ReconFile[100];          //!< Reconstructed Pictures
  char TraceFile[100];          //!< Trace Outputs
  int intra_period;             //!< Random Access period though intra

  // B pictures
  int successive_Bframe;        //!< number of B frames that will be used
  int qpB;                      //!< QP of B frames

  // SP Pictures
  int sp_periodicity;           //!< The periodicity of SP-pictures
  int qpsp;                     //!< SP Picture QP for prediction error
  int qpsp_pred;                //!< SP Picture QP for predicted block

  // Introduced by TOM
  int symbol_mode;              //!< Specifies the mode the symbols are mapped on bits
  int of_mode;                  //!< Specifies the mode of the output file
  int partition_mode;           //!< Specifies the mode of data partitioning

  int SequenceHeaderType;
  int TRModulus;
  int PicIdModulus;

  int InterSearch16x16;
  int InterSearch16x8;
  int InterSearch8x16;
  int InterSearch8x8;
  int InterSearch8x4;
  int InterSearch4x8;
  int InterSearch4x4;

  char PictureTypeSequence[MAXPICTURETYPESEQUENCELEN];

#ifdef _FULL_SEARCH_RANGE_
  int full_search;
#endif
#ifdef _ADAPT_LAST_GROUP_
  int last_frame;
#endif
#ifdef _CHANGE_QP_
  int qpN2, qpB2, qp2start;
#endif
  int rdopt;
#ifdef _ADDITIONAL_REFERENCE_FRAME_
  int add_ref_frame;
#endif
#ifdef _LEAKYBUCKET_
  int NumberLeakyBuckets;
  char LeakyBucketRateFile[100];
  char LeakyBucketParamFile[100];
#endif

  int InterlaceCodingOption;
  int Encapsulated_NAL_Payload;

  int LossRateA;
  int LossRateB;
  int LossRateC;
  int NoOfDecoders;
  int RestrictRef;
  int NumFramesInELSubSeq;
  int NumFrameIn2ndIGOP;

  int RandomIntraMBRefresh;     //! Number of pseudo-random intra-MBs per picture
  int FmoNumSliceGroups;        //! Number of slice groups, 0 == scan order 
                                //! slices, maximum value is 7
  int FmoType;                  //! Type pf FMO MBAmap description, see CD doc
  char FmoConfigFileName[100];  //! Filename for config info fot type 0, 2

} InputParameters;

//! ImageParameters
typedef struct
{
  int number;                  //!< current image number to be encoded
  int pn;                      //!< picture number
  int lindex;                  //!< next long term index to be used
  int max_lindex;              //!< max long term index 
  int nb_references;
  int current_mb_nr;
  int total_number_mb;
  int current_slice_nr;
  int type;
  int pstruct;                 //!< picture structure
  int types;                   /*!< This is for SP-Pictures, since all the syntax elements for SP-Pictures
                                    are the same as P-pictures, we keep the img->type as P_IMG but indicate
                                    SP-Pictures by img->types */
  int no_multpred;             /*!< 1: prediction from the last frame only.
                                    2: prediction from the last or second last frame etc. */
  int qp;                      //!< quant for the current frame
  int qpsp;                    //!< quant for the prediction frame of SP-frame
  int framerate;
  int width;                   //!< Number of pels
  int width_cr;                //!< Number of pels chroma
  int height;                  //!< Number of lines
  int height_cr;               //!< Number of lines  chroma
  int mb_y;                    //!< current MB vertical
  int mb_x;                    //!< current MB horizontal
  int mb_x_save;               //!< horizontal position of the last written MB
  int mb_y_save;               //!< vertical position of the last written MB
  int block_y;                 //!< current block vertical
  int block_x;                 //!< current block horizontal
  int subblock_y;              //!< current subblock vertical
  int subblock_x;              //!< current subblock horizontal
  int pix_y;                   //!< current pixel vertical
  int pix_x;                   //!< current pixel horizontal
  int mb_y_upd;
  int mb_y_intra;              //!< which GOB to intra code
  int pix_c_y;                 //!< current pixel chroma vertical
  int block_c_x;               //!< current block chroma vertical
  int pix_c_x;                 //!< current pixel chroma horizontal
  int **ipredmode;             //!< GH ipredmode[90][74];prediction mode for inter frames */ /* fix from ver 4.1
  int cod_counter;             //!< Current count of number of skipped macroblocks in a row

  // some temporal buffers
  int mprr[9][16][16];         //!< all 9 prediction modes? // enlarged from 4 to 16 for ABT (is that neccessary?)

  int mprr_2[5][16][16];       //!< all 4 new intra prediction modes
  int***** mv;                 //!< motion vectors for all block types and all reference frames
  int mpr[16][16];             //!< current best prediction mode
  int m7[16][16];              //!< the diff pixel values between orginal image and prediction

  int ****cofAC;               //!< AC coefficients [8x8block][4x4block][level/run][scan_pos]
  int ***cofDC;                //!< DC coefficients [yuv][level/run][scan_pos]

  Slice       *currentSlice;                                //!< pointer to current Slice data struct
  Macroblock    *mb_data;                                   //!< array containing all MBs of a whole frame
  SyntaxElement   MB_SyntaxElements[MAX_SYMBOLS_PER_MB];    //!< temporal storage for all chosen syntax elements of one MB

  int *quad;               //!< Array containing square values,used for snr computation  */                                         /* Values are limited to 5000 for pixel differences over 70 (sqr(5000)).
  int **intra_block;

  int tr;
  int refPicID;                        //!< temporal reference for reference frames (non-B frames)
  int refPicID_fld;
  int refPicID_frm; 
  int fld_type;                        //!< top or bottom field
  int fld_flag;
  int direct_intraP_ref[4][4];
  int pstruct_next_P;
  int imgtr_next_P_frm;
  int imgtr_last_P_frm;
  int imgtr_next_P_fld;
  int imgtr_last_P_fld;

  // B pictures
  int b_interval;
  int p_interval;
  int b_frame_to_code;
  int fw_mb_mode;
  int bw_mb_mode;
  int***** p_fwMV;       //!< for MVDFW
  int***** p_bwMV;       //!< for MVDBW

  int***** all_mv;       //!< replaces local all_mv
  int***** all_bmv;      //!< replaces local all_mv

  int buf_cycle;
  int i16offset;

  int layer;             //!< which layer this picture belonged to
  int old_layer;         //!< old layer number
} ImageParameters;

                                //!< statistics
typedef struct
{
  int   quant0;                 //!< quant for the first frame
  int   quant1;                 //!< average quant for the remaining frames
  float bitr;                   //!< bit rate for current frame, used only for output til terminal
  float bitr0;                  //!< stored bit rate for the first frame
  float bitrate;                //!< average bit rate for the sequence except first frame
  int   bit_ctr;                //!< counter for bit usage
  int   bit_ctr_0;              //!< stored bit use for the first frame
  int   bit_ctr_n;              //!< bit usage for the current frame
  int   bit_slice;              //!< number of bits in current slice
  int   bit_use_mode_inter[2][MAXMODE]; //!< statistics of bit usage
  int   bit_ctr_emulationprevention; //!< stored bits needed to prevent start code emulation
  int   mode_use_intra[25];     //!< Macroblock mode usage for Intra frames
  int   mode_use_inter[2][MAXMODE];

  // B pictures
  int   *mode_use_Bframe;
  int   *bit_use_mode_Bframe;
  int   bit_ctr_P;
  int   bit_ctr_B;
  float bitrate_P;
  float bitrate_B;

  int   bit_use_stuffingBits[3];
  int   bit_use_mb_type[3];
  int   bit_use_header[3];
  int   tmp_bit_use_cbp[3];
  int   bit_use_coeffY[3];
  int   bit_use_coeffC[3];
  int   bit_use_delta_quant[3];
} StatParameters;


extern InputParameters *input;
extern ImageParameters *img;
extern StatParameters *stat;
extern StatParameters stats_frame,stats_field;

extern SNRParameters *snr;

// files
FILE *p_dec,*p_dec_u,*p_dec_v;   //!< internal decoded image for debugging
FILE *p_stat;                    //!< status file for the last encoding session
FILE *p_log;                     //!< SNR file
FILE *p_in;                      //!< YUV
FILE *p_datpart;                 //!< file to write bitlength and id of all partitions
FILE *p_trace;                   //!< Trace file


/***********************************************************************
 * P r o t o t y p e s   f o r    T M L
 ***********************************************************************
 */


void intrapred_luma(int CurrPixX,int CurrPixY);
void intrapred_luma_ABT(int pos_x,int pos_y,int subblock_size_x,int subblock_size_y);
void init();
void find_snr();
void oneforthpix();
void oneforthpix_2();
int  encode_oneIorP_Frame();
int  encode_oneB_Frame();
int  find_sad(int hadamard, int m7[16][16]);
int  dct_luma(int pos_mb1,int pos_mb2,int *cnt_nonz,int);
int  dct_luma_sp(int pos_mb1,int pos_mb2,int *cnt_nonz);
void  copyblock_sp(int pos_mb1,int pos_mb2);
int  dct_chroma(int uv,int i11);
int  dct_chroma_sp(int uv,int i11);
int  motion_search(int isi);
void levrun_linfo_c2x2(int level,int run,int *len,int *info);
void levrun_linfo_intra(int level,int run,int *len,int *info);
void levrun_linfo_inter(int level,int run,int *len,int *info);
int  sign(int a,int b);
void intrapred_chroma(int,int,int uv);
void intrapred_luma_2();
int  find_sad2(int *intra_mode);

int dct_luma2(int);

void init_img();
// void init_stat();
void report();
void information_init();
void init_frame();
void select_picture_type(SyntaxElement *symbol);
void read_one_new_frame();
void write_reconstructed_image();
void DeblockFrame(ImageParameters *img, byte **, byte ***) ;


void  LumaPrediction4x4 (int, int, int, int, int);
int   SATD (int*, int);

pel_t* FastLineX (int, pel_t*, int, int);
pel_t* UMVLineX  (int, pel_t*, int, int);

void LumaResidualCoding ();
void ChromaResidualCoding (int*);
void SetRefFrameInfo (int, int);
int  writeMBHeader   ();

extern int*   refbits;
extern int*** motion_cost;

void  Get_Direct_Motion_Vectors ();
void  PartitionMotionSearch     (int, int, double, int);               // useABT added.
int   BIDPartitionCost          (int, int, int, int, int);             // useABT added.
int   LumaResidualCoding8x8     (int*, int*, int, int, int, int, int); // useABT added.
int   writeLumaCoeff8x8         (int, int, int);                       // useABT added.
int   writeMotionVector8x8      (int, int, int, int, int, int);
int   writeReferenceFrame       (int, int, int, int);
int   writeIntra4x4Modes        (int);

int Get_Direct_Cost8x8 (int, double, int);                             // abt_mode added.
int Get_Direct_CostMB  (double);
int B8Mode2Value (int b8mode, int b8pdir);
void writeCBP_BIT_CABAC (int b8, int bit, int cbp, Macroblock* currMB, int inter, EncodingEnvironmentPtr eep_dp);

int GetSkipCostMB (double lambda);
void FindSkipModeMotionVector ();


// dynamic mem allocation
int  init_global_buffers();
void free_global_buffers();
void no_mem_exit  (char *where);

int  get_mem_mv  (int******);
void free_mem_mv (int*****);
void free_img    ();

int  get_mem_ACcoeff  (int*****);
int  get_mem_DCcoeff  (int****);
void free_mem_ACcoeff (int****);
void free_mem_DCcoeff (int***);

void put_buffer_frame();
void init_field();
void split_field_top();
void split_field_bot();
void put_buffer_top();
void put_buffer_bot();
int  decide_fld_frame(float snr_frame_Y, float snr_field_Y, int bit_field, int bit_frame, double lambda_picture);
int  get_mem4global_buffers_field();
void free_mem4global_buffers_field();
void read_one_new_field();
void combine_field();
void find_distortion();

// Added for (re-) structuring the TML soft
int   encode_one_frame();
int   encode_one_slice(int SLiceGroupId);   //! returns the number of MBs in the slice
void  malloc_slice();
void  free_slice();
void  init_slice();
void  encode_one_macroblock();
void  start_macroblock();
void  set_MB_parameters (int mb);           //! sets up img-> according to input-> and currSlice->
int   write_ipred_modes(); // used for intra ABT modes

int   writeMotionInfo2NAL ();

void  terminate_macroblock(Boolean *end_of_slice, Boolean *recode_macroblock);
int   slice_too_big(int rlc_bits);
void  write_one_macroblock();
void  proceed2nextMacroblock();
void  LumaResidualCoding_P();
void  ChromaCoding_P(int *cr_cbp);
void  SetRefFrameInfo_P();
int   MakeIntraPrediction(int *intra_pred_mode_2);
void  CheckAvailabilityOfNeighbors();
int   writeMotionInfo2NAL_Pframe();
void  writeCBPandCoeffs2NAL();
int   writeSyntaxElement_UVLC(SyntaxElement *se, DataPartition *this_dataPart);
int   writeSyntaxElement_fixed(SyntaxElement *se, DataPartition *this_dataPart);
int   writeSyntaxElement2Buf_UVLC(SyntaxElement *se, Bitstream* this_streamBuffer );
void  writeUVLC2buffer(SyntaxElement *se, Bitstream *currStream);
int   symbol2uvlc(SyntaxElement *se);
void  n_linfo(int n, int *len,int *info);
void  n_linfo2(int n, int dummy, int *len,int *info);
void  intrapred_linfo(int ipred1, int ipred2, int *len,int *info);
void  mvd_linfo2(int mvd, int dummy, int *len,int *info);
void  dquant_linfo(int mvd, int dummy, int *len,int *info);
void  cbp_linfo_intra(int cbp, int dummy, int *len,int *info);
void  cbp_linfo_inter(int cbp, int dummy, int *len,int *info);
#if TRACE
void  trace2out(SyntaxElement *se);
#endif
Boolean dummy_slice_too_big(int bits_slice);

// CABAC
void arienco_start_encoding(EncodingEnvironmentPtr eep, unsigned char *code_buffer, int *code_len );
int  arienco_bits_written(EncodingEnvironmentPtr eep);
int  get_trailing_bits(EncodingEnvironmentPtr eep);
void arienco_done_encoding(EncodingEnvironmentPtr eep);
void biari_init_context( BiContextTypePtr ctx, int ini_count_0, int ini_count_1, int max_cum_freq );
void biari_copy_context( BiContextTypePtr ctx_orig, BiContextTypePtr ctx_dest );
void biari_print_context( BiContextTypePtr ctx );
void rescale_cum_freq(BiContextTypePtr bi_ct);
void biari_encode_symbol(EncodingEnvironmentPtr eep, signed short symbol, BiContextTypePtr bi_ct );
MotionInfoContexts* create_contexts_MotionInfo(void);
TextureInfoContexts* create_contexts_TextureInfo(void);
void init_contexts_MotionInfo(MotionInfoContexts *enco_ctx, int ini_flag);
void init_contexts_TextureInfo(TextureInfoContexts *enco_ctx, int ini_flag);
void delete_contexts_MotionInfo(MotionInfoContexts *enco_ctx);
void delete_contexts_TextureInfo(TextureInfoContexts *enco_ctx);
void writeHeaderToBuffer();
void writeEOSToBuffer();
int  writeSyntaxElement_CABAC(SyntaxElement *se, DataPartition *this_dataPart);
void writeMB_typeInfo2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeB8_typeInfo2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeIntraPredMode2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeRefFrame2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeMVD2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeCBP2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeDquant_inter_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeDquant_intra_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeRunLevel2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeBiDirBlkSize2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeBiMVD2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);

void error(char *text, int code);
int  start_sequence();
int  terminate_sequence();
int  start_slice();
int  terminate_slice(int write_out);

// B pictures
int  motion_search_Bframe(int tot_intra_sad);
int  get_fwMV(int *min_fw_sad, int tot_intra_sad);
void get_bwMV(int *min_bw_sad);
void get_bid(int *bid_sad, int fw_predframe_no);
void get_dir(int *dir_sad);
void compare_sad(int tot_intra_sad, int fw_sad, int bw_sad, int bid_sad, int dir_sad, int);
void LumaResidualCoding_B();
void ChromaCoding_B(int *cr_cbp);
void SetRefFrameInfo_B();
int  writeMotionInfo2NAL_Bframe();
int  BlkSize2CodeNumber(int blc_size_h, int blc_size_v);

// Introduced for 1/8-pel
void interpolate_frame();
void interpolate_frame_to_fb();
void oneeighthpix();

void InitRefbuf ();
void InitMotionVectorSearchModule();
void InitRefbuf_fld ();
void copy2mref_fld();
void store_field_MV(int frame_number);
void store_field_colB8mode();


void  SetRefFrameInfo (int, int);

void frame_picture(int *bits_frm, float *dis_frm_y, float *dis_frm_u, float *dis_frm_v);
void field_picture(int *bits_fld, float *dis_fld_y, float *dis_fld_u, float *dis_fld_v);
void top_field_picture(int *bits_fld);
void bottom_field_picture(int *bits_fld);
void distortion_fld(float *dis_fld_y, float *dis_fld_u, float *dis_fld_v);
void picture_structure_decision(int bit_frame, int bit_field, float snr_frame, float snr_field);
void field_mode_buffer(int bit_field, float snr_field_y, float snr_field_u, float snr_field_v);
void frame_mode_buffer(int bit_frame, float snr_frame_y, float snr_frame_u, float snr_frame_v);
void set_ref_field(int *k);
void rotate_buffer();

int   writeLumaCoeff4x4     (int, int, int);
int   writeCBPandLumaCoeff  ();
int   writeChromaCoeff      ();
int   writeMB_bits_for_4x4_luma   (int, int, int);
int   writeMB_bits_for_16x16_luma ();
int   writeMB_bits_for_luma       (int);
int   writeMB_bits_for_DC_chroma  (int);
int   writeMB_bits_for_AC_chroma  (int);
int   writeMB_bits_for_CBP        ();
//!TO Hack for the Dquant-Problem
int   writeMB_bits_for_Dquant_inter     ();
int   writeMB_bits_for_Dquant_intra     ();
//! End TO

int   SingleUnifiedMotionSearch   (int, int, int**, int***, int*****, int, int*****, double);

//============= rate-distortion optimization ===================
void  clear_rdopt      ();
void  init_rdopt       ();
void  RD_Mode_Decision ();
int   RDcode_intra_subblock(int transform_mode,int pred_ipred_mode,int offset_x,int offset_y); // ABT
int   RDcode_intrablockABT();
//============= rate-distortion opt with packet losses ===========
void decode_one_macroblock();
void decode_one_mb (int, Macroblock*);
void decode_one_b8block (int, int, int, int, int);
void Get_Reference_Block(byte **imY, int block_y, int block_x, int mvhor, int mvver, byte **out);
byte Get_Reference_Pixel(byte **imY, int y, int x);
int Half_Upsample(byte **imY, int j, int i);
void DecOneForthPix(byte **dY, byte ***dref);
void compute_residue(int mode);
void compute_residue_b8block (int, int);
void compute_residue_mb (int);
void UpdateDecoders();
void Build_Status_Map(byte **s_map);
void Error_Concealment(byte **inY, byte **s_map, byte ***refY);
void Conceal_Error(byte **inY, int mb_y, int mb_x, byte ***refY, byte **s_map);
//============= restriction of reference frames based on the latest intra-refreshes==========
int CheckReliabilityOfRefFrame(int ref_frame, int);
void UpdatePixelMap();

//============= fast full integer search =======================
#ifdef _FAST_FULL_ME_
void  ClearFastFullIntegerSearch    ();
void  ResetFastFullIntegerSearch    ();
#endif

void process_2nd_IGOP();
void SetImgType();

// Tian Dong: for IGOPs
extern Boolean In2ndIGOP;
extern int start_frame_no_in_this_IGOP;
extern int start_tr_in_this_IGOP;
extern int FirstFrameIn2ndIGOP;
#define IMG_NUMBER (img->number-start_frame_no_in_this_IGOP)
#define PAYLOAD_TYPE_IDERP 8
byte *NAL_Payload_buffer;
void SODBtoRBSP(Bitstream *currStream);
int RBSPtoEBSP(byte *streamBuffer, int begin_bytepos, int end_bytepos);
int Bytes_After_Header;
int startcodeprefix_len;

#endif
