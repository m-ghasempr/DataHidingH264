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
#include "nalucommon.h"
#include "parsetcommon.h"

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
  PAR_OF_ANNEXB,    //!< Annex B bytestream format
  PAR_OF_RTP,       //!< RTP packets in outfile
//  PAR_OF_IFF        //!< Interim File Format
} PAR_OF_TYPE;

//! Boolean Type
/*typedef enum {
  FALSE,
  TRUE
} Boolean;
*/
typedef enum {
  FRAME_CODING,
  ADAPTIVE_CODING,
  FIELD_CODING,
  MB_CODING
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


typedef enum {
  FRAME,
  TOP_FIELD,
  BOTTOM_FIELD
} PictureType;           //!< New enum for field processing


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
//  int           *Ecodestrm_laststartcode;
  // storage in case of recode MB
  unsigned int  ElowS, ErangeS;
  unsigned int  EbufferS;
  unsigned int  Ebits_to_goS;
  unsigned int  Ebits_to_followS;
  byte          *EcodestrmS;
  int           *Ecodestrm_lenS;
  int           C, CS;
  int           E, ES;
  int           B, BS;
} EncodingEnvironment;

typedef EncodingEnvironment *EncodingEnvironmentPtr;

//! struct for context management
typedef struct
{
  unsigned short state;         // index into state-table CP  
  unsigned char  MPS;           // Least Probable Symbol 0/1 CP

  unsigned long  count;

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
#define NUM_MB_AFF_CTX 4

typedef struct
{
  BiContextType mb_type_contexts [3][NUM_MB_TYPE_CTX];
  BiContextType b8_type_contexts [2][NUM_B8_TYPE_CTX];
  BiContextType mv_res_contexts  [2][NUM_MV_RES_CTX];
  BiContextType ref_no_contexts  [2][NUM_REF_NO_CTX];
  BiContextType delta_qp_contexts   [NUM_DELTA_QP_CTX];
  BiContextType mb_aff_contexts     [NUM_MB_AFF_CTX];

} MotionInfoContexts;


#define NUM_IPR_CTX    2
#define NUM_CIPR_CTX   4
#define NUM_CBP_CTX    4
#define NUM_BCBP_CTX   4
#define NUM_MAP_CTX   15
#define NUM_LAST_CTX  15
#define NUM_ONE_CTX    5
#define NUM_ABS_CTX    5


typedef struct
{
  BiContextType  ipr_contexts [NUM_IPR_CTX]; 
  BiContextType  cipr_contexts[NUM_CIPR_CTX]; 
  BiContextType  cbp_contexts [3][NUM_CBP_CTX];
  BiContextType  bcbp_contexts[NUM_BLOCK_TYPES][NUM_BCBP_CTX];
  BiContextType  map_contexts [NUM_BLOCK_TYPES][NUM_MAP_CTX];
  BiContextType  last_contexts[NUM_BLOCK_TYPES][NUM_LAST_CTX];
  BiContextType  one_contexts [NUM_BLOCK_TYPES][NUM_ONE_CTX];
  BiContextType  abs_contexts [NUM_BLOCK_TYPES][NUM_ABS_CTX];
  BiContextType  fld_map_contexts [NUM_BLOCK_TYPES][NUM_MAP_CTX];
  BiContextType  fld_last_contexts[NUM_BLOCK_TYPES][NUM_LAST_CTX];
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

typedef struct MMCObuffer_s
{
  int MMCO;
  int DPN;
  int LPIN;
  int MLIP1;
  struct MMCObuffer_s *Next;
} MMCObuffer_t;

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
  struct macroblock   *field_available[2];
  // some storage of macroblock syntax elements for global access
  int                 mb_type;
  int                 mvd[2][BLOCK_MULTIPLE][BLOCK_MULTIPLE][2];          //!< indices correspond to [forw,backw][block_y][block_x][x,y]
  int                 intra_pred_modes[BLOCK_MULTIPLE*BLOCK_MULTIPLE];
  int                 cbp ;
  int                 cbp_blk ;    //!< 1 bit set for every 4x4 block with coefs (not implemented for INTRA)
  int                 b8mode[4];
  int                 b8pdir[4];
  unsigned long       cbp_bits;

  int                 lf_disable;
  int                 lf_alpha_c0_offset;
  int                 lf_beta_offset;

  int                 c_ipred_mode;      //!< chroma intra prediction mode
  int                 IntraChromaPredModeFlag;
        int                                                                     mb_field;
} Macroblock;



//! Bitstream
typedef struct
{
  int             byte_pos;           //!< current position in bitstream;
  int             bits_to_go;         //!< current bitcounter
  byte            byte_buf;           //!< current buffer for last written byte
  int             stored_byte_pos;    //!< storage for position in bitstream;
  int             stored_bits_to_go;  //!< storage for bitcounter
  byte            stored_byte_buf;    //!< storage for buffer of last written byte

  byte            byte_buf_skip;      //!< current buffer for last written byte
  int             byte_pos_skip;      //!< storage for position in bitstream;
  int             bits_to_go_skip;    //!< storage for bitcounter

  byte            *streamBuffer;      //!< actual buffer for written bytes
  int             write_flag;         //!< Bitstream contains data and needs to be written

//  int             tmp_byte_pos;       //!< temp storage for position in bitstream

//  int             last_startcode;     //!< location of last valid startcode

} Bitstream;

//! DataPartition
typedef struct datapartition
{

  Bitstream           *bitstream;
//  Bitstream           *bitstream_frm;   //!< frame stream buffer
//  Bitstream           *bitstream_topfld;   //!< field stream buffer
//  Bitstream           *bitstream_botfld;   //!< field stream buffer
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

  // !KS: RMPNI buffer should be retired. just do some sore simple stuff
  RMPNIbuffer_t        *rmpni_buffer; //!< stores the slice temporary buffer remapping commands

  int                 ref_pic_list_reordering_flag_l0;
  int                 *remapping_of_pic_nums_idc_l0;
  int                 *abs_diff_pic_num_minus1_l0;
  int                 *long_term_pic_idx_l0;
  int                 ref_pic_list_reordering_flag_l1;
  int                 *remapping_of_pic_nums_idc_l1;
  int                 *abs_diff_pic_num_minus1_l1;
  int                 *long_term_pic_idx_l1;

  Boolean             (*slice_too_big)(int bits_slice); //!< for use of callback functions

  int                 write_is_real;//GB
  int                 field_ctx[3][2]; //GB

} Slice;


#define MAXSLICEPERPICTURE 100
typedef struct 
{
  int   no_slices;
  Slice *slices[MAXSLICEPERPICTURE];
  int bits_per_picture;
  float distortion_y;
  float distortion_u;
  float distortion_v;
} Picture;

Picture *top_pic;
Picture *bottom_pic;
Picture *frame_pic;


unsigned int toprefpoc[MAX_NO_POC_FRAMES];
unsigned int bottomrefpoc[MAX_NO_POC_FRAMES];


typedef struct
{
  // Size info
  int x_size, y_framesize, y_fieldsize;  
  char *yf, *uf, *vf;                    //!< frame representation
  char *yt, *ut, *vt;                    //!< top field
  char *yb, *ub, *vb;                    //!< bottom field
} Sourceframe;

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
byte  ***mref_w;               //!< 1/4 pix luma for weighted prediction ME
byte ****mcef;               //!< pix chroma
int    **img4Y_tmp;          //!< for quarter pel interpolation
int   ***tmp_mv;             //!< motion vector buffer

int   **moving_block;           //!< stationary block buffer
int   **moving_block_frm;       //!< stationary block buffer - frame
int   **moving_block_top;       //!< stationary block buffer - field
int   **moving_block_bot;       //!< stationary block buffer - field

int    **refFrArr;           //!< Array for reference frames of each block

// B pictures
// motion vector : forward, backward, direct
int  ***tmp_fwMV;
int  ***tmp_bwMV;
int  ***tmp_fwMV_top;   //!< For MB level field/frame coding tools
int  ***tmp_fwMV_bot;   //!< For MB level field/frame coding tools
int  ***tmp_bwMV_top;   //!< For MB level field/frame coding tools
int  ***tmp_bwMV_bot;   //!< For MB level field/frame coding tools
int  **field_mb;      //!< For MB level field/frame coding tools
int  mb_adaptive;     //!< For MB level field/frame coding tools
int  MBPairIsField;     //!< For MB level field/frame coding tools
int  WriteFrameFieldMBInHeader; //! For MB level field/frame coding tools

int  ***dfMV;
int  ***dbMV;
int   **fw_refFrArr;
int   **bw_refFrArr;
byte  **nextP_imgY;
byte ***nextP_imgUV;
pel_t **Refbuf11;            //!< 1/1th pel (full pel) reference frame buffer
pel_t **Refbuf11_w;           // for weighted reference frame buffer

//Weighted prediction
int ***wp_weight;  // weight in [list][index][component] order
int ***wp_offset;  // offset in [list][index][component] order
int ****wbp_weight;  // weight in [list][fwd_index][bwd_idx][component] order
int **weight;      // weight in [refbuf][component] order, because JM currently uses one list for fwd and bwd ref pics
int **offset;      // offset in [refbuf][component] order, because JM currently uses one list for fwd and bwd ref pics
int luma_log_weight_denom;
int chroma_log_weight_denom;
int wp_luma_round;
int wp_chroma_round;

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
pel_t **Refbuf11_fld_w;            //!< 1/1th pel (full pel) reference frame buffer
int    **refFrArr_top;           //!< Array for reference frames of each block
int    **refFrArr_bot;           //!< Array for reference frames of each block
// int    **refFrArr_top_save;      //!< For MB level field/frame coding tools
// int    **refFrArr_bot_save;      //!< For MB level field/frame coding tools


// global picture format dependend buffers, mem allocation in image.c (field picture)
byte  ***mref_fld;               //!< 1/4 pix luma
byte  ***mref_fld_w;               //!< 1/4 pix luma for wp
byte ****mcef_fld;               //!< pix chroma
byte ***mref_mbfld;        //!< For MB level field/frame coding tools 
byte ***mref_mbfld_w;        //!< For MB level field/frame coding tools for wp


// global picture format dependend buffers, mem allocation in image.c (frame buffer)
byte  ***mref_frm;               //!< 1/4 pix luma
byte  ***mref_frm_w;               //!< 1/4 pix luma
byte ****mcef_frm;               //!< pix chroma

// B pictures
// motion vector : forward, backward, direct
int   **fw_refFrArr_top;
int   **bw_refFrArr_top;
int   **fw_refFrArr_bot;
int   **bw_refFrArr_bot;
int   ***tmp_mv_top;             //!< motion vector buffer
int   ***tmp_mv_bot;             //!< motion vector buffer
// int   ***tmp_mv_top_save;    //!< motion vector buffer for MB level field/frame coding tools
// int   ***tmp_mv_bot_save;    //!< motion vector buffer for MB level field/frame coding tools

int   **fwdir_refFrArr;         //!< direct mode forward reference buffer
int   **bwdir_refFrArr;         //!< direct mode backward reference buffer



// global picture format dependend buffers, mem allocation in image.c (frame buffer)
byte   **imgY_org_frm;
byte  ***imgUV_org_frm;
byte   **imgY_frm;               //!< Encoded luma images
byte  ***imgUV_frm;              //!< Encoded croma images
pel_t **Refbuf11_frm;            //!< 1/1th pel (full pel) reference frame buffer
pel_t **Refbuf11_frm_w;            //!< 1/1th pel (full pel) reference frame buffer
int    **refFrArr_frm;           //!< Array for reference frames of each block
int   direct_mode;
int   TopFieldIsSkipped, TopFrameIsSkipped; //!< Flag to prevent bottom MB in a field MB pair to be skipped 
                //!< if the top MB is skipped

// B pictures
// motion vector : forward, backward, direct
int   **fw_refFrArr_frm;
int   **bw_refFrArr_frm;

// Buffers for rd optimization with packet losses, Dim. Kontopodis
byte **pixel_map;   //!< Shows the latest reference frame that is reliable for each pixel
byte **refresh_map; //!< Stores the new values for pixel_map  
int intras;         //!< Counts the intra updates in each frame.

int  Bframe_ctr, frame_no, nextP_tr_fld, nextP_tr_frm;
int  tot_time;

#define ET_SIZE 300      //!< size of error text buffer
char errortext[ET_SIZE]; //!< buffer for error message for exit with error()

int    **abp_type_FrArr;      //!< Array for abp_type of each block
int    **abp_type_FrArr_top;
int    **abp_type_FrArr_bot;

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
  int hadamard;                 /*!< 0: 'normal' SAD in 1/3 pixel search.  1: use 4x4 Haphazard transform and '
                                     Sum of absolute transform difference' in 1/3 pixel search                   */
  int search_range;             /*!< search range - integer pel search and 16x16 blocks.  The search window is
                                     generally around the predicted vector. Max vector is 2xmcrange.  For 8x8
                                     and 4x4 block sizes the search range is 1/2 of that for 16x16 blocks.       */
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
  int direct_type;              //!< Direct Mode type to be used (1: Temporal, 0: Spatial)

  // SP Pictures
  int sp_periodicity;           //!< The periodicity of SP-pictures
  int qpsp;                     //!< SP Picture QP for prediction error
  int qpsp_pred;                //!< SP Picture QP for predicted block

  int WeightedPrediction;        //!< Weighted prediciton for P frames (0: not used, 1: explicit)
  int WeightedBiprediction;      //!< Weighted prediciton for B frames (0: not used, 1: explicit, 2: implicit)
  int StoredBPictures;           //!< Stored (Reference) B pictures replace P pictures (0: not used, 1: used)
  int explicit_B_prediction;     //!< explicite B prediction  // Jill get rid of eventually
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

  int LossRateA;              //! assumed loss probablility of partition A (or full slice), in per cent, used for loss-aware R/D optimization
  int LossRateB;              //! assumed loss probablility of partition B, in per cent, used for loss-aware R/D 
  int LossRateC;              //! assumed loss probablility of partition C, in per cent, used for loss-aware R/D 
  int NoOfDecoders;
  int RestrictRef;
  int NumFramesInELSubSeq;
  int NumFrameIn2ndIGOP;

  int RandomIntraMBRefresh;     //! Number of pseudo-random intra-MBs per picture
  int FmoType;                  //! Type pf FMO MBAmap description, see CD doc
  char FmoConfigFileName[100];  //! Filename for config info fot type 0, 2

  int LFSendParameters;
  int LFDisableIdc;
  int LFAlphaC0Offset;
  int LFBetaOffset;

  int SparePictureOption;
  int SPDetectionThreshold;
  int SPPercentageThreshold;

  // JVT-D095, JVT-D097
  int num_slice_groups_minus1; //! "FmoNumSliceGroups" in encoder.cfg, same as FmoNumSliceGroups, which should be erased later
  int mb_allocation_map_type; //! "FmoType" in encoder.cfg, same as FmoType, FmoType should be erased latter
  int top_left_mb; //! "FmoTopLeftMB"  in encoder.cfg
  int bottom_right_mb; //! "FmoBottomRightMB" in encoder.cfg
  int slice_group_change_direction; //! "FmoChangeDirection" in encoder.cfg
  int slice_group_change_rate_minus1; //! "FmoChangeRate" in encoder.cfg
  // End JVT-D095, JVT-D097

  int redundant_slice_flag; //! whether redundant slices exist,  JVT-D101
  int pic_order_cnt_type;   // POC200301

  int context_init_method;
  int model_number;

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
  int is_intra_block;
  int is_v_block;
  int pix_y;                   //!< current pixel vertical
  int pix_x;                   //!< current pixel horizontal
  int mb_y_upd;
  int mb_y_intra;              //!< which GOB to intra code
  int pix_c_y;                 //!< current pixel chroma vertical
  int block_c_x;               //!< current block chroma vertical
  int pix_c_x;                 //!< current pixel chroma horizontal
  int **ipredmode;             //!< GH ipredmode[90][74];prediction mode for inter frames */ /* fix from ver 4.1
  int cod_counter;             //!< Current count of number of skipped macroblocks in a row
  int ****nz_coeff;            //!< number of coefficients per block (CAVLC)


  // some temporal buffers
  int mprr[9][16][16];         //!< all 9 prediction modes? // enlarged from 4 to 16 for ABT (is that neccessary?)

  int mprr_2[5][16][16];       //!< all 4 new intra prediction modes
  int mprr_c[2][4][8][8];      //!< new chroma 8x8 intra prediction modes
  int***** mv;                 //!< motion vectors for all block types and all reference frames
  int mpr[16][16];             //!< current best prediction mode
  int m7[16][16];              //!< the diff pixel values between orginal image and prediction

  int ****cofAC;               //!< AC coefficients [8x8block][4x4block][level/run][scan_pos]
  int ***cofDC;                //!< DC coefficients [yuv][level/run][scan_pos]

  Picture     *currentPicture; //!< The coded picture currently in the works (typically frame_pic, top_pic, or bottom_pic)
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
  unsigned int fld_flag;                                
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

  int***** abp_all_dmv;        //!< replaces local all_dmv for forward interpolative prediction
  int***** abp_all_dmv_top;
  int***** abp_all_dmv_bot;
  int num_ref_pic_active_fwd_minus1;
  int num_ref_pic_active_bwd_minus1;


  int field_mb_y;   // Macroblock number of a field MB
  int field_block_y;  // Vertical block number for the first block of a field MB
  int field_pix_y;    // Co-ordinates of current macroblock in terms of field pixels (luma)
  int field_pix_c_y;  // Co-ordinates of current macroblock in terms of field pixels (chroma)
  int *****mv_top;    //!< For MB level field/frame coding tools
  int *****mv_bot;    //!< For MB level field/frame coding tools
  int *****p_fwMV_top;    //!< For MB level field/frame coding tools
  int *****p_fwMV_bot;    //!< For MB level field/frame coding tools
  int *****p_bwMV_top;    //!< For MB level field/frame coding tools
  int *****p_bwMV_bot;    //!< For MB level field/frame coding tools
  int *****all_mv_top;    //!< For MB level field/frame coding tools
  int *****all_mv_bot;    //!< For MB level field/frame coding tools
  int *****all_bmv_top;   //!< For MB level field/frame coding tools
  int *****all_bmv_bot;   //!< For MB level field/frame coding tools
  int **ipredmode_top;    //!< For MB level field/frame coding tools
  int **ipredmode_bot;    //!< For MB level field/frame coding tools
//  int update_stats;     //!< For MB level field/frame -- update stats flag
  int field_mode;     //!< For MB level field/frame -- field mode on flag
  int top_field;      //!< For MB level field/frame -- top field flag

  int buf_cycle;
  int i16offset;

  int layer;             //!< which layer this picture belonged to
  int old_layer;         //!< old layer number
  int NoResidueDirect;

  int redundant_pic_cnt; // JVT-D101

  int **field_anchor;

  //the following should probably go in sequence parameters
  // unsigned int log2_max_frame_num_minus4;
  unsigned int pic_order_cnt_type;
  // for poc mode 0, POC200301
  // unsigned int log2_max_pic_order_cnt_lsb_minus4;  
  // for poc mode 1, POC200301
  unsigned int delta_pic_order_always_zero_flag;
           int offset_for_non_ref_pic;
           int offset_for_top_to_bottom_field;
  unsigned int num_ref_frames_in_pic_order_cnt_cycle;
           int offset_for_ref_frame[1];  // MAX_LENGTH_POC_CYCLE in decoder

  // POC200301
  //the following is for slice header syntax elements of poc
  // for poc mode 0.
  unsigned int pic_order_cnt_lsb;
           int delta_pic_order_cnt_bottom;
  // for poc mode 1.
           int delta_pic_order_cnt[2];




  // POC200301
  unsigned int field_picture;
    signed int toppoc;      //poc for this frame or field
    signed int bottompoc;   //for completeness - poc of bottom field of a frame (always = poc+1)
  unsigned int frame_num;   //frame_num for this frame
  

  //the following should probably go in picture parameters
  unsigned int pic_order_present_flag; // ????????

  //the following are sent in the slice header
//  int delta_pic_order_cnt[2];
  int idr_flag;
  int disposable_flag;
  
  int adaptive_ref_pic_buffering_flag;
  int no_output_of_prior_pics_flag;
  int long_term_reference_flag;

  MMCObuffer_t *mmco_buffer;

  int model_number;

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

#define NUM_PIC_TYPE 5
  int   bit_use_stuffingBits[NUM_PIC_TYPE];
  int   bit_use_mb_type[NUM_PIC_TYPE];
  int   bit_use_header[NUM_PIC_TYPE];
  int   tmp_bit_use_cbp[NUM_PIC_TYPE];
  int   bit_use_coeffY[NUM_PIC_TYPE];
  int   bit_use_coeffC[NUM_PIC_TYPE];
  int   bit_use_delta_quant[NUM_PIC_TYPE];

  int   em_prev_bits_frm;
  int   em_prev_bits_fld;
  int  *em_prev_bits;
  int   bit_ctr_parametersets;
} StatParameters;

//!< For MB level field/frame coding tools
//!< temporary structure to store MB data for field/frame coding
typedef struct
{
  double min_rdcost;
  int tmp_mv[2][4][4];        // to hold the motion vectors for each block
  int tmp_fwMV[2][4][4];      // to hold forward motion vectors for B MB's
  int tmp_bwMV[2][4][4];      // to hold backward motion vectors for B MBs
  int dfMV[2][4][4], dbMV[2][4][4]; // to hold direct motion vectors for B MB's
  int    rec_mbY[16][16];       // hold the Y component of reconstructed MB
  int    rec_mbU[8][8], rec_mbV[8][8]; 
  int    ****cofAC;
  int    ***cofDC;
  int    mb_type;
  int    b8mode[4], b8pdir[4];
  int    frefar[4][4], brefar[4][4];
  int    **ipredmode;
  int    intra_pred_modes[16];
  int    cbp, cbp_blk;
  int    mode;
  int    *****mv, *****p_fwMV, *****p_bwMV ;
  int    *****all_mv;
  int    *****all_bmv;
  int    i16offset;
  int    c_ipred_mode;
} RD_DATA;

RD_DATA *rdopt; 
RD_DATA rddata_top_frame_mb, rddata_bot_frame_mb; //!< For MB level field/frame coding tools
RD_DATA rddata_top_field_mb, rddata_bot_field_mb; //!< For MB level field/frame coding tools

extern InputParameters *input;
extern ImageParameters *img;
extern StatParameters *stat;

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
void init();
int  find_sad(int hadamard, int m7[16][16]);
int  dct_luma(int pos_mb1,int pos_mb2,int *cnt_nonz,int);
int  dct_luma_sp(int pos_mb1,int pos_mb2,int *cnt_nonz);
void  copyblock_sp(int pos_mb1,int pos_mb2);
int  dct_chroma(int uv,int i11);
int  dct_chroma_sp(int uv,int i11);
int  motion_search(int isi);
int  sign(int a,int b);
void intrapred_chroma(int,int,int uv);
void intrapred_luma_16x16();
int  find_sad_16x16(int *intra_mode);

int dct_luma_16x16(int);
void copy_mref();

void init_poc();
void push_poc(unsigned int topvalue, unsigned int bottomvalue, unsigned int ref_frame_ind );

void init_img();
void report();
void information_init();
int  get_picture_type();
void DeblockFrame(ImageParameters *img, byte **, byte ***) ;


void  LumaPrediction4x4 (int, int, int, int, int, int);
int   SATD (int*, int);

pel_t* FastLineX (int, pel_t*, int, int);
pel_t* UMVLineX  (int, pel_t*, int, int);

void LumaResidualCoding ();
void ChromaResidualCoding (int*);
void IntraChromaPrediction8x8 (int*, int*, int*);
void SetRefFrameInfo (int, int);
int  writeMBHeader   (int rdopt); 

extern int*   refbits;
extern int*** motion_cost;

void  Get_Direct_Motion_Vectors ();
void  PartitionMotionSearch     (int, int, double);
int   BIDPartitionCost          (int, int, int, int);
int   LumaResidualCoding8x8     (int*, int*, int, int, int, int, int);
int   writeLumaCoeff8x8         (int, int);
int   writeMotionVector8x8      (int, int, int, int, int, int, int, int);
int   writeReferenceFrame       (int, int, int, int, int);
int   ABIDPartitionCost         (int, int, int*, int*, int, int*);
int   BBIDPartitionCost         (int, int, int* , int* , int, int* , int);
int   writeAbpCoeffIndex        (int, int, int, int);
int   writeIntra4x4Modes        (int);
int   writeChromaIntraPredMode  ();

int Get_Direct_Cost8x8 (int, double);
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

void split_field_top();
void split_field_bot();
int  decide_fld_frame(float snr_frame_Y, float snr_field_Y, int bit_field, int bit_frame, double lambda_picture);
int  get_mem4global_buffers_field();
void free_mem4global_buffers_field();
// void read_one_new_field();
void combine_field();

// Added for (re-) structuring the TML soft
Picture *malloc_picture();
void free_picture (Picture *pic);
int   encode_one_slice(int SLiceGroupId, Picture *pic);   //! returns the number of MBs in the slice

void  encode_one_macroblock();
void  start_macroblock();
void  set_MB_parameters (int mb);           //! sets up img-> according to input-> and currSlice->

int   writeMotionInfo2NAL ();

void  terminate_macroblock(Boolean *end_of_slice, Boolean *recode_macroblock);
int   slice_too_big(int rlc_bits);
void  write_one_macroblock(int eos_bit);
void  proceed2nextMacroblock();
void  proceed2nextSuperMacroblock();    //!< For MB level field/frame coding tools
void  back2topmb();             //!< For MB level field/frame coding tools
void proceed2lowermb();           //!< For MB level field/frame coding tools

void  LumaResidualCoding_P();
void  ChromaCoding_P(int *cr_cbp);
void  SetRefFrameInfo_P();
int   MakeIntraPrediction(int *intra_pred_mode_2);
void  CheckAvailabilityOfNeighbors();

void free_slice_list(Picture *currPic);

#if TRACE
void  trace2out(SyntaxElement *se);
#endif


// CABAC
void arienco_start_encoding(EncodingEnvironmentPtr eep, unsigned char *code_buffer, int *code_len, /* int *last_startcode, */int slice_type);
int  arienco_bits_written(EncodingEnvironmentPtr eep);
void arienco_done_encoding(EncodingEnvironmentPtr eep);
void biari_init_context (BiContextTypePtr ctx, const int* ini);
void rescale_cum_freq(BiContextTypePtr bi_ct);
void biari_encode_symbol(EncodingEnvironmentPtr eep, signed short symbol, BiContextTypePtr bi_ct );
void biari_encode_symbol_eq_prob(EncodingEnvironmentPtr eep, signed short symbol);
void biari_encode_symbol_final(EncodingEnvironmentPtr eep, signed short symbol);
MotionInfoContexts* create_contexts_MotionInfo(void);
TextureInfoContexts* create_contexts_TextureInfo(void);
void init_contexts_MotionInfo (MotionInfoContexts  *enco_ctx);
void init_contexts_TextureInfo(TextureInfoContexts *enco_ctx);
void delete_contexts_MotionInfo(MotionInfoContexts *enco_ctx);
void delete_contexts_TextureInfo(TextureInfoContexts *enco_ctx);
void writeHeaderToBuffer();
int  writeSyntaxElement_CABAC(SyntaxElement *se, DataPartition *this_dataPart);
void writeMB_typeInfo2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeB8_typeInfo2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeIntraPredMode2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeRefFrame2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeBwdRefFrame2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeMVD2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeCBP2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeDquant_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeRunLevel2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeBiDirBlkSize2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeBiMVD2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeCIPredMode2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void print_ctx_TextureInfo(TextureInfoContexts *enco_ctx);
void writeMB_skip_flagInfo2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp);
void writeFieldModeInfo2Buffer_CABAC(SyntaxElement *se, EncodingEnvironmentPtr eep_dp); //GB
void CheckAvailabilityOfNeighborsForAff(); //


void error(char *text, int code);
int  start_sequence();
int  terminate_sequence();
int  start_slice();
int  terminate_slice();

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

void InitRefbuf ();
void InitMotionVectorSearchModule();
void InitRefbuf_fld ();
void copy2mref_fld();


void  SetRefFrameInfo (int, int);

void set_ref_field(int *k);
void set_mbaff_parameters();  // For MB AFF
void writeVlcByteAlign(Bitstream* currStream);


int   writeLumaCoeff4x4_CABAC     (int, int, int);
int   writeCBPandLumaCoeff  ();
int   writeChromaCoeff      ();
int   writeMB_bits_for_4x4_luma   (int, int, int);
int   writeMB_bits_for_16x16_luma ();
int   writeMB_bits_for_luma       (int);
int   writeMB_bits_for_DC_chroma  (int);
int   writeMB_bits_for_AC_chroma  (int);
int   writeMB_bits_for_CBP        ();

int   SingleUnifiedMotionSearch   (int, int, int**, int***, int*****, int, int*****, double);

//============= rate-distortion optimization ===================
void  clear_rdopt      ();
void  init_rdopt       ();
void  RD_Mode_Decision ();
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

void AllocNalPayloadBuffer();
void FreeNalPayloadBuffer();
void SODBtoRBSP(Bitstream *currStream);
int RBSPtoEBSP(byte *streamBuffer, int begin_bytepos, int end_bytepos, int min_num_bytes);
int Bytes_After_Header;

// JVT-D101: the bit for redundant_pic_cnt in slice header may be changed, 
// therefore the bit position in the bitstream must be stored.
int rpc_bytes_to_go;
int rpc_bits_to_go;
void modify_redundant_pic_cnt(unsigned char *streamBuffer);
// End JVT-D101

int poc_distance( int refa, int refb);

#endif

#include "context_ini.h"
