/*!
 ************************************************************************
 *  \file
 *     global.h
 *
 *  \brief
 *     global definitions for H.264 encoder.
 *
 *  \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *
 ************************************************************************
 */
#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include "win32.h"
#include "defines.h"
#include "typedefs.h"
#include "parsetcommon.h"
#include "ifunctions.h"
#include "frame.h"
#include "io_video.h"
#include "io_image.h"
#include "nalucommon.h"
#include "params.h"
#include "distortion.h"
#include "lagrangian.h"
#include "quant_params.h"
#include "rc_types.h"
#include "pred_struct_types.h"



/***********************************************************************
 * D a t a    t y p e s   f o r  C A B A C
 ***********************************************************************
 */
typedef struct inp_par_enc InputParameters;

// structures that will be declared somewhere else
struct storable_picture;
struct coding_state;


typedef struct image_structure
{  
  FrameFormat format;      //!< ImageStructure format Information
  imgpel **data[3];        //!< ImageStructure pixel data
} ImageStructure;

//! block 8x8 temporary RD info
typedef struct info_8x8
{
   char   mode;
   char   pdir;
   char   ref[2];
   char   bipred;
} Info8x8;

typedef struct search_window {
  int    min_x;
  int    max_x;
  int    min_y;
  int    max_y;
} SearchWindow;

//! bit counter for a macroblock. Note that it seems safe to change all to unsigned short for 16x16 MBs. May be an issue if MB > 32x32
typedef struct bit_counter
{
  int mb_total;
  unsigned short mb_mode;
  unsigned short mb_inter;
  unsigned short mb_cbp;
  unsigned short mb_delta_quant;
  int mb_y_coeff;
  int mb_uv_coeff;
  int mb_cb_coeff;
  int mb_cr_coeff;  
  int mb_stuffing;
} BitCounter;

//! struct to characterize the state of the arithmetic coding engine
typedef struct
{
  struct video_par *p_Vid;
  unsigned int  Elow, Erange;
  unsigned int  Ebuffer;
  unsigned int  Ebits_to_go;
  unsigned int  Echunks_outstanding;
  int           Epbuf;
  byte          *Ecodestrm;
  int           *Ecodestrm_len;
  int           C;
  int           E;
} EncodingEnvironment;

typedef EncodingEnvironment *EncodingEnvironmentPtr;

//! struct for context management
typedef struct
{
  unsigned long  count;
  uint16 state;         // index into state-table CP
  unsigned char  MPS;           // Least Probable Symbol 0/1 CP  
} BiContextType;

typedef BiContextType *BiContextTypePtr;


/**********************************************************************
 * C O N T E X T S   F O R   T M L   S Y N T A X   E L E M E N T S
 **********************************************************************
 */

typedef struct
{
  BiContextType mb_type_contexts [3][NUM_MB_TYPE_CTX];
  BiContextType b8_type_contexts [2][NUM_B8_TYPE_CTX];
  BiContextType mv_res_contexts  [2][NUM_MV_RES_CTX];
  BiContextType ref_no_contexts  [2][NUM_REF_NO_CTX];
} MotionInfoContexts;

typedef struct
{
  BiContextType  transform_size_contexts   [NUM_TRANSFORM_SIZE_CTX];
  BiContextType  ipr_contexts [NUM_IPR_CTX];
  BiContextType  cipr_contexts[NUM_CIPR_CTX];
  BiContextType  cbp_contexts [3][NUM_CBP_CTX];
  BiContextType  bcbp_contexts[NUM_BLOCK_TYPES][NUM_BCBP_CTX];
  BiContextType  delta_qp_contexts   [NUM_DELTA_QP_CTX];
  BiContextType  one_contexts [NUM_BLOCK_TYPES][NUM_ONE_CTX];
  BiContextType  abs_contexts [NUM_BLOCK_TYPES][NUM_ABS_CTX];
#if ENABLE_FIELD_CTX
  BiContextType  mb_aff_contexts [NUM_MB_AFF_CTX];
  BiContextType  map_contexts [2][NUM_BLOCK_TYPES][NUM_MAP_CTX];
  BiContextType  last_contexts[2][NUM_BLOCK_TYPES][NUM_LAST_CTX];
#else
  BiContextType  map_contexts [1][NUM_BLOCK_TYPES][NUM_MAP_CTX];
  BiContextType  last_contexts[1][NUM_BLOCK_TYPES][NUM_LAST_CTX];
#endif
} TextureInfoContexts;

//*********************** end of data type definition for CABAC *******************

//! Pixel position for checking neighbors
typedef struct pix_pos
{
  int   available;
  int   mb_addr;
  short x;
  short y;
  short pos_x;
  short pos_y;
} PixelPos;

//! Motion Vector structure
typedef struct
{
  short mv_x;
  short mv_y;
} MotionVector;

//! definition of pic motion parameters
typedef struct pic_motion_params
{
  int64 ***   ref_pic_id;    //!< reference picture identifier [list][subblock_y][subblock_x]
  int64 ***   ref_id;        //!< reference picture identifier [list][subblock_y][subblock_x]
  short ****  mv;            //!< motion vector       [list][subblock_x][subblock_y][component]
  char  ***   ref_idx;       //!< reference picture   [list][subblock_y][subblock_x]
  byte *      mb_field;      //!< field macroblock indicator
  byte **     field_frame;   //!< indicates if co_located is field or frame.
} PicMotionParams;


//! Set Explicit GOP Parameters.
//! Currently only supports Enhancement GOP but could be easily extended
typedef struct
{
  int slice_type;         //! Slice type
  int display_no;         //! GOP Display order
  int reference_idc;      //! Is reference?
  int slice_qp;           //! Assigned QP
  int hierarchy_layer;    //! Hierarchy layer (used with GOP Hierarchy option 2)
  int hierarchyPocDelta;  //! Currently unused
} GOP_DATA;

//! Buffer structure for decoded reference picture marking commands
typedef struct DecRefPicMarking_s
{
  int memory_management_control_operation;
  int difference_of_pic_nums_minus1;
  int long_term_pic_num;
  int long_term_frame_idx;
  int max_long_term_frame_idx_plus1;
  struct DecRefPicMarking_s *Next;
} DecRefPicMarking_t;

typedef struct me_block
{
  struct video_par  *p_Vid;

  short            pos_x_padded;    //!< position x in image
  short            pos_y_padded;    //!< position y in image

  short            blocktype;
  short            blocksize_x;     //!< block size x dimension
  short            blocksize_y;     //!< block size y dimension
  short            blocksize_cr_x;  //!< block size x dimension
  short            blocksize_cr_y;  //!< block size y dimension

  short            block_x;         //!< position x in MB (subpartitions)
  short            block_y;         //!< position y in MB

  short            pos_x;           //!< position x in image
  short            pos_y;           //!< position y in image
  short            pos_x2;          //!< position x >> 2 in image
  short            pos_y2;          //!< position y >> 2 in image

  short            pos_cr_x;        //!< position x in image
  short            pos_cr_y;        //!< position y in image

  char             list;             //!< current list (list 0 or list 1). This is needed for bipredictive ME
  char             ref_idx;          //!< current reference index
  MotionVector     mv[2];            //!< motion vectors (L0/L1)
  PixelPos         block[4];
  struct slice    *p_slice;
  struct search_window searchRange;
  int              cost;           //!< Rate Distortion cost
  imgpel         **orig_pic;      //!< Block Data
  Boolean          ChromaMEEnable;
  int              ChromaMEWeight;
  // use weighted prediction based ME
  int              apply_weights;
  int              apply_bi_weights;
  // 8x8 transform
  int              test8x8;
  // Single List Prediction weights
  short            weight_luma;
  short            weight_cr[2];
  short            offset_luma;
  short            offset_cr[2];

  // Biprediction weights
  short            weight1;
  short            weight2;
  short            offsetBi;
  short            weight1_cr[2];
  short            weight2_cr[2];
  short            offsetBi_cr[2];
  int              search_pos2;   // <--  search positions for    half-pel search  (default: 9)
  int              search_pos4;   // <--  search positions for quarter-pel search  (default: 9)
  int              iteration_no;  // <--  bi pred iteration number

  // functions
  distblk (*computePredFPel)    (struct storable_picture *, struct me_block *, distblk , MotionVector * );
  distblk (*computePredHPel)    (struct storable_picture *, struct me_block *, distblk , MotionVector * );
  distblk (*computePredQPel)    (struct storable_picture *, struct me_block *, distblk , MotionVector * );
  distblk (*computeBiPredFPel)  (struct storable_picture *, struct storable_picture *, struct me_block *, distblk , MotionVector *, MotionVector *);
  distblk (*computeBiPredHPel)  (struct storable_picture *, struct storable_picture *, struct me_block *, distblk , MotionVector *, MotionVector *);
  distblk (*computeBiPredQPel)  (struct storable_picture *, struct storable_picture *, struct me_block *, distblk , MotionVector *, MotionVector *);
} MEBlock;

//! Syntax Element
typedef struct syntaxelement
{
  int                 type;           //!< type of syntax element for data part.
  int                 value1;         //!< numerical value of syntax element
  int                 value2;         //!< for blocked symbols, e.g. run/level
  int                 len;            //!< length of code
  int                 inf;            //!< info part of UVLC code
  unsigned int        bitpattern;     //!< UVLC bitpattern
  int                 context;        //!< CABAC context

#if TRACE
  #define             TRACESTRING_SIZE 100            //!< size of trace string
  char                tracestring[TRACESTRING_SIZE];  //!< trace string
#endif

  //!< for mapping of syntaxElement to UVLC
  void    (*mapping)(int value1, int value2, int* len_ptr, int* info_ptr);

} SyntaxElement;

//! Macroblock
typedef struct macroblock
{
  struct slice       *p_slice;                    //!< pointer to the current slice
  struct video_par     *p_Vid;                      //!< pointer to VideoParameters
  InputParameters    *p_Inp;
  int                 mbAddrX;                    //!< current MB address
  short               mb_type;                    //!< current MB mode type
  short               slice_nr;                   //!< current MB slice id
  short               mb_x;                       //!< current MB horizontal
  short               mb_y;                       //!< current MB vertical
  short               block_x;                    //!< current block horizontal
  short               block_y;                    //!< current block vertical

  short               pix_x;                      //!< current pixel horizontal
  short               pix_y;                      //!< current pixel vertical
  short               pix_c_x;                    //!< current pixel chroma horizontal
  short               pix_c_y;                    //!< current pixel chroma vertical

  short               opix_y;                     //!< current original picture pixel vertical
  short               opix_c_y;                   //!< current original picture pixel chroma vertical

  short               subblock_x;                 //!< current subblock horizontal
  short               subblock_y;                 //!< current subblock vertical

  short               list_offset;
  Boolean             prev_recode_mb;

  int                 mbAddrA, mbAddrB, mbAddrC, mbAddrD;
  byte                mbAvailA, mbAvailB, mbAvailC, mbAvailD;

  short               qp;                         //!< QP luma  
  short               qpc[2];                     //!< QP chroma
  short               qp_scaled[MAX_PLANE];       //!< QP scaled for all comps.
  short               qpsp;
  int                 cbp;
  short               prev_qp;
  short               prev_dqp;
  int                 prev_cbp;
  int                 cr_cbp[3];        // 444. Should be added in an external structure
  int                 c_nzCbCr[3];

  short               i16offset;  

  BitCounter          bits;

  short               ar_mode; //!< mb type to store adaptive rounding parameters
  short               mvd[2][BLOCK_MULTIPLE][BLOCK_MULTIPLE][2];          //!< indices correspond to [list][block_y][block_x][x,y]
  char                c_ipred_mode;      //!< chroma intra prediction mode
  char                i16mode;

  Info8x8             b8x8[4];    
  //imgpel              intra4x4_pred[16][17]; //!< 4x4 Intra prediction samples
  imgpel              intra4x4_pred[3][17]; //!< 4x4 Intra prediction samples
  imgpel              intra8x8_pred[3][25]; //!< 8x8 Intra prediction samples

  char                IntraChromaPredModeFlag;
  byte                mb_field;
  byte                is_field_mode;
  byte                luma_transform_size_8x8_flag;
  byte                temp_transform_size_8x8_flag;
  byte                NoMbPartLessThan8x8Flag;
  byte                valid_8x8;
  byte                write_mb;
  byte                is_intra_block;

  char                DFDisableIdc;
  char                DFAlphaC0Offset;
  char                DFBetaOffset;

  int                 skip_flag;

  char                intra_pred_modes   [MB_BLOCK_PARTITIONS];
  char                intra_pred_modes8x8[MB_BLOCK_PARTITIONS];           //!< four 8x8 blocks in a macroblock

  int64               cbp_blk ;    //!< 1 bit set for every 4x4 block with coefs (not implemented for INTRA)
  int64               cbp_bits[3];
  int64               cbp_bits_8x8[3];

  distblk              min_rdcost;
  distblk              min_dcost;
  distblk              min_rate;

  short               best_mode;
  char                best_c_imode;
  short               best_i16offset;
  char                best_i16mode;
  int                 best_cbp;

  //For residual DPCM
  short               ipmode_DPCM;



  struct macroblock   *mb_up;   //!< pointer to neighboring MB (CABAC)
  struct macroblock   *mb_left; //!< pointer to neighboring MB (CABAC)
  struct macroblock   *PrevMB;

  // Functions
  void (*GetMVPredictor) (struct macroblock *currMB, PixelPos *block, short pmv[2], short ref_frame, 
    char **refPic, short ***tmp_mv, int mb_x, int mb_y, int blockshape_x, int blockshape_y);
  distblk (*IntPelME)       (struct macroblock *currMB, MotionVector *, struct me_block *mv_block, distblk, int);
  distblk (*SubPelBiPredME) (struct macroblock *currMB, struct me_block *, int list, 
    MotionVector *pred_mv1, MotionVector *pred_mv2, MotionVector *mv1, MotionVector *mv2, distblk min_mcost, int* lambda_factor);
  distblk (*SubPelME)       (struct macroblock *currMB, 
    MotionVector *pred_mv, struct me_block *mv_block, distblk min_mcost, int* lambda_factor);
  distblk (*BiPredME)       (struct macroblock *currMB, int, MotionVector *, MotionVector *, MotionVector *, MotionVector *, struct me_block *, int, distblk, int);
  int (*trans_4x4)       (struct macroblock *currMB, ColorPlane pl, int block_x, int block_y, int *coeff_cost, int intra);
  int (*trans_16x16)     (struct macroblock *currMB, ColorPlane pl);
  int (*trans_8x8)       (struct macroblock *currMB, ColorPlane pl, int b8, int *coeff_cost, int intra);
  int  (*trans_cr_4x4[2])(struct macroblock *currMB, int uv,int i11);

  void (*cbp_linfo_intra)(int cbp, int dummy, int *len,int *info);
  void (*cbp_linfo_inter)(int cbp, int dummy, int *len,int *info);    
} Macroblock;

//! Bitstream
typedef struct
{
  int     buffer_size;        //!< Buffer size      
  int     byte_pos;           //!< current position in bitstream;
  int     bits_to_go;         //!< current bitcounter
  
  int     stored_byte_pos;    //!< storage for position in bitstream;
  int     stored_bits_to_go;  //!< storage for bitcounter
  int     byte_pos_skip;      //!< storage for position in bitstream;
  int     bits_to_go_skip;    //!< storage for bitcounter
  int     write_flag;         //!< Bitstream contains data and needs to be written

  byte    byte_buf;           //!< current buffer for last written byte
  byte    stored_byte_buf;    //!< storage for buffer of last written byte
  byte    byte_buf_skip;      //!< current buffer for last written byte
  byte    *streamBuffer;      //!< actual buffer for written bytes

#if TRACE
  Boolean trace_enabled;
#endif
} Bitstream;


//! DataPartition
typedef struct datapartition
{
  struct slice        *p_Slice;
  struct video_par      *p_Vid;
  InputParameters     *p_Inp;   //!< pointer to the input parameters
  Bitstream           *bitstream;
  NALU_t              *nal_unit;
  EncodingEnvironment ee_cabac;
  EncodingEnvironment ee_recode;
} DataPartition;

/*! For MB level field/frame coding tools  
    temporary structure to store MB data for field/frame coding */
typedef struct rd_data
{
  double  min_rdcost;
  double  min_dcost;
  double  min_rate;

  imgpel  ***rec_mb;            //!< hold the components of reconstructed MB
  int     ****cofAC;
  int     ***cofDC;
  int     *****cofAC_new;
  short   mb_type;  
  int64   cbp_blk;
  int     cbp;
  int     prev_cbp;
  int     mode;
  short   i16offset;
  char    i16mode;
  Boolean NoMbPartLessThan8x8Flag;

  char    c_ipred_mode;
  short   qp;
  short   prev_qp;
  short   prev_dqp;
  byte    luma_transform_size_8x8_flag;
  Info8x8 block;
  Info8x8 b8x8[4];
  
  // These need to be changed to MotionVector parameters
  short   ******all_mv;               //!< all modes motion vectors
  short *******bipred_mv;             //!<Biprediction MVs  

  char    intra_pred_modes[16];
  char    intra_pred_modes8x8[16];
  char    **ipredmode;
  char    ***refar;                  //!< reference frame array [list][y][x]
} RD_DATA;

 //! Slice
typedef struct slice
{
  struct video_par      *p_Vid;   // pointer to the original video structure
  InputParameters     *p_Inp;   // pointer to the input parameters
  pic_parameter_set_rbsp_t *active_pps;
  seq_parameter_set_rbsp_t *active_sps;

  int                 picture_id;
  int                 qp;
  int                 qs;
  short               slice_type;   //!< picture type
  unsigned int        frame_num;
  unsigned int        max_frame_num;
  signed int          framepoc;     //!< min (toppoc, bottompoc)
  signed int          ThisPOC;      //!< current picture POC
  short               slice_nr;
  int                 model_number;
  char                colour_plane_id;   //!< colour plane id for 4:4:4 profile
  int                 lossless_qpprime_flag;
  int                 P444_joined;
  int                 disthres;
  int                 Transform8x8Mode;
  // RDOQ at the slice level (note this allows us to reduce passed parameters, 
  // but also could enable RDOQ control at the slice level.
  int                 UseRDOQuant;
  int                 RDOQ_QP_Num;

  char                num_ref_idx_active[2];
  int                 ref_pic_list_reordering_flag[2];
  int                 *reordering_of_pic_nums_idc[2];
  int                 *abs_diff_pic_num_minus1[2];
  int                 *long_term_pic_idx[2];

  int                 width_blk;               //!< Number of columns in blocks
  int                 height_blk;              //!< Number of lines in blocks

  PictureStructure    structure;
  Boolean             mb_aff_frame_flag;

  int                 start_mb_nr;
  int                 max_part_nr;  //!< number of different partitions
  int                 num_mb;       //!< number of MBs in the slice

  int                 cmp_cbp[3];
  int                 curr_cbp[2];
  char                symbol_mode;
  short               NoResidueDirect;
  short               partition_mode;
  short               idr_flag;
  int                 frame_no;
  unsigned int        PicSizeInMbs;
  int                 num_blk8x8_uv;
  int                 nal_reference_idc;                       //!< nal_reference_idc from NAL unit
  short               bitdepth_luma;
  short               bitdepth_chroma;


  DataPartition       *partArr;     //!< array of partitions
  MotionInfoContexts  *mot_ctx;     //!< pointer to struct of context models for use in CABAC
  TextureInfoContexts *tex_ctx;     //!< pointer to struct of context models for use in CABAC

  int                 mvscale[6][MAX_REFERENCE_PICTURES];
  char                direct_spatial_mv_pred_flag;              //!< Direct Mode type to be used (0: Temporal, 1: Spatial)
  // Deblocking filter parameters
  char                DFDisableIdc;                             //!< Deblocking Filter Disable indicator
  char                DFAlphaC0Offset;                          //!< Deblocking Filter Alpha offset
  char                DFBetaOffset;                             //!< Deblocking Filter Beta offset

  byte                weighted_prediction;                       //!< Use of weighted prediction 
  byte                weighted_bipred_idc;                      //!< Use of weighted biprediction (note that weighted_pred_flag is probably not needed here)

  short               luma_log_weight_denom;
  short               chroma_log_weight_denom;
  short               wp_luma_round;
  short               wp_chroma_round;

  short  max_num_references;      //!< maximum number of reference pictures that may occur
  // Motion vectors for a macroblock
  // These need to be changed to MotionVector parameters
  short ******all_mv;         //!< replaces local all_mv
  short *******bipred_mv;     //!< Biprediction MVs  
  //Weighted prediction
  short ***wp_weight;         //!< weight in [list][index][component] order
  short ***wp_offset;         //!< offset in [list][index][component] order
  short ****wbp_weight;       //!< weight in [list][fwd_index][bwd_idx][component] order

  int *****cofAC_new;          //!< AC coefficients [comp][8x8block][4x4block][level/run][scan_pos]

  int ****cofAC;               //!< AC coefficients [8x8block][4x4block][level/run][scan_pos]
  int *** cofDC;               //!< DC coefficients [yuv][level/run][scan_pos]

  // For rate control
  int diffy[16][16];
  
  int64 cur_cbp_blk[MAX_PLANE];
  int coeff_cost_cr[MAX_PLANE];


  int **tblk16x16;   //!< Transform related array
  int **tblk4x4;     //!< Transform related array
  int ****i16blk4x4;

  RD_DATA *rddata;
  // RD_DATA data. Moved here to enable parallelization at the slice level
  // of RDOQ
  RD_DATA rddata_trellis_best;
  RD_DATA rddata_trellis_curr;
  //!< For MB level field/frame coding tools
  RD_DATA rddata_top_frame_mb;
  RD_DATA rddata_bot_frame_mb; 
  RD_DATA rddata_top_field_mb;
  RD_DATA rddata_bot_field_mb;

  Boolean si_frame_indicator;
  Boolean sp2_frame_indicator;

  char    ***direct_ref_idx;           //!< direct mode reference index buffer
  char    **direct_pdir;               //!< direct mode direction buffer

  MotionVector ****tmp_mv8;
  MotionVector ****tmp_mv4;

  distblk    ***motion_cost8;
  distblk    ***motion_cost4;
  int     deltaQPTable[9]; 

  // RDOQ
  struct est_bits_cabac *estBitsCabac; // [NUM_BLOCK_TYPES]
  double norm_factor_4x4;
  double norm_factor_8x8;

  imgpel ****mpr_4x4;           //!< prediction samples for   4x4 intra prediction modes
  imgpel ****mpr_8x8;           //!< prediction samples for   8x8 intra prediction modes
  imgpel ****mpr_16x16;         //!< prediction samples for 16x16 intra prediction modes (and chroma)
  imgpel ***mb_pred;            //!< current best prediction mode
  int ***mb_rres;               //!< the diff pixel values between the original macroblock/block and its prediction (reconstructed)
  int ***mb_ores;               //!< the diff pixel values between the original macroblock/block and its prediction (original)

  // Residue Color Transform
  char b8_ipredmode8x8[4][4];
  char b8_intra_pred_modes8x8[16];
  char b4_ipredmode[16];
  char b4_intra_pred_modes[16];

  struct rdo_structure    *p_RDO;
  struct colocated_params *p_colocated;
  struct epzs_params      *p_EPZS;  

  // This should be the right location for this
  struct storable_picture **listX[6];
  char listXsize[6];

  // Some Cabac related parameters (could be put in a different structure so we can dynamically allocate them when needed)
  int  coeff[64];
  int  coeff_ctr;
  int  pos;


  // Function pointers
  int     (*Mode_Decision_for_4x4IntraBlocks)   (Macroblock *currMB, int  b8,  int  b4,  int  lambda, distblk*  min_cost);
  int     (*Mode_Decision_for_8x8IntraBlocks)   (Macroblock *currMB, int b8, int lambda, distblk *min_cost);
  distblk   (*rdcost_for_4x4_intra_blocks)        (Macroblock *currMB, int *c_nz, int b8, int b4, int ipmode, int lambda, int mostProbableMode, distblk min_rdcost);
  distblk   (*rdcost_for_8x8_intra_blocks)        (Macroblock *currMB, int *c_nz, int b8, int ipmode, int lambda, distblk min_rdcost, int mostProbableMode);

  void    (*Get_Direct_Motion_Vectors)          (Macroblock *currMB);
  Boolean (*slice_too_big)                      (int bits_slice); //!< for use of callback functions
  void    (*SetMotionVectorsMB)                 (Macroblock *currMB, struct pic_motion_params *motion);
  void    (*encode_one_macroblock)              (Macroblock *currMB);
  void    (*set_stored_mb_parameters)           (Macroblock *currMB);
  void    (*set_ref_and_motion_vectors)         (Macroblock *currMB, struct pic_motion_params *motion, Info8x8 *part, int block);

  void    (*compute_cost8x8)                    (struct video_par *p_Vid, imgpel **cur_img, imgpel **prd_img, int pic_opix_x, distblk *cost, distblk min_cost);
  void    (*compute_cost4x4)                    (struct video_par *p_Vid, imgpel **cur_img, imgpel **prd_img, int pic_opix_x, distblk *cost, distblk min_cost);
  distblk   (*find_sad_16x16)                   (Macroblock *currMB);
  distblk    (*distI16x16)                      (Macroblock *currMB, imgpel **img_org, imgpel **pred_img, distblk min_cost);

  void    (*SetLagrangianMultipliers)           (struct video_par *p_Vid, InputParameters *p_Inp);

  // Write MB info function pointers
  void (*writeMB_typeInfo)      (Macroblock *currMB, SyntaxElement *se, DataPartition *dP);
  void (*writeMB_Skip)          (Macroblock *currMB, SyntaxElement *se, DataPartition *dP);
  void (*writeIntraPredMode)    (SyntaxElement *se, DataPartition *dP);
  void (*writeB8_typeInfo)      (SyntaxElement *se, DataPartition *dP);
  void (*writeRefFrame[6])      (SyntaxElement *se, DataPartition *dP);
  int  (*writeMotionInfo2NAL)   (Macroblock* currMB);
  int  (*write_MB_layer)        (Macroblock *currMB, int rdopt, int *coeff_rate);
  void (*writeMVD)              (Macroblock *currMB, SyntaxElement *se, DataPartition *dP);
  void (*writeCBP)              (Macroblock* currMB, SyntaxElement *se, DataPartition *dP);
  void (*writeDquant)           (Macroblock* currMB, SyntaxElement *se, DataPartition *dP);
  void (*writeCIPredMode)       (Macroblock *currMB, SyntaxElement *se, DataPartition *dP);
  void (*writeFieldModeInfo)    (Macroblock *currMB, SyntaxElement *se, DataPartition *dP);
  void (*writeMB_transform_size)(Macroblock *currMB, SyntaxElement *se, DataPartition *dP);
  int  (*writeCoeff16x16)       (Macroblock* currMB, ColorPlane);
  int  (*writeCoeff4x4_CAVLC)   (Macroblock* currMB, int block_type, int b8, int b4, int param);
  void (*write_and_store_CBP_block_bit) (Macroblock* currMB, EncodingEnvironmentPtr eep_dp, int type, int cbp_bit, TextureInfoContexts*  tex_ctx);

  // Coding state
  void (*reset_coding_state) (Macroblock *currMB, struct coding_state *);
  void (*store_coding_state) (Macroblock *currMB, struct coding_state *);

  // Reordering
  void (*poc_ref_pic_reorder_frame) ( struct slice *currSlice, struct storable_picture **list, unsigned num_ref_idx_lX_active, 
    int *reordering_of_pic_nums_idc, int *abs_diff_pic_num_minus1, int *long_term_pic_idx, int list_no );

  // Quantization
  int (*quant_4x4)     (Macroblock *currMB, int **tblock, struct quant_methods *q_method);
  int (*quant_ac4x4cr) (Macroblock *currMB, int **tblock, struct quant_methods *q_method);
  int (*quant_dc4x4)   (Macroblock *currMB, int **tblock, int qp, int*  DCLevel, int*  DCRun, LevelQuantParams *q_params_4x4, const byte (*pos_scan)[2]);
  int (*quant_ac4x4)   (Macroblock *currMB, int **tblock, struct quant_methods *q_method);
  int (*quant_8x8)     (Macroblock *currMB, int **tblock, struct quant_methods *q_method);
  int (*quant_8x8cavlc)(Macroblock *currMB, int **tblock, struct quant_methods *q_method, int***  cofAC);

  int (*quant_dc_cr)     (Macroblock *currMB, int **tblock, int qp, int* DCLevel, int* DCRun, 
    LevelQuantParams *q_params_4x4, int **fadjust, const byte (*pos_scan)[2]);
  void (*rdoq_4x4)       (Macroblock *currMB, int **tblock, struct quant_methods *q_method, int levelTrellis[16]);

  void (*rdoq_dc)        (Macroblock *currMB, int **tblock, int qp_per, int qp_rem, LevelQuantParams *q_params_4x4, 
    const byte (*pos_scan)[2], int levelTrellis[16], int type);

  void (*rdoq_ac4x4)     (Macroblock *currMB, int **tblock , struct quant_methods *q_method, int levelTrellis[16]);

  void (*rdoq_dc_cr)     (Macroblock *currMB, int **tblock, int qp_per, int qp_rem, LevelQuantParams *q_params_4x4, 
    const byte (*pos_scan)[2], int levelTrellis[16], int type);

  void (*set_modes_and_reframe)    (Macroblock *currMB, int b8, short* p_dir, int list_mode[2], char *list_ref_idx);
  void (*store_8x8_motion_vectors) (struct slice *currSlice, int dir, int block8x8, Info8x8 *B8x8Info);

  distblk (*getDistortion)       ( Macroblock *currMB );
} Slice;

//! DistortionParams
typedef struct distortion_metric
{
  float value[3];                    //!< current frame distortion
  float average[3];                  //!< average frame distortion
  float avslice[NUM_SLICE_TYPES][3]; //!< average frame type distortion
} DistMetric;

typedef struct distortion_params
{
  int        frame_ctr;                     //!< number of coded frames
  DistMetric metric[TOTAL_DIST_TYPES];      //!< Distortion metrics
} DistortionParams;

typedef struct picture
{
  int   no_slices;
  int   bits_per_picture;
  struct slice *slices[MAXSLICEPERPICTURE];

  DistMetric distortion;
  byte  idr_flag;
} Picture;

//! block 8x8 temporary RD info
typedef struct block_8x8_info
{
  Info8x8 best[MAXMODE][4];
} Block8x8Info;
                             
//! VideoParameters
typedef struct video_par
{
  InputParameters          *p_Inp;
  pic_parameter_set_rbsp_t *active_pps;
  seq_parameter_set_rbsp_t *active_sps;
  struct sei_params        *p_SEI;
  struct decoders          *p_decs;

  int number;                  //!< current image number to be encoded (in first layer)  
  int LevelIndex;              //!< mapped level idc
  int MaxVmvR[6];              //!< maximum vertical motion vector
  int MaxHmvR[6];              //!< maximum horizontal motion vector
  int current_mb_nr;
  short current_slice_nr;
  short type;
  PictureStructure structure;  //!< picture structure
  int base_dist;
  int num_ref_frames;          //!< number of reference frames to be used
  int max_num_references;      //!< maximum number of reference pictures that may occur
  int masterQP;                //!< Master quantization parameter
  int qp;                      //!< quant for the current frame
  int qpsp;                    //!< quant for the prediction frame of SP-frame

  int prev_frame_no; // POC200301
  int consecutive_non_reference_pictures; // POC200301

  // prediction structure
  int frm_struct_buffer;      //!< length of the frame struct buffer (this may also be a define)
  int last_idr_code_order;
  int last_idr_disp_order;
  int last_mmco_5_code_order; //!< it is a good idea to re-initialize POCs after such a frame: while seemingly not required it is critical for good implicit WP performance (I *THINK*) since POCs are set internally to zero after an MMCO=5 command.
  int last_mmco_5_disp_order;
  int curr_frm_idx;           //!< frame we wish to code in coding order (points also to p_frm_struct)
  int FrameNumOffset;         //!< POC type 1
  int prevFrameNumOffset;     //!< POC type 1
  unsigned int prevFrameNum;  //!< POC type 1
  SeqStructure    *p_pred;
  FrameUnitStruct *p_curr_frm_struct;
  PicStructure    *p_curr_pic;
  SliceStructure  *p_curr_slice;

  struct search_window searchRange;
  ImageData imgData;           //!< Image data to be encoded
  ImageData imgData0;          //!< Input Image Data
  ImageData imgData1;
  ImageData imgData2;

  struct image_structure imgSRC;
  struct image_structure imgREF;
  struct image_structure imgRGB_src;
  struct image_structure imgRGB_ref;
  int       **imgY_sub_tmp;           //!< Y picture temporary component (Quarter pel)
  imgpel    **imgY_com;               //!< Encoded luma images
  imgpel   ***imgUV_com;              //!< Encoded croma images

  imgpel    **pCurImg;                //!< Reference image. Luma for other profiles, can be any component for 4:4:4
  // global picture format dependend buffers, mem allocation in image.c
  imgpel    **pImgOrg[MAX_PLANE];

  Picture *p_frame_pic;
  Picture **frame_pic;
  Picture **field_pic;
  Picture *frame_pic_si;

  byte *MapUnitToSliceGroupMap;

  byte *buf;
  byte *ibuf;

#ifdef _LEAKYBUCKET_
  long *Bit_Buffer;
  unsigned long total_frame_buffer;
#endif

  unsigned int log2_max_frame_num_minus4;
  unsigned int log2_max_pic_order_cnt_lsb_minus4;
  unsigned int max_frame_num;
  unsigned int max_pic_order_cnt_lsb;

  int64  me_tot_time;
  int64  tot_time;
  int64  me_time;

  byte mixedModeEdgeFlag;

  int *RefreshPattern;
  int *IntraMBs;
  int WalkAround;
  int NumberOfMBs;
  int NumberIntraPerPicture;

  short start_me_refinement_hp; //!< if set then recheck the center position when doing half-pel motion refinement
  short start_me_refinement_qp; //!< if set then recheck the center position when doing quarter-pel motion refinement

  // Motion Estimation
  struct umhex_struct *p_UMHex;
  struct umhex_smp_struct *p_UMHexSMP;
  struct me_full_fast *p_ffast_me;

  struct search_window *p_search_window;

  // EPZS
  struct epzs_struct *sdiamond;
  struct epzs_struct *square;
  struct epzs_struct *ediamond;
  struct epzs_struct *ldiamond;
  struct epzs_struct *sbdiamond;
  struct epzs_struct *pmvfast;


  // RDOQ
  int precalcUnaryLevelTab[128][MAX_PREC_COEFF];
  int AdaptRndWeight;
  int AdaptRndCrWeight;


  //////////////////////////////////////////////////////////////////////////
  // B pictures
  // motion vector : forward, backward, direct
  byte  MBPairIsField;     //!< For MB level field/frame coding tools

  // Buffers for rd optimization with packet losses, Dim. Kontopodis
  byte **pixel_map;   //!< Shows the latest reference frame that is reliable for each pixel
  byte **refresh_map; //!< Stores the new values for pixel_map
  int intras;         //!< Counts the intra updates in each frame.

  int RCMinQP;
  int RCMaxQP;

  float framerate;
  int width;                   //!< Number of pels
  int width_padded;            //!< Width in pels of padded picture
  int width_blk;               //!< Number of columns in blocks
  int width_cr;                //!< Number of pels chroma
  int height;                  //!< Number of lines
  int height_padded;           //!< Number in lines of padded picture
  int height_blk;              //!< Number of lines in blocks
  int height_cr;               //!< Number of lines  chroma
  int height_cr_frame;         //!< Number of lines  chroma frame
  int size;                    //!< Luma Picture size in pels
  int size_cr;                 //!< Chroma Picture size in pels

  int is_v_block;
  int mb_y_upd;
  int mb_y_intra;              //!< which GOB to intra code
  char **ipredmode;            //!< intra prediction mode
  char **ipredmode8x8;         //!< help storage for 8x8 modes, inserted by YV
  //fast intra prediction;
  char **ipredmode4x4_line;     //!< intra prediction mode
  char **ipredmode8x8_line;     //!< help storage for 8x8 modes, inserted by YV

  int cod_counter;             //!< Current count of number of skipped macroblocks in a row
  int ***nz_coeff;             //!< number of coefficients per block (CAVLC)
  int pix_x;                   //!< current pixel horizontal
  int pix_y;                   //!< current pixel vertical


  imgpel min_IPCM_value;
  // Cabac related parameters (we can put these in a structure)
  int pic_bin_count;

  // Adaptive rounding
  int ****ARCofAdj4x4;         //!< Transform coefficients for 4x4 luma/chroma. 
  int ****ARCofAdj8x8;         //!< Transform coefficients for 4x4 luma/chroma. 


  Picture       *currentPicture;         //!< The coded picture currently in the works (typically p_frame_pic, p_Vid->field_pic[0], or p_Vid->field_pic[1])
  struct slice  *currentSlice;           //!< pointer to current Slice data struct
  Macroblock    *mb_data;                //!< array containing all MBs of a whole frame
  Block8x8Info  *b8x8info;               //!< block 8x8 information for RDopt

  //FAST_REFPIC_DECISION
  int           mb_refpic_used; //<! [2][16] for fast reference decision;
  int           parent_part_refpic_used; //<! [2][16] for fast reference decision;
  int           submb_parent_part_refpic_used;//<! [2][16] for fast reference decision;
  //end;
  int *intra_block;

  int frame_no;
  int fld_type;                          //!< top or bottom field
  byte fld_flag;
  unsigned int rd_pass;

  int  redundant_coding;
  int  key_frame;
  int  redundant_ref_idx;
  int  frm_no_in_file;

  char DFDisableIdc;
  char DFAlphaC0Offset;
  char DFBetaOffset;

  char direct_spatial_mv_pred_flag;              //!< Direct Mode type to be used (0: Temporal, 1: Spatial)

  int pad_size_uv_x;
  int pad_size_uv_y;
  unsigned char chroma_mask_mv_y;
  unsigned char chroma_mask_mv_x;
  int chroma_shift_y, chroma_shift_x;
  int shift_cr_x, shift_cr_x2, shift_cr_y;
  int padded_size_x;  
  int padded_size_x_m8x8;
  int padded_size_x_m4x4;
  int cr_padded_size_x;
  int cr_padded_size_x_m8;
  int cr_padded_size_x2;
  int cr_padded_size_x4;

  int num_ref_idx_l0_active;
  int num_ref_idx_l1_active;

  Boolean field_mode;     //!< For MB level field/frame -- field mode on flag
  Boolean top_field;      //!< For MB level field/frame -- top field flag

  int layer;              //!< which layer this picture belonged to

  int AdaptiveRounding;   //!< Adaptive Rounding parameter based on JVT-N011

  int redundant_pic_cnt;

  Boolean mb_aff_frame_flag;    //!< indicates frame with mb aff coding

  //the following should probably go in sequence parameters
  unsigned int pic_order_cnt_type;

  // for poc mode 1
  Boolean      delta_pic_order_always_zero_flag;
  int          offset_for_non_ref_pic;
  int          offset_for_top_to_bottom_field;
  unsigned int num_ref_frames_in_pic_order_cnt_cycle;
  int          offset_for_ref_frame[1];

  //the following is for slice header syntax elements of poc
  // for poc mode 0.
  unsigned int pic_order_cnt_lsb;
  int          delta_pic_order_cnt_bottom;
  // for poc mode 1.
  int          delta_pic_order_cnt[2];

  int          frm_iter;   //!< frame variations to create (useful for multiple coding passes)

  unsigned int field_picture;
  signed int toppoc;       //!< poc for this frame or field
  signed int bottompoc;    //!< for completeness - poc of bottom field of a frame (always = poc+1)
  signed int framepoc;     //!< min (toppoc, bottompoc)
  signed int ThisPOC;      //!< current picture POC
  unsigned int frame_num;    //!< frame_num for this frame

  unsigned int PicWidthInMbs;
  unsigned int PicHeightInMapUnits;
  unsigned int FrameHeightInMbs;
  unsigned int PicSizeInMbs;
  unsigned int FrameSizeInMbs;

  //the following should probably go in picture parameters
  Boolean bottom_field_pic_order_in_frame_present_flag; // ????????

  //the following are sent in the slice header
  //  int delta_pic_order_cnt[2];
  int nal_reference_idc;

  int     adaptive_ref_pic_buffering_flag;
  int     no_output_of_prior_pics_flag;
  Boolean long_term_reference_flag;

  DecRefPicMarking_t *dec_ref_pic_marking_buffer;
  int* mvbits;

  int     max_mvd;  //for MVD overflow checking;
  int*    refbits;

  // rate control variables
  int NumberofCodedMacroBlocks;
  int BasicUnitQP;
  int NumberofMBTextureBits;
  int NumberofMBHeaderBits;
  unsigned int BasicUnit;
  byte write_macroblock;
  byte bot_MB;
  int write_mbaff_frame;

  int DeblockCall;

  int last_pic_bottom_field;
  int last_has_mmco_5;
  int pre_frame_num;

  int slice_group_change_cycle;

  short bitdepth_luma;
  short bitdepth_chroma;
  int bitdepth_scale[2];
  int bitdepth_luma_qp_scale;
  int bitdepth_chroma_qp_scale;
  int bitdepth_lambda_scale;
  int max_bitCount;
  int max_qp_delta;
  int min_qp_delta;
  // Lagrangian Parameters
  LambdaParams **lambda;  
  double  **lambda_md;     //!< Mode decision Lambda
  double ***lambda_me;     //!< Motion Estimation Lambda
  int    ***lambda_mf;     //!< Integer formatted Motion Estimation Lambda

  double **lambda_mf_factor; //!< Motion Estimation Lamda Scale Factor

  imgpel dc_pred_value_comp[MAX_PLANE]; //!< component value for DC prediction (depends on component pel bit depth)
  imgpel dc_pred_value;                 //!< DC prediction value for current component
  int max_pel_value_comp      [MAX_PLANE];       //!< max value that one picture element (pixel) can take (depends on pic_unit_bitdepth)
  int max_imgpel_value_comp_sq   [MAX_PLANE];       //!< max value that one picture element (pixel) can take (depends on pic_unit_bitdepth)
  short max_imgpel_value;              //!< max value that one picture element (pixel) can take (depends on pic_unit_bitdepth)

  int num_blk8x8_uv;
  int num_cdc_coeff;
  ColorFormat yuv_format;
  int P444_joined;
  int lossless_qpprime_flag;
  short mb_cr_size_x;
  short mb_cr_size_y;
  int mb_size[MAX_PLANE][2];
  struct edge_info *p_edge;

  int chroma_qp_offset[2];      //!< offset for qp for chroma [0-Cb, 1-Cr]

  int auto_crop_right;
  int auto_crop_bottom;

  short checkref;
  int last_valid_reference;
  int bytes_in_picture;

  int64 last_bit_ctr_n;

  int AverageFrameQP;
  int SumFrameQP;

  int ChromaArrayType;
  Macroblock    *mb_data_JV[MAX_PLANE];  //!< mb_data to be used during 4:4:4 independent mode encoding
  char colour_plane_id;    //!< colour_plane_id of the current coded slice (valid only when separate_colour_plane_flag is 1)

  int lastIntraNumber;
  int lastINTRA;
  int last_ref_idc;
  int idr_refresh;

  int p_dec;                      //!< decoded image file handle
  int frame_statistic_start;
  int initial_Bframes;
  int cabac_encoding;

  unsigned int primary_pic_type;

  int frameOffsetTotal[2][MAX_REFERENCE_PICTURES]; 
  int frameOffsetCount[2][MAX_REFERENCE_PICTURES]; 
  short frameOffset[2][MAX_REFERENCE_PICTURES];
  int frameOffsetAvail; 

  double *mb16x16_cost_frame;
  double mb16x16_cost;

  int **lrec;
  int ***lrec_uv;
  Boolean sp2_frame_indicator;
  int number_sp2_frames;

  Boolean giRDOpt_B8OnlyFlag;

  int  frameNuminGOP;
  // Redundant picture
  imgpel **imgY_tmp;
  imgpel **imgUV_tmp[2];

  // RDOQ related parameters
  int Motion_Selected;
  //int Intra_Selected; 

  int CbCr_predmode_8x8[4]; 
  distblk ****motion_cost;
  int*** initialized;
  int*** modelNumber;

  int num_mb_per_slice;
  int number_of_slices;

  int  imgpel_abs_range;
#if (JM_MEM_DISTORTION)
  int* imgpel_abs;
  int* imgpel_quad;
#endif

  GOP_DATA *gop_structure;
  byte *MBAmap;
  unsigned int PicSizeInMapUnits;
  int FirstMBInSlice[MAXSLICEGROUPIDS];

  MotionVector* spiral_search;
  MotionVector* spiral_hpel_search;
  MotionVector* spiral_qpel_search;

  // files
  FILE *p_log;                     //!< SNR file

  // ANNEX B output file
  FILE *f_annexb; 
  // RTP output file
  FILE *f_rtp;
  int CurrentRTPTimestamp;             //!< The RTP timestamp of the current packet,
  //! incremented with all P and I frames
  uint16 CurrentRTPSequenceNumber;     //!< The RTP sequence number of the current packet
  //!< incremented by one for each sent packet

  // This should be the right location for this
  struct storable_picture **listX[6];
  char listXsize[6];  

  DistortionParams *p_Dist;
  struct stat_parameters  *p_Stats;
  pic_parameter_set_rbsp_t *PicParSet[MAXPPS];
  struct decoded_picture_buffer *p_Dpb;
  struct frame_store            *out_buffer;
  struct storable_picture       *enc_picture;
  struct storable_picture       **enc_frame_picture;
  struct storable_picture       **enc_field_picture;
  struct storable_picture       *enc_frame_picture_JV[MAX_PLANE];  //!< enc_frame to be used during 4:4:4 independent mode encoding
  struct quant_params           *p_Quant;
  struct scaling_list           *p_QScale;

  // rate control
  RCGeneric   *p_rc_gen;
  RCGeneric   *p_rc_gen_init, *p_rc_gen_best;
  RCQuadratic *p_rc_quad;
  RCQuadratic *p_rc_quad_init, *p_rc_quad_best;

  double entropy[128];
  double enorm  [128];
  double probability[128];


  FILE       *expSFile;
  struct exp_seq_info *expSeq;
  // Weighted prediction
  struct wpx_object   *pWPX;

#ifdef BEST_NZ_COEFF
  int gaaiMBAFF_NZCoeff[4][12];
#endif

  int offset_y, offset_cr;
  int wka0, wka1, wka2, wka3, wka4;
  // Note that these function pointers definitely affect now parallelization since they are only
  // allocated once. We need to add such info at maybe within picture information or at a lower level
  // RC
  int  (*updateQP)                (struct video_par *p_Vid, InputParameters *p_Inp, RCQuadratic *p_quad, RCGeneric *p_gen, int topfield);
  void (*rc_update_pict_frame_ptr)(struct video_par *p_Vid, InputParameters *p_Inp, RCQuadratic *p_quad, RCGeneric *p_gen, int nbits);
  void (*rc_update_picture_ptr)   (struct video_par *p_Vid, InputParameters *p_Inp, int bits);
  void (*rc_init_pict_ptr)        (struct video_par *p_Vid, InputParameters *p_Inp, RCQuadratic *p_quad, RCGeneric *p_gen, int fieldpic, int topfield, int targetcomputation, float mult);
  
  //Various
  void (*buf2img)              (imgpel** imgX, unsigned char* buf, int size_x, int size_y, int o_size_x, int o_size_y, int symbol_size_in_bytes, int bitshift);
  void (*getNeighbour)         (Macroblock *currMB, int xN, int yN, int mb_size[2], PixelPos *pix);
  void (*get_mb_block_pos)     (int mb_addr, short *x, short *y);
  int  (*WriteNALU)            (struct video_par *p_Vid, NALU_t *n);     //! Hides the write function in Annex B or RTP
  void (*error_conceal_picture)(struct video_par *p_Vid, struct storable_picture *enc_pic, int decoder);
  
  // function pointer for different ways of obtaining chroma interpolation
  void (*OneComponentChromaPrediction4x4)   (Macroblock *currMB, imgpel* , int , int , short*** , struct storable_picture *listX, int );
  
  // deblocking
  void (*GetStrength)    (byte Strength[16], Macroblock *MbQ, int dir,int edge, int mvlimit);
  void (*EdgeLoopLuma)   (ColorPlane pl, imgpel** Img, byte Strength[16], Macroblock *MbQ, int dir, int edge, int width);
  void (*EdgeLoopChroma)(imgpel** Img, byte Strength[16], Macroblock *MbQ, int dir, int edge, int width, int uv);

  // We should move these at the slice level at some point.
  void (*EstimateWPBSlice) (struct slice *currSlice);
  void (*EstimateWPPSlice) (struct slice *currSlice, int offset);
  int  (*TestWPPSlice)     (struct video_par *p_Vid, int offset);
  int  (*TestWPBSlice)     (struct video_par *p_Vid, int method);
  distblk  (*distortion4x4)(int*, distblk);
  distblk  (*distortion8x8)(int*, distblk);

  // ME distortion Function pointers. We need to move this to the MB or slice level
  distblk (*computeUniPred[6])   (struct storable_picture *ref1, struct me_block *, distblk , MotionVector * );
  distblk (*computeBiPred1[3])   (struct storable_picture *ref1, struct storable_picture *ref2, struct me_block*, distblk , MotionVector *, MotionVector *);
  distblk (*computeBiPred2[3])   (struct storable_picture *ref1, struct storable_picture *ref2, struct me_block*, distblk , MotionVector *, MotionVector *);
} VideoParameters;


//! definition of pic motion parameters
typedef struct pic_motion_params2
{
  int64    ref_pic_id;    //!< reference picture identifier [list][subblock_y][subblock_x]
  int64    ref_id;        //!< reference picture identifier [list][subblock_y][subblock_x]
  short    mv[2];         //!< motion vector       [list][subblock_x][subblock_y][component]
  char     ref_idx;       //!< reference picture   [list][subblock_y][subblock_x]
  byte     mb_field;      //!< field macroblock indicator
  byte     field_frame;   //!< indicates if co_located is field or frame.
} PicMotionParams2;

typedef struct
{
  distblk  mb_p8x8_cost;  
  distblk  smb_p8x8_cost[4];  
  distblk  smb_p8x8_rdcost[4];
  Info8x8 part[4];

  int cbp8x8;
  int cbp_blk8x8;
  int cnt_nonz_8x8;    
  
  int **lrec; //!< transform and quantized coefficients will be stored here for SP frames
  imgpel **mpr8x8;
  imgpel ***mpr8x8CbCr;
  imgpel **rec_mbY8x8;
  imgpel ***rec_mb8x8_cr;
} RD_8x8DATA;

typedef struct
{
  double lambda_md;        //!< Mode decision Lambda
//#if JCOST_CALC_SCALEUP
  int    lambda_mdfp;       //!< Fixed point mode decision lambda;
//#endif
  double lambda_me[3];     //!< Motion Estimation Lambda
  int    lambda_mf[3];     //!< Integer formatted Motion Estimation Lambda
  int    best_mcost[2];

  short  valid[MAXMODE];
  short  curr_mb_field;
} RD_PARAMS;

typedef struct scaling_list      ScaleParameters;
typedef struct rdo_structure     RDOPTStructure;

typedef struct encoder_params
{
  InputParameters   *p_Inp;          //!< Input Parameters
  VideoParameters   *p_Vid;          //!< Image Parameters
  FILE              *p_trace;        //!< Trace file
  int64              bufferSize;     //!< buffer size for tiff reads (not currently supported)
} EncoderParams;

extern EncoderParams  *p_Enc;

/***********************************************************************
 * P r o t o t y p e s   f o r    T M L
 ***********************************************************************
 */
extern void intrapred_16x16 (Macroblock *currMB, ColorPlane pl);
// Transform function pointers

extern void copyblock_sp    (Macroblock *currMB, ColorPlane pl, int pos_mb1,int pos_mb2);
extern int  dct_chroma      (Macroblock *currMB, int uv, int cr_cbp);
extern int  dct_chroma_sp   (Macroblock *currMB, int uv,int i11);
extern int  dct_chroma_sp2  (Macroblock *currMB, int uv,int i11);



extern distblk   GetDirectCostMB           (Macroblock *currMB);
extern distblk   GetDirectCost8x8          (Macroblock *currMB, int, distblk*);

extern int   picture_coding_decision   (VideoParameters *p_Vid, Picture *picture1, Picture *picture2, int qp);

extern void no_mem_exit  (char *where);
extern int  get_mem_ACcoeff_new  (int****** cofAC, int chroma);
extern int  get_mem_ACcoeff      (VideoParameters *p_Vid, int*****);
extern int  get_mem_DCcoeff      (int****);
extern void free_mem_ACcoeff     (int****);
extern void free_mem_ACcoeff_new (int***** cofAC);
extern void free_mem_DCcoeff     (int***);

#if TRACE
extern void  trace2out(SyntaxElement *se);
extern void  trace2out_cabac(SyntaxElement *se);
#endif

extern void error(char *text, int code);
extern void set_mbaff_parameters(Macroblock  *currMB);

//============= restriction of reference frames based on the latest intra-refreshes==========
extern int  is_bipred_enabled          (InputParameters *p_Inp, int mode); 
extern void update_qp                  (Macroblock *currMB);
extern void select_plane               (VideoParameters *p_Vid, ColorPlane color_plane);
extern void select_dct                 (Macroblock *currMB);
extern void set_slice_type             (VideoParameters *p_Vid, InputParameters *p_Inp, int slice_type);
extern void free_encoder_memory        (VideoParameters *p_Vid, InputParameters *p_Inp);
extern void output_SP_coefficients     (VideoParameters *p_Vid, InputParameters *p_Inp);
extern void read_SP_coefficients       (VideoParameters *p_Vid, InputParameters *p_Inp);
extern void Init_redundant_frame       (VideoParameters *p_Vid, InputParameters *p_Inp);
extern void Set_redundant_frame        (VideoParameters *p_Vid, InputParameters *p_Inp);
extern void encode_one_redundant_frame (VideoParameters *p_Vid, InputParameters *p_Inp);


// struct with pointers to the sub-images
typedef struct 
{
  imgpel ****luma;    //!< component 0 (usually Y, X, or R)
  imgpel ****crcb[2]; //!< component 2 (usually U/V, Y/Z, or G/B)
} SubImageContainer;

// For 4:4:4 independent mode
extern void change_plane_JV      ( VideoParameters *p_Vid, int nplane );
extern void make_frame_picture_JV( VideoParameters *p_Vid );



short **PicPos;
char errortext[ET_SIZE]; //!< buffer for error message for exit with error()

#endif

