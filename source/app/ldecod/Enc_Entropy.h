#pragma once
#include "global.h"
#include "Enc_Input_param.h"
#include "nalucommon.h"


#define LARGE_FILE 4*1024
typedef struct bit_stream_enc Enc_Bitstream;
typedef struct coding_state Dec_CSobj;
typedef struct encoding_environment EncodingEnvironment;
typedef EncodingEnvironment *EncodingEnvironmentPtr;
typedef struct bit_counter BitCounter;
typedef struct inp_par_enc Enc_InputParameters;
typedef struct info_8x8 Info8x8;
typedef struct enc_bi_context_type Enc_BiContextType;
typedef Enc_BiContextType *Enc_BiContextTypePtr;

void Alloc_Init_CABAC();
void Alloc_Init_CAVLC();
unsigned char* Init_Output_Buffer(char *File_name, char *Input_String, int *size, unsigned char *large);
void Handle_Large_File(char *File_name, int *size, char *large, char *Output_string, unsigned int num);
int ReadExpGlomb(char *Bit_stream, int *Offset, int *Bit_offset_to_go);
void WriteExpGlomb(char *Bit_stream, int *Offset, int *Bit_offset_to_go, int data);
SliceType Find_Slice_type(FILE *Input_File, unsigned char *fileN);
boolean ReadFlag(char *Bit_stream, int *Offset, int *Bit_offset_to_go);
void WriteFlag(char *Bit_stream, int *Offset, int *Bit_offset_to_go, Boolean data);
void Copy2End(unsigned long int S_point, FILE *INP, FILE *OUTP);
//! struct for context management

struct enc_bi_context_type
{
	unsigned long  count;
	byte state; //uint16 state;         // index into state-table CP
	unsigned char  MPS;           // Least Probable Symbol 0/1 CP  
};

typedef struct
{
	Enc_BiContextType  transform_size_contexts[NUM_TRANSFORM_SIZE_CTX];
	Enc_BiContextType  ipr_contexts[NUM_IPR_CTX];
	Enc_BiContextType  cipr_contexts[NUM_CIPR_CTX];
	Enc_BiContextType  cbp_contexts[3][NUM_CBP_CTX];
	Enc_BiContextType  bcbp_contexts[NUM_BLOCK_TYPES][NUM_BCBP_CTX];
	Enc_BiContextType  map_contexts[2][NUM_BLOCK_TYPES][NUM_MAP_CTX];
	Enc_BiContextType  last_contexts[2][NUM_BLOCK_TYPES][NUM_LAST_CTX];
	Enc_BiContextType  one_contexts[NUM_BLOCK_TYPES][NUM_ONE_CTX];
	Enc_BiContextType  abs_contexts[NUM_BLOCK_TYPES][NUM_ABS_CTX];
} Enc_TextureInfoContexts;

struct bit_counter
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
};

struct info_8x8
{
	char   mode;
	char   pdir;
	char   ref[2];
	char   bipred;
};

struct encoding_environment
{
	struct enc_video_par *p_Vid;
	unsigned int  Elow, Erange;
	unsigned int  Ebuffer;
	unsigned int  Ebits_to_go;
	unsigned int  Echunks_outstanding;
	int           Epbuf;
	unsigned char  *Ecodestrm;
	int           *Ecodestrm_len;
	int           C;
	int           E;
};

struct bit_stream_enc
{
	int     buffer_size;        //!< Buffer size      
	int     byte_pos;           //!< current position in bitstream;
	int     bits_to_go;         //!< current bitcounter

	int     stored_byte_pos;    //!< storage for position in bitstream;
	int     stored_bits_to_go;  //!< storage for bitcounter
	int     byte_pos_skip;      //!< storage for position in bitstream;
	int     bits_to_go_skip;    //!< storage for bitcounter
	int     write_flag;         //!< Bitstream contains data and needs to be written

	unsigned char    byte_buf;           //!< current buffer for last written byte
	unsigned char    stored_byte_buf;    //!< storage for buffer of last written byte
	unsigned char    byte_buf_skip;      //!< current buffer for last written byte
	unsigned char    *streamBuffer;      //!< actual buffer for written bytes

};

struct coding_state {

	// important variables of data partition array
	int                  no_part;
	Bitstream            *bitstream;
	DecodingEnvironment  *decenv;

	// contexts for binary arithmetic coding
	MotionInfoContexts   *mot_ctx;
	TextureInfoContexts  *tex_ctx;

	// bit counter
	BitCounter            bits;

	// elements of current macroblock
	short                 mvd[2][BLOCK_MULTIPLE][BLOCK_MULTIPLE][2];
	int64                 cbp_bits[3];
	int64                 *cbp_bits_8x8;
};

typedef struct macroblock_enc
{
	struct enc_slice       *p_Slice;                    //!< pointer to the current slice
	struct enc_video_par   *p_Vid;                      //!< pointer to VideoParameters
	Enc_InputParameters    *p_Inp;
	int                    mbAddrX;                    //!< current MB address
	short                  mb_type;                    //!< current MB mode type
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
	int                 DeblockCall;

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
	BlockPos *PicPos;
	short               ar_mode; //!< mb type to store adaptive rounding parameters
	short               mvd[2][BLOCK_MULTIPLE][BLOCK_MULTIPLE][2];          //!< indices correspond to [list][block_y][block_x][x,y]
	char                c_ipred_mode;      //!< chroma intra prediction mode
	char                i16mode;

	Info8x8             b8x8[4];

	imgpel              **intra4x4_pred; //[3][17]; //!< 4x4 Intra prediction samples
	imgpel              **intra4x4_pred_buf[2]; //[3][17]; //!< 4x4 Intra prediction samples
	imgpel              **intra8x8_pred; //[3][25]; //!< 8x8 Intra prediction samples
	imgpel              **intra8x8_pred_buf[2]; //[3][25]; //!< 8x8 Intra prediction samples
	imgpel              **intra16x16_pred; //[3][33]; //!< 8x8 Intra prediction samples
	imgpel              **intra16x16_pred_buf[2]; //[3][33]; //!< 8x8 Intra prediction samples

	char                IntraChromaPredModeFlag;
	byte                mb_field;
	byte                is_field_mode;
	byte                luma_transform_size_8x8_flag;
	byte                temp_transform_size_8x8_flag;
	byte                NoMbPartLessThan8x8Flag;
	byte                valid_8x8;
	byte                valid_4x4;
	byte                write_mb;
	byte                is_intra_block;

	char                DFDisableIdc;
	char                DFAlphaC0Offset;
	char                DFBetaOffset;

	int                 skip_flag;

	char                intra_pred_modes[MB_BLOCK_PARTITIONS];
	char                intra_pred_modes8x8[MB_BLOCK_PARTITIONS];           //!< four 8x8 blocks in a macroblock

	int64               cbp_blk;    //!< 1 bit set for every 4x4 block with coefs (not implemented for INTRA)
	int64               cbp_bits[3];
	int64               cbp_bits_8x8[3];

	distblk              min_rdcost;
	distblk              min_dcost;
	distblk              min_rate;
	int                  min_bits;

	short               best_mode;
	short               best_i16offset;
	char                best_i16mode;
	char                best_c_imode;
	int                 best_cbp;

	//For residual DPCM
	short               ipmode_DPCM;
	int mby, mbx;
	struct macroblock_enc   *mb_up;   //!< pointer to neighboring MB (CABAC)
	struct macroblock_enc   *mb_left; //!< pointer to neighboring MB (CABAC)
	struct macroblock_enc   *PrevMB;
} Enc_Macroblock;

typedef struct enc_video_par
{
	Enc_InputParameters *p_Inp;
	pic_parameter_set_rbsp_t *active_pps;
	seq_parameter_set_rbsp_t *active_sps;
	int cabac_encoding;
	int num_blk8x8_uv;
	int num_cdc_coeff;
	int ***nz_coeff;             //!< number of coefficients per block (CAVLC)
	ColorFormat yuv_format;
	int mb_size[MAX_PLANE][2];
	short   *intra_block;
	BlockPos *PicPos;
	Enc_Macroblock    *mb_data;                                   //!< array containing all MBs of a whole frame
	int is_v_block;
	int width, height;
	PictureStructure structure;  //!< picture structure
} Enc_VideoParameters;

typedef struct datapartition_enc
{
	struct enc_slice        *p_Slice;
	struct enc_video_par      *p_Vid;
	Enc_InputParameters     *p_Inp;   //!< pointer to the input parameters
	Enc_Bitstream           *bitstream;
	NALU_t              *nal_unit;
	EncodingEnvironment ee_cabac;
	EncodingEnvironment ee_recode;
} Enc_DataPartition;

typedef struct enc_slice
{
	struct enc_video_par    *p_Vid;   // pointer to the original video structure
	Enc_InputParameters     *p_Inp;   // pointer to the input parameters
	pic_parameter_set_rbsp_t *active_pps;
	seq_parameter_set_rbsp_t *active_sps;

	struct decoded_picture_buffer *p_Dpb;

	int                 layer_id;
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
	int                 colour_plane_id;   //!< colour plane id for 4:4:4 profile
	int                 lossless_qpprime_flag;
	int                 P444_joined;
	int                 disthres;
	int                 Transform8x8Mode;
	int                 rdoq_motion_copy;
	// RDOQ at the slice level (note this allows us to reduce passed parameters, 
	// but also could enable RDOQ control at the slice level.
	int                 UseRDOQuant;
	int                 RDOQ_QP_Num;

	char                num_ref_idx_active[2];
	int                 ref_pic_list_reordering_flag[2];
	int                 *modification_of_pic_nums_idc[2];
	int                 *abs_diff_pic_num_minus1[2];
	int                 *long_term_pic_idx[2];
#if (MVC_EXTENSION_ENABLE)
	int                 *abs_diff_view_idx_minus1[2];
#endif
	int                 view_id;
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

																 // Bit depth related information
	short               bitdepth_luma;
	short               bitdepth_chroma;
	int                 bitdepth_scale[2];
	int                 bitdepth_luma_qp_scale;
	int                 bitdepth_chroma_qp_scale;
	int                 bitdepth_lambda_scale;

	Enc_DataPartition       *partArr;     //!< array of partitions
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
	MotionVector *****all_mv;         //!< replaces local all_mv
	MotionVector ******bipred_mv;     //!< Biprediction MVs  
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
	int    norm_shift_4x4;
	int    norm_shift_8x8;

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
	struct epzs_params      *p_EPZS;

	// This should be the right location for this
	struct storable_picture **listX[6];
	char listXsize[6];

	// Some Cabac related parameters (could be put in a different structure so we can dynamically allocate them when needed)
	int  coeff[64];
	int  coeff_ctr;
	int  pos;

	int default_pic_num[2][MAX_REFERENCE_PICTURES];
	int default_view_id[2][MAX_REFERENCE_PICTURES];
} Enc_Slice;

typedef struct enc_syntaxelement
{
	int                 type;           //!< type of syntax element for data part.
	int                 value1;         //!< numerical value of syntax element
	int                 value2;         //!< for blocked symbols, e.g. run/level
	int                 len;            //!< length of code
	int                 inf;            //!< info part of UVLC code
	unsigned int        bitpattern;     //!< UVLC bitpattern
	int                 context;        //!< CABAC context

										//!< for mapping of syntaxElement to UVLC
	void(*mapping)(int value1, int value2, int* len_ptr, int* info_ptr);

} Enc_SyntaxElement;

typedef struct cbp_enc
{
	short               mb_type;                    //!< current MB mode type
	int64               cbp_blk;    //!< 1 bit set for every 4x4 block with coefs (not implemented for INTRA)
	int64               cbp_bits[3];
	int64               cbp_bits_8x8[3];
} Enc_CBP;

Enc_Macroblock *Enc_MB;
Enc_Slice *Enc_currSlice;
Enc_VideoParameters *Enc_Vid;
EncodingEnvironment *Enc_Stream;
Enc_DataPartition *Enc_dp;
unsigned char Enc_Buffer[8 * 1024 * 1024];
int Enc_Buffer_len;
Enc_Bitstream *Enc_VStream;
typedef struct {
	unsigned char BuffFrame[1024 * 1024];
	unsigned char byteBuffFrame;
	int pBuffFrame;
	int bits2Go;
}OutputBuffer;
OutputBuffer outputBuffer;
int writepoint;
int FirstMB;
Enc_CBP Enc_MBs[512 * 256];
