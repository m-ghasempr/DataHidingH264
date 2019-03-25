
/*!
 ************************************************************************
 *
 * \file fast_me.h
 *
 * \brief
 *   Macro definitions and global variables for fast integer pel motion 
 *   estimation and fractional pel motion estimation
 *
 * \author
 *   Main contributors: (see contributors.h for copyright, address and affiliation details)
 *    - Zhibo Chen         <chenzhibo@tsinghua.org.cn>
 *    - JianFeng Xu        <fenax@video.mdc.tsinghua.edu.cn>  
 *    - Wenfang Fu         <fwf@video.mdc.tsinghua.edu.cn>
 *    - Xiaozhong Xu 	   <xxz@video.mdc.tsinghua.edu.cn>
 *
 * \date
 *   2006.1
 ************************************************************************
 */

#ifndef _FAST_ME_H_
#define _FAST_ME_H_

#include "mbuffer.h"

#define EARLY_TERMINATION  if ((min_mcost-pred_SAD)<pred_SAD*betaFourth_2)             \
  goto fourth_2_step;                                             \
  else if((min_mcost-pred_SAD)<pred_SAD*betaFourth_1)           \
  goto fourth_1_step;    

#define SEARCH_ONE_PIXEL  if(abs(cand_x - center_x) <=search_range && abs(cand_y - center_y)<= search_range) \
    { \
    if(!McostState[cand_y-center_y+search_range][cand_x-center_x+search_range]) \
    { \
    mcost = MV_COST (lambda_factor, mvshift, cand_x, cand_y, pred_x, pred_y); \
	if(mcost<min_mcost)					\
	{\
    mcost = PartCalMad(ref_pic, orig_pic, get_ref_line,blocksize_y,blocksize_x,blocksize_x4,mcost,min_mcost,cand_x,cand_y); \
    McostState[cand_y-center_y+search_range][cand_x-center_x+search_range] = 1; \
    if (mcost < min_mcost) \
    { \
    best_x = cand_x; \
    best_y = cand_y; \
    min_mcost = mcost; \
    } \
	}\
    } \
    }
#define SEARCH_ONE_PIXEL_BIPRED  if(abs(cand_x - center2_x) <=search_range && abs(cand_y - center2_y)<= search_range) \
    { \
    if(!McostState[cand_y-center2_y+search_range][cand_x-center2_x+search_range]) \
    { \
    mcost  = MV_COST (lambda_factor, mvshift, center1_x, center1_y, pred_x1, pred_y1); \
    mcost += MV_COST (lambda_factor, mvshift, cand_x,	 cand_y,	pred_x2, pred_y2); \
	if(mcost<min_mcost)					\
	{\
    mcost  = PartCalMadBiPred(cur_pic, blocksize_y, blocksize_x, blocksize_x4, mcost,min_mcost,center1_x, center1_y, cand_x, cand_y); \
    McostState[cand_y-center2_y+search_range][cand_x-center2_x+search_range] = 1; \
    if (mcost < min_mcost) \
    { \
    best_x = cand_x; \
    best_y = cand_y; \
    min_mcost = mcost; \
    } \
	}\
    } \
    }

byte **McostState;							//state for integer pel search
byte **SearchState;							//state for fractional pel search

int ****fastme_ref_cost;                    //store SAD information needed for forward ref-frame prediction
int ***fastme_l0_cost;                      //store SAD information needed for forward median and uplayer prediction
int ***fastme_l1_cost;                      //store SAD information needed for backward median and uplayer prediction
int ***fastme_l0_cost_bipred;               //store SAD information for bipred mode
int ***fastme_l1_cost_bipred;               //store SAD information for bipred mode
int bipred_flag;                            //flag for bipred
int **fastme_best_cost;						//for multi ref early termination threshold
int pred_SAD;								// SAD prediction in use. 
int pred_MV_ref[2], pred_MV_uplayer[2];     //pred motion vector by space or temporal correlation,Median is provided

int FME_blocktype;							//blocktype for FME SetMotionVectorPredictor
int predict_point[5][2];
int SAD_a,SAD_b,SAD_c,SAD_d;
int Threshold_DSR_MB[8];					// Threshold for usage of DSR. DSR refer to JVT-Q088
//for early termination
float  Bsize[8];
float AlphaFourth_1[8];
float AlphaFourth_2[8];
byte *flag_intra;
int  flag_intra_SAD;

void DefineThreshold(void);
void DefineThresholdMB(void);
int get_mem_FME(void);
void free_mem_FME(void);

void decide_intrabk_SAD(void);
void skip_intrabk_SAD(int best_mode, int ref_max);
void setup_FME(short ref, int list, int block_y, int block_x, int blocktype, short   ******all_mv);

int                                     //  ==> minimum motion cost after search
FastIntegerPelBlockMotionSearch  (pel_t**   orig_pic,      // <--  not used
                                  short     ref,           // <--  reference frame (0... or -1 (backward))
                                  int       list,
                                  int       pic_pix_x,     // <--  absolute x-coordinate of regarded AxB block
                                  int       pic_pix_y,     // <--  absolute y-coordinate of regarded AxB block
                                  int       blocktype,     // <--  block type (1-16x16 ... 7-4x4)
                                  short     pred_mv_x,     // <--  motion vector predictor (x) in sub-pel units
                                  short     pred_mv_y,     // <--  motion vector predictor (y) in sub-pel units
                                  short*    mv_x,          //  --> motion vector (x) - in pel units
                                  short*    mv_y,          //  --> motion vector (y) - in pel units
                                  int       search_range,  // <--  1-d search range in pel units                         
                                  int       min_mcost,     // <--  minimum motion cost (cost for center or huge value)
                                  int       lambda_factor);// <--  lagrangian parameter for determining motion cost

int AddUpSADQuarter(int pic_pix_x,int pic_pix_y,int blocksize_x,int blocksize_y,
                    int cand_mv_x,int cand_mv_y, StorablePicture *ref_picture, pel_t**   orig_pic, 
                    int Mvmcost, int min_mcost,int useABT,int blocktype);

int                                                   //  ==> minimum motion cost after search
FastSubPelBlockMotionSearch (pel_t**   orig_pic,      // <--  original pixel values for the AxB block
                             short       ref,           // <--  reference frame (0... or -1 (backward))
                             int       list,
                             int       pic_pix_x,     // <--  absolute x-coordinate of regarded AxB block
                             int       pic_pix_y,     // <--  absolute y-coordinate of regarded AxB block
                             int       blocktype,     // <--  block type (1-16x16 ... 7-4x4)
                             short     pred_mv_x,     // <--  motion vector predictor (x) in sub-pel units
                             short     pred_mv_y,     // <--  motion vector predictor (y) in sub-pel units
                             short*    mv_x,          // <--> in: search center (x) / out: motion vector (x) - in pel units
                             short*    mv_y,          // <--> in: search center (y) / out: motion vector (y) - in pel units
                             int       search_pos2,   // <--  search positions for    half-pel search  (default: 9)
                             int       search_pos4,   // <--  search positions for quarter-pel search  (default: 9)
                             int       min_mcost,     // <--  minimum motion cost (cost for center or huge value)
                             int       lambda_factor, // <--  lagrangian parameter for determining motion cost
                             int  useABT);

int                                               //  ==> minimum motion cost after search
SubPelBlockMotionSearch (pel_t**   orig_pic,      // <--  original pixel values for the AxB block
                         short     ref,           // <--  reference frame (0... or -1 (backward))
                         int       list,
                         int       pic_pix_x,     // <--  absolute x-coordinate of regarded AxB block
                         int       pic_pix_y,     // <--  absolute y-coordinate of regarded AxB block
                         int       blocktype,     // <--  block type (1-16x16 ... 7-4x4)
                         int       pred_mv_x,     // <--  motion vector predictor (x) in sub-pel units
                         int       pred_mv_y,     // <--  motion vector predictor (y) in sub-pel units
                         short*    mv_x,          // <--> in: search center (x) / out: motion vector (x) - in pel units
                         short*    mv_y,          // <--> in: search center (y) / out: motion vector (y) - in pel units
                         int       search_pos2,   // <--  search positions for    half-pel search  (default: 9)
                         int       search_pos4,   // <--  search positions for quarter-pel search  (default: 9)
                         int       min_mcost,     // <--  minimum motion cost (cost for center or huge value)
                         int      lambda_factor         // <--  lagrangian parameter for determining motion cost
                         );

int                                                //  ==> minimum motion cost after search
FastBipredIntegerPelBlockMotionSearch (pel_t**   orig_pic,      // <--  original pixel values for the AxB block
                          short       ref,         // <--  reference frame (0... or -1 (backward))
                          int       list,
                          int       pic_pix_x,     // <--  absolute x-coordinate of regarded AxB block
                          int       pic_pix_y,     // <--  absolute y-coordinate of regarded AxB block
                          int       blocktype,     // <--  block type (1-16x16 ... 7-4x4)
                          short     pred_mv_x1,    // <--  motion vector predictor (x) in sub-pel units
                          short     pred_mv_y1,    // <--  motion vector predictor (y) in sub-pel units
                          short     pred_mv_x2,    // <--  motion vector predictor (x) in sub-pel units
                          short     pred_mv_y2,    // <--  motion vector predictor (y) in sub-pel units
                          short*    mv_x,          // <--> in: search center (x) / out: motion vector (x) - in pel units
                          short*    mv_y,          // <--> in: search center (y) / out: motion vector (y) - in pel units
                          short*    s_mv_x,        // <--> in: search center (x) / out: motion vector (x) - in pel units
                          short*    s_mv_y,        // <--> in: search center (y) / out: motion vector (y) - in pel units
                          int       search_range,  // <--  1-d search range in pel units
                          int       min_mcost,     // <--  minimum motion cost (cost for center or huge value)
                          int       lambda_factor // <--  lagrangian parameter for determining motion cost
									);
int compute_BiPredSad1(pel_t** cur_pic,
                             int blocksize_y,
                             int blocksize_x, 
                             int blocksize_x4,
                             int mcost,
                             int min_mcost,
                             int cand_x1, int cand_y1, 
                             int cand_x2, int cand_y2);
int compute_BiPredSad2(pel_t** cur_pic,
                             int blocksize_y,
                             int blocksize_x, 
                             int blocksize_x4,
                             int mcost,
                             int min_mcost,
                             int cand_x1, int cand_y1, 
                             int cand_x2, int cand_y2);

#endif
