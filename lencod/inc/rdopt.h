/*!
 ***************************************************************************
 * \file
 *    rdopt.h
 *
 * \author
 *    Alexis Michael Tourapis
 *
 * \date
 *    2 January 2008
 *
 * \brief
 *    Headerfile for RDO
 **************************************************************************
 */

#ifndef _RDO_H_
#define _RDO_H_

extern int   **bestInterFAdjust4x4, **bestIntraFAdjust4x4;
extern int   **bestInterFAdjust8x8, **bestIntraFAdjust8x8;
extern int   ***bestInterFAdjust4x4Cr, ***bestIntraFAdjust4x4Cr;
extern int   ***bestInterFAdjust8x8Cr, ***bestIntraFAdjust8x8Cr;
extern int   **fadjust8x8, **fadjust4x4, ***fadjust4x4Cr, ***fadjust8x8Cr;

extern int   ****cofAC, ****cofAC8x8;        // [8x8block][4x4block][level/run][scan_pos]
extern int   ***cofDC;                       // [yuv][level/run][scan_pos]
extern int   **cofAC4x4, ****cofAC4x4intern; // [level/run][scan_pos]
extern int   cbp, cbp8x8, cnt_nonz_8x8;
extern int   cbp_blk8x8;
extern char  l0_refframe[4][4], l1_refframe[4][4];
extern short b8mode[4], b8pdir[4];
extern short best8x8mode [4];                // [block]
extern char  best8x8pdir  [MAXMODE][4];       // [mode][block]
extern char  best8x8l0ref [MAXMODE][4];       // [mode][block]
extern char  best8x8l1ref [MAXMODE][4];       // [mode][block]
extern char  best8x8ref[2][MAXMODE][4];       // [mode][block]

//CSptr cs_mb, cs_b8, cs_cm, cs_ib8, cs_ib4;
extern int   best_c_imode;
extern int   best_i16offset;
extern short best_mode;
extern short  bi_pred_me;

//mixed transform sizes definitions
extern int   luma_transform_size_8x8_flag;

extern short all_mv8x8[2][2][4][4][2];       //[8x8_data/temp_data][LIST][block_x][block_y][MVx/MVy]
extern short pred_mv8x8[2][2][4][4][2];

extern int   ****cofAC8x8ts[3];        // [plane][8x8block][4x4block][level/run][scan_pos]
extern int   ****cofAC8x8CbCr[2];
extern int   **cofAC4x4CbCr[2];
extern int   ****cofAC4x4CbCrintern[2];

extern int64    cbp_blk8_8x8ts;
extern int      cbp8_8x8ts;
extern int      cost8_8x8ts;
extern int      cnt_nonz8_8x8ts;

// adaptive langrangian parameters
extern double mb16x16_cost;
extern double lambda_mf_factor;

int (*Mode_Decision_for_4x4IntraBlocks) (Macroblock *currMB, int  b8,  int  b4,  double  lambda,  double*  min_cost, int cr_cbp[3]);
double RDCost_for_4x4IntraBlocks (Macroblock *currMB, int* nonzero, int b8, int b4, int ipmode, double lambda, int mostProbableMode, int c_nzCbCr[3]);
int valid_intra_mode(int ipmode);
void compute_comp_cost(imgpel **cur_img, imgpel prd_img[16][16], int pic_opix_x, int *cost);
void generate_pred_error(imgpel **cur_img, imgpel prd_img[16][16], imgpel cur_prd[16][16], 
                         int m7[16][16], int pic_opix_x, int block_x);


void store_adaptive_rounding (int block_y, int block_x);
void update_adaptive_rounding(int block_y, int block_x);

#endif

