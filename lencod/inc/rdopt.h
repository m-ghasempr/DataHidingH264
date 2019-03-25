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

extern int   *****cofAC8x8CbCr;
extern int   ****cofAC, ****cofAC8x8;        // [8x8block][4x4block][level/run][scan_pos]
extern int   ***cofDC;                       // [yuv][level/run][scan_pos]
extern int   **cofAC4x4; // [level/run][scan_pos]
extern int   cbp, cbp8x8, cnt_nonz_8x8;
extern int   cbp_blk8x8;
extern char  **l0_refframe, **l1_refframe;
extern short b8mode[4], b8pdir[4];

//CSptr cs_mb, cs_b8, cs_cm, cs_ib8, cs_ib4;
extern int   best_c_imode;
extern int   best_i16offset;
extern short best_mode;

//mixed transform sizes definitions
extern int   luma_transform_size_8x8_flag;

extern int   ****cofAC8x8ts[3];        // [plane][8x8block][4x4block][level/run][scan_pos]
extern int   **cofAC4x4CbCr[2];

extern int64    cbp_blk8_8x8ts;
extern int      cbp8_8x8ts;
extern int      cost8_8x8ts;
extern int      cnt_nonz8_8x8ts;

// adaptive langrangian parameters
extern double mb16x16_cost;
extern double lambda_mf_factor;

int (*Mode_Decision_for_4x4IntraBlocks) (Slice *currSlice, Macroblock *currMB, int  b8,  int  b4,  double  lambda,  double*  min_cost, int cr_cbp[3], int is_cavlc);
double RDCost_for_4x4IntraBlocks (Slice *currSlice, Macroblock *currMB, int* nonzero, int b8, int b4, int ipmode, double lambda, int mostProbableMode, int c_nzCbCr[3], int is_cavlc);
int valid_intra_mode(int ipmode);
void compute_comp_cost(imgpel **cur_img, imgpel **prd_img, int pic_opix_x, int *cost);
void generate_pred_error(imgpel **cur_img, imgpel **prd_img, imgpel **cur_prd, 
                         int **mb_rres, int pic_opix_x, int block_x);
                         
//============= rate-distortion optimization ===================
void  clear_rdopt (InputParameters *params);
void  init_rdopt  (InputParameters *params);

void copy_4x4block(imgpel **oblock, imgpel **iblock, int o_xoffset, int i_xoffset);

void   update_qp_cbp     (Macroblock *currMB, short best_mode);
void   update_qp_cbp_tmp (Macroblock *currMB, int cbp, int best_mode);


#endif

