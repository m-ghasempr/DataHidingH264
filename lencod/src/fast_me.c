          
/*!
 ************************************************************************
 *
 * \file fast_me.c
 *
 * \brief
 *   Fast integer pel motion estimation and fractional pel motion estimation
 *   algorithms are described in this file.
 *   1. get_mem_FME() and free_mem_FME() are functions for allocation and release
 *      of memories about motion estimation
 *   2. FME_BlockMotionSearch() is the function for fast integer pel motion 
 *      estimation and fractional pel motion estimation
 *   3. DefineThreshold() defined thresholds for early termination
 * \author 
 *    Main contributors: (see contributors.h for copyright, address and affiliation details)
 *    - Zhibo Chen         <chenzhibo@tsinghua.org.cn>
 *    - JianFeng Xu        <fenax@video.mdc.tsinghua.edu.cn>  
 *    - Wenfang Fu         <fwf@video.mdc.tsinghua.edu.cn>
 *	  - Xiaozhong Xu  	   <xxz@video.mdc.tsinghua.edu.cn>
 * \date    
 *    2006.1
 ************************************************************************
 */

#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "limits.h"
#include "memalloc.h"
#include "fast_me.h"
#include "refbuf.h"
#include "mb_access.h"
#include "image.h"

#define Q_BITS          15
#define MIN_IMG_WIDTH	176
extern  unsigned int*   byte_abs;
extern  int*   mvbits;
extern  short*   spiral_search_x;
extern  short*   spiral_search_y;


static pel_t *(*get_line) (pel_t**, int, int, int, int);
static pel_t*  ref_pic;
static pel_t *(*get_ref_line)(int, pel_t*, int, int, int, int);
static const int Diamond_x[4] = {-1, 0, 1, 0};
static const int Diamond_y[4] = {0, 1, 0, -1};
static const int Hexagon_x[6] = {2, 1, -1, -2, -1, 1};
static const int Hexagon_y[6] = {0, -2, -2, 0,  2, 2};
static const int Big_Hexagon_x[16] = {0,-2, -4,-4,-4, -4, -4, -2,  0,  2,  4,  4, 4, 4, 4, 2};
static const int Big_Hexagon_y[16] = {4, 3, 2,  1, 0, -1, -2, -3, -4, -3, -2, -1, 0, 1, 2, 3};
static const int   Quater_x[8] = {2, 1, 0, -1, -2, -1,  0,  1};
static const int   Quater_y[8] = {0, 1, 2,  1,  0, -1, -2, -1};

// for bipred mode
static int height,width;
static int pred_MV_ref_flag;  
static short weightSpic, weightRpic, offsetBi;
int (*PartCalMadBiPred)(pel_t **, int, int, int, int, int, int, int, int, int);
static pel_t *(*get_ref_line1)(int, pel_t *, int, int, int, int);
static pel_t *(*get_ref_line2)(int, pel_t *, int, int, int, int);
static pel_t *ref1_pic;
static pel_t *ref2_pic;

static const int   Multi_Ref_Thd[8]   = {0,  300,  120,  120,  60,  30,   30,  15};                        
static const int   Big_Hexagon_Thd[8] = {0, 3000, 1500, 1500, 800, 400,  400, 200};                        
static const int   Median_Pred_Thd[8] = {0,  750,  350,  350, 170,  80,   80,  40};
static const int   Threshold_DSR[8]   = {0, 2200, 1000, 1000, 500, 250,  250, 120};                                      

static int Median_Pred_Thd_MB[8];
static int Big_Hexagon_Thd_MB[8];
static int Multi_Ref_Thd_MB[8];


static const int quant_coef[6][4][4] = {
  {{13107, 8066,13107, 8066},{ 8066, 5243, 8066, 5243},{13107, 8066,13107, 8066},{ 8066, 5243, 8066, 5243}},
  {{11916, 7490,11916, 7490},{ 7490, 4660, 7490, 4660},{11916, 7490,11916, 7490},{ 7490, 4660, 7490, 4660}},
  {{10082, 6554,10082, 6554},{ 6554, 4194, 6554, 4194},{10082, 6554,10082, 6554},{ 6554, 4194, 6554, 4194}},
  {{ 9362, 5825, 9362, 5825},{ 5825, 3647, 5825, 3647},{ 9362, 5825, 9362, 5825},{ 5825, 3647, 5825, 3647}},
  {{ 8192, 5243, 8192, 5243},{ 5243, 3355, 5243, 3355},{ 8192, 5243, 8192, 5243},{ 5243, 3355, 5243, 3355}},
  {{ 7282, 4559, 7282, 4559},{ 4559, 2893, 4559, 2893},{ 7282, 4559, 7282, 4559},{ 4559, 2893, 4559, 2893}}
};


void DefineThreshold()
{
  AlphaFourth_1[1] = 0.01f;
  AlphaFourth_1[2] = 0.01f;
  AlphaFourth_1[3] = 0.01f;
  AlphaFourth_1[4] = 0.02f;
  AlphaFourth_1[5] = 0.03f;
  AlphaFourth_1[6] = 0.03f;
  AlphaFourth_1[7] = 0.04f;

  AlphaFourth_2[1] = 0.06f;
  AlphaFourth_2[2] = 0.07f;
  AlphaFourth_2[3] = 0.07f;
  AlphaFourth_2[4] = 0.08f;
  AlphaFourth_2[5] = 0.12f;
  AlphaFourth_2[6] = 0.11f;
  AlphaFourth_2[7] = 0.15f;

  DefineThresholdMB();
  return;
}

void DefineThresholdMB()
{
  int gb_qp_per    = (input->qpN-MIN_QP)/6;
  int gb_qp_rem    = (input->qpN-MIN_QP)%6;
  
  int gb_q_bits    = Q_BITS+gb_qp_per;
  int gb_qp_const,Thresh4x4;

  float Quantize_step;
  int i;
// scale factor: defined for different image sizes
  float scale_factor = (float)((1-input->FMEScale*0.1)+input->FMEScale*0.1*(img->width/MIN_IMG_WIDTH));
// QP factor: defined for different quantization steps
  float QP_factor = (float)((1.0-0.90*(input->qpN/51.0f)));

  gb_qp_const=(1<<gb_q_bits)/6;
  Thresh4x4 =   ((1<<gb_q_bits) - gb_qp_const)/quant_coef[gb_qp_rem][0][0];
  Quantize_step = Thresh4x4/(4*5.61f)*2.0f*scale_factor;
  Bsize[7]=(16*16)*Quantize_step;

  Bsize[6]=Bsize[7]*4;
  Bsize[5]=Bsize[7]*4;
  Bsize[4]=Bsize[5]*4;
  Bsize[3]=Bsize[4]*4;
  Bsize[2]=Bsize[4]*4;
  Bsize[1]=Bsize[2]*4;

  for(i=1;i<8;i++)
  {
//ET_Thd1: early termination after median prediction
	Median_Pred_Thd_MB[i]  = (int) (Median_Pred_Thd[i]* scale_factor*QP_factor); 
//ET_thd2: early termination after every circle of 16 points Big-Hex Search
    Big_Hexagon_Thd_MB[i]  = (int) (Big_Hexagon_Thd[i]* scale_factor*QP_factor); 
//threshold for multi ref case
	Multi_Ref_Thd_MB[i]    = (int) (Multi_Ref_Thd[i]  * scale_factor*QP_factor); 
//threshold for usage of DSR technique. DSR ref to JVT-R088
	Threshold_DSR_MB[i]    = (int) (Threshold_DSR[i]  * scale_factor*QP_factor); 
  }

}


int get_mem_FME()
{
  int memory_size = 0;
  if (NULL==(flag_intra = calloc ((img->width>>4)+1,sizeof(byte)))) no_mem_exit("get_mem_FME: flag_intra"); //fwf 20050330

  memory_size += get_mem2D(&McostState, 2*input->search_range+1, 2*input->search_range+1); 
  memory_size += get_mem4Dint(&(fastme_ref_cost), img->max_num_references, 9, 4, 4);
  memory_size += get_mem3Dint(&(fastme_l0_cost), 9, img->height/4, img->width/4);
  memory_size += get_mem3Dint(&(fastme_l1_cost), 9, img->height/4, img->width/4);
  memory_size += get_mem2D(&SearchState,7,7);
  memory_size += get_mem2Dint(&(fastme_best_cost), 7, img->width/4);
  if(input->BiPredMotionEstimation == 1)//memory allocation for bipred mode
  {
	  memory_size += get_mem3Dint(&(fastme_l0_cost_bipred), 9, img->height/4, img->width/4);//for bipred
	  memory_size += get_mem3Dint(&(fastme_l1_cost_bipred), 9, img->height/4, img->width/4);//for bipred
  }
  
  return memory_size;
}


void free_mem_FME()
{
  free_mem2D(McostState);
  free_mem4Dint(fastme_ref_cost, img->max_num_references, 9);
  free_mem3Dint(fastme_l0_cost, 9);
  free_mem3Dint(fastme_l1_cost, 9);
  free_mem2D(SearchState);
  free_mem2Dint(fastme_best_cost);
  free (flag_intra);
  if(input->BiPredMotionEstimation == 1)
  {
	  free_mem3Dint(fastme_l0_cost_bipred, 9);//for bipred
	  free_mem3Dint(fastme_l1_cost_bipred, 9);//for bipred
  }

}


int PartCalMad(pel_t *ref_pic,pel_t** orig_pic,pel_t *(*get_ref_line)(int, pel_t*, int, int, int, int), int blocksize_y,int blocksize_x, int blocksize_x4,int mcost,int min_mcost,int cand_x,int cand_y)
{
  int y,x4;
  int height=((img->MbaffFrameFlag)&&(img->mb_data[img->current_mb_nr].mb_field))?img->height/2:img->height;
  pel_t *orig_line, *ref_line;
  for (y=0; y<blocksize_y; y++)
    {
    ref_line  = get_ref_line (blocksize_x, ref_pic, cand_y+y, cand_x, height, img->width);//2004.3.3
    orig_line = orig_pic [y];
    
    for (x4=0; x4<blocksize_x4; x4++)
    {
      mcost += byte_abs[ *orig_line++ - *ref_line++ ];
      mcost += byte_abs[ *orig_line++ - *ref_line++ ];
      mcost += byte_abs[ *orig_line++ - *ref_line++ ];
      mcost += byte_abs[ *orig_line++ - *ref_line++ ];
    }
    if (mcost >= min_mcost)
    {
      break;
    }
    }
    return mcost;
}
/*!
************************************************************************
* \brief
*    BiPred Mode SAD computation (without weights)
************************************************************************
*/
int PartCalMadBiPred1(pel_t** cur_pic,
                             int blocksize_y,
                             int blocksize_x, 
                             int blocksize_x4,
                             int mcost,
                             int min_mcost,
                             int cand_x1, int cand_y1, 
                             int cand_x2, int cand_y2)
{
  pel_t *cur_line, *ref1_line, *ref2_line;
  int bi_diff; 
  int y,x4;  
  
  for (y = 0; y < blocksize_y; y++)
  {
    ref2_line = get_ref_line2 (blocksize_x, ref2_pic, cand_y2 + y, cand_x2, height, width);
    ref1_line = get_ref_line1 (blocksize_x, ref1_pic, cand_y1 + y, cand_x1, height, width);
    cur_line = cur_pic [y];
    
    for (x4 = 0; x4 < blocksize_x4; x4++)
    {         
      bi_diff = (*cur_line++) - ((*ref1_line++ + *ref2_line++)>>1);
      mcost += byte_abs[bi_diff];
      bi_diff = (*cur_line++) - ((*ref1_line++ + *ref2_line++)>>1);
      mcost += byte_abs[bi_diff];
      bi_diff = (*cur_line++) - ((*ref1_line++ + *ref2_line++)>>1);
      mcost += byte_abs[bi_diff];
      bi_diff = (*cur_line++) - ((*ref1_line++ + *ref2_line++)>>1);
      mcost += byte_abs[bi_diff];
    }        
    
    if (mcost >= min_mcost) break;
  }
  return mcost;
}


/*!
************************************************************************
* \brief
*    BiPred Mode SAD computation (with weights)
************************************************************************
*/
int PartCalMadBiPred2(pel_t** cur_pic,
                             int blocksize_y,
                             int blocksize_x, 
                             int blocksize_x4,
                             int mcost,
                             int min_mcost,
                             int cand_x1, int cand_y1, 
                             int cand_x2, int cand_y2)
{
  pel_t *cur_line, *ref1_line, *ref2_line;
  int bi_diff; 
  int denom = luma_log_weight_denom + 1;
  int lround = 2 * wp_luma_round;
  int y,x4;  
  int weightedpel, pixel1, pixel2;
  for (y=0; y<blocksize_y; y++)
  {
    ref2_line  = get_ref_line2 (blocksize_x, ref2_pic, cand_y2 + y, cand_x2, height, width);
    ref1_line  = get_ref_line1 (blocksize_x, ref1_pic, cand_y1 + y, cand_x1, height, width);
    cur_line = cur_pic [y];
    
    for (x4 = 0; x4 < blocksize_x4; x4++)
    { 
      pixel1 = weightSpic * (*ref1_line++);
      pixel2 = weightRpic * (*ref2_line++);
      weightedpel =  Clip3 (0, img->max_imgpel_value ,((pixel1 + pixel2 + lround) >> denom) + offsetBi);
      bi_diff = (*cur_line++)  - weightedpel;
      mcost += byte_abs[bi_diff];
      
      pixel1 = weightSpic * (*ref1_line++);
      pixel2 = weightRpic * (*ref2_line++);
      weightedpel =  Clip3 (0, img->max_imgpel_value ,((pixel1 + pixel2 + lround) >> denom) + offsetBi);
      bi_diff = (*cur_line++)  - weightedpel;
      mcost += byte_abs[bi_diff];
      
      pixel1 = weightSpic * (*ref1_line++);
      pixel2 = weightRpic * (*ref2_line++);
      weightedpel =  Clip3 (0, img->max_imgpel_value ,((pixel1 + pixel2 + lround) >> denom) + offsetBi);                     
      bi_diff = (*cur_line++)  - weightedpel;
      mcost += byte_abs[bi_diff];
      
      pixel1 = weightSpic * (*ref1_line++);
      pixel2 = weightRpic * (*ref2_line++);
      weightedpel =  Clip3 (0, img->max_imgpel_value ,((pixel1 + pixel2 + lround) >> denom) + offsetBi);
      bi_diff = (*cur_line++)  - weightedpel;
      mcost += byte_abs[bi_diff];
      if (mcost >= min_mcost) break;
    }    
    
    if (mcost >= min_mcost) break;
  }
  return mcost;
}

/*!
 ************************************************************************
 * \brief
 *    FastIntegerPelBlockMotionSearch: fast pixel block motion search 
 *    this algrithm is called UMHexagonS(see JVT-D016),which includes 
 *    four steps with different kinds of search patterns
 * \par Input:
 * pel_t**   orig_pic,     // <--  original picture
 * int       ref,          // <--  reference frame (0... or -1 (backward))
 * int       pic_pix_x,    // <--  absolute x-coordinate of regarded AxB block
 * int       pic_pix_y,    // <--  absolute y-coordinate of regarded AxB block
 * int       blocktype,    // <--  block type (1-16x16 ... 7-4x4)
 * int       pred_mv_x,    // <--  motion vector predictor (x) in sub-pel units
 * int       pred_mv_y,    // <--  motion vector predictor (y) in sub-pel units
 * int*      mv_x,         //  --> motion vector (x) - in pel units
 * int*      mv_y,         //  --> motion vector (y) - in pel units
 * int       search_range, // <--  1-d search range in pel units                         
 * int       min_mcost,    // <--  minimum motion cost (cost for center or huge value)
 * int       lambda_factor // <--  lagrangian parameter for determining motion cost
 * \par
 * Two macro definitions defined in this program:
 * 1. EARLY_TERMINATION: early termination algrithm, refer to JVT-D016.doc
 * 2. SEARCH_ONE_PIXEL: search one pixel in search range
 * \author
 *   Main contributors: (see contributors.h for copyright, address and affiliation details)
 *   - Zhibo Chen         <chenzhibo@tsinghua.org.cn>
 *   - JianFeng Xu        <fenax@video.mdc.tsinghua.edu.cn>
 *   - Xiaozhong Xu       <xxz@video.mdc.tsinghua.edu.cn>
 * \date   :
 *   2006.1
 ************************************************************************
 */
int                                     //  ==> minimum motion cost after search
FastIntegerPelBlockMotionSearch  (pel_t**   orig_pic,     // <--  not used
                                  short     ref,          // <--  reference frame (0... or -1 (backward))
                                  int       list,
                                  int       pic_pix_x,    // <--  absolute x-coordinate of regarded AxB block
                                  int       pic_pix_y,    // <--  absolute y-coordinate of regarded AxB block
                                  int       blocktype,    // <--  block type (1-16x16 ... 7-4x4)
                                  short     pred_mv_x,    // <--  motion vector predictor (x) in sub-pel units
                                  short     pred_mv_y,    // <--  motion vector predictor (y) in sub-pel units
                                  short*    mv_x,         //  --> motion vector (x) - in pel units
                                  short*    mv_y,         //  --> motion vector (y) - in pel units
                                  int       search_range, // <--  1-d search range in pel units                         
                                  int       min_mcost,    // <--  minimum motion cost (cost for center or huge value)
                                  int       lambda_factor)       // <--  lagrangian parameter for determining motion cost
{
  int   pos, cand_x, cand_y,  mcost;
  int   list_offset   = ((img->MbaffFrameFlag)&&(img->mb_data[img->current_mb_nr].mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;

  int   mvshift       = 2;                  // motion vector shift for getting sub-pel units
  int   blocksize_y   = input->blc_size[blocktype][1];            // vertical block size
  int   blocksize_x   = input->blc_size[blocktype][0];            // horizontal block size
  int   blocksize_x4  = blocksize_x >> 2;                         // horizontal block size in 4-pel units
  int   pred_x        = (pic_pix_x << mvshift) + pred_mv_x;       // predicted position x (in sub-pel units)
  int   pred_y        = (pic_pix_y << mvshift) + pred_mv_y;       // predicted position y (in sub-pel units)
  int   center_x      = pic_pix_x + *mv_x;                        // center position x (in pel units)
  int   center_y      = pic_pix_y + *mv_y;                        // center position y (in pel units)
  int   best_x = 0, best_y = 0;
  int   search_step,iYMinNow, iXMinNow;
  int   i,m,j; 
  float betaFourth_1,betaFourth_2;
  int	temp_Big_Hexagon_x[16];//  temp for Big_Hexagon_x;
  int	temp_Big_Hexagon_y[16];//  temp for Big_Hexagon_y; 
  short mb_x = pic_pix_x - img->opix_x; 
  short mb_y = pic_pix_y - img->opix_y;
  short pic_pix_x2 = pic_pix_x >> 2;
  short block_x = (mb_x >> 2);
  short block_y = (mb_y >> 2);
  int ET_Thred = Median_Pred_Thd_MB[blocktype];//ET threshold in use
  int   *SAD_prediction = fastme_best_cost[blocktype-1];//multi ref SAD prediction
  //===== Use weighted Reference for ME ====

  int  apply_weights = ( (active_pps->weighted_pred_flag && (img->type == P_SLICE || img->type == SP_SLICE)) ||
                         (active_pps->weighted_bipred_idc && (img->type == B_SLICE)));  
  height=((img->MbaffFrameFlag)&&(img->mb_data[img->current_mb_nr].mb_field))?img->height/2:img->height;
  ref_pic = (apply_weights && input->UseWeightedReferenceME) ? listX[list+list_offset][ref]->imgY_11_w : listX[list+list_offset][ref]->imgY_11;

  

  //===== set function for getting reference picture lines =====
  if ((center_x > search_range) && (center_x < img->width -1-search_range-blocksize_x) &&
    (center_y > search_range) && (center_y < height-1-search_range-blocksize_y)   )
  {
    get_ref_line = FastLineX;
  }
  else
  {
    get_ref_line = UMVLineX;
  }
  
  //////allocate memory for search state//////////////////////////
  memset(McostState[0],0,(2*input->search_range+1)*(2*input->search_range+1));


  //check the center median predictor
  cand_x = center_x ;
  cand_y = center_y ;
  mcost = MV_COST (lambda_factor, mvshift, cand_x, cand_y, pred_x, pred_y);
  mcost = PartCalMad(ref_pic, orig_pic, get_ref_line,blocksize_y,blocksize_x,blocksize_x4,mcost,min_mcost,cand_x,cand_y);
  McostState[search_range][search_range] = 1;
  if (mcost < min_mcost)
  {
    min_mcost = mcost;
    best_x = cand_x;
    best_y = cand_y;
  }

  iXMinNow = best_x;
  iYMinNow = best_y;
  for (m = 0; m < 4; m++)
  {   
    cand_x = iXMinNow + Diamond_x[m];
    cand_y = iYMinNow + Diamond_y[m];   
    SEARCH_ONE_PIXEL
  } 

  if(center_x != pic_pix_x || center_y != pic_pix_y)
  {
    cand_x = pic_pix_x ;
    cand_y = pic_pix_y ;
    SEARCH_ONE_PIXEL

    iXMinNow = best_x;
    iYMinNow = best_y;
    for (m = 0; m < 4; m++)
    {   
      cand_x = iXMinNow + Diamond_x[m];
      cand_y = iYMinNow + Diamond_y[m];   
      SEARCH_ONE_PIXEL
    } 
  }
/***********************************init process*************************/
//for multi ref
  if(ref>0 && img->structure == FRAME  && min_mcost > ET_Thred && SAD_prediction[pic_pix_x2]<Multi_Ref_Thd_MB[blocktype])
	goto terminate_step;

//ET_Thd1: early termination for low motion case
  if( min_mcost < ET_Thred)	
  {
    goto terminate_step;
  }
  else // hybrid search for main search loop
  {
/****************************(MV and SAD prediction)********************************/
    setup_FME(ref, list, block_y, block_x, blocktype, img->all_mv );
	ET_Thred = Big_Hexagon_Thd_MB[blocktype];  // ET_Thd2: early termination Threshold for strong motion



	// Threshold defined for EARLY_TERMINATION 
	  if (pred_SAD == 0) 
	  {
		betaFourth_1=0;
		betaFourth_2=0;
	  }
	  else
	  {
		betaFourth_1 = Bsize[blocktype]/(pred_SAD*pred_SAD)-AlphaFourth_1[blocktype];
		betaFourth_2 = Bsize[blocktype]/(pred_SAD*pred_SAD)-AlphaFourth_2[blocktype];
 
	  }  
/*********************************************end of init ***********************************************/
  }  
// first_step: initial start point prediction 

    if(blocktype>1)
  {
    cand_x = pic_pix_x + (pred_MV_uplayer[0]/4);
    cand_y = pic_pix_y + (pred_MV_uplayer[1]/4);
    SEARCH_ONE_PIXEL
  } 


  //prediciton using mV of last ref moiton vector
 if(pred_MV_ref_flag == 1)	  			//Notes: for interlace case, ref==1 should be added
  {
      cand_x = pic_pix_x + (pred_MV_ref[0]/4);
      cand_y = pic_pix_y + (pred_MV_ref[1]/4);
      SEARCH_ONE_PIXEL
  }
  //small local search
  iXMinNow = best_x;
  iYMinNow = best_y;
  for (m = 0; m < 4; m++)
  {   
    cand_x = iXMinNow + Diamond_x[m];
    cand_y = iYMinNow + Diamond_y[m];   
    SEARCH_ONE_PIXEL
  } 

  //early termination alogrithm, refer to JVT-G016
    EARLY_TERMINATION
  
  if(blocktype>6)
    goto fourth_1_step;
  else
    goto sec_step;
  
sec_step: //Unsymmetrical-cross search 
  iXMinNow = best_x;
  iYMinNow = best_y;
  
  for(i = 1; i < search_range; i+=2)
  {
    search_step = i;
    cand_x = iXMinNow + search_step;
    cand_y = iYMinNow ;
    SEARCH_ONE_PIXEL    
    cand_x = iXMinNow - search_step;
    cand_y = iYMinNow ;
    SEARCH_ONE_PIXEL
  }
  for(i = 1; i < (search_range/2);i+=2)
  {
    search_step = i;
    cand_x = iXMinNow ;
    cand_y = iYMinNow + search_step;
    SEARCH_ONE_PIXEL
    cand_x = iXMinNow ;
    cand_y = iYMinNow - search_step;
    SEARCH_ONE_PIXEL
  }


  //early termination alogrithm, refer to JVT-G016
    EARLY_TERMINATION
  
  iXMinNow = best_x;
  iYMinNow = best_y;

//third_step:    // Uneven Multi-Hexagon-grid Search 
//sub step 1: 5x5 squre search
  for(pos=1;pos<25;pos++)
  {
    cand_x = iXMinNow + spiral_search_x[pos];
    cand_y = iYMinNow + spiral_search_y[pos];
    SEARCH_ONE_PIXEL
  }

  //early termination alogrithm, refer to JVT-G016
  EARLY_TERMINATION

//sub step 2:  Multi-Hexagon-grid search
  memcpy(temp_Big_Hexagon_x,Big_Hexagon_x,64);
  memcpy(temp_Big_Hexagon_y,Big_Hexagon_y,64);      
  for(i=1;i<=(search_range/4); i++)
  {

    for (m = 0; m < 16; m++)
    {
      cand_x = iXMinNow + temp_Big_Hexagon_x[m];
      cand_y = iYMinNow + temp_Big_Hexagon_y[m];
	  temp_Big_Hexagon_x[m] += Big_Hexagon_x[m];
	  temp_Big_Hexagon_y[m] += Big_Hexagon_y[m];	

      SEARCH_ONE_PIXEL
    }
// ET_Thd2: early termination Threshold for strong motion
	if(min_mcost < ET_Thred)
	{
	  goto terminate_step;
	}
  }


//fourth_step:  //Extended Hexagon-based Search
// the fourth step with a small search pattern
fourth_1_step:  //sub step 1: small Hexagon search
      for(i=0; i < search_range; i++) //change into 1/2
      {
        iXMinNow = best_x;
        iYMinNow = best_y;
        for (m = 0; m < 6; m++)
        {   
          cand_x = iXMinNow + Hexagon_x[m];
          cand_y = iYMinNow + Hexagon_y[m];   
          SEARCH_ONE_PIXEL

        } 
        if(best_x == iXMinNow && best_y == iYMinNow)
            break;
      }
fourth_2_step: //sub step 2: small Diamond search

      for(i = 0; i < search_range; i++) //change into 1/2
      {
	    iXMinNow = best_x;
	    iYMinNow = best_y;
        for (m = 0; m < 4; m++)
        {   
          cand_x = iXMinNow + Diamond_x[m];
          cand_y = iYMinNow + Diamond_y[m];   
          SEARCH_ONE_PIXEL

        } 
        if(best_x == iXMinNow && best_y == iYMinNow)
            break;
      }

terminate_step:
	  
// store SAD infomation for prediction	  
    //FAST MOTION ESTIMATION. ZHIBO CHEN 2003.3
	  for (i=0; i < (blocksize_x>>2); i++)
	  {
	    for (j=0; j < (blocksize_y>>2); j++)
		{
		  if(list == 0) 
		  {
		    fastme_ref_cost[ref][blocktype][block_y+j][block_x+i] = min_mcost;
		    if (ref==0)
			  fastme_l0_cost[blocktype][(img->pix_y>>2)+block_y+j][(img->pix_x>>2)+block_x+i] = min_mcost;
		  }
		  else
		  {
		    fastme_l1_cost[blocktype][(img->pix_y>>2)+block_y+j][(img->pix_x>>2)+block_x+i] = min_mcost;
		  }
		}
	  }
//for multi ref SAD prediction
	  if ((ref==0) || (SAD_prediction[pic_pix_x2] > min_mcost))
		SAD_prediction[pic_pix_x2] = min_mcost;

      *mv_x = best_x - pic_pix_x;
      *mv_y = best_y - pic_pix_y; 
      return min_mcost;
}


  /*!
 ************************************************************************
 * \brief
 * Functions for fast fractional pel motion estimation.
 * 1. int AddUpSADQuarter() returns SADT of a fractiona pel MV
 * 2. int FastSubPelBlockMotionSearch () proceed the fast fractional pel ME
 * \authors  
 *    Zhibo Chen
 *    Dept.of EE, Tsinghua Univ.
 * \date 
 *    2003.4
 ************************************************************************
 */
int AddUpSADQuarter(int pic_pix_x,int pic_pix_y,int blocksize_x,int blocksize_y,
                    int cand_mv_x,int cand_mv_y, StorablePicture *ref_picture, pel_t**   orig_pic, 
                    int Mvmcost, int min_mcost,int useABT, int blocktype)
{

  int j, i, k;  
  int diff[16], *d; 
  int mcost = Mvmcost;
  int c_diff[MB_PIXELS];
  int y_offset, ypels =(128 - 64 * (blocktype == 3));
  int ry0, ry4, ry8, ry12;
  int y0, y1, y2, y3;
  int x0, x1, x2, x3;
  int abort_search, rx0; 
  int img_width  = ((ref_picture->size_x + 2*IMG_PAD_SIZE - 1)<<2);
  int img_height = ((ref_picture->size_y + 2*IMG_PAD_SIZE - 1)<<2);

  //===== Use weighted Reference for ME ====
  pel_t **ref_pic;      
  pel_t *ref_line;
  pel_t *orig_line;
  int  apply_weights = ( (active_pps->weighted_pred_flag && (img->type == P_SLICE || img->type == SP_SLICE)) ||
                         (active_pps->weighted_bipred_idc && (img->type == B_SLICE)));  
  
  if (apply_weights && input->UseWeightedReferenceME)
  {
    ref_pic = ref_picture->imgY_ups_w;
  }
  else
    ref_pic = ref_picture->imgY_ups;
  ///////////////////////////////////////////

  
  for (y0=0, abort_search=0; y0<blocksize_y && !abort_search; y0+=4)
  {
    y_offset = (y0>7)*ypels;
    ry0  = (y0<<2) + cand_mv_y;
    ry4  = ry0 + 4;
    ry8  = ry4 + 4;
    ry12 = ry8 + 4;
    y1 = y0 + 1;
    y2 = y1 + 1;
    y3 = y2 + 1;


    for (x0=0; x0<blocksize_x; x0+=4)
    {
        rx0 = (x0<<2) + cand_mv_x;
        x1  = x0 + 1;
        x2  = x1 + 1;
        x3  = x2 + 1;
        d   = diff;

        orig_line = orig_pic [y0];    
        ref_line  = get_line (ref_pic, ry0, rx0, img_height, img_width);
        *d++      = orig_line[x0] - *(ref_line     );
        *d++      = orig_line[x1] - *(ref_line + 4 );
        *d++      = orig_line[x2] - *(ref_line + 8 );
        *d++      = orig_line[x3] - *(ref_line + 12);

        orig_line = orig_pic [y1];    
        ref_line  = get_line (ref_pic, ry4, rx0, img_height, img_width);
        *d++      = orig_line[x0] - *(ref_line     );
        *d++      = orig_line[x1] - *(ref_line + 4 );
        *d++      = orig_line[x2] - *(ref_line + 8 );
        *d++      = orig_line[x3] - *(ref_line + 12);

        orig_line = orig_pic [y2];
        ref_line  = get_line (ref_pic, ry8, rx0, img_height, img_width);
        *d++      = orig_line[x0] - *(ref_line     );
        *d++      = orig_line[x1] - *(ref_line += 4 );
        *d++      = orig_line[x2] - *(ref_line += 4 );
        *d++      = orig_line[x3] - *(ref_line += 4);

        orig_line = orig_pic [y3];    
        ref_line  = get_line (ref_pic, ry12, rx0, img_height, img_width);
        *d++      = orig_line[x0] - *(ref_line     );
        *d++      = orig_line[x1] - *(ref_line += 4);
        *d++      = orig_line[x2] - *(ref_line += 4);
        *d        = orig_line[x3] - *(ref_line += 4);

      if (!useABT)
      {
        if ((mcost += SATD (diff, input->hadamard)) > min_mcost)
        {
          abort_search = 1;
          break;
        }
      }
      else  // copy diff to curr_diff for ABT SATD calculation
      {
          i = (x0&0x7) +  (x0>7) * 64 + y_offset;
          for(k=0, j=y0; j<BLOCK_SIZE + y0; j++, k+=BLOCK_SIZE)
            memcpy(&(c_diff[i + ((j&0x7)<<3)]), &diff[k], BLOCK_SIZE*sizeof(int));
      }
    }
  }
  
  if(useABT)
  {
    mcost += find_SATD (c_diff, blocktype);
  }

  return mcost;
}


int                                                   //  ==> minimum motion cost after search
FastSubPelBlockMotionSearch (pel_t**   orig_pic,      // <--  original pixel values for the AxB block
                             short     ref,           // <--  reference frame (0... or -1 (backward))
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
                             int       lambda_factor,
                             int       useABT)        // <--  lagrangian parameter for determining motion cost
{
  static int Diamond_x[4] = {-1, 0, 1, 0};
  static int Diamond_y[4] = {0, 1, 0, -1};
  int   mcost;
  int   cand_mv_x, cand_mv_y;
  
  int   list_offset   = ((img->MbaffFrameFlag)&&(img->mb_data[img->current_mb_nr].mb_field))? img->current_mb_nr%2 ? 4 : 2 : 0;
  StorablePicture *ref_picture = listX[list+list_offset][ref];
  
  int   mv_shift        = 0;
  int   blocksize_x     = input->blc_size[blocktype][0];
  int   blocksize_y     = input->blc_size[blocktype][1];
  int   pic4_pix_x      = ((pic_pix_x + IMG_PAD_SIZE)<< 2);
  int   pic4_pix_y      = ((pic_pix_y + IMG_PAD_SIZE)<< 2);
  short max_pos_x4      = ((ref_picture->size_x - blocksize_x + 2*IMG_PAD_SIZE)<<2);
  short max_pos_y4      = ((ref_picture->size_y - blocksize_y + 2*IMG_PAD_SIZE)<<2);
  
  int   search_range_dynamic,iXMinNow,iYMinNow,i;
  int   m,currmv_x = 0,currmv_y = 0;
  int   pred_frac_mv_x,pred_frac_mv_y,abort_search;
  int   mv_cost; 
  
  int   pred_frac_up_mv_x, pred_frac_up_mv_y;
  
  if ((pic4_pix_x + *mv_x > 1) && (pic4_pix_x + *mv_x < max_pos_x4 - 1) &&
      (pic4_pix_y + *mv_y > 1) && (pic4_pix_y + *mv_y < max_pos_y4 - 1)   )
  {
    get_line = FastLine4X;
  }
  else
  {
    get_line = UMVLine4X;    
  }

  search_range_dynamic = 3;
  pred_frac_mv_x = (pred_mv_x - *mv_x)%4;
  pred_frac_mv_y = (pred_mv_y - *mv_y)%4; 
  
  pred_frac_up_mv_x = (pred_MV_uplayer[0] - *mv_x)%4;
  pred_frac_up_mv_y = (pred_MV_uplayer[1] - *mv_y)%4;
  
  
  memset(SearchState[0],0,(2*search_range_dynamic+1)*(2*search_range_dynamic+1));
  
  if(input->hadamard)
  {
    cand_mv_x = *mv_x;    
    cand_mv_y = *mv_y;    
    mv_cost = MV_COST (lambda_factor, mv_shift, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);    
    mcost = AddUpSADQuarter(pic_pix_x,pic_pix_y,blocksize_x,blocksize_y,cand_mv_x + pic4_pix_x,cand_mv_y + pic4_pix_y,ref_picture,orig_pic,mv_cost,min_mcost,useABT, blocktype);
    SearchState[search_range_dynamic][search_range_dynamic] = 1;
    if (mcost < min_mcost)
    {
      min_mcost = mcost;
      currmv_x = cand_mv_x;
      currmv_y = cand_mv_y; 
    }
  }
  else
  {
    SearchState[search_range_dynamic][search_range_dynamic] = 1;
    currmv_x = *mv_x;
    currmv_y = *mv_y; 
  }
  
  if(pred_frac_mv_x!=0 || pred_frac_mv_y!=0)
  {
    cand_mv_x = *mv_x + pred_frac_mv_x;    
    cand_mv_y = *mv_y + pred_frac_mv_y;    
    mv_cost = MV_COST (lambda_factor, mv_shift, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);    
    mcost = AddUpSADQuarter(pic_pix_x,pic_pix_y,blocksize_x,blocksize_y,cand_mv_x + pic4_pix_x, cand_mv_y + pic4_pix_y,ref_picture,orig_pic,mv_cost,min_mcost,useABT, blocktype);
    SearchState[cand_mv_y -*mv_y + search_range_dynamic][cand_mv_x - *mv_x + search_range_dynamic] = 1;
    if (mcost < min_mcost)
    {
      min_mcost = mcost;
      currmv_x = cand_mv_x;
      currmv_y = cand_mv_y; 
    }
  }
  
  
  iXMinNow = currmv_x;
  iYMinNow = currmv_y;
  for(i=0;i<search_range_dynamic;i++) 
  {
    abort_search=1;
    for (m = 0; m < 4; m++)
    {
      cand_mv_x = iXMinNow + Diamond_x[m];    
      cand_mv_y = iYMinNow + Diamond_y[m]; 
      
      if(abs(cand_mv_x - *mv_x) <=search_range_dynamic && abs(cand_mv_y - *mv_y)<= search_range_dynamic)
      {
        if(!SearchState[cand_mv_y -*mv_y+ search_range_dynamic][cand_mv_x -*mv_x+ search_range_dynamic])
        {
          mv_cost = MV_COST (lambda_factor, mv_shift, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);    
          mcost = AddUpSADQuarter(pic_pix_x,pic_pix_y,blocksize_x,blocksize_y,cand_mv_x + pic4_pix_x, cand_mv_y + pic4_pix_y,ref_picture,orig_pic,mv_cost,min_mcost,useABT, blocktype);
          SearchState[cand_mv_y - *mv_y + search_range_dynamic][cand_mv_x - *mv_x + search_range_dynamic] = 1;
          if (mcost < min_mcost)
          {
            min_mcost = mcost;
            currmv_x = cand_mv_x;
            currmv_y = cand_mv_y; 
            abort_search = 0; 
            
          }
        }
      }
    }
    iXMinNow = currmv_x;
    iYMinNow = currmv_y;
    if(abort_search)
      break;
  }
  
  *mv_x = currmv_x;
  *mv_y = currmv_y;
  
  //===== return minimum motion cost =====
  return min_mcost;
}

/*!
 ************************************************************************
 * \brief
 * Functions for SAD prediction of intra block cases.
 * 1. void   decide_intrabk_SAD() judges the block coding type(intra/inter) 
 *    of neibouring blocks
 * 2. void skip_intrabk_SAD() set the SAD to zero if neigouring block coding 
 *    type is intra
 * \date
 *    2003.4
 ************************************************************************
 */
void   decide_intrabk_SAD()
{
  if (img->type != I_SLICE)
  {
    if (img->pix_x == 0 && img->pix_y == 0)
    {
      flag_intra_SAD = 0;
    }
    else if (img->pix_x == 0)
    {
      flag_intra_SAD = flag_intra[(img->pix_x)>>4];
    }
    else if (img->pix_y == 0)
    {
      flag_intra_SAD = flag_intra[((img->pix_x)>>4)-1];
    }
    else 
    {
      flag_intra_SAD = ((flag_intra[(img->pix_x)>>4])||(flag_intra[((img->pix_x)>>4)-1])||(flag_intra[((img->pix_x)>>4)+1])) ;
    }
  }
  return;
}

void skip_intrabk_SAD(int best_mode, int ref_max)
{
  int i,j,k, ref;
  if (img->number > 0) 
    flag_intra[(img->pix_x)>>4] = (best_mode == 9 || best_mode == 10) ? 1:0;
  if (img->type != I_SLICE  && (best_mode == 9 || best_mode == 10))
  {
    for (i=0; i < 4; i++)
    {
      for (j=0; j < 4; j++)
      {
        for (k=0; k < 9;k++)
        {
	  fastme_l0_cost[k][j][i] = 0;
	  fastme_l1_cost[k][j][i] = 0;
          for (ref=0; ref<ref_max;ref++)
          {
            fastme_ref_cost[ref][k][j][i] = 0;
          }
        }
      }
    }
  
  }
  return;
}


void setup_FME(short ref, int list, int block_y, int block_x, int blocktype, short   ******all_mv)
{
  int N_Bframe=0;
  int n_Bframe=0;
  int temp_blocktype = 0;
  int indication_blocktype[8]={0,0,1,1,2,4,4,5};
  N_Bframe = input->successive_Bframe;
  n_Bframe =(N_Bframe) ? (frame_ctr[B_SLICE]%(N_Bframe+1)): 0;


  /**************************** MV prediction **********************/ 
  //MV uplayer prediction
  if (blocktype>1) 
  {
    temp_blocktype = indication_blocktype[blocktype];
    pred_MV_uplayer[0] = all_mv[block_y][block_x][list][ref][temp_blocktype][0];
    pred_MV_uplayer[1] = all_mv[block_y][block_x][list][ref][temp_blocktype][1];

  }


  //MV ref-frame prediction
  pred_MV_ref_flag = 0;
  if(list==0)
  {
    if (img->field_picture) 
    {
      if ( ref > 1)
      {
        pred_MV_ref[0] = all_mv[block_y][block_x][0][ref-2][blocktype][0];
        pred_MV_ref[0] = (int)(pred_MV_ref[0]*((ref>>1)+1)/(float)((ref>>1)));
        pred_MV_ref[1] = all_mv[block_y][block_x][0][ref-2][blocktype][1];
        pred_MV_ref[1] = (int)(pred_MV_ref[1]*((ref>>1)+1)/(float)((ref>>1)));
        pred_MV_ref_flag = 1;
      }
      if (img->type == B_SLICE &&  (ref==0 || ref==1) )
      {
        pred_MV_ref[0] =(int) (all_mv[block_y][block_x][1][0][blocktype][0]*(-n_Bframe)/(N_Bframe-n_Bframe+1.0f));
        pred_MV_ref[1] =(int) (all_mv[block_y][block_x][1][0][blocktype][1]*(-n_Bframe)/(N_Bframe-n_Bframe+1.0f));
        pred_MV_ref_flag = 1;
      }
    }
    else //frame case
    {
      if ( ref > 0)
      {
        pred_MV_ref[0] = all_mv[block_y][block_x][0][ref-1][blocktype][0];
        pred_MV_ref[0] = (int)(pred_MV_ref[0]*(ref+1)/(float)(ref));
        pred_MV_ref[1] = all_mv[block_y][block_x][0][ref-1][blocktype][1];
        pred_MV_ref[1] = (int)(pred_MV_ref[1]*(ref+1)/(float)(ref));
        pred_MV_ref_flag = 1;
      }
      if (img->type == B_SLICE && (ref==0)) //B frame forward prediction, first ref
      {
        pred_MV_ref[0] =(int) (all_mv[block_y][block_x][1][0][blocktype][0]*(-n_Bframe)/(N_Bframe-n_Bframe+1.0f));
        pred_MV_ref[1] =(int) (all_mv[block_y][block_x][1][0][blocktype][1]*(-n_Bframe)/(N_Bframe-n_Bframe+1.0f));
        pred_MV_ref_flag = 1;
      }
    }
  }
  /******************************SAD prediction**********************************/
  if (list==0 && ref>0)  //pred_SAD_ref
  {

    if (flag_intra_SAD) //add this for irregular motion
    {
      pred_SAD = 0;
    }
    else 
    {
      if (img->field_picture)
      {
        if (ref > 1)
        {
          pred_SAD = fastme_ref_cost[ref-2][blocktype][block_y][block_x];
        }
        else
        {
          pred_SAD = fastme_ref_cost[0][blocktype][block_y][block_x];
        }
      }
      else
      {
        pred_SAD = fastme_ref_cost[ref-1][blocktype][block_y][block_x];
      }

    }
  }
  else if (blocktype>1)  // pred_SAD_uplayer
  {
    if (flag_intra_SAD) 
    {
      pred_SAD = 0;
    }
    else
    {
      pred_SAD = (list==1) ? (fastme_l1_cost[temp_blocktype][(img->pix_y>>2)+block_y][(img->pix_x>>2)+block_x]) : (fastme_l0_cost[temp_blocktype][(img->pix_y>>2)+block_y][(img->pix_x>>2)+block_x]);
      pred_SAD /= 2; 
    }
  }
  else pred_SAD = 0 ;  // pred_SAD_space

}

/*!
 ************************************************************************
 * \brief
 *    FastBipredIntegerPelBlockMotionSearch: fast pixel block motion search for bipred mode
 *    this algrithm is called UMHexagonS(see JVT-D016),which includes 
 *    four steps with different kinds of search patterns
 * \author
 *   Main contributors: (see contributors.h for copyright, address and affiliation details)
 *   - Zhibo Chen         <chenzhibo@tsinghua.org.cn>
 *   - JianFeng Xu        <fenax@video.mdc.tsinghua.edu.cn>
 *   - Xiaozhong Xu       <xxz@video.mdc.tsinghua.edu.cn>
 * \date   :
 *   2006.1
 ************************************************************************
 */
int                                                //  ==> minimum motion cost after search
FastBipredIntegerPelBlockMotionSearch (pel_t**   cur_pic,      // <--  original pixel values for the AxB block
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
                          int       lambda_factor) // <--  lagrangian parameter for determining motion cost

{
  int	temp_Big_Hexagon_x[16];// = Big_Hexagon_x;
  int	temp_Big_Hexagon_y[16];// = Big_Hexagon_y; 
  int   mvshift       = 2;                  // motion vector shift for getting sub-pel units

  int   search_step,iYMinNow, iXMinNow;
  int   i,m,j; 
  float betaFourth_1,betaFourth_2;
  int   pos, cand_x, cand_y,mcost;
  int   list_offset   = img->mb_data[img->current_mb_nr].list_offset; 
  int   blocksize_y   = input->blc_size[blocktype][1];            // vertical block size
  int   blocksize_x   = input->blc_size[blocktype][0];            // horizontal block size
  int   blocksize_x4  = blocksize_x >> 2;                         // horizontal block size in 4-pel units
  int   pred_x1        = (pic_pix_x << 2) + pred_mv_x1;       // predicted position x (in sub-pel units)
  int   pred_y1        = (pic_pix_y << 2) + pred_mv_y1;       // predicted position y (in sub-pel units)
  int   pred_x2        = (pic_pix_x << 2) + pred_mv_x2;       // predicted position x (in sub-pel units)
  int   pred_y2        = (pic_pix_y << 2) + pred_mv_y2;       // predicted position y (in sub-pel units)
  short center2_x      = pic_pix_x + *mv_x;                      // center position x (in pel units)
  short center2_y      = pic_pix_y + *mv_y;                      // center position y (in pel units)
  short center1_x	   = pic_pix_x + *s_mv_x;                      // mvx of second pred (in pel units)
  short center1_y	   = pic_pix_y + *s_mv_y;                      // mvy of second pred (in pel units)
  short apply_weights   = (active_pps->weighted_bipred_idc>0);  
  short offsetSpic = (apply_weights ? (list == 0?  wp_offset[list_offset    ][ref]     [0]:  wp_offset[list_offset + 1][0  ]     [0]) : 0);
  short offsetRpic = (apply_weights ? (list == 0?  wp_offset[list_offset + 1][ref]     [0]:  wp_offset[list_offset    ][0  ]     [0]) : 0);
  short mb_x = pic_pix_x - img->opix_x; 
  short mb_y = pic_pix_y - img->opix_y;
  short block_x = (mb_x >> 2);
  short block_y = (mb_y >> 2); 
  int   best_x = center2_x;
  int   best_y = center2_y;
  int ET_Thred = Median_Pred_Thd_MB[blocktype];
  offsetBi		= (offsetRpic + offsetSpic + 1)>>1;
  weightSpic	= (apply_weights ? (list == 0? wbp_weight[list_offset    ][ref][0  ][0]: wbp_weight[list_offset + 1][0  ][ref][0]) : 1<<luma_log_weight_denom);
  weightRpic	= (apply_weights ? (list == 0? wbp_weight[list_offset + 1][ref][0  ][0]: wbp_weight[list_offset    ][0  ][ref][0]) : 1<<luma_log_weight_denom);
  ref1_pic		= listX[list + list_offset          ][ref]->imgY_11;
  ref2_pic		= listX[list == 0 ? 1 + list_offset: list_offset][ 0 ]->imgY_11;  
  width			= listX[list+list_offset            ][ref]->size_x;
  height		= listX[list+list_offset            ][ref]->size_y;
  PartCalMadBiPred = apply_weights ? PartCalMadBiPred2 : PartCalMadBiPred1;
  
  //===== set function for getting reference picture lines =====
  if ((center2_x > search_range) && (center2_x < width -1-search_range-blocksize_x) &&
    (center2_y > search_range) && (center2_y < height-1-search_range-blocksize_y)   )
  {
    get_ref_line2 = FastLineX;
  }
  else
  {
    get_ref_line2 = UMVLineX2;
  }
  
  //===== set function for getting reference picture lines =====
  if ((center1_y > search_range) && (center1_y < height-1-search_range-blocksize_y)   )
  {
    get_ref_line1 = FastLineX;
  }
  else
  {
    get_ref_line1 = UMVLineX;
  }

//////////////////////////////////////////////////////////////////////////
   
  //////allocate memory for search state//////////////////////////
  memset(McostState[0],0,(2*search_range+1)*(2*search_range+1));


  //check the center median predictor
  cand_x = center2_x ;
  cand_y = center2_y ;
  mcost  = MV_COST (lambda_factor, mvshift, center1_x, center1_y, pred_x1, pred_y1);
  mcost += MV_COST (lambda_factor, mvshift, cand_x,	   cand_y,	  pred_x2, pred_y2);
  mcost  = PartCalMadBiPred(cur_pic, blocksize_y, blocksize_x, blocksize_x4, mcost,INT_MAX,center1_x, center1_y, cand_x, cand_y);
  McostState[search_range][search_range] = 1;
  if (mcost < min_mcost)
  {
    min_mcost = mcost;
    best_x = cand_x;
    best_y = cand_y;
  }

  iXMinNow = best_x;
  iYMinNow = best_y;
  for (m = 0; m < 4; m++)
  {   
    cand_x = iXMinNow + Diamond_x[m];
    cand_y = iYMinNow + Diamond_y[m];   
    SEARCH_ONE_PIXEL_BIPRED
  } 

  if(center2_x != pic_pix_x || center2_y != pic_pix_y)
  {
    cand_x = pic_pix_x ;
    cand_y = pic_pix_y ;
    SEARCH_ONE_PIXEL_BIPRED

    iXMinNow = best_x;
    iYMinNow = best_y;
    for (m = 0; m < 4; m++)
    {   
      cand_x = iXMinNow + Diamond_x[m];
      cand_y = iYMinNow + Diamond_y[m];   
      SEARCH_ONE_PIXEL_BIPRED
    } 
  }  
/***********************************init process*************************/

  if( min_mcost < ET_Thred)	
  {
    goto terminate_step;
  }
  else
  {
    int  N_Bframe=0;
	int  n_Bframe=0;
    short****** bipred_mv = list ? img->bipred_mv1 : img->bipred_mv2;
    N_Bframe = input->successive_Bframe;
    n_Bframe = frame_ctr[B_SLICE]%(N_Bframe+1);


  /**************************** MV prediction **********************/ 
  //MV uplayer prediction
  // non for bipred mode
  
  //MV ref-frame prediction

		if(list==0)
		{
		  if (img->field_picture) 
		  {
			  pred_MV_ref[0] =(int) (bipred_mv[block_y][block_x][1][0][blocktype][0]*(-n_Bframe)/(N_Bframe-n_Bframe+1.0f));
			  pred_MV_ref[1] =(int) (bipred_mv[block_y][block_x][1][0][blocktype][1]*(-n_Bframe)/(N_Bframe-n_Bframe+1.0f));
		  }
		  else //frame case
		  {
			  pred_MV_ref[0] =(int) (bipred_mv[block_y][block_x][1][0][blocktype][0]*(-n_Bframe)/(N_Bframe-n_Bframe+1.0f));
			  pred_MV_ref[1] =(int) (bipred_mv[block_y][block_x][1][0][blocktype][1]*(-n_Bframe)/(N_Bframe-n_Bframe+1.0f));
		  }
		}
  /******************************SAD prediction**********************************/

        pred_SAD =min(min(SAD_a,SAD_b),SAD_c);  // pred_SAD_space
		ET_Thred = Big_Hexagon_Thd_MB[blocktype];

 ///////Threshold defined for early termination///////////////////  
	  if (pred_SAD == 0) 
	  {
		betaFourth_1=0;
		betaFourth_2=0;
	  }
	  else
	  {
		betaFourth_1 = Bsize[blocktype]/(pred_SAD*pred_SAD)-AlphaFourth_1[blocktype];
		betaFourth_2 = Bsize[blocktype]/(pred_SAD*pred_SAD)-AlphaFourth_2[blocktype];
 
	  }  

  }

/***********************************end of init *************************/

 

// first_step: initial start point prediction 
  
//prediciton using mV of last ref moiton vector
 if(list == 0)	  		
  {
      cand_x = pic_pix_x + (pred_MV_ref[0]/4);
      cand_y = pic_pix_y + (pred_MV_ref[1]/4);
      SEARCH_ONE_PIXEL_BIPRED
  }


  //small local search
  iXMinNow = best_x;
  iYMinNow = best_y;
  for (m = 0; m < 4; m++)
  {   
    cand_x = iXMinNow + Diamond_x[m];
    cand_y = iYMinNow + Diamond_y[m];   
    SEARCH_ONE_PIXEL_BIPRED
  } 

  //early termination alogrithm, refer to JVT-G016
    EARLY_TERMINATION
  
  
//sec_step: //Unsymmetrical-cross search 
  iXMinNow = best_x;
  iYMinNow = best_y;
  
  for(i = 1; i < search_range; i+=2)
  {
    search_step = i;
    cand_x = iXMinNow + search_step;
    cand_y = iYMinNow ;
    SEARCH_ONE_PIXEL_BIPRED    
    cand_x = iXMinNow - search_step;
    cand_y = iYMinNow ;
    SEARCH_ONE_PIXEL_BIPRED
  }
 
  for(i = 1; i < (search_range/2);i+=2)
  {
    search_step = i;
    cand_x = iXMinNow ;
    cand_y = iYMinNow + search_step;
    SEARCH_ONE_PIXEL_BIPRED
    cand_x = iXMinNow ;
    cand_y = iYMinNow - search_step;
    SEARCH_ONE_PIXEL_BIPRED
  }
  //early termination alogrithm, refer to JVT-G016
    EARLY_TERMINATION

		
//third_step:     // Uneven Multi-Hexagon-grid Search 
  iXMinNow = best_x;
  iYMinNow = best_y;
//sub step1: 5x5 square search
  for(pos=1;pos<25;pos++)
  {
    cand_x = iXMinNow + spiral_search_x[pos];
    cand_y = iYMinNow + spiral_search_y[pos];
    SEARCH_ONE_PIXEL_BIPRED
  }

  //early termination alogrithm, refer to JVT-G016
	EARLY_TERMINATION			//added back by xxz

//sub step2: multi-grid-hexagon-search
  memcpy(temp_Big_Hexagon_x,Big_Hexagon_x,64);
  memcpy(temp_Big_Hexagon_y,Big_Hexagon_y,64);    		
  for(i=1;i<=(input->search_range>>2); i++)
  {

    for (m = 0; m < 16; m++)
    {
      cand_x = iXMinNow + temp_Big_Hexagon_x[m];
      cand_y = iYMinNow + temp_Big_Hexagon_y[m];
	  temp_Big_Hexagon_x[m] += Big_Hexagon_x[m];
	  temp_Big_Hexagon_y[m] += Big_Hexagon_y[m];	

      SEARCH_ONE_PIXEL_BIPRED
    }
	if(min_mcost < ET_Thred)
	{
	   		  goto terminate_step;

	}
  }
//fourth step: Local Refinement: Extended Hexagon-based Search
fourth_1_step:  

      for(i=0; i < search_range; i++) 
      {
        iXMinNow = best_x;
        iYMinNow = best_y;
        for (m = 0; m < 6; m++)
        {   
          cand_x = iXMinNow + Hexagon_x[m];
          cand_y = iYMinNow + Hexagon_y[m];   
          SEARCH_ONE_PIXEL_BIPRED
        } 
        if(best_x == iXMinNow && best_y == iYMinNow)
            break;
      }
fourth_2_step: 

      for(i = 0; i < search_range; i++) 
      {
	    iXMinNow = best_x;
	    iYMinNow = best_y;
        for (m = 0; m < 4; m++)
        {   
          cand_x = iXMinNow + Diamond_x[m];
          cand_y = iYMinNow + Diamond_y[m];   
          SEARCH_ONE_PIXEL_BIPRED
        } 
        if(best_x == iXMinNow && best_y == iYMinNow)
            break;
      }

terminate_step:
		  for (i=0; i < (blocksize_x>>2); i++)
		  {
			for (j=0; j < (blocksize_y>>2); j++)
			{
				if(list == 0) 
				{
				  fastme_l0_cost_bipred[blocktype][(img->pix_y>>2)+block_y+j][(img->pix_x>>2)+block_x+i] = min_mcost;
				}
				else
				{
				  fastme_l1_cost_bipred[blocktype][(img->pix_y>>2)+block_y+j][(img->pix_x>>2)+block_x+i] = min_mcost;
				}
			}
		  }

      *mv_x = best_x - pic_pix_x; 
      *mv_y = best_y - pic_pix_y; 


      return min_mcost;
}
