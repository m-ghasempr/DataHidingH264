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
 * \brief
 * Fast integer pel motion estimation and fractional pel motion estimation
 * algorithms are described in this file.
 * 1. get_mem_FME() and free_mem_FME() are functions for allocation and release
 *    of memories about motion estimation
 * 2. FME_BlockMotionSearch() is the function for fast integer pel motion 
 *    estimation and fractional pel motion estimation
 * 3. DefineThreshold() defined thresholds for early termination
 * \ Main contributors: (see contributors.h for copyright, address and affiliation details)
 *   Zhibo Chen         <chenzhibo@tsinghua.org.cn>
 *   JianFeng Xu        <fenax@video.mdc.tsinghua.edu.cn>  
 * \date   : 2003.8
 ************************************************************************
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <memory.h>
#include <assert.h>



#include "memalloc.h"
#include "fast_me.h"
#include "refbuf.h"
#include "mbuffer.h"
#include "image.h"

#define Q_BITS          15

extern  int*   byte_abs;
extern  int*   mvbits;
extern  int*   spiral_search_x;
extern  int*   spiral_search_y;


static pel_t (*PelY_14) (pel_t**, int, int);
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
  static float ThresholdFac[8] = {0,8,4,4,2.5,1.5,1.5,1}; 
  static int ThreshUp[8] = {0, 1024,512,512,448,384,384,384};

  AlphaSec[1] = 0.01f;
  AlphaSec[2] = 0.01f;
  AlphaSec[3] = 0.01f;
  AlphaSec[4] = 0.02f;
  AlphaSec[5] = 0.03f;
  AlphaSec[6] = 0.03f;
  AlphaSec[7] = 0.04f;

  AlphaThird[1] = 0.06f;
  AlphaThird[2] = 0.07f;
  AlphaThird[3] = 0.07f;
  AlphaThird[4] = 0.08f;
  AlphaThird[5] = 0.12f;
  AlphaThird[6] = 0.11f;
  AlphaThird[7] = 0.15f;

  DefineThresholdMB();
  return;
}

void DefineThresholdMB()
{
  int gb_qp_per    = (input->qpN-MIN_QP)/6;
  int gb_qp_rem    = (input->qpN-MIN_QP)%6;
  
  int gb_q_bits    = Q_BITS+gb_qp_per;
  int gb_qp_const,Thresh4x4;

  if (img->type == I_SLICE)
    gb_qp_const=(1<<gb_q_bits)/3;    // intra
  else
    gb_qp_const=(1<<gb_q_bits)/6;    // inter
  
  Thresh4x4 =   ((1<<gb_q_bits) - gb_qp_const)/quant_coef[gb_qp_rem][0][0];
  Quantize_step = Thresh4x4/(4*5.61f);
  Bsize[7]=(16*16)*Quantize_step;

  Bsize[6]=Bsize[7]*4;
  Bsize[5]=Bsize[7]*4;
  Bsize[4]=Bsize[5]*4;
  Bsize[3]=Bsize[4]*4;
  Bsize[2]=Bsize[4]*4;
  Bsize[1]=Bsize[2]*4;
}

/*!
 ************************************************************************
 * \brief
 *    Dynamic memory allocation of all infomation needed for Fast ME
 * \par Input:
 * \return Number of allocated bytes
 * \Date: 2003/3
 ************************************************************************
 */

int get_mem_mincost (int****** mv)
{
  int i, j, k, l;

  if(input->InterlaceCodingOption != FRAME_CODING)   img->buf_cycle *= 2;	

  if ((*mv = (int*****)calloc(input->img_width/4,sizeof(int****))) == NULL)
    no_mem_exit ("get_mem_mv: mv");
  for (i=0; i<input->img_width/4; i++)
  {
    if (((*mv)[i] = (int****)calloc(input->img_height/4,sizeof(int***))) == NULL)
      no_mem_exit ("get_mem_mv: mv");
    for (j=0; j<input->img_height/4; j++)
    {
      if (((*mv)[i][j] = (int***)calloc(img->buf_cycle,sizeof(int**))) == NULL)
        no_mem_exit ("get_mem_mv: mv");
      for (k=0; k<img->buf_cycle; k++)
      {
        if (((*mv)[i][j][k] = (int**)calloc(9,sizeof(int*))) == NULL)
          no_mem_exit ("get_mem_mv: mv");
        for (l=0; l<9; l++)
          if (((*mv)[i][j][k][l] = (int*)calloc(3,sizeof(int))) == NULL)
            no_mem_exit ("get_mem_mv: mv");
      }
    }
  }
  if(input->InterlaceCodingOption != FRAME_CODING)   img->buf_cycle /= 2;	

  return input->img_width/4*input->img_height/4*img->buf_cycle*9*3*sizeof(int);
}
/*!
 *******************************************************************************
 * \brief
 *    Dynamic memory allocation of all infomation needed for backward prediction
 * \par Input:
 * \return Number of allocated bytes
 * \Date: 2003/3
 *******************************************************************************
 */
int get_mem_bwmincost (int****** mv)
{
  int i, j, k, l;

  if(input->InterlaceCodingOption != FRAME_CODING)   img->buf_cycle *= 2;	

  if ((*mv = (int*****)calloc(input->img_width/4,sizeof(int****))) == NULL)
    no_mem_exit ("get_mem_mv: mv");
  for (i=0; i<input->img_width/4; i++)
  {
    if (((*mv)[i] = (int****)calloc(input->img_height/4,sizeof(int***))) == NULL)
      no_mem_exit ("get_mem_mv: mv");
    for (j=0; j<input->img_height/4; j++)
    {
      if (((*mv)[i][j] = (int***)calloc(img->buf_cycle,sizeof(int**))) == NULL)
        no_mem_exit ("get_mem_mv: mv");
      for (k=0; k<img->buf_cycle; k++)
      {
        if (((*mv)[i][j][k] = (int**)calloc(9,sizeof(int*))) == NULL)
          no_mem_exit ("get_mem_mv: mv");
        for (l=0; l<9; l++)
          if (((*mv)[i][j][k][l] = (int*)calloc(3,sizeof(int))) == NULL)
            no_mem_exit ("get_mem_mv: mv");
      }
    }
  }
 // if(input->InterlaceCodingOption != FRAME_CODING)   img->buf_cycle /= 2;	

  return input->img_width/4*input->img_height/4*1*9*3*sizeof(int);
}

int get_mem_FME()
{
	int memory_size = 0;
	memory_size += get_mem2Dint(&McostState, 2*input->search_range+1, 2*input->search_range+1);
	memory_size += get_mem_mincost (&(all_mincost));
	memory_size += get_mem_bwmincost(&(all_bwmincost));
	memory_size += get_mem2D(&SearchState,7,7);
	
	return memory_size;
}
/*!
 ************************************************************************
 * \brief
 *    free the memory allocated for of all infomation needed for Fast ME
 * \par Input:
 * \Date: 2003/3
 ************************************************************************
 */
void free_mem_mincost (int***** mv)
{
  int i, j, k, l;
  if(input->InterlaceCodingOption != FRAME_CODING)   img->buf_cycle *= 2;	

  for (i=0; i<input->img_width/4; i++)
  {
    for (j=0; j<input->img_height/4; j++)
    {
      for (k=0; k<img->buf_cycle; k++)
      {
        for (l=0; l<9; l++)
          free (mv[i][j][k][l]);
        free (mv[i][j][k]);
      }
      free (mv[i][j]);
    }
    free (mv[i]);
  }
  free (mv);
  if(input->InterlaceCodingOption != FRAME_CODING)   img->buf_cycle /= 2;	
}

/*!
 ***********************************************************************************
 * \brief
 *    free the memory allocated for of all infomation needed for backward prediction
 * \par Input:
 * \Date: 2003/3
 ***********************************************************************************
 */
void free_mem_bwmincost (int***** mv)
{
  int i, j, k, l;

  for (i=0; i<input->img_width/4; i++)
  {
    for (j=0; j<input->img_height/4; j++)
    {
      for (k=0; k<1; k++)
      {
        for (l=0; l<9; l++)
          free (mv[i][j][k][l]);
        free (mv[i][j][k]);
      }
      free (mv[i][j]);
    }
    free (mv[i]);
  }
  free (mv);
}

void free_mem_FME()
{
  free_mem2Dint(McostState);
  free_mem_mincost (all_mincost);
  free_mem_bwmincost(all_bwmincost);

  free_mem2D(SearchState);
}

void
FME_SetMotionVectorPredictor (int  pmv[2],
                          int  **refFrArray,
                          int  ***tmp_mv,
                          int  ref_frame,
						  int  list,
                          int  mb_x,
                          int  mb_y,
                          int  blockshape_x,
                          int  blockshape_y,
						  int  blocktype,
						  int  ref)
{
  int pic_block_x          = img->block_x + (mb_x>>2);
  int pic_block_y          = img->block_y + (mb_y>>2);
  int mb_nr                = img->current_mb_nr;
  int mb_width             = img->width/16;
  int mb_available_up      = (img->mb_y == 0          ) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width  ].slice_nr);
  int mb_available_left    = (img->mb_x == 0          ) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-1         ].slice_nr);
  int mb_available_upleft  = (img->mb_x == 0 ||
                              img->mb_y == 0          ) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width-1].slice_nr);
  int mb_available_upright = (img->mb_x >= mb_width-1 ||
                              img->mb_y == 0          ) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width+1].slice_nr);
  int block_available_up, block_available_left, block_available_upright, block_available_upleft;
  int mv_a, mv_b, mv_c, mv_d, pred_vec=0;
  int mvPredType, rFrameL, rFrameU, rFrameUR;
  int hv;

  //FAST MOTION ESTIMATION. ZHIBO CHEN 2003.3
  int SAD_a, SAD_b, SAD_c, SAD_d;
  int temp_pred_SAD[2];
  pred_SAD_space = 0;

  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    pic_block_y          = img->field_block_y + (mb_y>>2);
    mb_available_up      = (img->field_mb_y == 0          ) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width  ].slice_nr);
    mb_available_left    = (img->mb_x == 0          ) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-1         ].slice_nr);
    mb_available_upleft  = (img->mb_x == 0 ||
                            img->field_mb_y == 0          ) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width-1].slice_nr);
    mb_available_upright = (img->mb_x >= mb_width-1 ||
                            img->field_mb_y == 0          ) ? 0 : (img->mb_data[mb_nr].slice_nr == img->mb_data[mb_nr-mb_width+1].slice_nr);
  }
  
  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
  {
    if(img->field_mode && !img->top_field)
      mb_available_upright=0; // set mb_available_upright to 0 for bottom MBs in a MB pair
    else if(!img->field_mode && img->mb_y%2)
      mb_available_upright=0; // set mb_available_upright to 0 for bottom MBs in a MB pair

  }

  /* D B C */
  /* A X   */
  
  /* 1 A, B, D are set to 0 if unavailable       */
  /* 2 If C is not available it is replaced by D */
  block_available_up   = mb_available_up   || (mb_y > 0);
  block_available_left = mb_available_left || (mb_x > 0);

  if (mb_y > 0)
  {
    if (mb_x < 8)  // first column of 8x8 blocks
    {
      if (mb_y==8)
      {
        if (blockshape_x == 16)      block_available_upright = 0;
        else                         block_available_upright = 1;
      }
      else
      {
        if (mb_x+blockshape_x != 8)  block_available_upright = 1;
        else                         block_available_upright = 0;
      }
    }
    else
    {
      if (mb_x+blockshape_x != 16)   block_available_upright = 1;
      else                           block_available_upright = 0;
    }
  }
  else if (mb_x+blockshape_x != MB_BLOCK_SIZE)
  {
    block_available_upright = block_available_up;
  }
  else
  {
    block_available_upright = mb_available_upright;
  }
  
  if (mb_x > 0)
  {
    block_available_upleft = (mb_y > 0 ? 1 : block_available_up);
  }
  else if (mb_y > 0)
  {
    block_available_upleft = block_available_left;
  }
  else
  {
    block_available_upleft = mb_available_upleft;
  }
  
  	//write refFrArr[][] for current subMB partition
		refFrArray[pic_block_y][pic_block_x]= ref;	//needed for other subMB partitions in this MB partition

//  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive)
//    block_available_upright = 0;      // temp fix for MB level field/frame coding

  mvPredType = MVPRED_MEDIAN;
  rFrameL    = block_available_left    ? refFrArray[pic_block_y]  [pic_block_x-1] : -1;
  rFrameU    = block_available_up      ? refFrArray[pic_block_y-1][pic_block_x]   : -1;
  rFrameUR   = block_available_upright ? refFrArray[pic_block_y-1][pic_block_x+blockshape_x/4] :
               block_available_upleft  ? refFrArray[pic_block_y-1][pic_block_x-1] : -1;
  
  
  
  /* Prediction if only one of the neighbors uses the reference frame
  * we are checking
  */
  
  if(rFrameL == ref_frame && rFrameU != ref_frame && rFrameUR != ref_frame)
    mvPredType = MVPRED_L;
  else if(rFrameL != ref_frame && rFrameU == ref_frame && rFrameUR != ref_frame)
    mvPredType = MVPRED_U;
  else if(rFrameL != ref_frame && rFrameU != ref_frame && rFrameUR == ref_frame)
    mvPredType = MVPRED_UR;
  // Directional predictions 
  if(blockshape_x == 8 && blockshape_y == 16)
  {
    if(mb_x == 0)
    {
      if(rFrameL == ref_frame)
        mvPredType = MVPRED_L;
    }
    else
    {
      //if( block_available_upright && refFrArray[pic_block_y-1][pic_block_x+blockshape_x/4] == ref_frame)
      if(rFrameUR == ref_frame) 
        mvPredType = MVPRED_UR;
    }
  }
  else if(blockshape_x == 16 && blockshape_y == 8)
  {
    if(mb_y == 0)
    {
      if(rFrameU == ref_frame)
        mvPredType = MVPRED_U;
    }
    else
    {
      if(rFrameL == ref_frame)
        mvPredType = MVPRED_L;
    }
  }
  
#define MEDIAN(a,b,c)  (a>b?a>c?b>c?b:c:a:b>c?a>c?a:c:b)
  
  for (hv=0; hv < 2; hv++)
  {
    mv_a = block_available_left    ? tmp_mv[hv][pic_block_y  ][4+pic_block_x-1]              : 0;
    mv_b = block_available_up      ? tmp_mv[hv][pic_block_y-1][4+pic_block_x]                : 0;
    mv_d = block_available_upleft  ? tmp_mv[hv][pic_block_y-1][4+pic_block_x-1]              : 0;
    mv_c = block_available_upright ? tmp_mv[hv][pic_block_y-1][4+pic_block_x+blockshape_x/4] : mv_d;
    
	SAD_a = block_available_left    ? ((list==1) ? all_bwmincost[((img->pix_x+mb_x)>>2) -1][((img->pix_y+mb_y)>>2)][0][blocktype][0] : all_mincost[((img->pix_x+mb_x)>>2) -1][((img->pix_y+mb_y)>>2)][ref_frame][blocktype][0]) : 0;
    SAD_b = block_available_up      ? ((list==1) ? all_bwmincost[((img->pix_x+mb_x)>>2)][((img->pix_y+mb_y)>>2) -1][0][blocktype][0] : all_mincost[((img->pix_x+mb_x)>>2)][((img->pix_y+mb_y)>>2) -1][ref_frame][blocktype][0]) : 0;
    SAD_d = block_available_upleft  ? ((list==1) ? all_bwmincost[((img->pix_x+mb_x)>>2) -1][((img->pix_y+mb_y)>>2) -1][0][blocktype][0] : all_mincost[((img->pix_x+mb_x)>>2) -1][((img->pix_y+mb_y)>>2) -1][ref_frame][blocktype][0]) : 0;
    SAD_c = block_available_upright ? ((list==1) ? all_bwmincost[((img->pix_x+mb_x)>>2) +1][((img->pix_y+mb_y)>>2) -1][0][blocktype][0] : all_mincost[((img->pix_x+mb_x)>>2) +1][((img->pix_y+mb_y)>>2) -1][ref_frame][blocktype][0]) : SAD_d;
	
    switch (mvPredType)
    {
    case MVPRED_MEDIAN:
      if(!(block_available_upleft || block_available_up || block_available_upright))
	  {
        pred_vec = mv_a;
//FAST MOTION ESTIMATION. ZHIBO CHEN 2003.3
	  	  temp_pred_SAD[hv] = SAD_a;
		}
      else
        pred_vec = MEDIAN (mv_a, mv_b, mv_c);
	  if (pred_vec == mv_a && SAD_a != 0) temp_pred_SAD[hv] = SAD_a;
	  else if (pred_vec == mv_b && SAD_b!=0) temp_pred_SAD[hv] = SAD_b;
	  else temp_pred_SAD[hv] = SAD_c;
      break;
      
    case MVPRED_L:
      pred_vec = mv_a;
	  temp_pred_SAD[hv] = SAD_a;
      break;
    case MVPRED_U:
      pred_vec = mv_b;
	  temp_pred_SAD[hv] = SAD_b;
      break;
    case MVPRED_UR:
      pred_vec = mv_c;
	  temp_pred_SAD[hv] = SAD_c;
      break;
    default:
      break;
    }

    pmv[hv] = pred_vec;
  }
//FAST MOTION ESTIMATION. ZHIBO CHEN 2003.3
  pred_SAD_space = temp_pred_SAD[0]>temp_pred_SAD[1]?temp_pred_SAD[1]:temp_pred_SAD[0];
#undef MEDIAN
}

int                                         //  ==> minimum motion cost after search
FME_BlockMotionSearch (int       ref,           // <--  reference frame (0... )
									 int       list,
                   int       pic_pix_x,     // <--  absolute x-coordinate of regarded AxB block
                   int       pic_pix_y,     // <--  absolute y-coordinate of regarded AxB block
                   int       blocktype,     // <--  block type (1-16x16 ... 7-4x4)
                   int       search_range,  // <--  1-d search range for integer-position search
                   double    lambda         // <--  lagrangian parameter for determining motion cost
                   )
{
  static pel_t   orig_val [256];
  static pel_t  *orig_pic  [16] = {orig_val,     orig_val+ 16, orig_val+ 32, orig_val+ 48,
                                   orig_val+ 64, orig_val+ 80, orig_val+ 96, orig_val+112,
                                   orig_val+128, orig_val+144, orig_val+160, orig_val+176,
                                   orig_val+192, orig_val+208, orig_val+224, orig_val+240};

  int       pred_mv_x, pred_mv_y, mv_x, mv_y, i, j;

  int       max_value = (1<<20);
  int       min_mcost = max_value;
  int       mb_x      = pic_pix_x-img->pix_x;
  int       mb_y      = pic_pix_y-img->pix_y;
  int       block_x   = (mb_x>>2);
  int       block_y   = (mb_y>>2);
  int       bsx       = input->blc_size[blocktype][0];
  int       bsy       = input->blc_size[blocktype][1];
  int       refframe  = ref;		//why chage name?
  int*      pred_mv;
  //int**     refFrArray = ((img->type!=B_SLICE) ? refFrArr : list==0 ? fw_refFrArr : bw_refFrArr);
  //int***    tmp_mv_array  = ((img->type!=B_SLICE) ? tmp_mv   : list==0 ? tmp_fwMV    : tmp_bwMV);
  
  //sw 10.23
  int**     refFrArray    = enc_picture->ref_idx[list];
  int***    tmp_mv_array  = enc_picture->mv[list];
  
  int*****  all_bmv   = img->all_bmv;		//probably redundant g050
  int*****  all_mv    = (list ? img->all_bmv : img->all_mv);	
  byte**    imgY_org_pic = imgY_org;

//FAST MOTION ESTIMATION. ZHIBO CHEN 2003.3
  int       N_Bframe = input->successive_Bframe, n_Bframe =(N_Bframe) ? ((Bframe_ctr%N_Bframe)+1) : 0 ;
assert(ref>=0);
  if(input->InterlaceCodingOption >= MB_CODING && mb_adaptive && img->field_mode)
  {
    mb_y    = pic_pix_y - img->field_pix_y;
    block_y = mb_y >> 2;
    if(img->top_field)
    {
      pred_mv   = ((img->type!=B_SLICE) ? img->pred_mv_top  : list==0 ? img->p_fwMV_top : img->p_bwMV_top)[mb_x>>2][mb_y>>2][refframe][blocktype];
      refFrArray = ((img->type!=B_SLICE) ? refFrArr_top : list==0 ? fw_refFrArr_top : bw_refFrArr_top);
      tmp_mv_array  = ((img->type!=B_SLICE) ? tmp_mv_top   : list==0 ? tmp_fwMV_top    : tmp_bwMV_top);
      all_bmv   = img->all_bmv_top;
      all_mv    = (list ? img->all_bmv_top    : img->all_mv_top);
      imgY_org_pic = imgY_org_top;
    }
    else
    {
      pred_mv   = ((img->type!=B_SLICE) ? img->pred_mv_bot   : list==0 ? img->p_fwMV_bot : img->p_bwMV_bot)[mb_x>>2][mb_y>>2][refframe][blocktype];
      refFrArray = ((img->type!=B_SLICE) ? refFrArr_bot  : list==0 ? fw_refFrArr_bot : bw_refFrArr_bot);
      tmp_mv_array  = ((img->type!=B_SLICE) ? tmp_mv_bot    : list==0 ? tmp_fwMV_bot    : tmp_bwMV_bot);
      all_bmv   = img->all_bmv_bot;
      all_mv    = (list ? img->all_bmv_bot    : img->all_mv_bot);
      imgY_org_pic = imgY_org_bot;
    }
  }
  else
	  pred_mv = ((img->type!=B_SLICE) ? img->pred_mv  : list==0 ? img->p_fwMV : img->p_bwMV)[mb_x>>2][mb_y>>2][refframe][blocktype];


											//NB  pred_mv is only 2D here

  //==================================
  //=====   GET ORIGINAL BLOCK   =====
  //==================================
  for (j = 0; j < bsy; j++)
  {
    for (i = 0; i < bsx; i++)
    {
      orig_pic[j][i] = imgY_org_pic[pic_pix_y+j][pic_pix_x+i];
    }
  }

//FAST MOTION ESTIMATION. ZHIBO CHEN 2003.3
  if(blocktype>6)
  {
	  pred_MV_uplayer[0] = all_mv[block_x][block_y][refframe][5][0];
	  pred_MV_uplayer[1] = all_mv[block_x][block_y][refframe][5][1];
	  pred_SAD_uplayer    = (list==1) ? (all_bwmincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][0][5][0]) : (all_mincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][refframe][5][0]);
	  pred_SAD_uplayer   /= 2; 
	  
  }
  else if(blocktype>4)
  {
	  pred_MV_uplayer[0] = all_mv[block_x][block_y][refframe][4][0];
	  pred_MV_uplayer[1] = all_mv[block_x][block_y][refframe][4][1];
	  pred_SAD_uplayer    = (list==1) ? (all_bwmincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][0][4][0]) : (all_mincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][refframe][4][0]);
	  pred_SAD_uplayer   /= 2; 
	  
  }
  else if(blocktype == 4)
  {
	  pred_MV_uplayer[0] = all_mv[block_x][block_y][refframe][2][0];
	  pred_MV_uplayer[1] = all_mv[block_x][block_y][refframe][2][1];
	  pred_SAD_uplayer    = (list==1) ? (all_bwmincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][0][2][0]) : (all_mincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][refframe][2][0]);
	  pred_SAD_uplayer   /= 2; 
  }
  else if(blocktype > 1)
  {
	  pred_MV_uplayer[0] = all_mv[block_x][block_y][refframe][1][0];
	  pred_MV_uplayer[1] = all_mv[block_x][block_y][refframe][1][1];
	  pred_SAD_uplayer    = (list==1) ? (all_bwmincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][0][1][0]) : (all_mincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][refframe][1][0]);
	  pred_SAD_uplayer   /= 2; 
  }
  
  if ((img->type==B_SLICE)&& (img->nal_reference_idc>0))
  {
	  if(blocktype>6)
	  {
		  pred_MV_uplayer[0] = all_bmv[block_x][block_y][refframe][5][0];
		  pred_MV_uplayer[1] = all_bmv[block_x][block_y][refframe][5][1];
		  pred_SAD_uplayer    = (list==1) ? (all_bwmincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][0][5][0]) : (all_mincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][refframe][5][0]);
		  pred_SAD_uplayer   /= 2; 
	  }
	  else if(blocktype>4)
	  {
		  pred_MV_uplayer[0] = all_bmv[block_x][block_y][refframe][4][0];
		  pred_MV_uplayer[1] = all_bmv[block_x][block_y][refframe][4][1];
		  pred_SAD_uplayer    = (list==1) ? (all_bwmincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][0][4][0]) : (all_mincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][refframe][4][0]);
		  pred_SAD_uplayer   /= 2; 
	  }
	  else if(blocktype == 4)
	  {
		  pred_MV_uplayer[0] = all_bmv[block_x][block_y][refframe][2][0];
		  pred_MV_uplayer[1] = all_bmv[block_x][block_y][refframe][2][1];
		  pred_SAD_uplayer    = (list==1) ? (all_bwmincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][0][2][0]) : (all_mincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][refframe][2][0]);
		  pred_SAD_uplayer   /= 2; 
	  }
	  else if(blocktype > 1)
	  {
		  pred_MV_uplayer[0] = all_bmv[block_x][block_y][refframe][1][0];
		  pred_MV_uplayer[1] = all_bmv[block_x][block_y][refframe][1][1];
		  pred_SAD_uplayer    = (list==1) ? (all_bwmincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][0][1][0]) : (all_mincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][refframe][1][0]);
		  pred_SAD_uplayer   /= 2; 
	  }
  }

  pred_SAD_uplayer = flag_intra_SAD ? 0 : pred_SAD_uplayer;// for irregular motion
 
  //Coordinate prediction
  if (img->number > refframe+1)
  {
	  pred_SAD_time = all_mincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][refframe][blocktype][0];
	  pred_MV_time[0] = all_mincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][refframe][blocktype][1];
	  pred_MV_time[1] = all_mincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][refframe][blocktype][2];
  }
 
  if(list==1 && (Bframe_ctr%N_Bframe) > 1) 
  {
	  pred_SAD_time = all_bwmincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][refframe][blocktype][0];
	  pred_MV_time[0] = (int)(all_bwmincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][0][blocktype][1] * ((n_Bframe==1) ? (N_Bframe) : (N_Bframe-n_Bframe+1.0)/(N_Bframe-n_Bframe+2.0)) );//should add a factor
	  pred_MV_time[1] = (int)(all_bwmincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][0][blocktype][2] *((n_Bframe==1) ? (N_Bframe) : (N_Bframe-n_Bframe+1.0)/(N_Bframe-n_Bframe+2.0)) );//should add a factor
  }

  if (input->InterlaceCodingOption == FIELD_CODING) 
  {
	  if (img->type == P_SLICE && refframe > 1)
	  {
		  pred_SAD_ref = all_mincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][(refframe-2)][blocktype][0];
		  pred_SAD_ref = flag_intra_SAD ? 0 : pred_SAD_ref;//add this for irregular motion
		  pred_MV_ref[0] = all_mincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][(refframe-2)][blocktype][1];
		  pred_MV_ref[0] = (int)(pred_MV_ref[0]*((refframe>>1)+1)/(float)((refframe>>1)));
		  pred_MV_ref[1] = all_mincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][(refframe-2)][blocktype][2];
		  pred_MV_ref[1] = (int)(pred_MV_ref[1]*((refframe>>1)+1)/(float)((refframe>>1)));
	  }
	  if (img->type == B_SLICE && list==0 && (refframe==0 || refframe==1) )
	  {
		  pred_SAD_ref = all_bwmincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][0][blocktype][0];
		  pred_SAD_ref = flag_intra_SAD ? 0 : pred_SAD_ref;//add this for irregular motion
		  pred_MV_ref[0] =(int) (all_bwmincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][refframe][blocktype][1]*(-n_Bframe)/(N_Bframe-n_Bframe+1.0f)); //should add a factor
		  pred_MV_ref[1] =(int) ( all_bwmincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][refframe][blocktype][2]*(-n_Bframe)/(N_Bframe-n_Bframe+1.0f)); 
	  }
  }
 else //frame case
  {
	  if (refframe > 0)
	  {//field_mode top_field
		  pred_SAD_ref = all_mincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][(refframe-1)][blocktype][0];
		  pred_SAD_ref = flag_intra_SAD ? 0 : pred_SAD_ref;//add this for irregular motion
		  pred_MV_ref[0] = all_mincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][(refframe-1)][blocktype][1];
		  pred_MV_ref[0] = (int)(pred_MV_ref[0]*(refframe+1)/(float)(refframe));
		  pred_MV_ref[1] = all_mincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][(refframe-1)][blocktype][2];
		  pred_MV_ref[1] = (int)(pred_MV_ref[1]*(refframe+1)/(float)(refframe));
	  }
	  if (img->type == B_SLICE && (list==0 && refframe==0)) //B frame forward prediction, first ref
	  {
		  pred_SAD_ref = all_bwmincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][0][blocktype][0];
		  pred_SAD_ref = flag_intra_SAD ? 0 : pred_SAD_ref;//add this for irregular motion
		  pred_MV_ref[0] =(int) (all_bwmincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][refframe][blocktype][1]*(-n_Bframe)/(N_Bframe-n_Bframe+1.0f)); //should add a factor
		  pred_MV_ref[1] =(int) ( all_bwmincost[(img->pix_x>>2)+block_x][(img->pix_y>>2)+block_y][refframe][blocktype][2]*(-n_Bframe)/(N_Bframe-n_Bframe+1.0f)); 
	  }
  }
  //===========================================
  //=====   GET MOTION VECTOR PREDICTOR   =====
  //===========================================
											//NB  pred_mv is only 2D here
  FME_SetMotionVectorPredictor (pred_mv, refFrArray, tmp_mv_array, refframe, list, mb_x, mb_y, bsx, bsy, blocktype, ref);
  pred_mv_x = pred_mv[0];
  pred_mv_y = pred_mv[1];


  //==================================
  //=====   INTEGER-PEL SEARCH   =====
  //==================================

//FAST MOTION ESTIMATION. ZHIBO CHEN 2003.3
  mv_x = pred_mv_x / 4;
  mv_y = pred_mv_y / 4;
  if (!input->rdopt)
  {
    //--- adjust search center so that the (0,0)-vector is inside ---
    mv_x = max (-search_range, min (search_range, mv_x));
    mv_y = max (-search_range, min (search_range, mv_y));
  }

  min_mcost = FastIntegerPelBlockMotionSearch(orig_pic, ref, list, pic_pix_x, pic_pix_y, blocktype,
                                            pred_mv_x, pred_mv_y, &mv_x, &mv_y, search_range,
                                            min_mcost, lambda);



//FAST MOTION ESTIMATION. ZHIBO CHEN 2003.3
  for (i=0; i < (bsx>>2); i++)
  {
	  for (j=0; j < (bsy>>2); j++)
	  {
		  if(list == 0) 
			  all_mincost[(img->pix_x>>2)+block_x+i][(img->pix_y>>2)+block_y+j][refframe][blocktype][0] = min_mcost;
		  else
			  all_bwmincost[(img->pix_x>>2)+block_x+i][(img->pix_y>>2)+block_y+j][refframe][blocktype][0] = min_mcost; 
	  }
  }


 //==============================
  //=====   SUB-PEL SEARCH   =====
  //==============================
  if (input->hadamard)
  {
    min_mcost = max_value;
  }
  
if(blocktype >3)
{
  min_mcost =  FastSubPelBlockMotionSearch (orig_pic, ref, list, pic_pix_x, pic_pix_y, blocktype,
                                        pred_mv_x, pred_mv_y, &mv_x, &mv_y, 9, 9,
                                        min_mcost, lambda, /*useABT*/0);
}
else
{
  min_mcost =  SubPelBlockMotionSearch (orig_pic, ref, list, pic_pix_x, pic_pix_y, blocktype,
                                        pred_mv_x, pred_mv_y, &mv_x, &mv_y, 9, 9,
                                        min_mcost, lambda);
}


  for (i=0; i < (bsx>>2); i++)
  {
	  for (j=0; j < (bsy>>2); j++)
	  {
		  if(list == 0)
		  {
			  all_mincost[(img->pix_x>>2)+block_x+i][(img->pix_y>>2)+block_y+j][refframe][blocktype][1] = mv_x;
			  all_mincost[(img->pix_x>>2)+block_x+i][(img->pix_y>>2)+block_y+j][refframe][blocktype][2] = mv_y;
		  }
		  else
		  {
			  all_bwmincost[(img->pix_x>>2)+block_x+i][(img->pix_y>>2)+block_y+j][refframe][blocktype][1] = mv_x;
			  all_bwmincost[(img->pix_x>>2)+block_x+i][(img->pix_y>>2)+block_y+j][refframe][blocktype][2] = mv_y;

		  }
	  }
  }


  if (!input->rdopt)
  {
    // Get the skip mode cost
    if (blocktype == 1 && img->type == P_SLICE)
    {
      int cost;

      FindSkipModeMotionVector ();

      cost  = GetSkipCostMB (lambda);
      cost -= (int)floor(8*lambda+0.4999);

      if (cost < min_mcost)
      {
        min_mcost = cost;
        mv_x      = img->all_mv [0][0][0][0][0];
        mv_y      = img->all_mv [0][0][0][0][1];
      }
    }
  }

  //===============================================
  //=====   SET MV'S AND RETURN MOTION COST   =====
  //===============================================
  {int h4x4blkno = 4+(pic_pix_x>>2);  int v4x4blkno = pic_pix_y>>2;
  for (i=0; i < (bsx>>2); i++){
    for (j=0; j < (bsy>>2); j++){
      all_mv[block_x+i][block_y+j][refframe][blocktype][0] = mv_x;
      all_mv[block_x+i][block_y+j][refframe][blocktype][1] = mv_y;
			tmp_mv_array[0][v4x4blkno+j][h4x4blkno+i] = mv_x;	//adding 4 is weird, also y,x instead of x,y????
			tmp_mv_array[1][v4x4blkno+j][h4x4blkno+i] = mv_y;
    }
  }}
  if (img->type==B_SLICE && img->nal_reference_idc>0)	//probably wrong, executed even for list = 0. g050
  {
    for (i=0; i < (bsx>>2); i++)
    for (j=0; j < (bsy>>2); j++)
    {
      //  Backward
      all_bmv[block_x+i][block_y+j][ref][blocktype][0] = mv_x;
      all_bmv[block_x+i][block_y+j][ref][blocktype][1] = mv_y;
    }
  }
  return min_mcost;
}

_inline int PartCalMad(pel_t *ref_pic,pel_t** orig_pic,pel_t *(*get_ref_line)(int, pel_t*, int, int), int blocksize_y,int blocksize_x, int blocksize_x4,int mcost,int min_mcost,int cand_x,int cand_y)
{
	int y,x4;
	pel_t *orig_line, *ref_line;
	for (y=0; y<blocksize_y; y++)
    {
		ref_line  = get_ref_line (blocksize_x, ref_pic, cand_y+y, cand_x);
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
 * double    lambda        // <--  lagrangian parameter for determining motion cost
 * \par
 * Three macro definitions defined in this program:
 * 1. EARLY_TERMINATION: early termination algrithm, refer to JVT-D016.doc
 * 2. SEARCH_ONE_PIXEL: search one pixel in search range
 * 3. SEARCH_ONE_PIXEL1(value_iAbort): search one pixel in search range,
 *                                 but give a parameter to show if mincost refeshed
 * \ Main contributors: (see contributors.h for copyright, address and affiliation details)
 *   Zhibo Chen         <chenzhibo@tsinghua.org.cn>
 *   JianFeng Xu        <fenax@video.mdc.tsinghua.edu.cn>  
 * \date   : 2003.8
 ************************************************************************
 */
int                                     //  ==> minimum motion cost after search
FastIntegerPelBlockMotionSearch  (pel_t**   orig_pic,     // <--  not used
								  int       ref,          // <--  reference frame (0... or -1 (backward))
									int       list,
								  int       pic_pix_x,    // <--  absolute x-coordinate of regarded AxB block
								  int       pic_pix_y,    // <--  absolute y-coordinate of regarded AxB block
								  int       blocktype,    // <--  block type (1-16x16 ... 7-4x4)
								  int       pred_mv_x,    // <--  motion vector predictor (x) in sub-pel units
								  int       pred_mv_y,    // <--  motion vector predictor (y) in sub-pel units
								  int*      mv_x,         //  --> motion vector (x) - in pel units
								  int*      mv_y,         //  --> motion vector (y) - in pel units
								  int       search_range, // <--  1-d search range in pel units                         
								  int       min_mcost,    // <--  minimum motion cost (cost for center or huge value)
								  double    lambda)       // <--  lagrangian parameter for determining motion cost
{
	static int Diamond_x[4] = {-1, 0, 1, 0};
	static int Diamond_y[4] = {0, 1, 0, -1};
	static int Hexagon_x[6] = {2, 1, -1, -2, -1, 1};
	static int Hexagon_y[6] = {0, -2, -2, 0,  2, 2};
	static int Big_Hexagon_x[16] = {0,-2, -4,-4,-4, -4, -4, -2,  0,  2,  4,  4, 4, 4, 4, 2};
	static int Big_Hexagon_y[16] = {4, 3, 2,  1, 0, -1, -2, -3, -4, -3, -2, -1, 0, 1, 2, 3};

	int   pos, cand_x, cand_y,  mcost;
	pel_t *(*get_ref_line)(int, pel_t*, int, int);
	pel_t*  ref_pic       = listX[list][ref]->imgY_11;//img->type==B_IMG? Refbuf11 [ref+((mref==mref_fld)) +1] : Refbuf11[ref];
	int   best_pos      = 0;                                        // position with minimum motion cost
	int   max_pos       = (2*search_range+1)*(2*search_range+1);    // number of search positions
	int   lambda_factor = LAMBDA_FACTOR (lambda);                   // factor for determining lagragian motion cost
	int   mvshift       = 2;                  // motion vector shift for getting sub-pel units
	int   blocksize_y   = input->blc_size[blocktype][1];            // vertical block size
	int   blocksize_x   = input->blc_size[blocktype][0];            // horizontal block size
	int   blocksize_x4  = blocksize_x >> 2;                         // horizontal block size in 4-pel units
	int   pred_x        = (pic_pix_x << mvshift) + pred_mv_x;       // predicted position x (in sub-pel units)
	int   pred_y        = (pic_pix_y << mvshift) + pred_mv_y;       // predicted position y (in sub-pel units)
	int   center_x      = pic_pix_x + *mv_x;                        // center position x (in pel units)
	int   center_y      = pic_pix_y + *mv_y;                        // center position y (in pel units)
	int    best_x, best_y;
	int   check_for_00  = (blocktype==1 && !input->rdopt && img->type!=B_SLICE && ref==0);
	int   search_step,iYMinNow, iXMinNow;
	int   i,m, iSADLayer; 
	int   iAbort;
	int       N_Bframe = input->successive_Bframe;
	float betaSec,betaThird;
	 

	//===== set function for getting reference picture lines =====
	if ((center_x > search_range) && (center_x < img->width -1-search_range-blocksize_x) &&
		(center_y > search_range) && (center_y < img->height-1-search_range-blocksize_y)   )
	{
		get_ref_line = FastLineX;
	}
	else
	{
		get_ref_line = UMVLineX;
	}
	
	//////allocate memory for search state//////////////////////////
	memset(McostState[0],0,(2*search_range+1)*(2*search_range+1)*4);
	
   ///////Threshold defined for early termination///////////////////	
	if(ref>0) 
	{
		if(pred_SAD_ref!=0)
		{
			betaSec = Bsize[blocktype]/(pred_SAD_ref*pred_SAD_ref)-AlphaSec[blocktype];
			betaThird = Bsize[blocktype]/(pred_SAD_ref*pred_SAD_ref)-AlphaThird[blocktype];
		}
		else
		{
			betaSec = 0;
			betaThird = 0;
		}
	}
	else 
	{
		if(blocktype==1)
		{
			if(pred_SAD_space !=0)
			{
				betaSec = Bsize[blocktype]/(pred_SAD_space*pred_SAD_space)-AlphaSec[blocktype];
				betaThird = Bsize[blocktype]/(pred_SAD_space*pred_SAD_space)-AlphaThird[blocktype];
			}
			else
			{
				betaSec = 0;
				betaThird = 0;
			}
		}
		else
		{
			if(pred_SAD_uplayer !=0)
			{
				betaSec = Bsize[blocktype]/(pred_SAD_uplayer*pred_SAD_uplayer)-AlphaSec[blocktype];
				betaThird = Bsize[blocktype]/(pred_SAD_uplayer*pred_SAD_uplayer)-AlphaThird[blocktype];
			}
			else
			{
				betaSec = 0;
				betaThird = 0;
			}
		}
	}
	/*****************************/

	//check the center median predictor
	cand_x = center_x ;
	cand_y = center_y ;
	mcost = MV_COST (lambda_factor, mvshift, cand_x, cand_y, pred_x, pred_y);
	mcost = PartCalMad(ref_pic, orig_pic, get_ref_line,blocksize_y,blocksize_x,blocksize_x4,mcost,min_mcost,cand_x,cand_y);
	McostState[search_range][search_range] = mcost;
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
	
    if(blocktype>1)
	{
		cand_x = pic_pix_x + (pred_MV_uplayer[0]/4);
		cand_y = pic_pix_y + (pred_MV_uplayer[1]/4);
		SEARCH_ONE_PIXEL
		if ((min_mcost-pred_SAD_uplayer)<pred_SAD_uplayer*betaThird)
			goto third_step;
		else if((min_mcost-pred_SAD_uplayer)<pred_SAD_uplayer*betaSec)
			goto sec_step;
	} 

	//coordinate position prediction
	if ((img->number > 1 + ref && ref!=-1) || (list == 1 && (Bframe_ctr%N_Bframe) > 1))  //for debug
	{
		cand_x = pic_pix_x + pred_MV_time[0]/4;
		cand_y = pic_pix_y + pred_MV_time[1]/4;
		SEARCH_ONE_PIXEL
  }

	//prediciton using mV of last ref moiton vector
	if (input->InterlaceCodingOption == FIELD_CODING)
	{
		if ((list==0 && ref > 0) || (img->type == B_SLICE && list == 0 && (ref==0 ||ref==2 ) )) 
			//Notes: for interlace case, ref==1 should be added
		{
			cand_x = pic_pix_x + pred_MV_ref[0]/4;
			cand_y = pic_pix_y + pred_MV_ref[1]/4;
			SEARCH_ONE_PIXEL
		}
	}
	else
	{
		if ((list==0 && ref > 0) || (img->type == B_SLICE && list == 0 && ref==0 )) 
			//Notes: for interlace case, ref==1 should be added
		{
			cand_x = pic_pix_x + pred_MV_ref[0]/4;
			cand_y = pic_pix_y + pred_MV_ref[1]/4;
			SEARCH_ONE_PIXEL
		}
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

	//early termination algrithm, refer to JVT-D016
    EARLY_TERMINATION
	
	if(blocktype>6)
		goto sec_step;
	else
		goto first_step;
	
first_step: //Unsymmetrical-cross search 
	iXMinNow = best_x;
	iYMinNow = best_y;
	
	for(i=1;i<=search_range/2;i++)
	{
		search_step = 2*i - 1;
		cand_x = iXMinNow + search_step;
		cand_y = iYMinNow ;
		SEARCH_ONE_PIXEL		
		cand_x = iXMinNow - search_step;
		cand_y = iYMinNow ;
		SEARCH_ONE_PIXEL
	}
	
	for(i=1;i<=search_range/4;i++)
	{
		search_step = 2*i - 1;
		cand_x = iXMinNow ;
		cand_y = iYMinNow + search_step;
		SEARCH_ONE_PIXEL
		cand_x = iXMinNow ;
		cand_y = iYMinNow - search_step;
		SEARCH_ONE_PIXEL
	}
	//early termination algrithm, refer to JVT-D016
    EARLY_TERMINATION
	
	iXMinNow = best_x;
	iYMinNow = best_y;
    // Uneven Multi-Hexagon-grid Search	
	for(pos=1;pos<25;pos++)
	{
		cand_x = iXMinNow + spiral_search_x[pos];
		cand_y = iYMinNow + spiral_search_y[pos];
		SEARCH_ONE_PIXEL
	}
	//early termination algrithm, refer to JVT-D016
    EARLY_TERMINATION
	
	for(i=1;i<=search_range/4; i++)
	{
		iAbort = 0;   
		for (m = 0; m < 16; m++)
		{
			cand_x = iXMinNow + Big_Hexagon_x[m]*i;
			cand_y = iYMinNow + Big_Hexagon_y[m]*i; 
			SEARCH_ONE_PIXEL1(1)
		}
		if (iAbort)
		{	
			//early termination algrithm, refer to JVT-D016
			EARLY_TERMINATION
		}
	}
sec_step:  //Extended Hexagon-based Search
		  iXMinNow = best_x;
		  iYMinNow = best_y;
		  for(i=0;i<search_range;i++) 
		  {
			  iAbort = 1;   
			  for (m = 0; m < 6; m++)
			  {		
				  cand_x = iXMinNow + Hexagon_x[m];
				  cand_y = iYMinNow + Hexagon_y[m];   
				  SEARCH_ONE_PIXEL1(0)
			  } 
			  if(iAbort)
				  break;
			  iXMinNow = best_x;
			  iYMinNow = best_y;
		  }
third_step: // the third step with a small search pattern
		  iXMinNow = best_x;
		  iYMinNow = best_y;
		  for(i=0;i<search_range;i++) 
		  {
			  iSADLayer = 65536;
			  iAbort = 1;   
			  for (m = 0; m < 4; m++)
			  {		
				  cand_x = iXMinNow + Diamond_x[m];
				  cand_y = iYMinNow + Diamond_y[m];   
				  SEARCH_ONE_PIXEL1(0)
			  } 
			  if(iAbort)
				  break;
			  iXMinNow = best_x;
			  iYMinNow = best_y;
		  }

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
 * \authors: Zhibo Chen
 *           Dept.of EE, Tsinghua Univ.
 * \date   : 2003.4
 ************************************************************************
 */
  int AddUpSADQuarter(int pic_pix_x,int pic_pix_y,int blocksize_x,int blocksize_y,
    int cand_mv_x,int cand_mv_y, pel_t **ref_pic, pel_t**   orig_pic, int Mvmcost, int min_mcost,int useABT)
  {
    int abort_search, y0, x0, rx0, ry0, ry; 
    pel_t *orig_line;
    int   diff[16], *d; 
    int  mcost = Mvmcost;
    int yy,kk,xx;
    int   curr_diff[MB_BLOCK_SIZE][MB_BLOCK_SIZE]; // for ABT SATD calculation
    
    for (y0=0, abort_search=0; y0<blocksize_y && !abort_search; y0+=4)
    {
      ry0 = ((pic_pix_y+y0)<<2) + cand_mv_y;
      
      for (x0=0; x0<blocksize_x; x0+=4)
      {
        rx0 = ((pic_pix_x+x0)<<2) + cand_mv_x;
        d   = diff;
        
        orig_line = orig_pic [y0  ];    ry=ry0;
        *d++      = orig_line[x0  ]  -  PelY_14 (ref_pic, ry, rx0   );
        *d++      = orig_line[x0+1]  -  PelY_14 (ref_pic, ry, rx0+ 4);
        *d++      = orig_line[x0+2]  -  PelY_14 (ref_pic, ry, rx0+ 8);
        *d++      = orig_line[x0+3]  -  PelY_14 (ref_pic, ry, rx0+ 12);
        
        orig_line = orig_pic [y0+1];    ry=ry0+4;
        *d++      = orig_line[x0  ]  -  PelY_14 (ref_pic, ry, rx0   );
        *d++      = orig_line[x0+1]  -  PelY_14 (ref_pic, ry, rx0+ 4);
        *d++      = orig_line[x0+2]  -  PelY_14 (ref_pic, ry, rx0+ 8);
        *d++      = orig_line[x0+3]  -  PelY_14 (ref_pic, ry, rx0+ 12);
        
        orig_line = orig_pic [y0+2];    ry=ry0+8;
        *d++      = orig_line[x0  ]  -  PelY_14 (ref_pic, ry, rx0   );
        *d++      = orig_line[x0+1]  -  PelY_14 (ref_pic, ry, rx0+ 4);
        *d++      = orig_line[x0+2]  -  PelY_14 (ref_pic, ry, rx0+ 8);
        *d++      = orig_line[x0+3]  -  PelY_14 (ref_pic, ry, rx0+ 12);
        
        orig_line = orig_pic [y0+3];    ry=ry0+12;
        *d++      = orig_line[x0  ]  -  PelY_14 (ref_pic, ry, rx0   );
        *d++      = orig_line[x0+1]  -  PelY_14 (ref_pic, ry, rx0+ 4);
        *d++      = orig_line[x0+2]  -  PelY_14 (ref_pic, ry, rx0+ 8);
        *d        = orig_line[x0+3]  -  PelY_14 (ref_pic, ry, rx0+ 12);
        
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
          for (yy=y0,kk=0; yy<y0+4; yy++)
            for (xx=x0; xx<x0+4; xx++, kk++)
              curr_diff[yy][xx] = diff[kk];
        }
      }
    }
    
    return mcost;
  }
  

  int                                                   //  ==> minimum motion cost after search
    FastSubPelBlockMotionSearch (pel_t**   orig_pic,      // <--  original pixel values for the AxB block
    int       ref,           // <--  reference frame (0... or -1 (backward))
	int       list,
    int       pic_pix_x,     // <--  absolute x-coordinate of regarded AxB block
    int       pic_pix_y,     // <--  absolute y-coordinate of regarded AxB block
    int       blocktype,     // <--  block type (1-16x16 ... 7-4x4)
    int       pred_mv_x,     // <--  motion vector predictor (x) in sub-pel units
    int       pred_mv_y,     // <--  motion vector predictor (y) in sub-pel units
    int*      mv_x,          // <--> in: search center (x) / out: motion vector (x) - in pel units
    int*      mv_y,          // <--> in: search center (y) / out: motion vector (y) - in pel units
    int       search_pos2,   // <--  search positions for    half-pel search  (default: 9)
    int       search_pos4,   // <--  search positions for quarter-pel search  (default: 9)
    int       min_mcost,     // <--  minimum motion cost (cost for center or huge value)
    double    lambda,
    int	useABT)        // <--  lagrangian parameter for determining motion cost
  {
	static int Diamond_x[4] = {-1, 0, 1, 0};
	static int Diamond_y[4] = {0, 1, 0, -1};
    int   mcost;
    int   cand_mv_x, cand_mv_y;
    
	int   incr            = list==1 ? ((!img->fld_type)&&(enc_picture!=enc_frame_picture)&&(img->type==B_SLICE)) : (enc_picture==enc_frame_picture)&&(img->type==B_SLICE) ;
	pel_t **ref_pic = listX[list][ref]->imgY_ups;
    
    int   lambda_factor   = LAMBDA_FACTOR (lambda);
    int   mv_shift        = 0;
	int   check_position0 = (blocktype==1 && *mv_x==0 && *mv_y==0 && input->hadamard && !input->rdopt && img->type!=B_SLICE && ref==0);
    int   blocksize_x     = input->blc_size[blocktype][0];
    int   blocksize_y     = input->blc_size[blocktype][1];
    int   pic4_pix_x      = (pic_pix_x << 2);
    int   pic4_pix_y      = (pic_pix_y << 2);
    int   max_pos_x4      = ((img->width -blocksize_x+1)<<2);
    int   max_pos_y4      = ((img->height-blocksize_y+1)<<2);
    
    int   min_pos2        = (input->hadamard ? 0 : 1);
    int   max_pos2        = (input->hadamard ? max(1,search_pos2) : search_pos2);
    int   search_range_dynamic,iXMinNow,iYMinNow,i;
    int   iSADLayer,m,currmv_x,currmv_y,iCurrSearchRange;
    int   search_range = input->search_range;
    int   pred_frac_mv_x,pred_frac_mv_y,abort_search;
    int   mv_cost; 
    
    int   pred_frac_up_mv_x, pred_frac_up_mv_y;
    
    *mv_x <<= 2;
    *mv_y <<= 2;
    if ((pic4_pix_x + *mv_x > 1) && (pic4_pix_x + *mv_x < max_pos_x4 - 2) &&
      (pic4_pix_y + *mv_y > 1) && (pic4_pix_y + *mv_y < max_pos_y4 - 2)   )
    {
      PelY_14 = FastPelY_14;
    }
    else
    {
      PelY_14 = UMVPelY_14;
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
      mcost = AddUpSADQuarter(pic_pix_x,pic_pix_y,blocksize_x,blocksize_y,cand_mv_x,cand_mv_y,ref_pic,orig_pic,mv_cost,min_mcost,useABT);
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
      mcost = AddUpSADQuarter(pic_pix_x,pic_pix_y,blocksize_x,blocksize_y,cand_mv_x,cand_mv_y,ref_pic,orig_pic,mv_cost,min_mcost,useABT);
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
    iCurrSearchRange = 2*search_range_dynamic+1; 
    for(i=0;i<iCurrSearchRange;i++) 
    {
      abort_search=1;
      iSADLayer = 65536;
      for (m = 0; m < 4; m++)
      {
        cand_mv_x = iXMinNow + Diamond_x[m];    
        cand_mv_y = iYMinNow + Diamond_y[m]; 
        
        if(abs(cand_mv_x - *mv_x) <=search_range_dynamic && abs(cand_mv_y - *mv_y)<= search_range_dynamic)
        {
          if(!SearchState[cand_mv_y -*mv_y+ search_range_dynamic][cand_mv_x -*mv_x+ search_range_dynamic])
          {
            mv_cost = MV_COST (lambda_factor, mv_shift, cand_mv_x, cand_mv_y, pred_mv_x, pred_mv_y);		
            mcost = AddUpSADQuarter(pic_pix_x,pic_pix_y,blocksize_x,blocksize_y,cand_mv_x,cand_mv_y,ref_pic,orig_pic,mv_cost,min_mcost,useABT);
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
  * \date   : 2003.4
 ************************************************************************
 */
void   decide_intrabk_SAD()
{
  if (img->type != 0)
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
  if (img->type!=0  && (best_mode == 9 || best_mode == 10))
  {
	  for (i=0; i < 4; i++)
	  {
		  for (j=0; j < 4; j++)
		  {
			  for (k=1; k < 8;k++)
			  {
				  for (ref=0; ref<ref_max;ref++)
				  {
					  all_mincost[(img->pix_x>>2)+i][(img->pix_y>>2)+j][ref][k][0] = 0;   
				  }
			  }
		  }
	  }
  
  }
  return;
}