/*!
 *************************************************************************************
 * \file mc_prediction.c
 *
 * \brief
 *    Motion Compensation
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Alexis Michael Tourapis         <alexismt@ieee.org>
 *
 *************************************************************************************
 */

#include "contributors.h"

#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <math.h>

#include "block.h"
#include "global.h"

#include "macroblock.h"
#include "mc_prediction.h"
#include "refbuf.h"
#include "image.h"
#include "mb_access.h"
#include "me_distortion.h"

/*!
 ************************************************************************
 * \brief
 *    Weighted BiPrediction
 ************************************************************************
 */
static inline void MCWeightedBiPrediction(imgpel** mb_pred, imgpel* l0pred, imgpel *l1pred, 
                                          int block_size_y, int block_x, int block_size_x,
                                          int max_imgpel_value,
                                          short wbp0, short wbp1, short offset, short wp_round, short weight_denom)
{
  int i, j;
  int block_x4 = block_x + block_size_x;

  for   (j = 0; j< block_size_y; j++)
  {
    for (i=block_x; i<block_x4; i++)  
      mb_pred[j][i] = (imgpel) iClip1( max_imgpel_value, 
      ((wbp0 * *l0pred++ + wbp1 * *l1pred++ + wp_round) >> (weight_denom)) + offset); 
  }
}

/*!
 ************************************************************************
 * \brief
 *    Weighted Prediction
 ************************************************************************
 */
static inline void MCWeightedPrediction(imgpel** mb_pred, imgpel* lpred, 
                                        int block_size_y, int block_x, int block_size_x,
                                        int max_imgpel_value,
                                        short wp, short offset, short wp_round, short weight_denom)
{
  int i, j;
  int block_x4 = block_x + block_size_x;

  for   (j = 0; j < block_size_y; j++)
  {
    for (i=block_x; i<block_x4; i++)
      mb_pred[j][i] = (imgpel) iClip1( max_imgpel_value, 
      ((wp * *lpred++  + wp_round) >> weight_denom) + offset);
  }
}

/*!
 ************************************************************************
 * \brief
 *    BiPrediction
 ************************************************************************
 */
static inline void MCBiPrediction(imgpel** mb_pred, imgpel* l0pred, imgpel *l1pred, 
                                  int block_size_y, int block_x, int block_size_x)
{
  int i, j;
  int block_x4 = block_x + block_size_x;

  for   (j = 0; j < block_size_y; j++)
  {
    for (i=block_x; i<block_x4; i++)
      mb_pred[j][i] = (*l0pred++ + *l1pred++ + 1) >> 1;
  }
}

/*!
 ************************************************************************
 * \brief
 *    BiPrediction
 ************************************************************************
 */
static inline void MCPrediction(imgpel** mb_pred, imgpel* lpred, int block_size_y, int block_x, int block_size_x)
{
  int j;
  for (j = 0; j < block_size_y; j++)
  {
    memcpy(&(mb_pred[j][block_x]), lpred, block_size_x * sizeof(imgpel));
    lpred += block_size_x;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Predict Luma block
 ************************************************************************
 */
static inline void OneComponentLumaPrediction ( VideoParameters *p_Vid, //!< video encoding parameters for current picture
                                               imgpel*   mpred,       //!< array of prediction values (row by row)
                                               int    pic_pix_x,      //!< motion shifted horizontal coordinate of block
                                               int    pic_pix_y,      //!< motion shifted vertical   coordinate of block
                                               int    block_size_x,   //!< horizontal block size
                                               int    block_size_y,   //!< vertical block size
                                               StorablePicture *list //!< reference picture list
                                               )
{
  int     j;
  imgpel *ref_line = UMVLine4X (list, pic_pix_y, pic_pix_x);

  for (j = 0; j < block_size_y; j++) 
  {
    memcpy(mpred, ref_line, block_size_x * sizeof(imgpel));
    ref_line += p_Vid->padded_size_x;
    mpred += block_size_x;
  }  
}


/*!
 ************************************************************************
 * \brief
 *    Predict one Luma block
 ************************************************************************
 */
void luma_prediction (Macroblock* currMB, //!< Current Macroblock
                     int   block_x,      //!< relative horizontal block coordinate of block
                     int   block_y,      //!< relative vertical   block coordinate of block
                     int   block_size_x, //!< relative horizontal block coordinate of block
                     int   block_size_y, //!< relative vertical   block coordinate of block
                     int   p_dir,        //!< prediction direction (0=list0, 1=list1, 2=bipred)
                     int   list_mode[2], //!< list prediction mode (1-7, 0=DIRECT)
                     char  *ref_idx,     //!< reference pictures
                     short bipred_me     //!< use bi prediction mv (0=no bipred, 1 = use set 1, 2 = use set 2)
                     )
{
  VideoParameters *p_Vid     = currMB->p_Vid;
  InputParameters *p_Inp     = currMB->p_Inp;
  Slice           *currSlice = currMB->p_slice;
  imgpel l0_pred[MB_PIXELS];
  imgpel l1_pred[MB_PIXELS];

#if (PAD_AFTER)
  int  pic_opix_x   = ((currMB->pix_x + block_x) << 2);
  int  pic_opix_y   = ((currMB->opix_y + block_y) << 2);
#else
  int  pic_opix_x   = ((currMB->pix_x + block_x) << 2) + IMG_PAD_SIZE_TIMES4;
  int  pic_opix_y   = ((currMB->opix_y + block_y) << 2) + IMG_PAD_SIZE_TIMES4;
#endif
  int  bx           = block_x >> 2;
  int  by           = block_y >> 2;
  short****** mv_array = currSlice->all_mv;
  short   *curr_mv = NULL;
  imgpel **mb_pred = currSlice->mb_pred[0];

  int  apply_weights = ( currSlice->weighted_prediction != 0 );

  if (bipred_me && ref_idx[0] == 0 && ref_idx[1] == 0 && p_dir == 2 && is_bipred_enabled(p_Inp, list_mode[0]) && is_bipred_enabled(p_Inp, list_mode[1]))
    mv_array = currSlice->bipred_mv[bipred_me - 1]; 

  switch (p_dir)
  {
  case 0:
    curr_mv = mv_array[LIST_0][(short) ref_idx[0]][list_mode[0]][by][bx];
    OneComponentLumaPrediction (p_Vid, l0_pred, pic_opix_x + curr_mv[0], pic_opix_y + curr_mv[1], block_size_x, block_size_y, p_Vid->listX[LIST_0 + currMB->list_offset][(short) ref_idx[0]]);
    break;
  case 1:
    curr_mv = mv_array[LIST_1][(short) ref_idx[1]][list_mode[1]][by][bx];
    OneComponentLumaPrediction (p_Vid, l1_pred, pic_opix_x + curr_mv[0], pic_opix_y + curr_mv[1], block_size_x, block_size_y, p_Vid->listX[LIST_1 + currMB->list_offset][(short)ref_idx[1]]);
    break;
  case 2:
    curr_mv = mv_array[LIST_0][(short) ref_idx[0]][list_mode[0]][by][bx];
    OneComponentLumaPrediction (p_Vid, l0_pred, pic_opix_x + curr_mv[0], pic_opix_y + curr_mv[1], block_size_x, block_size_y, p_Vid->listX[LIST_0 + currMB->list_offset][(short)ref_idx[0]]);
    curr_mv = mv_array[LIST_1][(short) ref_idx[1]][list_mode[1]][by][bx];
    OneComponentLumaPrediction (p_Vid, l1_pred, pic_opix_x + curr_mv[0], pic_opix_y + curr_mv[1], block_size_x, block_size_y, p_Vid->listX[LIST_1 + currMB->list_offset][(short)ref_idx[1]]);
    break;
  default:
    break;
  }

  if (apply_weights)
  {
    if (p_dir==2)
    {
      MCWeightedBiPrediction(&mb_pred[block_y], l0_pred, l1_pred, block_size_y, block_x, block_size_x, 
        p_Vid->max_imgpel_value,
        currSlice->wbp_weight[0][(short)ref_idx[0]][(short)ref_idx[1]][0], currSlice->wbp_weight[1][(short)ref_idx[0]][(short)ref_idx[1]][0],
        (currSlice->wp_offset[0][(short)ref_idx[0]][0] + currSlice->wp_offset[1][(short)ref_idx[1]][0] + 1)>>1, 
         (currSlice->wp_luma_round << 1), currSlice->luma_log_weight_denom + 1);
    }
    else if (p_dir==0)
    {
      MCWeightedPrediction(&mb_pred[block_y], l0_pred, block_size_y, block_x, block_size_x,
        p_Vid->max_imgpel_value,
        currSlice->wp_weight[0][(short)ref_idx[0]][0], currSlice->wp_offset[0][(short)ref_idx[0]][0], currSlice->wp_luma_round, currSlice->luma_log_weight_denom);
    }
    else // (p_dir==1)
    {
      MCWeightedPrediction(&mb_pred[block_y], l1_pred, block_size_y, block_x, block_size_x,
        p_Vid->max_imgpel_value,
        currSlice->wp_weight[1][(short)ref_idx[1]][0], currSlice->wp_offset[1][(short)ref_idx[1]][0], currSlice->wp_luma_round, currSlice->luma_log_weight_denom);
    }
  }
  else
  {
    if (p_dir==2)
    {
      MCBiPrediction(&mb_pred[block_y], l0_pred, l1_pred, block_size_y, block_x, block_size_x);    
    }
    else if (p_dir==0)
    {
      MCPrediction(&mb_pred[block_y], l0_pred, block_size_y, block_x, block_size_x);
    }
    else // (p_dir==1)
    {
      MCPrediction(&mb_pred[block_y], l1_pred, block_size_y, block_x, block_size_x);
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Predict one Luma block
 ************************************************************************
 */
void luma_prediction_bi (Macroblock* currMB, //!< Current Macroblock
                       int   block_x,      //!< relative horizontal block coordinate of 4x4 block
                       int   block_y,      //!< relative vertical   block coordinate of 4x4 block
                       int   block_size_x, //!< horizontal block size
                       int   block_size_y, //!< vertical   block size
                       int   l0_mode,      //!< list0 prediction mode (1-7, 0=DIRECT if l1_mode=0)
                       int   l1_mode,      //!< list1 prediction mode (1-7, 0=DIRECT if l0_mode=0)
                       short l0_ref_idx,   //!< reference frame for list0 prediction (-1: Intra4x4 pred. with l0_mode)
                       short l1_ref_idx,   //!< reference frame for list1 prediction 
                       int   list          //!< current list for prediction.
                       )
{
  VideoParameters *p_Vid     = currMB->p_Vid;
  Slice           *currSlice = currMB->p_slice;

  imgpel l0_pred[MB_PIXELS];
  imgpel l1_pred[MB_PIXELS];
#if (PAD_AFTER)
  int  pic_opix_x   = ((currMB->pix_x + block_x) << 2);
  int  pic_opix_y   = ((currMB->opix_y + block_y) << 2);
#else
  int  pic_opix_x   = ((currMB->pix_x + block_x) << 2) + IMG_PAD_SIZE_TIMES4;
  int  pic_opix_y   = ((currMB->opix_y + block_y) << 2) + IMG_PAD_SIZE_TIMES4;
#endif
  int  bx        = block_x >> 2;
  int  by        = block_y >> 2;

  int  apply_weights = ( currSlice->weighted_prediction != 0 );

  short   ******mv_array = currSlice->bipred_mv[list]; 
  short   *mv_arrayl0 = mv_array[LIST_0][l0_ref_idx][l0_mode][by][bx];
  short   *mv_arrayl1 = mv_array[LIST_1][l1_ref_idx][l1_mode][by][bx];
  imgpel **mb_pred = currSlice->mb_pred[0];

  OneComponentLumaPrediction (p_Vid, l0_pred, pic_opix_x + mv_arrayl0[0], pic_opix_y + mv_arrayl0[1], block_size_x, block_size_y, p_Vid->listX[0+currMB->list_offset][l0_ref_idx]);
  OneComponentLumaPrediction (p_Vid, l1_pred, pic_opix_x + mv_arrayl1[0], pic_opix_y + mv_arrayl1[1], block_size_x, block_size_y, p_Vid->listX[1+currMB->list_offset][l1_ref_idx]);

  if (apply_weights)
  {
    MCWeightedBiPrediction(&mb_pred[block_y], l0_pred, l1_pred, block_size_y, block_x, block_size_x, 
      p_Vid->max_imgpel_value,
      currSlice->wbp_weight[0][l0_ref_idx][l1_ref_idx][0], currSlice->wbp_weight[1][l0_ref_idx][l1_ref_idx][0],
      (currSlice->wp_offset[0][l0_ref_idx][0] + currSlice->wp_offset[1][l1_ref_idx][0] + 1)>>1, 
      (currSlice->wp_luma_round << 1), currSlice->luma_log_weight_denom + 1);
  }
  else
  {
    MCBiPrediction(&mb_pred[block_y], l0_pred, l1_pred, block_size_y, block_x, block_size_x);
  }
}


/*!
 ************************************************************************
 * \brief
 *    Predict (on-the-fly) one component of a chroma 4x4 block
 ************************************************************************
 */
void OneComponentChromaPrediction4x4_regenerate (
                                 Macroblock  *currMB,    //!< Current Macroblock
                                 imgpel*     mpred,      //!< array to store prediction values
                                 int         block_c_x,  //!< horizontal pixel coordinate of 4x4 block
                                 int         block_c_y,  //!< vertical   pixel coordinate of 4x4 block
                                 short***    mv,         //!< motion vector array
                                 StorablePicture *list,  //!< image components (color planes)
                                 int         uv)         //!< chroma component
{
  VideoParameters *p_Vid = currMB->p_Vid;
  int     i, j, ii, jj, ii0, jj0, ii1, jj1, if0, if1, jf0, jf1;
  short*  mvb;

  int     f1_x = 64/p_Vid->mb_cr_size_x;
  int     f2_x=f1_x-1;

  int     f1_y = 64/p_Vid->mb_cr_size_y;
  int     f2_y=f1_y-1;

  int     f3=f1_x*f1_y, f4=f3>>1;
  int     list_offset = p_Vid->mb_data[currMB->mbAddrX].list_offset;
  int     max_y_cr = (int) (list_offset ? (p_Vid->height_cr >> 1) - 1 : p_Vid->height_cr - 1);
  int     max_x_cr = (int) (p_Vid->width_cr - 1);
  int     jjx, iix;
  int     mb_cr_y_div4 = p_Vid->mb_cr_size_y>>2;
  int     mb_cr_x_div4 = p_Vid->mb_cr_size_x>>2;
  int     jpos;

  imgpel** refimage = list->imgUV[uv];

  for (j=block_c_y; j < block_c_y + BLOCK_SIZE; j++)
  {
    jjx = j/mb_cr_y_div4;
    jpos = (j + currMB->opix_c_y)*f1_y;

    for (i=block_c_x; i < block_c_x + BLOCK_SIZE; i++)
    {
      iix = i/mb_cr_x_div4;
      mvb  = mv [jjx][iix];

      ii   = (i + currMB->pix_c_x)*f1_x + mvb[0];
      jj   = jpos + mvb[1];

      if (p_Vid->active_sps->chroma_format_idc == 1)
        jj  += list->chroma_vector_adjustment;

      ii0  = iClip3 (0, max_x_cr, ii/f1_x);
      jj0  = iClip3 (0, max_y_cr, jj/f1_y);
      ii1  = iClip3 (0, max_x_cr, (ii+f2_x)/f1_x);
      jj1  = iClip3 (0, max_y_cr, (jj+f2_y)/f1_y);

      if1  = (ii&f2_x);  if0 = f1_x-if1;
      jf1  = (jj&f2_y);  jf0 = f1_y-jf1;

      *mpred++ = (imgpel) (
       (if0 * jf0 * refimage[jj0][ii0] +
        if1 * jf0 * refimage[jj0][ii1] +
        if0 * jf1 * refimage[jj1][ii0] +
        if1 * jf1 * refimage[jj1][ii1] + f4) / f3);
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Retrieve one component of a chroma 4x4 block from the buffer
 ************************************************************************
 */
void OneComponentChromaPrediction4x4_retrieve (
                                 Macroblock  *currMB,    //!< Current Macroblock
                                 imgpel*     mpred,      //!< array to store prediction values
                                 int         block_c_x,  //!< horizontal pixel coordinate of 4x4 block
                                 int         block_c_y,  //!< vertical   pixel coordinate of 4x4 block
                                 short***    mv,         //!< motion vector array
                                 StorablePicture *list,  //!< image components (color planes)
                                 int         uv)         //!< chroma component
{
  VideoParameters *p_Vid = currMB->p_Vid;
  int     j, ii, jj;
  short*  mvb;

  int     jjx;
  int     right_shift_x = 4 - p_Vid->chroma_shift_x;
  int     right_shift_y = 4 - p_Vid->chroma_shift_y;
  int     jpos;

  int     pos_x1 = block_c_x >> right_shift_x;
  int     pos_x2 = (block_c_x + 2) >> right_shift_x;

#if (PAD_AFTER)
  int     ipos1 = ((block_c_x + currMB->pix_c_x    ) << p_Vid->chroma_shift_x);
  int     ipos2 = ((block_c_x + currMB->pix_c_x + 2) << p_Vid->chroma_shift_x);
  int jj_chroma = ((p_Vid->active_sps->chroma_format_idc == 1) ? list->chroma_vector_adjustment : 0);
#else
  int     ipos1 = ((block_c_x + currMB->pix_c_x    ) << p_Vid->chroma_shift_x) + IMG_PAD_SIZE_TIMES4;
  int     ipos2 = ((block_c_x + currMB->pix_c_x + 2) << p_Vid->chroma_shift_x) + IMG_PAD_SIZE_TIMES4;
  int jj_chroma = ((p_Vid->active_sps->chroma_format_idc == 1) ? list->chroma_vector_adjustment : 0) + IMG_PAD_SIZE_TIMES4;
#endif

  imgpel *line_ptr;

  for (j=block_c_y; j < block_c_y + BLOCK_SIZE; j++)
  {
    jjx = j >> right_shift_y; // translate into absolute block (luma) coordinates

    jpos = ( (j + currMB->opix_c_y) << p_Vid->chroma_shift_y ) + jj_chroma;

    mvb  = mv [jjx][pos_x1];

    ii   = ipos1 + mvb[0];
    jj   = jpos  + mvb[1];

    line_ptr = UMVLine8X_chroma ( list, uv + 1, jj, ii);
    *mpred++ = *line_ptr++;
    *mpred++ = *line_ptr;

    mvb  = mv [jjx][pos_x2];

    ii   = ipos2 + mvb[0];
    jj   = jpos  + mvb[1];

    line_ptr = UMVLine8X_chroma ( list, uv + 1, jj, ii);
    *mpred++ = *line_ptr++;
    *mpred++ = *line_ptr;
  }
}


/*!
 ************************************************************************
 * \brief
 *    Retrieve one component of a chroma block from the buffer
 ************************************************************************
 */
static inline 
void OneComponentChromaPrediction (VideoParameters *p_Vid, //!< video encoding parameters for current picture
                                   imgpel* mpred,      //!< array to store prediction values
                                   int    pic_pix_x,      //!< motion shifted horizontal coordinate of block
                                   int    pic_pix_y,      //!< motion shifted vertical  block
                                   int    block_size_x,   //!< horizontal block size
                                   int    block_size_y,   //!< vertical block size                                      
                                   StorablePicture *list, //!< reference picture list
                                   int    uv)         //!< chroma component
{
  int     j;
  imgpel *ref_line = UMVLine4Xcr (list, uv + 1, pic_pix_y, pic_pix_x);

  for (j = 0; j < block_size_y; j++) 
  {
    memcpy(mpred, ref_line, block_size_x * sizeof(imgpel));
    ref_line += p_Vid->cr_padded_size_x;
    mpred += block_size_x;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Predict an intra chroma 4x4 block
 ************************************************************************
 */
static inline
void IntraChromaPrediction4x4 (Macroblock* currMB, //!< Current Macroblock
                               int  uv,            //!< color component
                               int  block_x,       //!< relative horizontal block coordinate of 4x4 block
                               int  block_y)       //!< relative vertical   block coordinate of 4x4 block
{
  int j;
  Slice *currSlice = currMB->p_slice;
  imgpel **mb_pred        = currSlice->mb_pred[ uv ];
  imgpel **curr_mpr_16x16 = currSlice->mpr_16x16[uv][ (short) currMB->c_ipred_mode];

  //===== prediction =====
  for (j=block_y; j<block_y + BLOCK_SIZE; j++)
    memcpy(&mb_pred[j][block_x],&curr_mpr_16x16[j][block_x], BLOCK_SIZE * sizeof(imgpel));
}

/*!
 ************************************************************************
 * \brief
 *    Predict one chroma block
 ************************************************************************
 */
void chroma_prediction (Macroblock* currMB, // <-- Current Macroblock
                       int   uv,            // <-- colour component
                       int   block_x,       // <-- relative horizontal block coordinate of block
                       int   block_y,       // <-- relative vertical   block coordinate of block
                       int   block_size_x,  // <-- relative horizontal block coordinate of block
                       int   block_size_y,  // <-- relative vertical   block coordinate of block                        
                       int   p_dir,         // <-- prediction direction (0=list0, 1=list1, 2=bipred)
                       int   l0_mode,       // <-- list0  prediction mode (1-7, 0=DIRECT if l1_mode=0)
                       int   l1_mode,       // <-- list1 prediction mode (1-7, 0=DIRECT if l0_mode=0)
                       short l0_ref_idx,    // <-- reference frame for list0 prediction (if (<0) -> intra prediction)
                       short l1_ref_idx,    // <-- reference frame for list1 prediction 
                       short bipred_me      // <-- use bi prediction mv (0=no bipred, 1 = use set 1, 2 = use set 2)
                       )    
{
  VideoParameters *p_Vid     = currMB->p_Vid;
  InputParameters *p_Inp     = currMB->p_Inp;
  Slice           *currSlice = currMB->p_slice;

  imgpel l0_pred[MB_PIXELS];
  imgpel l1_pred[MB_PIXELS];

#if (PAD_AFTER)
  int  pic_opix_x   = ((currMB->pix_c_x + block_x) << 2);
  int  pic_opix_y   = ((currMB->opix_c_y + block_y) << 2);
#else
  int  pic_opix_x   = ((currMB->pix_c_x + block_x) << 2) + IMG_PAD_SIZE_TIMES4;
  int  pic_opix_y   = ((currMB->opix_c_y + block_y) << 2) + IMG_PAD_SIZE_TIMES4;
#endif
  int  bx           = block_x >> 2;
  int  by           = block_y >> 2;
  short****** mv_array = currSlice->all_mv;    
  int uv_comp = uv + 1;
  imgpel **mb_pred = currSlice->mb_pred[ uv_comp];

  int  apply_weights = ( currSlice->weighted_prediction != 0 );

  if (bipred_me && l0_ref_idx == 0 && l1_ref_idx == 0 && p_dir == 2 && is_bipred_enabled(p_Inp, l0_mode)  && is_bipred_enabled(p_Inp, l1_mode))
    mv_array = currSlice->bipred_mv[bipred_me - 1]; 

  //===== INTRA PREDICTION =====
  if (p_dir==-1)
  {
    IntraChromaPrediction4x4 (currMB, uv_comp, block_x, block_y);
    return;
  }

  //===== INTER PREDICTION =====
  switch (p_dir)
  {
  case 0:
    OneComponentChromaPrediction (p_Vid, l0_pred, pic_opix_x + mv_array[LIST_0][l0_ref_idx][l0_mode][by][bx][0], pic_opix_y + mv_array[LIST_0][l0_ref_idx][l0_mode][by][bx][1], block_size_x, block_size_y, p_Vid->listX[0+currMB->list_offset][l0_ref_idx], uv);
    break;
  case 1: 
    OneComponentChromaPrediction (p_Vid, l1_pred, pic_opix_x + mv_array[LIST_1][l1_ref_idx][l1_mode][by][bx][0], pic_opix_y + mv_array[LIST_1][l1_ref_idx][l1_mode][by][bx][1], block_size_x, block_size_y, p_Vid->listX[1+currMB->list_offset][l1_ref_idx], uv);
    break;
  case 2:
    OneComponentChromaPrediction (p_Vid, l0_pred, pic_opix_x + mv_array[LIST_0][l0_ref_idx][l0_mode][by][bx][0], pic_opix_y + mv_array[LIST_0][l0_ref_idx][l0_mode][by][bx][1], block_size_x, block_size_y, p_Vid->listX[0+currMB->list_offset][l0_ref_idx], uv);
    OneComponentChromaPrediction (p_Vid, l1_pred, pic_opix_x + mv_array[LIST_1][l1_ref_idx][l1_mode][by][bx][0], pic_opix_y + mv_array[LIST_1][l1_ref_idx][l1_mode][by][bx][1], block_size_x, block_size_y, p_Vid->listX[1+currMB->list_offset][l1_ref_idx], uv);
    break;
  default:
    break;
  }

  if (apply_weights)
  {
    if (p_dir==2)
    {
      MCWeightedBiPrediction(&mb_pred[block_y], l0_pred, l1_pred, block_size_y, block_x, block_size_x, 
      p_Vid->max_pel_value_comp[1],
        currSlice->wbp_weight[0][l0_ref_idx][l1_ref_idx][uv_comp], currSlice->wbp_weight[1][l0_ref_idx][l1_ref_idx][uv_comp],
        (currSlice->wp_offset[0][l0_ref_idx][uv_comp] + currSlice->wp_offset[1][l1_ref_idx][uv_comp] + 1)>>1, 
         (currSlice->wp_chroma_round << 1), currSlice->chroma_log_weight_denom + 1);
    }
    else if (p_dir==0)
    {
      MCWeightedPrediction(&mb_pred[block_y], l0_pred, block_size_y, block_x, block_size_x,
      p_Vid->max_pel_value_comp[1],
        currSlice->wp_weight[0][l0_ref_idx][uv_comp], currSlice->wp_offset[0][l0_ref_idx][uv_comp], currSlice->wp_chroma_round, currSlice->chroma_log_weight_denom );
    }
    else // (p_dir==1)
    {
      MCWeightedPrediction(&mb_pred[block_y], l1_pred, block_size_y, block_x, block_size_x,
      p_Vid->max_pel_value_comp[1],
        currSlice->wp_weight[1][l1_ref_idx][uv_comp], currSlice->wp_offset[1][l1_ref_idx][uv_comp], currSlice->wp_chroma_round, currSlice->chroma_log_weight_denom );
    }
  }
  else
  {
    if (p_dir==2)
    {
      MCBiPrediction(&mb_pred[block_y], l0_pred, l1_pred, block_size_y, block_x, block_size_x);
    }
    else if (p_dir==0)
    {
      MCPrediction(&mb_pred[block_y], l0_pred, block_size_y, block_x, block_size_x);
    }
    else // (p_dir==1)
    {
      MCPrediction(&mb_pred[block_y], l1_pred, block_size_y, block_x, block_size_x);
    }
  }
}



/*!
 ************************************************************************
 * \brief
 *    Predict one chroma 4x4 block
 ************************************************************************
 */
void chroma_prediction_4x4 (Macroblock* currMB,  // <-- Current Macroblock
                           int   uv,           // <-- colour component
                           int   block_x,      // <-- relative horizontal block coordinate of 4x4 block
                           int   block_y,      // <-- relative vertical   block coordinate of 4x4 block
                           int   p_dir,        // <-- prediction direction (0=list0, 1=list1, 2=bipred)
                           int   l0_mode,      // <-- list0  prediction mode (1-7, 0=DIRECT if l1_mode=0)
                           int   l1_mode,      // <-- list1 prediction mode (1-7, 0=DIRECT if l0_mode=0)
                           short l0_ref_idx,   // <-- reference frame for list0 prediction (if (<0) -> intra prediction)
                           short l1_ref_idx,   // <-- reference frame for list1 prediction 
                           short bipred_me     // <-- use bi prediction mv (0=no bipred, 1 = use set 1, 2 = use set 2)
                           )   
{
  VideoParameters *p_Vid     = currMB->p_Vid;
  InputParameters *p_Inp     = currMB->p_Inp;
  Slice           *currSlice = currMB->p_slice;

  imgpel l0_pred[MB_PIXELS];
  imgpel l1_pred[MB_PIXELS];

  short****** mv_array = currSlice->all_mv;
  int uv_comp = uv + 1;
  imgpel **mb_pred = currSlice->mb_pred[uv_comp];
  int     list_offset = currMB->list_offset;
  
  int  apply_weights = ( currSlice->weighted_prediction != 0 );

  if (bipred_me && l0_ref_idx == 0 && l1_ref_idx == 0 && p_dir == 2 && is_bipred_enabled(p_Inp, l0_mode)  && is_bipred_enabled(p_Inp, l1_mode) )
    mv_array = currSlice->bipred_mv[bipred_me - 1]; 
  //===== INTRA PREDICTION =====
  if (p_dir==-1)
  {
    IntraChromaPrediction4x4 (currMB, uv_comp, block_x, block_y);
    return;
  }

  //===== INTER PREDICTION =====
  switch (p_dir)
  {
  case 0: // LIST_0
    p_Vid->OneComponentChromaPrediction4x4 (currMB, l0_pred, block_x, block_y, mv_array[LIST_0][l0_ref_idx][l0_mode], p_Vid->listX[LIST_0 + list_offset][l0_ref_idx], uv);
    break;
  case 1: // LIST_1
    p_Vid->OneComponentChromaPrediction4x4 (currMB, l1_pred, block_x, block_y, mv_array[LIST_1][l1_ref_idx][l1_mode], p_Vid->listX[LIST_1 + list_offset][l1_ref_idx], uv);
    break;
  case 2: // BI_PRED
    p_Vid->OneComponentChromaPrediction4x4 (currMB, l0_pred, block_x, block_y, mv_array[LIST_0][l0_ref_idx][l0_mode], p_Vid->listX[LIST_0 + list_offset][l0_ref_idx], uv);
    p_Vid->OneComponentChromaPrediction4x4 (currMB, l1_pred, block_x, block_y, mv_array[LIST_1][l1_ref_idx][l1_mode], p_Vid->listX[LIST_1 + list_offset][l1_ref_idx], uv);
    break;
  default:
    break;
  }

  if (apply_weights)
  {
    if (p_dir==2)
    {
      MCWeightedBiPrediction(&mb_pred[block_y], l0_pred, l1_pred, BLOCK_SIZE, block_x, BLOCK_SIZE, 
        p_Vid->max_pel_value_comp[1],
        currSlice->wbp_weight[0][l0_ref_idx][l1_ref_idx][uv_comp], currSlice->wbp_weight[1][l0_ref_idx][l1_ref_idx][uv_comp],
        (currSlice->wp_offset[0][l0_ref_idx][uv_comp] + currSlice->wp_offset[1][l1_ref_idx][uv_comp] + 1)>>1, 
        (currSlice->wp_chroma_round << 1), currSlice->chroma_log_weight_denom + 1);

    }
    else if (p_dir==0)
    {
      MCWeightedPrediction(&mb_pred[block_y], l0_pred, BLOCK_SIZE, block_x, BLOCK_SIZE,
      p_Vid->max_pel_value_comp[1],
        currSlice->wp_weight[0][l0_ref_idx][uv_comp], currSlice->wp_offset[0][l0_ref_idx][uv_comp], currSlice->wp_chroma_round, currSlice->chroma_log_weight_denom );
    }
    else // (p_dir==1)
    {
      MCWeightedPrediction(&mb_pred[block_y], l1_pred, BLOCK_SIZE, block_x, BLOCK_SIZE,
      p_Vid->max_pel_value_comp[1],
        currSlice->wp_weight[1][l1_ref_idx][uv_comp], currSlice->wp_offset[1][l1_ref_idx][uv_comp], currSlice->wp_chroma_round, currSlice->chroma_log_weight_denom );
    }
  }
  else
  {
    if (p_dir==2)
    {
      MCBiPrediction(&mb_pred[block_y], l0_pred, l1_pred, BLOCK_SIZE, block_x, BLOCK_SIZE);
    }
    else if (p_dir==0)
    {
      MCPrediction(&mb_pred[block_y], l0_pred, BLOCK_SIZE, block_x, BLOCK_SIZE);
    }
    else // (p_dir==1)
    {
      MCPrediction(&mb_pred[block_y], l1_pred, BLOCK_SIZE, block_x, BLOCK_SIZE);
    }
  }
}


/*!
 ************************************************************************
 * \brief
 *    Intra prediction of the chrminance layers of one macroblock
 ************************************************************************
 */
void intra_chroma_prediction (Macroblock *currMB, int *mb_up, int *mb_left, int*mb_up_left)
{
  int s, s0, s1, s2, s3, i, j, k;
  int ih,iv, ib, ic, iaa;
  int      uv;
  int      blk_x, blk_y;
  int      b8,b4;
  imgpel**  image;
  imgpel   vline[16];
  int      block_x, block_y;
  int      mb_available_up;
  int      mb_available_left[2];
  int      mb_available_up_left;

  int      mode;
  int      best_mode = DC_PRED_8;  //just an initilaization here, should always be overwritten
  distblk      cost;
  distblk      min_cost;
  PixelPos up;        //!< pixel position  p(0,-1)
  PixelPos left[17];  //!< pixel positions p(-1, -1..15)
  int diff  [16];

  Slice *currSlice = currMB->p_slice;
  VideoParameters *p_Vid = currSlice->p_Vid;
  InputParameters *p_Inp = currSlice->p_Inp;
  int      cr_MB_x = p_Vid->mb_cr_size_x;
  int      cr_MB_y = p_Vid->mb_cr_size_y;
  imgpel **cur_pred = NULL;
  imgpel ***curr_mpr_16x16 = NULL;
  imgpel *hline = NULL;
  imgpel *img_org = NULL, *img_prd = NULL;

  int      yuv = p_Vid->yuv_format - 1;
  int      dc_pred_value_chroma = p_Vid->dc_pred_value_comp[1];
  int      max_imgpel_value_uv  = p_Vid->max_pel_value_comp[1];

  static const int block_pos[3][4][4]= //[yuv][b8][b4]
  {
    { {0, 1, 2, 3},{0, 0, 0, 0},{0, 0, 0, 0},{0, 0, 0, 0}},
    { {0, 1, 2, 3},{2, 3, 2, 3},{0, 0, 0, 0},{0, 0, 0, 0}},
    { {0, 1, 2, 3},{1, 1, 3, 3},{2, 3, 2, 3},{3, 3, 3, 3}}
  };

  for (i=0;i<cr_MB_y+1;i++)
  {
    p_Vid->getNeighbour(currMB, -1 , i-1 , p_Vid->mb_size[IS_CHROMA], &left[i]);
  }
  p_Vid->getNeighbour(currMB, 0 , -1 , p_Vid->mb_size[IS_CHROMA], &up);

  mb_available_up                             = up.available;
  mb_available_up_left                        = left[0].available;
  mb_available_left[0] = mb_available_left[1] = left[1].available;

  if(p_Inp->UseConstrainedIntraPred)
  {
    mb_available_up = up.available ? p_Vid->intra_block[up.mb_addr] : 0;
    for (i=0, mb_available_left[0]=1; i<(cr_MB_y>>1);i++)
      mb_available_left[0]  &= left[i+1].available ? p_Vid->intra_block[left[i+1].mb_addr]: 0;
    for (i=(cr_MB_y>>1), mb_available_left[1]=1; i<cr_MB_y;i++)
      mb_available_left[1] &= left[i+1].available ? p_Vid->intra_block[left[i+1].mb_addr]: 0;
    mb_available_up_left = left[0].available ? p_Vid->intra_block[left[0].mb_addr]: 0;
  }

  if (mb_up)
    *mb_up = mb_available_up;
  if (mb_left)
    *mb_left = mb_available_left[0] && mb_available_left[1];
  if (mb_up_left)
    *mb_up_left = mb_available_up_left;


  // compute all chroma intra prediction modes for both U and V
  for (uv=0; uv<2; uv++)
  {
    image          = p_Vid->enc_picture->imgUV[uv];
    curr_mpr_16x16 = currSlice->mpr_16x16[uv + 1];

    // DC prediction
    for(b8=0; b8<p_Vid->num_blk8x8_uv >> 1;b8++)
    {
      for (b4 = 0; b4 < 4; b4++)
      {
        block_y = subblk_offset_y[yuv][b8][b4];
        block_x = subblk_offset_x[yuv][b8][b4];
        blk_x = block_x;
        blk_y = block_y + 1;

        s = dc_pred_value_chroma;
        s0 = s1 = s2 = s3 = 0;

        //===== get prediction value =====
        switch (block_pos[yuv][b8][b4])
        {
        case 0:  //===== TOP LEFT =====
          if      (mb_available_up)       
            for (i = blk_x; i < (blk_x + 4); i++)  
              s0 += image[up.pos_y][up.pos_x + i];
          if      (mb_available_left[0])  
            for (i = blk_y; i < (blk_y + 4);i++)  
              s2 += image[left[i].pos_y][left[i].pos_x];
          if      (mb_available_up && mb_available_left[0])  
            s  = (s0 + s2 + 4) >> 3;
          else if (mb_available_up)                          
            s  = (s0   + 2) >> 2;
          else if (mb_available_left[0])                     
            s  = (s2   +2) >> 2;
          break;
        case 1: //===== TOP RIGHT =====
          if      (mb_available_up)       
            for (i=blk_x;i<(blk_x+4);i++)  
              s1 += image[up.pos_y][up.pos_x + i];
          else if (mb_available_left[0])  
            for (i=blk_y;i<(blk_y+4);i++) 
              s2 += image[left[i].pos_y][left[i].pos_x];
          if      (mb_available_up)       
            s  = (s1   +2) >> 2;
          else if (mb_available_left[0])                    
            s  = (s2   +2) >> 2;
          break;
        case 2: //===== BOTTOM LEFT =====
          if      (mb_available_left[1])  
            for (i=blk_y;i<(blk_y+4);i++)  
              s3 += image[left[i].pos_y][left[i].pos_x];
          else if (mb_available_up)       
            for (i=blk_x;i<(blk_x+4);i++)  
              s0 += image[up.pos_y][up.pos_x + i];
          if      (mb_available_left[1])                     
            s  = (s3   +2) >> 2;
          else if (mb_available_up)                          
            s  = (s0   +2) >> 2;
          break;
        case 3: //===== BOTTOM RIGHT =====
          if      (mb_available_up)       
            for (i=blk_x;i<(blk_x+4);i++)  
              s1 += image[up.pos_y][up.pos_x + i];
          if      (mb_available_left[1])  
            for (i=blk_y;i<(blk_y+4);i++)  
              s3 += image[left[i].pos_y][left[i].pos_x];
          if      (mb_available_up && mb_available_left[1])  
            s  = (s1+s3+4) >> 3;
          else if (mb_available_up)                          
            s  = (s1   +2) >> 2;
          else if (mb_available_left[1])                     
            s  = (s3   +2) >> 2;
          break;
        }

        //===== prediction =====
        cur_pred = curr_mpr_16x16[DC_PRED_8];
        for (j=block_y; j<block_y+4; j++)
          for (i=block_x; i<block_x+4; i++)
          {
            cur_pred[j][i] = (imgpel) s;
          }
      }
    }

    // vertical prediction    
    if (mb_available_up)
    {
      cur_pred = curr_mpr_16x16[VERT_PRED_8];
      //memcpy(hline,&image[up.pos_y][up.pos_x], cr_MB_x * sizeof(imgpel));
      hline = &image[up.pos_y][up.pos_x];
      for (j=0; j<cr_MB_y; j++)
        memcpy(cur_pred[j], hline, cr_MB_x * sizeof(imgpel));
    }

    // horizontal prediction
    if (mb_available_left[0] && mb_available_left[1])
    {
      cur_pred = curr_mpr_16x16[HOR_PRED_8];
      for (i=0; i<cr_MB_y; i++)
        vline[i] = image[left[i+1].pos_y][left[i+1].pos_x];
      for (j=0; j<cr_MB_y; j++)
      {
        int predictor = vline[j];
        for (i=0; i<cr_MB_x; i++)        
          cur_pred[j][i] = (imgpel) predictor;
      }
    }

    // plane prediction
    if (mb_available_left[0] && mb_available_left[1] && mb_available_up && mb_available_up_left)
    {
      ih = (cr_MB_x>>1)*(hline[cr_MB_x-1] - image[left[0].pos_y][left[0].pos_x]);
      for (i=0;i<(cr_MB_x>>1)-1;i++)
        ih += (i+1)*(hline[(cr_MB_x>>1)+i] - hline[(cr_MB_x>>1)-2-i]);

      iv = (cr_MB_y>>1)*(vline[cr_MB_y-1] - image[left[0].pos_y][left[0].pos_x]);
      for (i=0;i<(cr_MB_y>>1)-1;i++)
        iv += (i+1)*(vline[(cr_MB_y>>1)+i] - vline[(cr_MB_y>>1)-2-i]);

      ib= ((cr_MB_x == 8?17:5)*ih+2*cr_MB_x)>>(cr_MB_x == 8?5:6);
      ic= ((cr_MB_y == 8?17:5)*iv+2*cr_MB_y)>>(cr_MB_y == 8?5:6);

      iaa=16*(hline[cr_MB_x-1] + vline[cr_MB_y-1]);
      cur_pred = curr_mpr_16x16[PLANE_8];
      for (j=0; j<cr_MB_y; j++)
        for (i=0; i<cr_MB_x; i++)
          cur_pred[j][i]= (imgpel) iClip1( max_imgpel_value_uv, (iaa+(i-(cr_MB_x>>1)+1)*ib+(j-(cr_MB_y>>1)+1)*ic+16)>>5);
    }
  }

  if (!p_Inp->rdopt)      // the rd-opt part does not work correctly (see encode_one_macroblock)
  {                       // since ipredmodes could be overwritten => encoder-decoder-mismatches
    // pick lowest cost prediction mode
    min_cost = DISTBLK_MAX;
    for (i=0;i<cr_MB_y;i++)
    {
      p_Vid->getNeighbour(currMB, 0 , i, p_Vid->mb_size[IS_CHROMA], &left[i]);
    }

    if ( p_Vid->mb_aff_frame_flag && p_Vid->field_mode )
    {
      for (i=0;i<cr_MB_y;i++)
      {
        left[i].pos_y = left[i].pos_y >> 1;
      }
    }

    for (mode=DC_PRED_8; mode<=PLANE_8; mode++)
    {
      if (((currSlice->slice_type != I_SLICE && currSlice->slice_type != SI_SLICE) || !p_Inp->IntraDisableInterOnly) && p_Inp->ChromaIntraDisable == 1 && mode!=DC_PRED_8)
        continue;

      if ((mode==VERT_PRED_8 && !mb_available_up) ||
        (mode==HOR_PRED_8 && (!mb_available_left[0] || !mb_available_left[1])) ||
        (mode==PLANE_8 && (!mb_available_left[0] || !mb_available_left[1] || !mb_available_up || !mb_available_up_left)))
        continue;

      cost = 0;
      for (uv = 1; uv < 3; uv++)
      {
        image = p_Vid->pImgOrg[uv];
        curr_mpr_16x16 = currSlice->mpr_16x16[uv];
        for (block_y=0; block_y<cr_MB_y; block_y+=4)
          for (block_x = 0; block_x < cr_MB_x; block_x += 4)
          {
            for (k=0, j = block_y; j < block_y + 4; j++)
            {
              img_prd = curr_mpr_16x16[mode][j];
              img_org = &image[left[j].pos_y][left[j].pos_x];
              for (i = block_x; i < block_x + 4; i++)
                diff[k++] = img_org[i] - img_prd[i];
            }
            cost += p_Vid->distortion4x4(diff, min_cost);
          }
      }
      if (cost < min_cost)
      {
        best_mode = mode;
        min_cost = cost;
      }
    }
    currMB->c_ipred_mode = (char) best_mode;
  }
}


