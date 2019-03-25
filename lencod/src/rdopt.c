
/*!
 ***************************************************************************
 * \file rdopt.c
 *
 * \brief
 *    Rate-Distortion optimized mode decision
 *
 * \author
 *    - Heiko Schwarz              <hschwarz@hhi.de>
 *    - Valeri George              <george@hhi.de>
 *    - Lowell Winger              <lwinger@lsil.com>
 *    - Alexis Michael Tourapis    <alexismt@ieee.org>
 * \date
 *    12. April 2001
 **************************************************************************
 */

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <memory.h>
#include <string.h>

#include "global.h"

#include "rdopt_coding_state.h"
#include "memalloc.h"
#include "mb_access.h"
#include "elements.h"
#include "intrarefresh.h"
#include "image.h"
#include "transform8x8.h"
#include "cabac.h"
#include "vlc.h"
#include "me_umhex.h"
#include "ratectl.h"            // head file for rate control
#include "mode_decision.h"
#include "fmo.h"
#include "macroblock.h"
#include "symbol.h"


imgpel pred[16][16];

#define FASTMODE 1
//#define RESET_STATE

extern const int LEVELMVLIMIT[17][6];
extern int   QP2QUANT[40];

const int AdaptRndCrPos[2][5] =
{
  //  P,   B,   I,  SP,  SI
  {   4,   7,   1,   4,   1}, // Intra MB
  {  10,  13,  10,  10,  10}  // Inter MB
};

const int AdaptRndPos[4][5] =
{
  //  P,   B,   I,  SP,  SI
  {   3,   6,   0,   3,   0}, // 4x4 Intra MB
  {   1,   2,   0,   1,   2}, // 8x8 Intra MB
  {   9,  12,   9,   9,   9}, // 4x4 Inter MB
  {   3,   4,   3,   3,   3}, // 8x8 Inter MB
};

imgpel   rec_mbY[16][16], rec_mbU[16][16], rec_mbV[16][16];    // reconstruction values

int lrec_rec[16][16],lrec_rec_U[16][16],lrec_rec_V[16][16]; // store the transf. and quantized coefficients for SP frames

static int diff[16];
static int diff4x4[64];
static int diff8x8[64];
RD_8x8DATA tr4x4, tr8x8;

int   **bestInterFAdjust4x4=NULL, **bestIntraFAdjust4x4=NULL;
int   **bestInterFAdjust8x8=NULL, **bestIntraFAdjust8x8=NULL;
int   ***bestInterFAdjust4x4Cr=NULL, ***bestIntraFAdjust4x4Cr=NULL;
int   ***bestInterFAdjust8x8Cr=NULL, ***bestIntraFAdjust8x8Cr=NULL;
int   **fadjust8x8=NULL, **fadjust4x4=NULL, ***fadjust4x4Cr=NULL, ***fadjust8x8Cr=NULL;


int   ****cofAC=NULL, ****cofAC8x8=NULL;        // [8x8block][4x4block][level/run][scan_pos]
int   ***cofDC=NULL;                       // [yuv][level/run][scan_pos]
int   **cofAC4x4=NULL, ****cofAC4x4intern=NULL; // [level/run][scan_pos]
int   cbp, cbp8x8, cnt_nonz_8x8;
int   cbp_blk8x8;
char  l0_refframe[4][4], l1_refframe[4][4];
int   b8mode[4], b8pdir[4];
short best8x8mode [4];                // [block]
char  best8x8pdir  [MAXMODE][4];       // [mode][block]
char  best8x8l0ref [MAXMODE][4];       // [mode][block]
char  best8x8l1ref [MAXMODE][4];       // [mode][block]


CSptr cs_mb=NULL, cs_b8=NULL, cs_cm=NULL, cs_ib8=NULL, cs_ib4=NULL;
int   best_c_imode;
int   best_i16offset;
short best_mode;
short  bi_pred_me;

//mixed transform sizes definitions
int   luma_transform_size_8x8_flag;

short all_mv8x8[2][2][4][4][2];       //[8x8_data/temp_data][LIST][block_x][block_y][MVx/MVy]
short pred_mv8x8[2][2][4][4][2];

int   ****cofAC8x8ts[3] = {NULL, NULL, NULL};        // [plane][8x8block][4x4block][level/run][scan_pos]
int   ****cofAC8x8CbCr[2];
int   **cofAC4x4CbCr[2];
int   ****cofAC4x4CbCrintern[2];

int64    cbp_blk8_8x8ts;
int      cbp8_8x8ts;
int      cost8_8x8ts;
int      cnt_nonz8_8x8ts;

// adaptive langrangian parameters
double mb16x16_cost;
double lambda_mf_factor;

void StoreMV8x8(int dir);
void RestoreMV8x8(int dir);
// end of mixed transform sizes definitions

//Adaptive Rounding update function
void update_offset_params(Macroblock *currMB, int mode, int luma_transform_size_8x8_flag);

char  b4_ipredmode[16], b4_intra_pred_modes[16];

/*!
 ************************************************************************
 * \brief
 *    delete structure for RD-optimized mode decision
 ************************************************************************
 */
void clear_rdopt ()
{
  free_mem_DCcoeff (cofDC);
  free_mem_ACcoeff (cofAC);
  free_mem_ACcoeff (cofAC8x8);
  free_mem_ACcoeff (cofAC4x4intern);

  if (input->Transform8x8Mode)
  {
    free_mem_ACcoeff (cofAC8x8ts[0]);
    //if (img->yuv_format==YUV444 && !IS_INDEPENDENT(input))
    {
      free_mem_ACcoeff (cofAC8x8ts[1]);
      free_mem_ACcoeff (cofAC8x8ts[2]);
    }
  }
  //if (img->yuv_format==YUV444 && !IS_INDEPENDENT(input))
  {
    free_mem_ACcoeff (cofAC8x8CbCr[0]);
    free_mem_ACcoeff (cofAC8x8CbCr[1]);
    free_mem_ACcoeff (cofAC4x4CbCrintern[0]);
    free_mem_ACcoeff (cofAC4x4CbCrintern[1]);
    
  }

  if (input->AdaptiveRounding)
  {
    free_mem2Dint(bestInterFAdjust4x4);
    free_mem2Dint(bestIntraFAdjust4x4);
    free_mem2Dint(bestInterFAdjust8x8);
    free_mem2Dint(bestIntraFAdjust8x8);
    free_mem2Dint(fadjust8x8);
    free_mem2Dint(fadjust4x4);
    if (input->yuv_format != 0)
    {
      if (input->yuv_format == YUV444)
      {
        free_mem3Dint(bestInterFAdjust8x8Cr, 2);
        free_mem3Dint(bestIntraFAdjust8x8Cr, 2);
      }
      free_mem3Dint(bestInterFAdjust4x4Cr, 2);
      free_mem3Dint(bestIntraFAdjust4x4Cr, 2);
      free_mem3Dint(fadjust4x4Cr, 2);
      free_mem3Dint(fadjust8x8Cr, 2);
    }
  }

  // structure for saving the coding state
  delete_coding_state (cs_mb);
  delete_coding_state (cs_b8);
  delete_coding_state (cs_cm);
  delete_coding_state (cs_ib8);
  delete_coding_state (cs_ib4);
}


/*!
 ************************************************************************
 * \brief
 *    create structure for RD-optimized mode decision
 ************************************************************************
 */
void init_rdopt ()
{
  rdopt = NULL;

  get_mem_DCcoeff (&cofDC);
  get_mem_ACcoeff (&cofAC);
  get_mem_ACcoeff (&cofAC8x8);
  get_mem_ACcoeff (&cofAC4x4intern);
  cofAC4x4 = cofAC4x4intern[0][0];

  if (input->Transform8x8Mode)
  {
    get_mem_ACcoeff (&cofAC8x8ts[0]);
    //if (img->yuv_format==YUV444 && !IS_INDEPENDENT(input))
    {
      get_mem_ACcoeff (&cofAC8x8ts[1]);
      get_mem_ACcoeff (&cofAC8x8ts[2]);
    }
  }

  //if (img->yuv_format==YUV444 && !IS_INDEPENDENT(input))
  {
    get_mem_ACcoeff (&cofAC8x8CbCr[0]);
    get_mem_ACcoeff (&cofAC8x8CbCr[1]);
    
    get_mem_ACcoeff (&cofAC4x4CbCrintern[0]);
    get_mem_ACcoeff (&cofAC4x4CbCrintern[1]);
    cofAC4x4CbCr[0] = cofAC4x4CbCrintern[0][0][0];
    cofAC4x4CbCr[1] = cofAC4x4CbCrintern[1][0][0];    
  }
  
  switch (input->rdopt)
  {
  case 0:
    encode_one_macroblock = encode_one_macroblock_low;
    break;
  case 1:
  default:
    encode_one_macroblock = encode_one_macroblock_high;
    break;
  case 2:
    encode_one_macroblock = encode_one_macroblock_highfast;
    break;
  case 3:
    encode_one_macroblock = encode_one_macroblock_highloss;
    break;
  }
  if (input->AdaptiveRounding)
  {
    get_mem2Dint(&bestInterFAdjust4x4, 16, 16);
    get_mem2Dint(&bestIntraFAdjust4x4, 16, 16);
    get_mem2Dint(&bestInterFAdjust8x8, 16, 16);
    get_mem2Dint(&bestIntraFAdjust8x8, 16, 16);
    get_mem2Dint(&fadjust8x8, 16, 16);
    get_mem2Dint(&fadjust4x4, 16, 16);
    if (input->yuv_format != 0 )
    {
      if  (input->yuv_format == YUV444)
      {
        get_mem3Dint(&bestInterFAdjust8x8Cr, 2, 16, 16);
        get_mem3Dint(&bestIntraFAdjust8x8Cr, 2, 16, 16);
      }
      get_mem3Dint(&bestInterFAdjust4x4Cr, 2, img->mb_cr_size_y, img->mb_cr_size_x);
      get_mem3Dint(&bestIntraFAdjust4x4Cr, 2, img->mb_cr_size_y, img->mb_cr_size_x);
      get_mem3Dint(&fadjust4x4Cr, 2, img->mb_cr_size_y, img->mb_cr_size_x);
      get_mem3Dint(&fadjust8x8Cr, 2, img->mb_cr_size_y, img->mb_cr_size_x);
    }
  }

  // structure for saving the coding state
  cs_mb  = create_coding_state ();
  cs_b8  = create_coding_state ();
  cs_cm  = create_coding_state ();
  cs_ib8 = create_coding_state ();
  cs_ib4 = create_coding_state ();

  if (input->CtxAdptLagrangeMult == 1)
  {
    mb16x16_cost = CALM_MF_FACTOR_THRESHOLD;
    lambda_mf_factor = 1.0;
  }

  getDistortion = distortionSSE;
}



/*!
 *************************************************************************************
 * \brief
 *    Updates the pixel map that shows, which reference frames are reliable for
 *    each MB-area of the picture.
 *
 * \note
 *    The new values of the pixel_map are taken from the temporary buffer refresh_map
 *
 *************************************************************************************
 */
void UpdatePixelMap()
{
  int mx,my,y,x,i,j;
  if (img->type==I_SLICE)
  {
    for (y=0; y<img->height; y++)
      for (x=0; x<img->width; x++)
      {
        pixel_map[y][x]=1;
      }
  }
  else
  {
    for (my=0; my<img->height >> 3; my++)
      for (mx=0; mx<img->width >> 3;  mx++)
      {
        j = my*8 + 8;
        i = mx*8 + 8;
        if (refresh_map[my][mx])
        {
          for (y=my*8; y<j; y++)
            for (x=mx*8; x<i; x++)
              pixel_map[y][x] = 1;
        }
        else
        {
          for (y=my*8; y<j; y++)
            for (x=mx*8; x<i; x++)
            {
              pixel_map[y][x] = imin(pixel_map[y][x] + 1, input->num_ref_frames+1);
            }
        }
     }
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Checks if a given reference frame is reliable for the current
 *    macroblock, given the motion vectors that the motion search has
 *    returned.
 *
 * \return
 *    If the return value is 1, the reference frame is reliable. If it
 *    is 0, then it is not reliable.
 *
 * \note
 *    A specific area in each reference frame is assumed to be unreliable
 *    if the same area has been intra-refreshed in a subsequent frame.
 *    The information about intra-refreshed areas is kept in the pixel_map.
 *
 *************************************************************************************
 */
int CheckReliabilityOfRef (int block, int list_idx, int ref, int mode)
{
  int y,x, block_y, block_x, dy, dx, y_pos, x_pos, yy, xx, pres_x, pres_y;
  int maxold_x  = img->width-1;
  int maxold_y  = img->height-1;
  int ref_frame = ref + 1;

  int by0 = (mode>=4?2*(block >> 1):mode==2?2*block:0);
  int by1 = by0 + (mode>=4||mode==2?2:4);
  int bx0 = (mode>=4?2*(block & 0x01):mode==3?2*block:0);
  int bx1 = bx0 + (mode>=4||mode==3?2:4);

  for (block_y=by0; block_y<by1; block_y++)
  {
    for (block_x=bx0; block_x<bx1; block_x++)
    {
      y_pos  = img->all_mv[block_y][block_x][list_idx][ref][mode][1];
      y_pos += (img->block_y + block_y) * BLOCK_SIZE * 4;
      x_pos  = img->all_mv[block_y][block_x][list_idx][ref][mode][0];
      x_pos += (img->block_x + block_x) * BLOCK_SIZE * 4;

      /* Here we specify which pixels of the reference frame influence
      the reference values and check their reliability. This is
      based on the function Get_Reference_Pixel */

      dy = y_pos & 3;
      dx = x_pos & 3;

      y_pos = (y_pos - dy) >> 2;
      x_pos = (x_pos - dx) >> 2;

      if (dy==0 && dx==0) //full-pel
      {
        for (y=y_pos ; y < y_pos + BLOCK_SIZE ; y++)
          for (x=x_pos ; x < x_pos + BLOCK_SIZE ; x++)
            if (pixel_map[iClip3(0,maxold_y,y)][iClip3(0,maxold_x,x)] < ref_frame)
              return 0;
      }
      else  /* other positions */
      {
        if (dy == 0)
        {
          for (y = y_pos ; y < y_pos + BLOCK_SIZE ; y++)
          {
            pres_y = iClip3(0, maxold_y, y);
            for (x = x_pos ; x < x_pos + BLOCK_SIZE ; x++)
            {
              for(xx = -2 ; xx < 4 ; xx++) {
                pres_x = iClip3(0, maxold_x, x + xx);
                if (pixel_map[pres_y][pres_x] < ref_frame)
                  return 0;
              }
            }
          }
        }
        else if (dx == 0)
        {
          for (y = y_pos ; y < y_pos + BLOCK_SIZE ; y++)
            for (x=x_pos ; x < x_pos + BLOCK_SIZE ; x++)
            {
              pres_x = iClip3(0,maxold_x,x);
              for(yy=-2;yy<4;yy++) {
                pres_y = iClip3(0,maxold_y, yy + y);
                if (pixel_map[pres_y][pres_x] < ref_frame)
                  return 0;
              }
            }
        }
        else if (dx == 2)
        {
          for (y = y_pos ; y < y_pos + BLOCK_SIZE ; y++)
            for (x = x_pos ; x < x_pos + BLOCK_SIZE ; x++)
            {
              for(yy=-2;yy<4;yy++) {
                pres_y = iClip3(0,maxold_y, yy + y);
                for(xx=-2;xx<4;xx++) {
                  pres_x = iClip3(0,maxold_x, xx + x);
                  if (pixel_map[pres_y][pres_x] < ref_frame)
                    return 0;
                }
              }
            }
        }
        else if (dy == 2)
        {
          for (y = y_pos ; y < y_pos + BLOCK_SIZE ; y++)
            for (x = x_pos ; x < x_pos + BLOCK_SIZE ; x++)
            {
              for(xx=-2;xx<4;xx++) {
                pres_x = iClip3(0,maxold_x, xx + x);
                for(yy=-2;yy<4;yy++) {
                  pres_y = iClip3(0,maxold_y, yy + y);
                  if (pixel_map[pres_y][pres_x] < ref_frame)
                    return 0;
                }
              }
            }
        }
        else
        {
          for (y = y_pos ; y < y_pos + BLOCK_SIZE ; y++)
          {
            for (x = x_pos ; x < x_pos + BLOCK_SIZE ; x++)
            {
              pres_y = dy == 1 ? y : y + 1;
              pres_y = iClip3(0,maxold_y,pres_y);

              for(xx=-2;xx<4;xx++)
              {
                pres_x = iClip3(0,maxold_x,xx + x);
                if (pixel_map[pres_y][pres_x] < ref_frame)
                  return 0;
              }

              pres_x = dx == 1 ? x : x + 1;
              pres_x = iClip3(0,maxold_x,pres_x);

              for(yy=-2;yy<4;yy++)
              {
                pres_y = iClip3(0,maxold_y, yy + y);
                if (pixel_map[pres_y][pres_x] < ref_frame)
                  return 0;
              }
            }
          }
        }
      }
    }
  }
  return 1;
}



/*!
 *************************************************************************************
 * \brief
 *    R-D Cost for an 4x4 Intra block
 *************************************************************************************
 */
double RDCost_for_4x4IntraBlocks (Macroblock    *currMB,
                                  int*    nonzero,
                                  int     b8,
                                  int     b4,
                                  int     ipmode,
                                  double  lambda,
                                  int mostProbableMode,
                                  int c_nzCbCr[3])
{
  double  rdcost;
  int     dummy = 0, x, y, rate;
  int64   distortion  = 0;
  int     block_x     = ((b8 & 0x01) << 3) + ((b4 & 0x01) << 2);
  int     block_y     = ((b8 >> 1) << 3) + ((b4 >> 1) << 2);
  int     pic_pix_x   = img->pix_x  + block_x;
  int     pic_pix_y   = img->pix_y  + block_y;
  int     pic_opix_y  = img->opix_y + block_y;
  imgpel  *img_org, *img_enc;

  Slice          *currSlice = img->currentSlice;
  SyntaxElement  se;
  const int      *partMap   = assignSE2partition[input->partition_mode];
  DataPartition  *dataPart;

  //===== perform DCT, Q, IQ, IDCT, Reconstruction =====
  select_dct(currMB);
  *nonzero = pDCT_4x4 (currMB, PLANE_Y, block_x, block_y, &dummy, 1);

  //===== get distortion (SSD) of 4x4 block =====
  for (y=0; y<4; y++)
  {
    img_org = pCurImg[pic_opix_y+y];
    img_enc = enc_picture->imgY[pic_pix_y+y];
    for (x=pic_pix_x; x<pic_pix_x+4; x++)
    {
      distortion += iabs2( img_org[x] - img_enc[x] );
    }
  }
  if((active_sps->chroma_format_idc == YUV444) && (IS_INDEPENDENT(input)==0))
  {
    ColorPlane k;
    for (k = PLANE_U; k <= PLANE_V; k++)
    {
      select_plane(k);
      c_nzCbCr[k] = pDCT_4x4(currMB, k, block_x, block_y, &dummy, 1);

      for( y = 0; y < 4; y++ )
      {
        for (x=pic_pix_x; x<pic_pix_x+4; x++)
        {
          distortion +=iabs2( pCurImg[pic_opix_y+y][x] - enc_picture->p_curr_img[pic_pix_y+y][x] );
        }
      }
    }
    ipmode_DPCM=NO_INTRA_PMODE;  //For residual DPCM
    select_plane(PLANE_Y);
  }
  else if((active_sps->chroma_format_idc == YUV444) && IS_INDEPENDENT(input))   //For residual DPCM
  {
    ipmode_DPCM=NO_INTRA_PMODE;  
  }


  //===== RATE for INTRA PREDICTION MODE  (SYMBOL MODE MUST BE SET TO CAVLC) =====
  se.value1 = (mostProbableMode == ipmode) ? -1 : ipmode < mostProbableMode ? ipmode : ipmode - 1;

  //--- set position and type ---
  se.context = (b8 << 2) + b4;
  se.type    = SE_INTRAPREDMODE;

  //--- choose data partition ---
  dataPart = &(currSlice->partArr[partMap[SE_INTRAPREDMODE]]);
  //--- encode and update rate ---
  writeIntraPredMode (&se, dataPart);
  rate = se.len;

  //===== RATE for LUMINANCE COEFFICIENTS =====
  if (input->symbol_mode == CAVLC)
  {
    rate  += writeCoeff4x4_CAVLC (currMB, LUMA, b8, b4, 0);
    if((active_sps->chroma_format_idc == YUV444) && (IS_INDEPENDENT(input)==0) )  
    {
      rate  += writeCoeff4x4_CAVLC (currMB, CB, b8, b4, 0);
      rate  += writeCoeff4x4_CAVLC (currMB, CR, b8, b4, 0);
    }
  }
  else
  {
    rate  += writeCoeff4x4_CABAC (currMB, PLANE_Y, b8, b4, 1);
    if((active_sps->chroma_format_idc == YUV444) && (IS_INDEPENDENT(input)==0) )  
    {
      rate  += writeCoeff4x4_CABAC (currMB, PLANE_U, b8, b4, 1);
      rate  += writeCoeff4x4_CABAC (currMB, PLANE_V, b8, b4, 1);
    }
  }
  //reset_coding_state (currMB, cs_cm);
  rdcost = (double)distortion + lambda * (double) rate;

  return rdcost;
}


/*!
 *************************************************************************************
 * \brief
 *    Mode Decision for an 4x4 Intra block
 *************************************************************************************
 */
int Mode_Decision_for_4x4IntraBlocks (Macroblock *currMB, int  b8,  int  b4,  double  lambda,  double*  min_cost, int cr_cbp[3])
{
  int     ipmode, best_ipmode = 0, i, j, y, cost, dummy;
  int     c_nz, nonzero = 0;
  int*    ACLevel = img->cofAC[b8][b4][0];
  int*    ACRun   = img->cofAC[b8][b4][1];
  static imgpel *org_img, *prd_img;
  static int *m7;
  int    c_nzCbCr[3]= {999,999, 999};
  static imgpel  rec4x4[4][4];
  static imgpel  rec4x4CbCr[2][4][4];
  int    uv;
  
  double rdcost;
  int    block_x     = ((b8 & 0x01) << 3) + ((b4 & 0x01) << 2);
  int    block_y     = ((b8 >> 1) << 3)  + ((b4 >> 1) << 2);
  int    pic_pix_x   = img->pix_x  + block_x;
  int    pic_pix_y   = img->pix_y  + block_y;
  int    pic_opix_x  = img->opix_x + block_x;
  int    pic_opix_y  = img->opix_y + block_y;
  int    pic_block_x = pic_pix_x >> 2;
  int    pic_block_y = pic_pix_y >> 2;
  double min_rdcost  = 1e30;
  int *d;

  int left_available, up_available, all_available;
  imgpel (*curr_mpr)[16] = img->mpr[0];

  char   upMode;
  char   leftMode;
  int    mostProbableMode;

  PixelPos left_block;
  PixelPos top_block;

  int  lrec4x4[4][4];
  int  fixedcost = (int) floor(4 * lambda );

  //For residual DPCM
  Boolean lossless_qpprime = ((currMB->qp + img->bitdepth_luma_qp_scale)==0 && img->lossless_qpprime_flag==1);  

#ifdef BEST_NZ_COEFF
  int best_nz_coeff = 0;
  int best_coded_block_flag = 0;
  int bit_pos = 1 + ((((b8>>1)<<1)+(b4>>1))<<2) + (((b8&1)<<1)+(b4&1));
  static int64 cbp_bits;

  if (b8==0 && b4==0)
   cbp_bits = 0;
#endif

  getLuma4x4Neighbour(currMB, block_x - 1, block_y,     &left_block);
  getLuma4x4Neighbour(currMB, block_x,     block_y - 1, &top_block);

  // constrained intra pred
  if (input->UseConstrainedIntraPred)
  {
    left_block.available = left_block.available ? img->intra_block[left_block.mb_addr] : 0;
    top_block.available  = top_block.available  ? img->intra_block[top_block.mb_addr]  : 0;
  }

  upMode            =  top_block.available ? img->ipredmode[top_block.pos_y ][top_block.pos_x ] : -1;
  leftMode          = left_block.available ? img->ipredmode[left_block.pos_y][left_block.pos_x] : -1;

  mostProbableMode  = (upMode < 0 || leftMode < 0) ? DC_PRED : upMode < leftMode ? upMode : leftMode;

  *min_cost = INT_MAX;

  ipmode_DPCM = NO_INTRA_PMODE; ////For residual DPCM

  //===== INTRA PREDICTION FOR 4x4 BLOCK =====
  intrapred_4x4 (currMB, PLANE_Y, pic_pix_x, pic_pix_y, &left_available, &up_available, &all_available);

  if ((img->yuv_format == YUV444) && !IS_INDEPENDENT(input))
  {
    select_plane(PLANE_U);
    intrapred_4x4 (currMB, PLANE_U, pic_pix_x, pic_pix_y, &left_available, &up_available, &all_available);
    select_plane(PLANE_V);
    intrapred_4x4 (currMB, PLANE_V, pic_pix_x, pic_pix_y, &left_available, &up_available, &all_available);
    select_plane(PLANE_Y);
  }

  //===== LOOP OVER ALL 4x4 INTRA PREDICTION MODES =====
  for (ipmode = 0; ipmode < NO_INTRA_PMODE; ipmode++)
  {
    int available_mode =  (all_available) || (ipmode==DC_PRED) ||
      (up_available && (ipmode==VERT_PRED||ipmode==VERT_LEFT_PRED||ipmode==DIAG_DOWN_LEFT_PRED)) ||
      (left_available && (ipmode==HOR_PRED||ipmode==HOR_UP_PRED));

    if (input->IntraDisableInterOnly==0 || img->type != I_SLICE)
    {
      if (input->Intra4x4ParDisable && (ipmode==VERT_PRED||ipmode==HOR_PRED))
        continue;

      if (input->Intra4x4DiagDisable && (ipmode==DIAG_DOWN_LEFT_PRED||ipmode==DIAG_DOWN_RIGHT_PRED))
        continue;

      if (input->Intra4x4DirDisable && ipmode>=VERT_RIGHT_PRED)
        continue;
    }

    if( available_mode)
    {
      if (!input->rdopt)
      {
        d = &diff[0];
        for (j = 0; j<4; j++)
        {
          org_img = &pCurImg[pic_opix_y+j][pic_opix_x];
          prd_img = img->mpr_4x4[0][ipmode][j];
          
          for (i=0; i<4; i++)
          {
            *d++ = *org_img++ - *prd_img++;
          }
        }
        cost  = (ipmode == mostProbableMode) ? 0 : fixedcost;
        cost += distortion4x4 (diff);
        if ((img->yuv_format == YUV444) && !IS_INDEPENDENT(input))
        {
          int k;
          for (k=j=0; j<4; j++)
          {
            for (i=0; i<4; i++, k++)
            {
              diff[k] = imgUV_org[0][pic_opix_y+j][pic_opix_x+i] - img->mpr_4x4[1][ipmode][j][i]; 
            }
          }
          
          cost += distortion4x4 (diff);
          for (k=j=0; j<4; j++)
          {
            for (i=0; i<4; i++, k++)
            {
              diff[k] = imgUV_org[1][pic_opix_y+j][pic_opix_x+i] - img->mpr_4x4[2][ipmode][j][i]; 
            }
          }
          cost += distortion4x4 (diff);
        }
        if (cost < *min_cost)
        {
          best_ipmode = ipmode;
          *min_cost   = cost;
        }
      }
      else
      {
        // get prediction and prediction error
        for (j=0; j<4; j++)
        {
          org_img = &pCurImg[pic_opix_y+j][pic_opix_x];
          prd_img = img->mpr_4x4[0][ipmode][j];
          m7 = &img->m7[0][block_y+j][block_x];
          
          memcpy(&curr_mpr[block_y+j][block_x], prd_img, BLOCK_SIZE * sizeof(imgpel));
          for (i=0; i<4; i++)
          {
            *m7++ = (int) (*org_img++ - *prd_img++);
          }
        }
        if ((img->yuv_format == YUV444) && !IS_INDEPENDENT(input))
        {
          for (uv = 1; uv < 3; uv++)
          {
            for (j=0; j<4; j++)
            {
              org_img = &imgUV_org[uv - 1][pic_opix_y+j][pic_opix_x];
              prd_img = img->mpr_4x4[uv][ipmode][j];
              m7 = &img->m7[uv][block_y+j][block_x];

              memcpy(&curr_mpr[block_y+j][block_x], prd_img, BLOCK_SIZE * sizeof(imgpel));
              for (i=0; i<4; i++)
              {
                *m7++ = (int) (*org_img++ - *prd_img++);
              }
            }
          }
        }

          if ((img->yuv_format == YUV444) && !IS_INDEPENDENT(input))  //For residual DPCM
          {
            if((lossless_qpprime)&&(ipmode<2))  
            {
              Residual_DPCM_4x4(ipmode, 0, block_y, block_x);
              Residual_DPCM_4x4(ipmode, 1, block_y, block_x);
              Residual_DPCM_4x4(ipmode, 2, block_y, block_x);
              ipmode_DPCM=ipmode;
            }
          }
          else if ((img->yuv_format == YUV444) && IS_INDEPENDENT(input))
          {
            if((lossless_qpprime)&&(ipmode<2))  
            {
              Residual_DPCM_4x4(ipmode, 0, block_y, block_x);
              ipmode_DPCM=ipmode;
            }
          }

        //===== store the coding state =====
        //store_coding_state (currMB, cs_cm);
        // get and check rate-distortion cost
#ifdef BEST_NZ_COEFF
        currMB->cbp_bits[0] = cbp_bits;
#endif
        if ((rdcost = RDCost_for_4x4IntraBlocks (currMB, &c_nz, b8, b4, ipmode, lambda, mostProbableMode, c_nzCbCr)) < min_rdcost)
        {
          //--- set coefficients ---
          memcpy(cofAC4x4[0],ACLevel, 18 * sizeof(int));
          memcpy(cofAC4x4[1],ACRun, 18 * sizeof(int));

          //--- set reconstruction ---
          for (y=0; y<4; y++)
          {
            memcpy(rec4x4[y],&enc_picture->imgY[pic_pix_y+y][pic_pix_x], BLOCK_SIZE * sizeof(imgpel));
          }

          // SP/SI reconstruction
          if(img->type==SP_SLICE &&(!si_frame_indicator && !sp2_frame_indicator))
          {
            for (y=0; y<4; y++)
            {
              memcpy(lrec4x4[y],&lrec[pic_pix_y+y][pic_pix_x], BLOCK_SIZE * sizeof(int));// stores the mode coefficients
            }
          }
            if( img->yuv_format==YUV444 && !IS_INDEPENDENT(input) ) 
            { 
              //--- set coefficients ---
              for (uv=0; uv < 2; uv++)
              {
                memcpy(cofAC4x4CbCr[uv][0],img->cofAC[b8+4+uv*4][b4][0], 18 * sizeof(int));
                memcpy(cofAC4x4CbCr[uv][1],img->cofAC[b8+4+uv*4][b4][1], 18 * sizeof(int));
                cr_cbp[uv + 1] = c_nzCbCr[uv + 1];

                //--- set reconstruction ---
                for (y=0; y<4; y++)
                {
                  memcpy(rec4x4CbCr[uv][y],&enc_picture->imgUV[uv][pic_pix_y+y][pic_pix_x], BLOCK_SIZE * sizeof(imgpel));
                } 
              }
            }
          //--- flag if dct-coefficients must be coded ---
          nonzero = c_nz;

          //--- set best mode update minimum cost ---
          *min_cost     = rdcost;
          min_rdcost    = rdcost;
          best_ipmode   = ipmode;
#ifdef BEST_NZ_COEFF
          best_nz_coeff = img->nz_coeff [img->current_mb_nr][block_x4][block_y4];
          best_coded_block_flag = (int)((currMB->cbp_bits[0] >> bit_pos)&(int64)(1));
#endif
          //store_coding_state (currMB, cs_ib4);
          if (img->AdaptiveRounding)
          {
            for (j = block_y; j < block_y + 4; j++)
              memcpy(&fadjust4x4[j][block_x],&img->fadjust4x4[1][j][block_x], BLOCK_SIZE * sizeof(int));
              if(img->yuv_format == YUV444 && !IS_INDEPENDENT(input))
              {
                for (j=0; j<4; j++)
                {
                  memcpy(&fadjust4x4Cr[0][block_y+j][block_x],&img->fadjust4x4Cr[0][1][block_y+j][block_x], BLOCK_SIZE * sizeof(int));
                  memcpy(&fadjust4x4Cr[1][block_y+j][block_x],&img->fadjust4x4Cr[1][1][block_y+j][block_x], BLOCK_SIZE * sizeof(int));
                }
              }
          }
        }

#ifndef RESET_STATE
        reset_coding_state (currMB, cs_cm);
#endif
      }
    }
  }

#ifdef BEST_NZ_COEFF
  img->nz_coeff [img->current_mb_nr][block_x4][block_y4] = best_nz_coeff;
  cbp_bits &= (~(int64)(1<<bit_pos));
  cbp_bits |= (int64)(best_coded_block_flag<<bit_pos);
#endif
  //===== set intra mode prediction =====
  img->ipredmode[pic_block_y][pic_block_x] = (char) best_ipmode;
  currMB->intra_pred_modes[4*b8+b4] =
    (char) (mostProbableMode == best_ipmode ? -1 : (best_ipmode < mostProbableMode ? best_ipmode : best_ipmode-1));

  if (!input->rdopt)
  {
    // get prediction and prediction error
    for (j=0; j < BLOCK_SIZE; j++)
    {
      m7 = &img->m7[0][block_y+j][block_x];
      org_img = &pCurImg[pic_opix_y+j][pic_opix_x];
      prd_img = img->mpr_4x4[0][best_ipmode][j];      
      memcpy(&curr_mpr[block_y+j][block_x], prd_img, BLOCK_SIZE * sizeof(imgpel));
      for (i=0; i<BLOCK_SIZE; i++)
      {
        m7[i] = org_img[i] - prd_img[i];
      }
    }

    if(lossless_qpprime)   //For residual DPCM
    {
      ipmode_DPCM=best_ipmode;  
      if((best_ipmode<2))
      {
        Residual_DPCM_4x4(best_ipmode, 0, block_y, block_x);
      }
    }

    select_dct(currMB);
    nonzero = pDCT_4x4 (currMB, PLANE_Y, block_x, block_y, &dummy, 1);

    if ((img->yuv_format == YUV444) && !IS_INDEPENDENT(input))
    {
      ColorPlane k;
      for (k = PLANE_U; k <= PLANE_V; k++)
      {
        select_plane(k);
        for (j=0; j<4; j++)
        {
          for (i=0; i<4; i++)
          {
            img->mpr[k][block_y+j][block_x+i]  = img->mpr_4x4[k][best_ipmode][j][i];    
            img->m7[k][block_y+j][block_x+i]   = pImgOrg[k][pic_opix_y+j][pic_opix_x+i] - img->mpr_4x4[k][best_ipmode][j][i]; 
          }
        }

        if(lossless_qpprime)   //For residual DPCM
        { 
          if((best_ipmode<2))
          {
            Residual_DPCM_4x4(best_ipmode, k, block_y, block_x);
          }
        }

        cr_cbp[k] = pDCT_4x4(currMB, k, block_x,block_y,&dummy,1);
      }
      select_plane(PLANE_Y);
    }
  }
  else
  {
    if( (img->yuv_format==YUV444) && !IS_INDEPENDENT(input) )
    {
      select_plane(PLANE_U);
      for (j=0; j<4; j++)
      {
        for (i=0; i<4; i++)
        {
          img->mpr[1][block_y+j][block_x+i]  = img->mpr_4x4[1][best_ipmode][j][i];
          img->m7[1][block_y+j][block_x+i]   = imgUV_org[0][img->pix_y+block_y+j][img->pix_x+block_x+i] - img->mpr_4x4[1][best_ipmode][j][i];
        }
      }
      cr_cbp[1] = pDCT_4x4(currMB, PLANE_U, block_x,block_y,&dummy,1);
      select_plane(PLANE_V);
      for (j=0; j<4; j++)
      {
        for (i=0; i<4; i++)
        {
          img->mpr[2][block_y+j][block_x+i]  = img->mpr_4x4[2][best_ipmode][j][i];
          img->m7[2][block_y+j][block_x+i]   = imgUV_org[1][img->pix_y+block_y+j][img->pix_x+block_x+i] - img->mpr_4x4[2][best_ipmode][j][i];
        }
      }
      cr_cbp[2]  = pDCT_4x4(currMB, PLANE_V,block_x,block_y,&dummy,1);
      select_plane(PLANE_Y);
    }
    //===== restore coefficients =====
    memcpy (ACLevel,cofAC4x4[0], 18 * sizeof(int));
    memcpy (ACRun,cofAC4x4[1], 18 * sizeof(int));

    //===== restore reconstruction and prediction (needed if single coeffs are removed) =====
    for (y=0; y<BLOCK_SIZE; y++)
    {
      memcpy (&enc_picture->imgY[pic_pix_y + y][pic_pix_x],rec4x4[y],    BLOCK_SIZE * sizeof(imgpel));
      memcpy (&curr_mpr[block_y + y][block_x],img->mpr_4x4[0][best_ipmode][y], BLOCK_SIZE * sizeof(imgpel));
    }
    
    // SP/SI reconstuction
    if(img->type==SP_SLICE &&(!si_frame_indicator && !sp2_frame_indicator))
    {
      for (y=0; y<BLOCK_SIZE; y++)
      {
        memcpy (&lrec[pic_pix_y+y][pic_pix_x],lrec4x4[y], BLOCK_SIZE * sizeof(int));//restore coefficients when encoding primary SP frame
      }
    }
    if( img->yuv_format == YUV444 && !IS_INDEPENDENT(input) ) 
    {
      for (uv=0; uv < 2; uv++ )
      {
        //===== restore coefficients =====
        memcpy(img->cofAC[b8+4+uv*4][b4][0], cofAC4x4CbCr[uv][0], 18 * sizeof(int));
        memcpy(img->cofAC[b8+4+uv*4][b4][1], cofAC4x4CbCr[uv][1], 18 * sizeof(int));
        //===== restore reconstruction and prediction (needed if single coeffs are removed) =====
        for (y=0; y<BLOCK_SIZE; y++)
        {
          memcpy(&enc_picture->imgUV[uv][pic_pix_y+y][pic_pix_x],rec4x4CbCr[uv][y], BLOCK_SIZE * sizeof(imgpel));
          memcpy(&img->mpr[uv + 1][block_y+y][block_x], img->mpr_4x4[uv + 1][best_ipmode][y], BLOCK_SIZE * sizeof(imgpel));
        } 
      }
    }

    if (img->AdaptiveRounding)
    {
      for (j=block_y; j<block_y + BLOCK_SIZE; j++)
        memcpy (&img->fadjust4x4[1][j][block_x],&fadjust4x4[j][block_x], BLOCK_SIZE * sizeof(int));
      if (img->yuv_format == YUV444 && !IS_INDEPENDENT(input))
      {
        for (j=0; j<BLOCK_SIZE; j++)
        {
          memcpy (&img->fadjust4x4Cr[0][1][block_y+j][block_x],&fadjust4x4Cr[0][block_y+j][block_x], BLOCK_SIZE * sizeof(int));
          memcpy (&img->fadjust4x4Cr[1][1][block_y+j][block_x],&fadjust4x4Cr[1][block_y+j][block_x], BLOCK_SIZE * sizeof(int));
        }
      }
    }
   
  }  

  return nonzero;
}


/*!
 *************************************************************************************
 * \brief
 *    Mode Decision for an 8x8 Intra block
 *************************************************************************************
 */
int Mode_Decision_for_8x8IntraBlocks(Macroblock *currMB, int b8,double lambda,double *cost, int non_zero[2])
{
  int  nonzero = 0, b4;
  double  cost4x4;
  int CbCr_cbp[3]={0, 0, 0};
  non_zero[0] = 0;
  non_zero[1] = 0;

  *cost = (int)floor(6.0 * lambda + 0.4999);

  for (b4=0; b4<4; b4++)
  {
    if (Mode_Decision_for_4x4IntraBlocks (currMB, b8, b4, lambda, &cost4x4, CbCr_cbp))
    {
      nonzero = 1;
    }
    *cost += cost4x4;

    non_zero[0] = (CbCr_cbp[1] != 0);
    non_zero[1] = (CbCr_cbp[2] != 0);
  }
#ifdef RESET_STATE
  //reset_coding_state (currMB, cs_cm);
#endif

  return nonzero;
}

/*!
 *************************************************************************************
 * \brief
 *    4x4 Intra mode decision for an macroblock
 *************************************************************************************
 */
int Mode_Decision_for_Intra4x4Macroblock (Macroblock *currMB, double lambda,  double* cost)
{
  int  cbp=0, b8;
  double cost8x8;
  int non_zero[2] = {0, 0};
  
  cmp_cbp[1] = cmp_cbp[2] = 0;

  for (*cost=0, b8=0; b8<4; b8++)
  {
    if (Mode_Decision_for_8x8IntraBlocks (currMB, b8, lambda, &cost8x8, non_zero))
    {
      cbp |= (1<<b8);
    }
    *cost += cost8x8;
    if (non_zero[0])
    {
      cmp_cbp[1] |= (1<<b8);
      cbp |= cmp_cbp[1];
      cmp_cbp[1] = cbp;
      cmp_cbp[2] = cbp;
    }
    if (non_zero[1])
    {
      cmp_cbp[2] |= (1<<b8);
      cbp |= cmp_cbp[2];
      cmp_cbp[1] = cbp;
      cmp_cbp[2] = cbp;
    }
  }
  return cbp;
}


/*!
 *************************************************************************************
 * \brief
 *    R-D Cost for an 8x8 Partition
 *************************************************************************************
 */
double RDCost_for_8x8blocks (Macroblock *currMB, // --> Current macroblock to code
                             int*    cnt_nonz,   // --> number of nonzero coefficients
                             int64*  cbp_blk,    // --> cbp blk
                             double  lambda,     // <-- lagrange multiplier
                             int     block,      // <-- 8x8 block number
                             int     mode,       // <-- partitioning mode
                             short   pdir,       // <-- prediction direction
                             short   l0_ref,     // <-- L0 reference picture
                             short   l1_ref)     // <-- L1 reference picture
{
  int  i, j, k;
  int  rate=0;
  int64 distortion=0;
  int  dummy = 0, mrate;
  int  fw_mode, bw_mode;
  int  cbp     = 0;
  int  pax     = 8*(block & 0x01);
  int  pay     = 8*(block >> 1);
  int  i0      = pax >> 2;
  int  j0      = pay >> 2;
  int  bframe  = (img->type==B_SLICE);
  int  direct  = (bframe && mode==0);
  int  b8value = B8Mode2Value (mode, pdir);

  SyntaxElement se;
  Slice         *currSlice = img->currentSlice;
  DataPartition *dataPart;
  const int     *partMap   = assignSE2partition[input->partition_mode];
  imgpel *img_enc, *img_org;

  EncodingEnvironmentPtr eep_dp;

  //=====
  //=====  GET COEFFICIENTS, RECONSTRUCTIONS, CBP
  //=====
  currMB->bi_pred_me=0;

  if (direct)
  {
    if (direct_pdir[img->block_y+j0][img->block_x+i0]<0) // mode not allowed
      return (1e20);
    else
      *cnt_nonz = LumaResidualCoding8x8 (currMB, &cbp, cbp_blk, block, direct_pdir[img->block_y+j0][img->block_x+i0], 0, 0,
      (short)imax(0,direct_ref_idx[LIST_0][img->block_y+j0][img->block_x+i0]),
      direct_ref_idx[LIST_1][img->block_y+j0][img->block_x+i0]);
  }
  else
  {
    if (pdir == 2 && active_pps->weighted_bipred_idc == 1)
    {
      int weight_sum = (active_pps->weighted_bipred_idc == 1)? wbp_weight[0][l0_ref][l1_ref][0] + wbp_weight[1][l0_ref][l1_ref][0] : 0;
      if (weight_sum < -128 ||  weight_sum > 127)
      {
        return (1e20);
      }
    }

    fw_mode   = (pdir==0||pdir==2 ? mode : 0);
    bw_mode   = (pdir==1||pdir==2 ? mode : 0);
    *cnt_nonz = LumaResidualCoding8x8 (currMB, &cbp, cbp_blk, block, pdir, fw_mode, bw_mode, l0_ref, l1_ref);
  }

  if( img->yuv_format==YUV444 && !IS_INDEPENDENT(input) ) 
  {
    *cnt_nonz += ( coeff_cost_cr[1] + coeff_cost_cr[2] );
  }
  //===== get residue =====
  if (input->rdopt==3 && img->type!=B_SLICE)
  {
    // We need the reconstructed prediction residue for the simulated decoders.
    compute_residue_b8block (block, -1);
  }

  //=====
  //=====   GET DISTORTION
  //=====
  if (input->rdopt==3 && img->type!=B_SLICE)
  {
    for (k=0; k<input->NoOfDecoders ;k++)
    {
      decode_one_b8block (k, P8x8, block, mode, l0_ref);
      for (j=img->opix_y+pay; j<img->opix_y+pay+8; j++)
      {
        img_org = pCurImg[j];
        img_enc =  decs->decY[k][j];
        for (i=img->opix_x+pax; i<img->opix_x+pax+8; i++)
        {
          distortion += iabs2( img_org[i] - img_enc[i]);
        }
      }
    }
    distortion /= input->NoOfDecoders;
  }
  else
  {
    for (j=pay; j<pay+8; j++)
    {
      img_org = pCurImg[img->opix_y + j];
      img_enc =  enc_picture->imgY[img->pix_y + j];
      for (i=img->pix_x+pax; i<img->pix_x+pax+8; i++)
      {
        distortion += iabs2( img_org[i] -img_enc[i]);
      }
    }
    if (img->yuv_format==YUV444 && !IS_INDEPENDENT(input))
    {
      for (j=pay; j<pay+8; j++)
      {
        for (i=img->pix_x+pax; i<img->pix_x+pax+8; i++)
        {
          distortion += iabs2( imgUV_org[0][img->opix_y+j][i] - enc_picture->imgUV[0][img->pix_y+j][i]);
          distortion += iabs2( imgUV_org[1][img->opix_y+j][i] - enc_picture->imgUV[1][img->pix_y+j][i]);
        }
      }
    }
    
  }
  
  if( img->yuv_format==YUV444 && !IS_INDEPENDENT(input) ) 
  {   
    cbp |= cmp_cbp[1];
    cbp |= cmp_cbp[2];
    
    cmp_cbp[1] = cbp;
    cmp_cbp[2] = cbp;
  }

  //=====
  //=====   GET RATE
  //=====
  //----- block 8x8 mode -----
  if (input->symbol_mode == CAVLC)
  {
    ue_linfo (b8value, dummy, &mrate, &dummy);
    rate += mrate;
  }
  else
  {
    se.value1 = b8value;
    se.type   = SE_MBTYPE;
    dataPart  = &(currSlice->partArr[partMap[se.type]]);
    writeB8_typeInfo(&se, dataPart);
    rate += se.len;
  }

  //----- motion information -----
  if (!direct)
  {
    if ((img->num_ref_idx_l0_active > 1 ) && (pdir==0 || pdir==2))
      rate  += writeReferenceFrame (currMB, mode, i0, j0, 1, l0_ref);
    
    if ((img->num_ref_idx_l1_active > 1 && img->type== B_SLICE ) && (pdir==1 || pdir==2))
    {
      rate  += writeReferenceFrame (currMB, mode, i0, j0, 0, l1_ref);
    }
    
    if (pdir==0 || pdir==2)
    {
      rate  += writeMotionVector8x8 (currMB, i0, j0, i0 + 2, j0 + 2, l0_ref, LIST_0, mode);
    }
    if (pdir==1 || pdir==2)
    {
      rate  += writeMotionVector8x8 (currMB, i0, j0, i0 + 2, j0 + 2, l1_ref, LIST_1, mode);
    }
  }

  //----- coded block pattern (for CABAC only) -----
  if (input->symbol_mode == CABAC)
  {
    dataPart = &(currSlice->partArr[partMap[SE_CBP]]);
    eep_dp   = &(dataPart->ee_cabac);
    mrate    = arienco_bits_written (eep_dp);
    writeCBP_BIT_CABAC (block, ((*cnt_nonz>0)?1:0), cbp8x8, currMB, 1, eep_dp);
    mrate    = arienco_bits_written (eep_dp) - mrate;
    rate    += mrate;
  }

  //----- luminance coefficients -----
  if (*cnt_nonz)
  {
    rate += writeCoeff8x8 (currMB, PLANE_Y, block, mode, currMB->luma_transform_size_8x8_flag);
  }
  if( active_sps->chroma_format_idc == YUV444 && !IS_INDEPENDENT(input) )
  {
    rate += writeCoeff8x8( currMB, PLANE_U, block, mode, currMB->luma_transform_size_8x8_flag );
    rate += writeCoeff8x8( currMB, PLANE_V, block, mode, currMB->luma_transform_size_8x8_flag );
  }

  return (double)distortion + lambda * (double)rate;
}


/*!
 *************************************************************************************
 * \brief
 *    Gets mode offset for intra16x16 mode
 *************************************************************************************
 */
int I16Offset (int cbp, int i16mode)
{
  return (cbp&15?13:1) + i16mode + ((cbp&0x30)>>2);
}


/*!
 *************************************************************************************
 * \brief
 *    Sets modes and reference frames for a macroblock
 *************************************************************************************
 */
void SetModesAndRefframeForBlocks (Macroblock *currMB, int mode)
{
  int i,j,k,l;
  int  bframe  = (img->type==B_SLICE);
  int  block_x, block_y;
  int  cur_ref;
  int  clist;
  char refl0, refl1;

  //--- macroblock type ---
  currMB->mb_type = mode;
  currMB->bi_pred_me= (mode == 1 ? img->bi_pred_me[mode] : 0);

  //--- block 8x8 mode and prediction direction ---
  switch (mode)
  {
  case 0:
    for(i=0;i<4;i++)
    {
      currMB->b8mode[i] = 0;
      currMB->b8pdir[i] = (bframe ? direct_pdir[img->block_y + ((i >> 1)<<1)][img->block_x + ((i & 0x01)<<1)] : 0);
    }
    break;
  case 1:
  case 2:
  case 3:
    for(i=0;i<4;i++)
    {
      currMB->b8mode[i] = mode;
      currMB->b8pdir[i] = best8x8pdir[mode][i];
    }
    break;
  case P8x8:
    for(i=0;i<4;i++)
    {
      currMB->b8mode[i]   = best8x8mode[i];
      currMB->b8pdir[i]   = best8x8pdir[mode][i];
    }
    break;
  case I4MB:
    for(i=0;i<4;i++)
    {
      currMB->b8mode[i] = IBLOCK;
      currMB->b8pdir[i] = -1;
    }
    break;
  case I16MB:
    for(i=0;i<4;i++)
    {
      currMB->b8mode[i] =  0;
      currMB->b8pdir[i] = -1;
    }
    break;
  case I8MB:
    for(i=0;i<4;i++)
    {
      currMB->b8mode[i] = I8MB;
      currMB->b8pdir[i] = -1;
    }
    //switch to 8x8 transform
    currMB->luma_transform_size_8x8_flag = 1;
    break;
  case IPCM:
    for(i=0;i<4;i++)
    {
      currMB->b8mode[i] = IPCM;
      currMB->b8pdir[i] = -1;
    }
    currMB->luma_transform_size_8x8_flag = 0;
    break;
  default:
    printf ("Unsupported mode in SetModesAndRefframeForBlocks!\n");
    exit (1);
  }

#define IS_FW ((best8x8pdir[mode][k]==0 || best8x8pdir[mode][k]==2) && (mode!=P8x8 || best8x8mode[k]!=0 || !bframe))
#define IS_BW ((best8x8pdir[mode][k]==1 || best8x8pdir[mode][k]==2) && (mode!=P8x8 || best8x8mode[k]!=0))
  //--- reference frame arrays ---
  if (mode==0 || mode==I4MB || mode==I16MB || mode==I8MB)
  {
    if (bframe)
    {
      if (!mode) // Direct
      {
        for (clist = LIST_0; clist <= LIST_1; clist++)
        {
          for (j = img->block_y; j < img->block_y + 4; j++)
          {
            memcpy(&enc_picture->ref_idx[clist][j][img->block_x], &direct_ref_idx[clist][j][img->block_x], 4 * sizeof(char));
          }
        }
      }
      else // Intra
      {
        for (clist = LIST_0; clist <= LIST_1; clist++)
        {
          for (j = img->block_y; j < img->block_y + 4; j++)
          {
            memset(&enc_picture->ref_idx[clist][j][img->block_x],-1, 4 * sizeof(char));
          }
        }
      }
    }
    else
    {
      if (!mode) // Skip
      {
        for (j = img->block_y; j < img->block_y + 4; j++)
          memset(&enc_picture->ref_idx[LIST_0][j][img->block_x],0, 4 * sizeof(char));
      }
      else // Intra
      {
        for (j = img->block_y; j < img->block_y + 4; j++)
          memset(&enc_picture->ref_idx[LIST_0][j][img->block_x],-1, 4 * sizeof(char));
      }
    }
  }
  else
  {
    if (bframe)
    {
      if (mode == 1)
      {
        if (currMB->bi_pred_me && best8x8pdir[mode][0] == 2)
        {
          for (j = img->block_y; j < img->block_y + 4;j++)
          {          
            memset(&enc_picture->ref_idx[LIST_0][j][img->block_x], 0, 4 * sizeof(char));
            memset(&enc_picture->ref_idx[LIST_1][j][img->block_x], 0, 4 * sizeof(char));
          }
        }
        else
        {
          refl0 = (best8x8pdir[mode][0] == 0 || best8x8pdir[mode][0] == 2) ? best8x8l0ref[mode][0] : -1;
          refl1 = (best8x8pdir[mode][0] == 1 || best8x8pdir[mode][0] == 2) ? best8x8l1ref[mode][0] : -1;
          for (j = img->block_y; j < img->block_y + 4;j++)
          {          
            memset(&enc_picture->ref_idx[LIST_0][j][img->block_x], refl0, 4 * sizeof(char));
            memset(&enc_picture->ref_idx[LIST_1][j][img->block_x], refl1, 4 * sizeof(char));
          }
        }
      }
      else if (mode == 2)
      {
        j = img->block_y;
        refl0 = (best8x8pdir[mode][0] == 0 || best8x8pdir[mode][0] == 2) ? best8x8l0ref[mode][0] : -1;
        refl1 = (best8x8pdir[mode][0] == 1 || best8x8pdir[mode][0] == 2) ? best8x8l1ref[mode][0] : -1;        
        memset(&enc_picture->ref_idx[LIST_0][j  ][img->block_x], refl0, 4 * sizeof(char));
        memset(&enc_picture->ref_idx[LIST_1][j++][img->block_x], refl1, 4 * sizeof(char));
        memset(&enc_picture->ref_idx[LIST_0][j  ][img->block_x], refl0, 4 * sizeof(char));
        memset(&enc_picture->ref_idx[LIST_1][j++][img->block_x], refl1, 4 * sizeof(char));
        
        refl0 = (best8x8pdir[mode][2] == 0 || best8x8pdir[mode][2] == 2) ? best8x8l0ref[mode][2] : -1;
        refl1 = (best8x8pdir[mode][2] == 1 || best8x8pdir[mode][2] == 2) ? best8x8l1ref[mode][2] : -1;
        memset(&enc_picture->ref_idx[LIST_0][j  ][img->block_x], refl0, 4 * sizeof(char));
        memset(&enc_picture->ref_idx[LIST_1][j++][img->block_x], refl1, 4 * sizeof(char));
        memset(&enc_picture->ref_idx[LIST_0][j  ][img->block_x], refl0, 4 * sizeof(char));
        memset(&enc_picture->ref_idx[LIST_1][j++][img->block_x], refl1, 4 * sizeof(char));
      }
      else if (mode == 3)
      {
        j = img->block_y;
        i = img->block_x;
        refl0 = (best8x8pdir[mode][0] == 0 || best8x8pdir[mode][0] == 2) ? best8x8l0ref[mode][0] : -1;
        refl1 = (best8x8pdir[mode][0] == 1 || best8x8pdir[mode][0] == 2) ? best8x8l1ref[mode][0] : -1;
        enc_picture->ref_idx[LIST_0][j  ][i  ] = refl0;
        enc_picture->ref_idx[LIST_1][j  ][i++] = refl1;
        enc_picture->ref_idx[LIST_0][j  ][i  ] = refl0;
        enc_picture->ref_idx[LIST_1][j  ][i++] = refl1;
        refl0 = (best8x8pdir[mode][1] == 0 || best8x8pdir[mode][1] == 2) ? best8x8l0ref[mode][1] : -1;
        refl1 = (best8x8pdir[mode][1] == 1 || best8x8pdir[mode][1] == 2) ? best8x8l1ref[mode][1] : -1;
        enc_picture->ref_idx[LIST_0][j  ][i  ] = refl0;
        enc_picture->ref_idx[LIST_1][j  ][i++] = refl1;
        enc_picture->ref_idx[LIST_0][j  ][i  ] = refl0;
        enc_picture->ref_idx[LIST_1][j++][i  ] = refl1;
        memcpy(&enc_picture->ref_idx[LIST_0][j  ][img->block_x], &enc_picture->ref_idx[LIST_0][img->block_y][img->block_x], 4 * sizeof(char));
        memcpy(&enc_picture->ref_idx[LIST_1][j++][img->block_x], &enc_picture->ref_idx[LIST_1][img->block_y][img->block_x], 4 * sizeof(char));
        memcpy(&enc_picture->ref_idx[LIST_0][j  ][img->block_x], &enc_picture->ref_idx[LIST_0][img->block_y][img->block_x], 4 * sizeof(char));
        memcpy(&enc_picture->ref_idx[LIST_1][j++][img->block_x], &enc_picture->ref_idx[LIST_1][img->block_y][img->block_x], 4 * sizeof(char));
        memcpy(&enc_picture->ref_idx[LIST_0][j  ][img->block_x], &enc_picture->ref_idx[LIST_0][img->block_y][img->block_x], 4 * sizeof(char));
        memcpy(&enc_picture->ref_idx[LIST_1][j++][img->block_x], &enc_picture->ref_idx[LIST_1][img->block_y][img->block_x], 4 * sizeof(char));
      }      
      else
      {
        for (j=0;j<4;j++)
        {
          block_y = img->block_y + j;
          for (i=0;i<4;i++)
          {
            block_x = img->block_x + i;
            k = 2*(j >> 1) + (i >> 1);
            l = 2*(j & 0x01) + (i & 0x01);

            if(mode == P8x8 && best8x8mode[k]==0)
            {
              enc_picture->ref_idx[LIST_0][block_y][block_x] = direct_ref_idx[LIST_0][block_y][block_x];
              enc_picture->ref_idx[LIST_1][block_y][block_x] = direct_ref_idx[LIST_1][block_y][block_x];
            }
            else
            {
              enc_picture->ref_idx[LIST_0][block_y][block_x] = (IS_FW ? best8x8l0ref[mode][k] : -1);
              enc_picture->ref_idx[LIST_1][block_y][block_x] = (IS_BW ? best8x8l1ref[mode][k] : -1);
            }
          }        
        }
      }
    }
    else
    {
      if (mode == 1)
      {
        char refl0 = best8x8pdir[mode][0] == 0 ? best8x8l0ref[mode][0] : -1;
        j = img->block_y;
        memset(&enc_picture->ref_idx[LIST_0][j++][img->block_x], refl0, 4 * sizeof(char));
        memset(&enc_picture->ref_idx[LIST_0][j++][img->block_x], refl0, 4 * sizeof(char));
        memset(&enc_picture->ref_idx[LIST_0][j++][img->block_x], refl0, 4 * sizeof(char));
        memset(&enc_picture->ref_idx[LIST_0][j  ][img->block_x], refl0, 4 * sizeof(char));
      }
      else if (mode == 2)
      {
        refl0 = best8x8pdir[mode][0] == 0 ? best8x8l0ref[mode][0] : -1;
        j = img->block_y;
        memset(&enc_picture->ref_idx[LIST_0][j++][img->block_x], refl0, 4 * sizeof(char));
        memset(&enc_picture->ref_idx[LIST_0][j++][img->block_x], refl0, 4 * sizeof(char));
        refl0 = best8x8pdir[mode][2] == 0 ? best8x8l0ref[mode][2] : -1;
        memset(&enc_picture->ref_idx[LIST_0][j++][img->block_x], refl0, 4 * sizeof(char));
        memset(&enc_picture->ref_idx[LIST_0][j  ][img->block_x], refl0, 4 * sizeof(char));
      }      
      else if (mode == 3)
      {
        j = img->block_y;
        i = img->block_x;
        refl0 = (best8x8pdir[mode][0] == 0) ? best8x8l0ref[mode][0] : -1;
        enc_picture->ref_idx[LIST_0][j  ][i++] = refl0;
        enc_picture->ref_idx[LIST_0][j  ][i++] = refl0;
        refl0 = (best8x8pdir[mode][1] == 0) ? best8x8l0ref[mode][1] : -1;
        enc_picture->ref_idx[LIST_0][j  ][i++] = refl0;
        enc_picture->ref_idx[LIST_0][j++][i  ] = refl0;
        memcpy(&enc_picture->ref_idx[LIST_0][j++][img->block_x], &enc_picture->ref_idx[LIST_0][img->block_y][img->block_x], 4 * sizeof(char));
        memcpy(&enc_picture->ref_idx[LIST_0][j++][img->block_x], &enc_picture->ref_idx[LIST_0][img->block_y][img->block_x], 4 * sizeof(char));
        memcpy(&enc_picture->ref_idx[LIST_0][j  ][img->block_x], &enc_picture->ref_idx[LIST_0][img->block_y][img->block_x], 4 * sizeof(char));
      }      
      else
      {
        for (j=0;j<4;j++)
        {
          block_y = img->block_y + j;
          for (i=0;i<4;i++)
          {
            block_x = img->block_x + i;
            k = 2*(j >> 1) + (i >> 1);
            l = 2*(j & 0x01) + (i & 0x01);
            enc_picture->ref_idx[LIST_0][block_y][block_x] = (IS_FW ? best8x8l0ref[mode][k] : -1);
          }
        }
      }
    }
  }

  if (bframe)
  {
    for (clist = LIST_0; clist <= LIST_1; clist++)
    {
      for (j = img->block_y; j < img->block_y + 4; j++)
        for (i = img->block_x; i < img->block_x + 4;i++)
        {
          cur_ref = (int) enc_picture->ref_idx[clist][j][i];
          enc_picture->ref_pic_id [clist][j][i] = (cur_ref>=0
            ? enc_picture->ref_pic_num[clist + currMB->list_offset][cur_ref]
            : -1);
        }
    }
  }
  else
  {
    for (j = img->block_y; j < img->block_y + 4; j++)
      for (i = img->block_x; i < img->block_x + 4;i++)
      {
        cur_ref = (int) enc_picture->ref_idx[LIST_0][j][i];
        enc_picture->ref_pic_id [LIST_0][j][i] = (cur_ref>=0
          ? enc_picture->ref_pic_num[LIST_0 + currMB->list_offset][cur_ref]
          : -1);
      }
  }

#undef IS_FW
#undef IS_BW
}


/*!
 *************************************************************************************
 * \brief
 *    Intra 16x16 mode decision
 *************************************************************************************
 */
void Intra16x16_Mode_Decision (Macroblock* currMB, int* i16mode)
{
  /* generate intra prediction samples for all 4 16x16 modes */
  intrapred_16x16 (currMB, PLANE_Y);
  if ((active_sps->chroma_format_idc == YUV444)&&(IS_INDEPENDENT(input)==0))
  {
    select_plane(PLANE_U);
    intrapred_16x16 (currMB, PLANE_U);
    select_plane(PLANE_V);
    intrapred_16x16 (currMB, PLANE_V);
    select_plane(PLANE_Y);
  }
  find_sad_16x16 (currMB, i16mode);   /* get best new intra mode */
  
  currMB->cbp = dct_16x16 (currMB, PLANE_Y, *i16mode);
  if ((active_sps->chroma_format_idc == YUV444)&&(IS_INDEPENDENT(input)==0))
  {
    select_plane(PLANE_U);
    cmp_cbp[1] = dct_16x16 (currMB, PLANE_U, *i16mode);
    select_plane(PLANE_V);
    cmp_cbp[2] = dct_16x16 (currMB, PLANE_V, *i16mode);
    select_plane(PLANE_Y);
    currMB->cbp |= cmp_cbp[1];
    currMB->cbp |= cmp_cbp[2];
    cmp_cbp[1] = currMB->cbp;
    cmp_cbp[2] = currMB->cbp;
  }
}



/*!
 *************************************************************************************
 * \brief
 *    Sets Coefficients and reconstruction for an 8x8 block
 *************************************************************************************
 */
void SetCoeffAndReconstruction8x8 (Macroblock* currMB)
{
  int block, k, j, i, uv;
  int cur_ref;

  //============= MIXED TRANSFORM SIZES FOR 8x8 PARTITION ==============
  //--------------------------------------------------------------------
  int l;
  int bframe = img->type==B_SLICE;

  if (currMB->luma_transform_size_8x8_flag)
  {

    //============= set mode and ref. frames ==============
    for(i = 0;i<4;i++)
    {
      currMB->b8mode[i]   = tr8x8.part8x8mode[i];
      currMB->b8pdir[i]   = tr8x8.part8x8pdir[i];
    }

    if (bframe)
    {
      for (j = 0;j<4;j++)
        for (i = 0;i<4;i++)
        {
          k = 2*(j >> 1)+(i >> 1);
          l = 2*(j & 0x01)+(i & 0x01);
          enc_picture->ref_idx[LIST_0][img->block_y+j][img->block_x+i] = ((currMB->b8pdir[k] & 0x01) == 0) ? tr8x8.part8x8l0ref[k] : - 1;
          enc_picture->ref_idx[LIST_1][img->block_y+j][img->block_x+i] = (currMB->b8pdir[k] > 0) ? tr8x8.part8x8l1ref[k] : - 1;
        }
    }
    else
    {
      for (j = 0;j<4;j++)
        for (i = 0;i<4;i++)
        {
          k = 2*(j >> 1)+(i >> 1);
          l = 2*(j & 0x01)+(i & 0x01);
          enc_picture->ref_idx[LIST_0][img->block_y+j][img->block_x+i] = tr8x8.part8x8l0ref[k];
        }
    }


    for (j = img->block_y;j<img->block_y + BLOCK_MULTIPLE;j++)
    {
      for (i = img->block_x;i<img->block_x + BLOCK_MULTIPLE;i++)
      {
        cur_ref = (int) enc_picture->ref_idx[LIST_0][j][i];

        enc_picture->ref_pic_id [LIST_0][j][i] =(cur_ref>=0
          ? enc_picture->ref_pic_num[LIST_0 + currMB->list_offset][cur_ref]
          : -1);
      }
    }

    if (bframe)
    {
      for (j = img->block_y; j < img->block_y + BLOCK_MULTIPLE; j++)
      {
        for (i = img->block_x;i<img->block_x + BLOCK_MULTIPLE;i++)
        {
          cur_ref = (int) enc_picture->ref_idx[LIST_1][j][i];

          enc_picture->ref_pic_id [LIST_1][j][i] = (cur_ref>=0
            ? enc_picture->ref_pic_num[LIST_1 + currMB->list_offset][cur_ref]
            : -1);
        }

      }
    }

    //====== set the mv's for 8x8 partition with transform size 8x8 ======
    //save the mv data for 4x4 transform

    StoreMV8x8(1);
    //set new mv data for 8x8 transform
    RestoreMV8x8(0);

    //============= get pre-calculated data ==============
    //restore coefficients from 8x8 transform

    for (block = 0; block<4; block++)
    {
      for (k = 0; k<4; k++)
        for (j = 0; j<2; j++)
          memcpy (img->cofAC[block][k][j],cofAC8x8ts[0][block][k][j], 65 * sizeof(int));
    }
    
    if( img->yuv_format == YUV444 && !IS_INDEPENDENT(input) )
    {
      for (uv=0; uv<2; uv++)
      {
        for (block = 0; block<4; block++)
        {
          for (k = 0; k<4; k++)
            for (j = 0; j<2; j++)
              memcpy (img->cofAC[4+block+uv*4][k][j],cofAC8x8ts[uv][block][k][j], 65 * sizeof(int));
        }
      }
    }
    //restore reconstruction
    if (cnt_nonz8_8x8ts <= _LUMA_8x8_COEFF_COST_ &&
      ((currMB->qp_scaled[0])!=0 || img->lossless_qpprime_flag==0) &&
      (img->type!=SP_SLICE))// modif ES added last condition (we probably never go there so is the next modification useful ? check)
    {
      currMB->cbp     = 0;
      currMB->cbp_blk = 0;
      for (j = 0; j < MB_BLOCK_SIZE; j++)
      {
        memcpy(&enc_picture->imgY[img->pix_y+j][img->pix_x], tr8x8.mpr8x8[j], MB_BLOCK_SIZE * sizeof(imgpel));
        if(img->type==SP_SLICE &&(!si_frame_indicator && !sp2_frame_indicator ))
          memcpy(&lrec[img->pix_y+j][img->pix_x],tr8x8.lrec[j], MB_BLOCK_SIZE * sizeof(int));
      }
      for(block=0;block<4;block++)
      {
        for( i = 0; i < 4; i++ )
          for( j = 0; j < 2; j++ )
            memset( img->cofAC[block][i][j], 0, 65 * sizeof(int));
      }
      if( img->yuv_format==YUV444 && !IS_INDEPENDENT(input) )
      {
        for (j = 0; j < MB_BLOCK_SIZE; j++)
        {
          memcpy(&enc_picture->imgUV[0][img->pix_y+j][img->pix_x], tr8x8.mpr8x8CbCr[0][j], MB_BLOCK_SIZE * sizeof(imgpel));
          memcpy(&enc_picture->imgUV[1][img->pix_y+j][img->pix_x], tr8x8.mpr8x8CbCr[1][j], MB_BLOCK_SIZE * sizeof(imgpel));
        }
        for (uv=0; uv<2; uv++)
        {
          for (block = 0; block<4; block++)
          {
            for (k = 0; k<4; k++)
              for (j = 0; j<2; j++)
                memset( img->cofAC[4+block+uv*4][k][j], 0, 65 * sizeof(int));
          }
        }
      }
    }
    else
    {
      currMB->cbp     = cbp8_8x8ts;
      currMB->cbp_blk = cbp_blk8_8x8ts;
      for (j = 0; j < MB_BLOCK_SIZE; j++)
      {
        memcpy (&enc_picture->imgY[img->pix_y+j][img->pix_x],tr8x8.rec_mbY8x8[j], MB_BLOCK_SIZE * sizeof(imgpel));
        if(img->type==SP_SLICE &&(!si_frame_indicator && !sp2_frame_indicator))
          memcpy (&lrec[img->pix_y+j][img->pix_x],tr8x8.lrec[j], MB_BLOCK_SIZE * sizeof(int));
      }
      
      if( img->yuv_format == YUV444 && !IS_INDEPENDENT(input) ) 
      {
        cmp_cbp[1] = cmp_cbp[2] = cbp8_8x8ts;
        for (j = 0; j < MB_BLOCK_SIZE; j++)
        {
          memcpy (&enc_picture->imgUV[0][img->pix_y+j][img->pix_x],tr8x8.rec_mbU8x8[j], MB_BLOCK_SIZE * sizeof(imgpel)); 
          memcpy (&enc_picture->imgUV[1][img->pix_y+j][img->pix_x],tr8x8.rec_mbV8x8[j], MB_BLOCK_SIZE * sizeof(imgpel)); 
        }
      }
    }
  }
  else
  {
    //============= get pre-calculated data ==============
    //---------------------------------------------------
    //--- restore coefficients ---
    for (block = 0; block<4+img->num_blk8x8_uv; block++)
    {
      for (k = 0; k<4; k++)
        for (j = 0; j<2; j++)
          memcpy (img->cofAC[block][k][j],cofAC8x8[block][k][j], 65 * sizeof(int));
    }
    if( img->yuv_format == YUV444 && !IS_INDEPENDENT(input) ) 
    {
      for (block = 0; block<4; block++)
      {
        for (k = 0; k<4; k++)
        {
          for (j = 0; j<2; j++)
          {
            memcpy (img->cofAC[block+4][k][j],cofAC8x8CbCr[0][block][k][j], 65 * sizeof(int));     
            memcpy (img->cofAC[block+8][k][j],cofAC8x8CbCr[1][block][k][j], 65 * sizeof(int));   
          }
        }
      }
    }

    if (cnt_nonz_8x8<=5 && img->type!=SP_SLICE &&
      ((currMB->qp_scaled[0])!=0 || img->lossless_qpprime_flag==0))
    {
      currMB->cbp     = 0;
      currMB->cbp_blk = 0;
      for (j = 0; j < MB_BLOCK_SIZE; j++)
      {
        memcpy (&enc_picture->imgY[img->pix_y+j][img->pix_x],tr4x4.mpr8x8[j], MB_BLOCK_SIZE * sizeof(imgpel));
        if(img->type ==SP_SLICE &&(!si_frame_indicator && !sp2_frame_indicator))
          memcpy (&lrec[img->pix_y+j][img->pix_x],tr4x4.lrec[j], MB_BLOCK_SIZE * sizeof(int)); // restore coeff. SP frame
      }
      for(block=0;block<4;block++)
      {
        for( i = 0; i < 4; i++ )
          for( j = 0; j < 2; j++ )
            memset( img->cofAC[block][i][j], 0, 65 * sizeof(int));
      }
      if (img->yuv_format ==YUV444 && !IS_INDEPENDENT(input))
      {
        for (j = 0; j < MB_BLOCK_SIZE; j++)
        {
          memcpy (&enc_picture->imgUV[0][img->pix_y+j][img->pix_x],tr4x4.mpr8x8CbCr[0][j], MB_BLOCK_SIZE * sizeof(imgpel));    
          memcpy (&enc_picture->imgUV[1][img->pix_y+j][img->pix_x],tr4x4.mpr8x8CbCr[1][j], MB_BLOCK_SIZE * sizeof(imgpel));  
        }
        for (uv=0; uv<2; uv++)
        {
          for (block = 0; block<4; block++)
          {
            for (k = 0; k<4; k++)
              for (j = 0; j<2; j++)
                memset( img->cofAC[4+block+uv*4][k][j], 0, 65 * sizeof(int));
          }
        }
      }
    }
    else
    {
      currMB->cbp     = cbp8x8;
      currMB->cbp_blk = cbp_blk8x8;
      for (j = 0; j < MB_BLOCK_SIZE; j++)
      {
        memcpy (&enc_picture->imgY[img->pix_y+j][img->pix_x],tr4x4.rec_mbY8x8[j], MB_BLOCK_SIZE * sizeof(imgpel));
        if(img->type==SP_SLICE &&(!si_frame_indicator && !sp2_frame_indicator))
          memcpy (&lrec[img->pix_y+j][img->pix_x],tr4x4.lrec[j], MB_BLOCK_SIZE * sizeof(int));
      }
      if( img->yuv_format == YUV444 && !IS_INDEPENDENT(input)) 
      {
        cmp_cbp[1] = cmp_cbp[2] = cbp8x8;
        for (j = 0; j < MB_BLOCK_SIZE; j++)
        {
          memcpy (&enc_picture->imgUV[0][img->pix_y+j][img->pix_x],tr4x4.rec_mbU8x8[j], MB_BLOCK_SIZE * sizeof(imgpel));
          memcpy (&enc_picture->imgUV[1][img->pix_y+j][img->pix_x],tr4x4.rec_mbV8x8[j], MB_BLOCK_SIZE * sizeof(imgpel));
        }
      }
    }
  }
}


/*!
 *************************************************************************************
 * \brief
 *    Sets motion vectors for a macroblock
 *************************************************************************************
 */
void SetMotionVectorsMB (Macroblock* currMB, int bframe)
{
  int i, j, k, l, m, mode8, pdir8, ref, by, bx;
  short ******all_mv  = img->all_mv;
  short ******pred_mv = img->pred_mv;
  int  l1_ref;
  int jdiv, jmod;

  if (!bframe)
  {
    for (j = 0; j<4; j++)
    {
      jmod = j & 0x01;
      jdiv = j >>   1;
      by    = img->block_y+j;
      for (i = 0; i<4; i++)
      {
        mode8 = currMB->b8mode[k=2*jdiv+(i>>1)];
        l     = 2*jmod + (i & 0x01);

        bx   = img->block_x+i;

        pdir8 = currMB->b8pdir[k];
        ref    = enc_picture->ref_idx[LIST_0][by][bx];

        if (pdir8>=0)
        {
          enc_picture->mv[LIST_0][by][bx][0] = all_mv [j][i][LIST_0][ ref][mode8][0];
          enc_picture->mv[LIST_0][by][bx][1] = all_mv [j][i][LIST_0][ ref][mode8][1];
        }
        else
        {
          enc_picture->mv[LIST_0][by][bx][0] = 0;
          enc_picture->mv[LIST_0][by][bx][1] = 0;
        }
      }
    }
  }
  else
  {
    for (j = 0; j<4; j++)
    {
      jmod = j & 0x01;
      jdiv = j >>   1;
      by    = img->block_y+j;
      for (i = 0; i<4; i++)
      {
        mode8 = currMB->b8mode[k=2*jdiv+(i>>1)];
        l     = 2*jmod + (i & 0x01);

        bx    = img->block_x+i;

        pdir8 = currMB->b8pdir[k];
        ref    = enc_picture->ref_idx[LIST_0][by][bx];
        l1_ref = enc_picture->ref_idx[LIST_1][by][bx];

        if (currMB->bi_pred_me && (pdir8 == 2) && currMB->mb_type==1)
        {
          all_mv  = currMB->bi_pred_me == 1 ? img->bipred_mv1 : img->bipred_mv2;
          ref = 0;
          l1_ref = 0;
        }

        if (pdir8==-1) // intra
        {
          enc_picture->mv[LIST_0][by][bx][0] = 0;
          enc_picture->mv[LIST_0][by][bx][1] = 0;
          enc_picture->mv[LIST_1][by][bx][0] = 0;
          enc_picture->mv[LIST_1][by][bx][1] = 0;
        }
        else if (pdir8==0) // list 0
        {
          enc_picture->mv[LIST_0][by][bx][0] = all_mv [j][i][LIST_0][ ref][mode8][0];
          enc_picture->mv[LIST_0][by][bx][1] = all_mv [j][i][LIST_0][ ref][mode8][1];
          enc_picture->mv[LIST_1][by][bx][0] = 0;
          enc_picture->mv[LIST_1][by][bx][1] = 0;
          enc_picture->ref_idx[LIST_1][by][bx] = -1;
        }
        else if (pdir8==1) // list 1
        {
          enc_picture->mv[LIST_0][by][bx][0] = 0;
          enc_picture->mv[LIST_0][by][bx][1] = 0;
          enc_picture->ref_idx[LIST_0][by][bx] = -1;
          enc_picture->mv[LIST_1][by][bx][0] = all_mv [j][i][LIST_1][l1_ref][mode8][0];
          enc_picture->mv[LIST_1][by][bx][1] = all_mv [j][i][LIST_1][l1_ref][mode8][1];
        }
        else if (pdir8==2) // bipredictive
        {
          enc_picture->mv[LIST_0][by][bx][0] = all_mv [j][i][LIST_0][ ref][mode8][0];
          enc_picture->mv[LIST_0][by][bx][1] = all_mv [j][i][LIST_0][ ref][mode8][1];
          enc_picture->mv[LIST_1][by][bx][0] = all_mv [j][i][LIST_1][l1_ref][mode8][0];
          enc_picture->mv[LIST_1][by][bx][1] = all_mv [j][i][LIST_1][l1_ref][mode8][1];
        }
        else
        {
          error("invalid direction mode", 255);
        }
      }
    }
  }

  // copy all the motion vectors into rdopt structure
  // Can simplify this by copying the MV's of the best mode (TBD)
  if(img->MbaffFrameFlag)
  {
    for(i = 0;i<4;i++)
    {
      for(j = 0;j<4;j++)
      {
        for (k = 0;k<2;k++)
        {
          for(l = 0;l<img->max_num_references;l++)
          {
            for(m = 0;m<9;m++)
            {
              rdopt->all_mv [j][i][k][l][m][0]  = all_mv [j][i][k][l][m][0];
              rdopt->pred_mv[j][i][k][l][m][0]  = pred_mv[j][i][k][l][m][0];

              rdopt->all_mv [j][i][k][l][m][1]  = all_mv [j][i][k][l][m][1];
              rdopt->pred_mv[j][i][k][l][m][1]  = pred_mv[j][i][k][l][m][1];
            }
          }
        }
      }
    }
  }
}



/*!
 *************************************************************************************
 * \brief
 *    R-D Cost for a macroblock
 *************************************************************************************
 */
int RDCost_for_macroblocks (Macroblock  *currMB,   // <-- Current Macroblock to code
                            double   lambda,       // <-- lagrange multiplier
                            int      mode,         // <-- modus (0-COPY/DIRECT, 1-16x16, 2-16x8, 3-8x16, 4-8x8(+), 5-Intra4x4, 6-Intra16x16)
                            double*  min_rdcost,   // <-> minimum rate-distortion cost
                            double*  min_dcost,   // <-> distortion of mode which has minimum rate-distortion cost.
                            double*  min_rate,     // --> bitrate of mode which has minimum rate-distortion cost.
                            int i16mode )
{
  int         i, j, k; //, k, ****ip4;
  int         j1, j2;
  int         rate = 0, coeff_rate = 0;
  int64       distortion = 0;
  double      rdcost;
  int         prev_mb_nr  = FmoGetPreviousMBNr(img->current_mb_nr);  
  Macroblock  *prevMB   = (prev_mb_nr >= 0) ? &img->mb_data[prev_mb_nr] : NULL;
  int         bframe    = (img->type==B_SLICE);
  int         tmp_cc;
  int         use_of_cc =  (img->type!=I_SLICE &&  input->symbol_mode!=CABAC);
  int         cc_rate, dummy;
  double      dummy_d;
  imgpel     (*curr_mpr)[16] = img->mpr[0];
  imgpel     (*curr_mpr_16x16)[16][16] = img->mpr_16x16[0];

  //=====
  //=====  SET REFERENCE FRAMES AND BLOCK MODES
  //=====
  SetModesAndRefframeForBlocks (currMB, mode);

  //=====
  //=====  GET COEFFICIENTS, RECONSTRUCTIONS, CBP
  //=====
  if (bframe && mode==0)
  {
    int block_x = (img->pix_x >> 2);
    int block_y = (img->pix_y >> 2);
    for (j = block_y; j < block_y + 4;j++)
      for (i = block_x; i < block_x + 4;i++)
        if (direct_pdir[j][i] < 0)
          return 0;
  }

  // Test MV limits for Skip Mode. This could be necessary for MBAFF case Frame MBs.
  if ((img->MbaffFrameFlag) && (!currMB->mb_field) && (img->type==P_SLICE) && (mode==0) )
  {
    if ( img->all_mv[0][0][0][0][0][0] < -8192
      || img->all_mv[0][0][0][0][0][0] > 8191
      || img->all_mv[0][0][0][0][0][1] < LEVELMVLIMIT[img->LevelIndex][4]
      || img->all_mv[0][0][0][0][0][1] > LEVELMVLIMIT[img->LevelIndex][5])
        return 0;
  }

  if (img->AdaptiveRounding)
  {
    memset(&(img->fadjust4x4[0][0][0]), 0, MB_PIXELS * sizeof(int));
    memset(&(img->fadjust8x8[0][0][0]), 0, MB_PIXELS * sizeof(int));
    if (img->yuv_format != 0)
    {
      memset(&(img->fadjust4x4Cr[0][0][0][0]), 0, img->mb_cr_size_y * img->mb_cr_size_x * sizeof(int));
      memset(&(img->fadjust4x4Cr[1][0][0][0]), 0, img->mb_cr_size_y * img->mb_cr_size_x  * sizeof(int));
      memset(&(img->fadjust8x8Cr[0][0][0][0]), 0, img->mb_cr_size_y * img->mb_cr_size_x  * sizeof(int));
      memset(&(img->fadjust8x8Cr[1][0][0][0]), 0, img->mb_cr_size_y * img->mb_cr_size_x  * sizeof(int));
    }
  }
  
  if (mode<P8x8)
  {
    LumaResidualCoding (currMB);

    // This code seems unnecessary 
    if(mode==0 && currMB->cbp!=0 && (img->type != B_SLICE || img->NoResidueDirect==1))
      return 0;
    if (img->yuv_format==YUV444 && !IS_INDEPENDENT(input))      
    {
      if(mode==0 && (currMB->cbp==0 && cmp_cbp[1] == 0 && cmp_cbp[2] == 0)&& currMB->luma_transform_size_8x8_flag == 1) //for B_skip, luma_transform_size_8x8_flag=0 only
        return 0;
    }
    else
    {
      if(mode==0 && currMB->cbp==0 && currMB->luma_transform_size_8x8_flag == 1) //for B_skip, luma_transform_size_8x8_flag=0 only        
        return 0;
    }

  }
  else if (mode==P8x8)
  {
    SetCoeffAndReconstruction8x8 (currMB);
  }
  else if (mode==I4MB)
  {
    currMB->cbp = Mode_Decision_for_Intra4x4Macroblock (currMB, lambda, &dummy_d);
  }
  else if (mode==I16MB)
  {
    Intra16x16_Mode_Decision  (currMB, &i16mode);
  }
  else if(mode==I8MB)
  {
    currMB->cbp = Mode_Decision_for_new_Intra8x8Macroblock(currMB, lambda, &dummy_d);
  }
  else if(mode==IPCM)
  {
    for (j = 0; j < MB_BLOCK_SIZE; j++)
    {
      memcpy(&enc_picture->imgY[j + img->pix_y][img->opix_x], &pCurImg[j + img->opix_y][img->opix_x], MB_BLOCK_SIZE * sizeof(imgpel));
    }
    if ((img->yuv_format != YUV400) && !IS_INDEPENDENT(input))
    {
      // CHROMA
      for (j = 0; j<img->mb_cr_size_y; j++)
      {
        j1 = j + img->opix_c_y;
        j2 = j + img->pix_c_y;
        memcpy(&enc_picture->imgUV[0][j2][img->opix_c_x], &imgUV_org[0][j1][img->opix_c_x], img->mb_cr_size_x * sizeof(imgpel));
        memcpy(&enc_picture->imgUV[1][j2][img->opix_c_x], &imgUV_org[1][j1][img->opix_c_x], img->mb_cr_size_x * sizeof(imgpel));
      }
    }
    for (j=0;j<4;j++)
      for (i=0; i<(4+img->num_blk8x8_uv); i++)
        img->nz_coeff[img->current_mb_nr][j][i] = 16;

  }

  if (input->rdopt==3 && img->type!=B_SLICE)
  {
    // We need the reconstructed prediction residue for the simulated decoders.
    compute_residue_mb (mode == I16MB ? i16mode : -1);
  }

  //Rate control
  if (input->RCEnable)
  {
    if (mode == I16MB)
      memcpy(pred, curr_mpr_16x16[i16mode], MB_PIXELS * sizeof(imgpel));
    else
      memcpy(pred, curr_mpr, MB_PIXELS * sizeof(imgpel));
  }

  img->i16offset = 0;
  dummy = 0;

  if (((img->yuv_format!=YUV400) && (active_sps->chroma_format_idc != YUV444)) && (mode != IPCM))
    ChromaResidualCoding (currMB);

  if (mode==I16MB)
    img->i16offset = I16Offset  (currMB->cbp, i16mode);

  //=====
  //=====   GET DISTORTION
  //=====
  // LUMA
  if (input->rdopt==3 && img->type!=B_SLICE)
  {
    for (k = 0; k<input->NoOfDecoders ;k++)
    {
      decode_one_mb (k, currMB);
      for (j = 0; j<MB_BLOCK_SIZE; j++)
      {
        for (i=img->opix_x; i<img->opix_x+MB_BLOCK_SIZE; i++)
          distortion += iabs2( pCurImg[img->opix_y+j][i] - decs->decY[k][img->opix_y+j][i] );
      }
    }
    distortion /= input->NoOfDecoders;

    if ((img->yuv_format != YUV400) && (active_sps->chroma_format_idc != YUV444))
    {
      // CHROMA
      for (j = 0; j<img->mb_cr_size_y; j++)
      {
        j1 = j + img->opix_c_y;
        j2 = j + img->pix_c_y;
        for (i = img->opix_c_x; i < img->opix_c_x + img->mb_cr_size_x; i++)
        {
          distortion += iabs2( imgUV_org[0][j1][i] - enc_picture->imgUV[0][j2][i] );
          distortion += iabs2( imgUV_org[1][j1][i] - enc_picture->imgUV[1][j2][i] );
        }
      }
    }
  }
  else
  {
    distortion = getDistortion(currMB);
  }

  //=====   S T O R E   C O D I N G   S T A T E   =====
  //---------------------------------------------------
  store_coding_state (currMB, cs_cm);

  //=====
  //=====   GET RATE
  //=====
  //----- macroblock header -----
  if (use_of_cc)
  {
    if (currMB->mb_type!=0 || (bframe && currMB->cbp!=0))
    {
      // cod counter and macroblock mode are written ==> do not consider code counter
      tmp_cc = img->cod_counter;
      rate   = writeMBLayer (currMB, 1, &coeff_rate);
      ue_linfo (tmp_cc, dummy, &cc_rate, &dummy);
      rate  -= cc_rate;
      img->cod_counter = tmp_cc;
    }
    else
    {
      // cod counter is just increased  ==> get additional rate
      ue_linfo (img->cod_counter + 1, dummy, &rate,    &dummy);
      ue_linfo (img->cod_counter    , dummy, &cc_rate, &dummy);
      rate -= cc_rate;
    }
  }
  else
  {
    rate = writeMBLayer (currMB, 1, &coeff_rate);
  }

  //=====   R E S T O R E   C O D I N G   S T A T E   =====
  //-------------------------------------------------------
  reset_coding_state (currMB, cs_cm);

  rdcost = (double)distortion + lambda * dmax(0.5,(double)rate);

  if (rdcost >= *min_rdcost ||
    ((currMB->qp_scaled[0]) == 0 && img->lossless_qpprime_flag == 1 && distortion != 0))
  {
#if FASTMODE
    // Reordering RDCost comparison order of mode 0 and mode 1 in P_SLICE
    // if RDcost of mode 0 and mode 1 is same, we choose best_mode is 0
    // This might not always be good since mode 0 is more biased towards rate than quality.
    if((img->type!=P_SLICE || mode != 0 || rdcost != *min_rdcost) || IS_FREXT_PROFILE(input->ProfileIDC))
#endif
      return 0;
  }


  if ((img->MbaffFrameFlag) && (mode ? 0: ((img->type == B_SLICE) ? !currMB->cbp:1)))  // AFF and current is skip
  {
    if (img->current_mb_nr & 0x01) //bottom
    {
      if (prevMB->mb_type ? 0:((img->type == B_SLICE) ? !prevMB->cbp:1)) //top is skip
      {
        if (!(field_flag_inference(currMB) == currMB->mb_field)) //skip only allowed when correct inference
          return 0;
      }
    }
  }

  //=====   U P D A T E   M I N I M U M   C O S T   =====
  //-----------------------------------------------------
  *min_rdcost = rdcost;
  *min_dcost = (double) distortion;
  *min_rate = lambda * (double)coeff_rate;

#ifdef BEST_NZ_COEFF
  for (j=0;j<4;j++)
    memcpy(&gaaiMBAFF_NZCoeff[j][0], &img->nz_coeff[img->current_mb_nr][j][0], (4+img->num_blk8x8_uv) * sizeof(int));
#endif

  return 1;
}

/*!
 *************************************************************************************
 * \brief
 *    Store adaptive rounding parameters
 *************************************************************************************
 */
void store_adaptive_rounding_parameters (Macroblock *currMB, int mode)
{
  int j;
  int is_inter = (mode != I4MB)&&(mode != I16MB)&&(mode != I8MB);

  if (currMB->luma_transform_size_8x8_flag)
  {
    if ((mode == P8x8))
      memcpy(&(bestInterFAdjust8x8[0][0]),&(img->fadjust8x8[2][0][0]),MB_PIXELS * sizeof(int));
    else if (is_inter)
      memcpy(&(bestInterFAdjust8x8[0][0]),&(img->fadjust8x8[0][0][0]),MB_PIXELS * sizeof(int));
    else
      memcpy(&(bestIntraFAdjust8x8[0][0]),&(img->fadjust8x8[1][0][0]),MB_PIXELS * sizeof(int));
  }
  else
  {
    if ((mode == P8x8))
      memcpy(&(bestInterFAdjust4x4[0][0]),&(img->fadjust4x4[3][0][0]),MB_PIXELS * sizeof(int));
    else if (is_inter)
      memcpy(&(bestInterFAdjust4x4[0][0]),&(img->fadjust4x4[0][0][0]),MB_PIXELS * sizeof(int));
    else
      memcpy(&(bestIntraFAdjust4x4[0][0]),&(img->fadjust4x4[1 + mode == I16MB][0][0]),MB_PIXELS * sizeof(int));
  }

  if ((active_sps->chroma_format_idc == YUV444)&& !IS_INDEPENDENT(input))
  {
    if (currMB->luma_transform_size_8x8_flag)
    {
      if ((mode == P8x8))
      {
        memcpy(&(bestInterFAdjust8x8Cr[0][0][0]),&(img->fadjust8x8Cr[0][2][0][0]),MB_PIXELS * sizeof(int));
        memcpy(&(bestInterFAdjust8x8Cr[1][0][0]),&(img->fadjust8x8Cr[1][2][0][0]),MB_PIXELS * sizeof(int));
      }
      else if (is_inter)
      {
        memcpy(&(bestInterFAdjust8x8Cr[0][0][0]),&(img->fadjust8x8Cr[0][0][0][0]),MB_PIXELS * sizeof(int));
        memcpy(&(bestInterFAdjust8x8Cr[1][0][0]),&(img->fadjust8x8Cr[1][0][0][0]),MB_PIXELS * sizeof(int));
      }
      else
      {
        memcpy(&(bestIntraFAdjust8x8Cr[0][0][0]),&(img->fadjust8x8Cr[0][1][0][0]),MB_PIXELS * sizeof(int));
        memcpy(&(bestIntraFAdjust8x8Cr[1][0][0]),&(img->fadjust8x8Cr[1][1][0][0]),MB_PIXELS * sizeof(int));
      }
    }
    else
    {
      if ((mode == P8x8))
      {
        memcpy(&(bestInterFAdjust4x4Cr[0][0][0]),&(img->fadjust4x4Cr[0][3][0][0]),MB_PIXELS * sizeof(int));
        memcpy(&(bestInterFAdjust4x4Cr[1][0][0]),&(img->fadjust4x4Cr[1][3][0][0]),MB_PIXELS * sizeof(int));
      }
      else if (is_inter)
      {
        memcpy(&(bestInterFAdjust4x4Cr[0][0][0]),&(img->fadjust4x4Cr[0][0][0][0]),MB_PIXELS * sizeof(int));
        memcpy(&(bestInterFAdjust4x4Cr[1][0][0]),&(img->fadjust4x4Cr[1][0][0][0]),MB_PIXELS * sizeof(int));
      }
      else
      {
        memcpy(&(bestIntraFAdjust4x4Cr[0][0][0]),&(img->fadjust4x4Cr[0][1 + mode == I16MB][0][0]),MB_PIXELS * sizeof(int));
        memcpy(&(bestIntraFAdjust4x4Cr[1][0][0]),&(img->fadjust4x4Cr[1][1 + mode == I16MB][0][0]),MB_PIXELS * sizeof(int));
      }
    }
  }
    if ( (input->yuv_format == YUV420 || input->yuv_format == YUV422) && (input->AdaptRndChroma))  
  {
    if (currMB->luma_transform_size_8x8_flag && mode == P8x8)
    {
      for (j = 0; j < img->mb_cr_size_y; j++)
      {
        memcpy(bestInterFAdjust4x4Cr[0][j],img->fadjust8x8Cr[0][0][j],img->mb_cr_size_x * sizeof(int));
        memcpy(bestInterFAdjust4x4Cr[1][j],img->fadjust8x8Cr[1][0][j],img->mb_cr_size_x * sizeof(int));
      }
    }
    else if (mode == P8x8)
    {
      for (j = 0; j < img->mb_cr_size_y; j++)
      {
        memcpy(bestInterFAdjust4x4Cr[0][j],img->fadjust4x4Cr[0][2][j],img->mb_cr_size_x * sizeof(int));
        memcpy(bestInterFAdjust4x4Cr[1][j],img->fadjust4x4Cr[1][2][j],img->mb_cr_size_x * sizeof(int));
      }
    }
    else if (is_inter)
    {
      for (j = 0; j < img->mb_cr_size_y; j++)
      {
        memcpy(bestInterFAdjust4x4Cr[0][j],img->fadjust4x4Cr[0][0][j],img->mb_cr_size_x * sizeof(int));
        memcpy(bestInterFAdjust4x4Cr[1][j],img->fadjust4x4Cr[1][0][j],img->mb_cr_size_x * sizeof(int));
      }
    }
    else
    {
      for (j = 0; j < img->mb_cr_size_y; j++)
      {
        memcpy(bestIntraFAdjust4x4Cr[0][j],img->fadjust4x4Cr[0][1][j],img->mb_cr_size_x * sizeof(int));
        memcpy(bestIntraFAdjust4x4Cr[1][j],img->fadjust4x4Cr[1][1][j],img->mb_cr_size_x * sizeof(int));
      }
    }
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Store macroblock parameters
 *************************************************************************************
 */
void store_macroblock_parameters (Macroblock *currMB, int mode)
{
  int  i, j, k, ****i4p, ***i3p;
  int        bframe   = (img->type==B_SLICE);

  //--- store best mode ---
  best_mode = mode;
  best_c_imode = currMB->c_ipred_mode;
  best_i16offset = img->i16offset;

  // If condition is not really necessary.
  bi_pred_me = (mode == 1) ? currMB->bi_pred_me : 0;

  memcpy(b8mode, currMB->b8mode, BLOCK_MULTIPLE * sizeof(int));
  memcpy(b8pdir, currMB->b8pdir, BLOCK_MULTIPLE * sizeof(int));
  memcpy(b4_intra_pred_modes,currMB->intra_pred_modes, MB_BLOCK_PARTITIONS * sizeof(char));
  memcpy(b8_intra_pred_modes8x8,currMB->intra_pred_modes8x8, MB_BLOCK_PARTITIONS * sizeof(char));

  for (j = 0 ; j < BLOCK_MULTIPLE; j++)
  {
    memcpy(&b4_ipredmode[j * BLOCK_MULTIPLE],&img->ipredmode[img->block_y + j][img->block_x],BLOCK_MULTIPLE * sizeof(char));
    memcpy(b8_ipredmode8x8[j],&img->ipredmode8x8[img->block_y + j][img->block_x],BLOCK_MULTIPLE * sizeof(char));
  }
  //--- reconstructed blocks ----
  for (j = 0; j < MB_BLOCK_SIZE; j++)
  {
    memcpy(rec_mbY[j],&enc_picture->imgY[img->pix_y+j][img->pix_x], MB_BLOCK_SIZE * sizeof(imgpel));
  }
  if((img->type==SP_SLICE) && (si_frame_indicator==0 && sp2_frame_indicator==0))
  {
    for (j = 0; j < MB_BLOCK_SIZE; j++)
    {
      memcpy(lrec_rec[j],&lrec[img->pix_y+j][img->pix_x], MB_BLOCK_SIZE * sizeof(int));//store coefficients SP frame
    }
  }

  if (img->AdaptiveRounding)
    store_adaptive_rounding_parameters (currMB, mode);

  if (img->yuv_format != YUV400)
  {
    for (j = 0; j<img->mb_cr_size_y; j++)
    {
      memcpy(rec_mbU[j],&enc_picture->imgUV[0][img->pix_c_y+j][img->pix_c_x], img->mb_cr_size_x * sizeof(imgpel));
      memcpy(rec_mbV[j],&enc_picture->imgUV[1][img->pix_c_y+j][img->pix_c_x], img->mb_cr_size_x * sizeof(imgpel));
    }
    if((img->type==SP_SLICE) && (si_frame_indicator==0 && sp2_frame_indicator==0))
    {
      //store uv coefficients SP frame
      for (j = 0; j<img->mb_cr_size_y; j++)
      {
        memcpy(lrec_rec_U[j],&lrec_uv[0][img->pix_c_y+j][img->pix_c_x], img->mb_cr_size_x * sizeof(int));
        memcpy(lrec_rec_V[j],&lrec_uv[1][img->pix_c_y+j][img->pix_c_x], img->mb_cr_size_x * sizeof(int));
      }
    }
  }

  //--- store results of decoders ---
  if (input->rdopt==3 && img->type!=B_SLICE)
  {
    for (k = 0; k<input->NoOfDecoders; k++)
    {
      for (j=img->pix_y; j<img->pix_y+16; j++)
        for (i=img->pix_x; i<img->pix_x+16; i++)
        {
          // Keep the decoded values of each MB for updating the ref frames
          decs->decY_best[k][j][i] = decs->decY[k][j][i];
        }
    }
  }

  //--- coeff, cbp, kac ---
  if (mode || bframe)
  {
    i4p=cofAC; cofAC=img->cofAC; img->cofAC=i4p;
    i3p=cofDC; cofDC=img->cofDC; img->cofDC=i3p;
    cbp     = currMB->cbp;
    curr_cbp[0] = cmp_cbp[1];  
    curr_cbp[1] = cmp_cbp[2]; 

    cur_cbp_blk[0] = currMB->cbp_blk;
  }
  else
  {
    cur_cbp_blk[0] = cbp = 0;
    cmp_cbp[1] = cmp_cbp[2] = 0; 
  }

  //--- store transform size ---
  luma_transform_size_8x8_flag = currMB->luma_transform_size_8x8_flag;


  for (j = 0; j<4; j++)
    memcpy(l0_refframe[j],&enc_picture->ref_idx[LIST_0][img->block_y+j][img->block_x], BLOCK_MULTIPLE * sizeof(char));

  if (bframe)
  {
    for (j = 0; j<4; j++)
      memcpy(l1_refframe[j],&enc_picture->ref_idx[LIST_1][img->block_y+j][img->block_x], BLOCK_MULTIPLE * sizeof(char));
  }
}


/*!
 *************************************************************************************
 * \brief
 *    Set stored macroblock parameters
 *************************************************************************************
 */
void set_stored_macroblock_parameters (Macroblock *currMB)
{
  imgpel     **imgY  = enc_picture->imgY;
  imgpel    ***imgUV = enc_picture->imgUV;

  int         mode   = best_mode;
  int         bframe = (img->type==B_SLICE);
  int         i, j, k, ****i4p, ***i3p;
  int         block_x, block_y;
  char    **ipredmodes = img->ipredmode;
  short   *cur_mv;

  //===== reconstruction values =====
  for (j = 0; j < MB_BLOCK_SIZE; j++)
  {
    memcpy(&imgY[img->pix_y+j][img->pix_x],rec_mbY[j], MB_BLOCK_SIZE * sizeof(imgpel));
  }
  if((img->type==SP_SLICE) &&(si_frame_indicator==0 && sp2_frame_indicator==0 ))
  {
    for (j = 0; j < MB_BLOCK_SIZE; j++)
      memcpy(&lrec[img->pix_y+j][img->pix_x],lrec_rec[j], MB_BLOCK_SIZE * sizeof(int)); //restore coeff SP frame
  }
  if(img->MbaffFrameFlag)
  {
    for (j = 0; j < MB_BLOCK_SIZE; j++)
      memcpy(rdopt->rec_mbY[j],rec_mbY[j], MB_BLOCK_SIZE * sizeof(imgpel));
  }

  if (img->AdaptiveRounding)
  {
    update_offset_params(currMB, mode,luma_transform_size_8x8_flag);
  }

  if (img->yuv_format != YUV400)
  {
    for (j = 0; j<img->mb_cr_size_y; j++)
    {
      memcpy(&imgUV[0][img->pix_c_y+j][img->pix_c_x],rec_mbU[j], img->mb_cr_size_x * sizeof(imgpel));
      memcpy(&imgUV[1][img->pix_c_y+j][img->pix_c_x],rec_mbV[j], img->mb_cr_size_x * sizeof(imgpel));
      if((img->type==SP_SLICE) &&(!si_frame_indicator && !sp2_frame_indicator))
      {
        memcpy(&lrec_uv[0][img->pix_c_y+j][img->pix_c_x],lrec_rec_U[j], img->mb_cr_size_x * sizeof(int));
        memcpy(&lrec_uv[1][img->pix_c_y+j][img->pix_c_x],lrec_rec_V[j], img->mb_cr_size_x * sizeof(int));
      }
      if(img->MbaffFrameFlag)
      {
        memcpy(rdopt->rec_mbU[j],rec_mbU[j], img->mb_cr_size_x * sizeof(imgpel));
        memcpy(rdopt->rec_mbV[j],rec_mbV[j], img->mb_cr_size_x * sizeof(imgpel));
      }
    }

    if((img->type==SP_SLICE) &&(!si_frame_indicator && !sp2_frame_indicator))
    {
      for (j = 0; j<img->mb_cr_size_y; j++)
      {
        memcpy(&lrec_uv[0][img->pix_c_y+j][img->pix_c_x],lrec_rec_U[j], img->mb_cr_size_x * sizeof(int));
        memcpy(&lrec_uv[1][img->pix_c_y+j][img->pix_c_x],lrec_rec_V[j], img->mb_cr_size_x * sizeof(int));
      }
    }

    if(img->MbaffFrameFlag)
    {
      for (j = 0; j<img->mb_cr_size_y; j++)
      {

        memcpy(rdopt->rec_mbU[j],rec_mbU[j], img->mb_cr_size_x * sizeof(imgpel));
        memcpy(rdopt->rec_mbV[j],rec_mbV[j], img->mb_cr_size_x * sizeof(imgpel));
      }
    }
  }

  //===== coefficients and cbp =====
  i4p=cofAC; cofAC=img->cofAC; img->cofAC=i4p;
  i3p=cofDC; cofDC=img->cofDC; img->cofDC=i3p;
  currMB->cbp      = cbp;
  currMB->cbp_blk = cur_cbp_blk[0];
  cmp_cbp[1] = curr_cbp[0]; 
  cmp_cbp[2] = curr_cbp[1]; 
  currMB->cbp |= cmp_cbp[1];
  currMB->cbp |= cmp_cbp[2];
  cmp_cbp[1] = currMB->cbp; 
  cmp_cbp[2] = currMB->cbp;

  //==== macroblock type ====
  currMB->mb_type = mode;

  if(img->MbaffFrameFlag)
  {
    rdopt->mode = mode;
    rdopt->i16offset = img->i16offset;
    rdopt->cbp = cbp;
    rdopt->cbp_blk = cur_cbp_blk[0];
    rdopt->mb_type  = mode;

    rdopt->prev_qp  = currMB->prev_qp;
    rdopt->prev_dqp = currMB->prev_dqp;
    rdopt->delta_qp = currMB->delta_qp;
    rdopt->qp       = currMB->qp;
    memcpy(rdopt->qp_scaled, currMB->qp_scaled, 3 * sizeof(int));
    rdopt->prev_cbp = currMB->prev_cbp;

    for(i = 0;i<4+img->num_blk8x8_uv;i++)
    {
      for(j = 0;j<4;j++)
        for(k = 0;k<2;k++)
          memcpy(rdopt->cofAC[i][j][k], img->cofAC[i][j][k], 65 * sizeof(int));
    }
    for(i = 0;i<3;i++)
      for(k = 0;k<2;k++)
        memcpy(rdopt->cofDC[i][k], img->cofDC[i][k], 18 * sizeof(int));
  }


  memcpy(currMB->b8mode,b8mode, BLOCK_MULTIPLE * sizeof(int));
  memcpy(currMB->b8pdir,b8pdir, BLOCK_MULTIPLE * sizeof(int));
  if(img->MbaffFrameFlag)
  {
    memcpy(rdopt->b8mode,b8mode, BLOCK_MULTIPLE * sizeof(int));
    memcpy(rdopt->b8pdir,b8pdir, BLOCK_MULTIPLE * sizeof(int));
  }

  currMB->bi_pred_me = currMB->mb_type == 1 ? bi_pred_me : 0;


  //if P8x8 mode and transform size 4x4 choosen, restore motion vector data for this transform size
  if (mode == P8x8 && !luma_transform_size_8x8_flag && input->Transform8x8Mode)
    RestoreMV8x8(1);

  //==== transform size flag ====
  if (img->yuv_format==YUV444 && !IS_INDEPENDENT(input))
  {
    if (((currMB->cbp == 0) && cmp_cbp[1] == 0 && cmp_cbp[2] == 0) && !(IS_OLDINTRA(currMB) || currMB->mb_type == I8MB))
      currMB->luma_transform_size_8x8_flag = 0;
    else
      currMB->luma_transform_size_8x8_flag = luma_transform_size_8x8_flag;

  }
  else
  {

    if (((currMB->cbp & 15) == 0) && !(IS_OLDINTRA(currMB) || currMB->mb_type == I8MB))
      currMB->luma_transform_size_8x8_flag = 0;
    else
      currMB->luma_transform_size_8x8_flag = luma_transform_size_8x8_flag;
  }

  rdopt->luma_transform_size_8x8_flag  = currMB->luma_transform_size_8x8_flag;

  if (input->rdopt==3 && img->type!=B_SLICE)
  {
    //! save the MB Mode of every macroblock
    decs->dec_mb_mode[img->mb_y][img->mb_x] = mode;
  }

  //==== reference frames =====
  for (j = 0; j < 4; j++)
  {
    block_y = img->block_y + j;
    for (i = 0; i < 4; i++)
    {
      block_x = img->block_x + i;
      k = 2*(j >> 1)+(i >> 1);

      // backward prediction or intra
      if ((currMB->b8pdir[k] == 1) || IS_INTRA(currMB))
      {
        enc_picture->ref_idx    [LIST_0][block_y][block_x]    = -1;
        enc_picture->ref_pic_id [LIST_0][block_y][block_x]    = -1;
        enc_picture->mv         [LIST_0][block_y][block_x][0] = 0;
        enc_picture->mv         [LIST_0][block_y][block_x][1] = 0;
        if(img->MbaffFrameFlag)
          rdopt->refar[LIST_0][j][i] = -1;
      }
      else
      {
        if (currMB->bi_pred_me && (currMB->b8pdir[k] == 2) && currMB->mb_type==1)
        {
          cur_mv = currMB->bi_pred_me == 1
            ? img->bipred_mv1[j][i][LIST_0][0][currMB->b8mode[k]]
            : img->bipred_mv2[j][i][LIST_0][0][currMB->b8mode[k]];

          enc_picture->ref_idx    [LIST_0][block_y][block_x] = 0;
          enc_picture->ref_pic_id [LIST_0][block_y][block_x] = enc_picture->ref_pic_num[LIST_0 + currMB->list_offset][0];
          enc_picture->mv         [LIST_0][block_y][block_x][0] = cur_mv[0];
          enc_picture->mv         [LIST_0][block_y][block_x][1] = cur_mv[1];
          if(img->MbaffFrameFlag)
            rdopt->refar[LIST_0][j][i] = 0;
        }
        else
        {
          cur_mv = img->all_mv[j][i][LIST_0][(short)l0_refframe[j][i]][currMB->b8mode[k]];

          enc_picture->ref_idx    [LIST_0][block_y][block_x]    = l0_refframe[j][i];
          enc_picture->ref_pic_id [LIST_0][block_y][block_x]    = enc_picture->ref_pic_num[LIST_0 + currMB->list_offset][(short)l0_refframe[j][i]];
          enc_picture->mv         [LIST_0][block_y][block_x][0] = cur_mv[0];
          enc_picture->mv         [LIST_0][block_y][block_x][1] = cur_mv[1];
          if(img->MbaffFrameFlag)
            rdopt->refar[LIST_0][j][i] = l0_refframe[j][i];
        }
      }

      // forward prediction or intra
      if ((currMB->b8pdir[k] == 0) || IS_INTRA(currMB))
      {
        enc_picture->ref_idx    [LIST_1][block_y][block_x]    = -1;
        enc_picture->ref_pic_id [LIST_1][block_y][block_x]    = -1;
        enc_picture->mv         [LIST_1][block_y][block_x][0] = 0;
        enc_picture->mv         [LIST_1][block_y][block_x][1] = 0;
        if(img->MbaffFrameFlag)
          rdopt->refar[LIST_1][j][i] = -1;
      }
    }
  }

  if (bframe)
  {
    for (j=0; j<4; j++)
    {
      block_y = img->block_y + j;
      for (i=0; i<4; i++)
      {
        block_x = img->block_x + i;
        k = 2*(j >> 1)+(i >> 1);

        // forward
        if (IS_INTRA(currMB)||(currMB->b8pdir[k] == 0))
        {
          enc_picture->ref_idx    [LIST_1][block_y][block_x]    = -1;
          enc_picture->ref_pic_id [LIST_1][block_y][block_x]    = -1;
          enc_picture->mv         [LIST_1][block_y][block_x][0] = 0;
          enc_picture->mv         [LIST_1][block_y][block_x][1] = 0;
          if(img->MbaffFrameFlag)
            rdopt->refar[LIST_1][j][i] = -1;
        }
        else
        {
          if (currMB->bi_pred_me && (currMB->b8pdir[k] == 2) && currMB->mb_type==1)
          {
            cur_mv = currMB->bi_pred_me == 1
              ? img->bipred_mv1[j][i][LIST_1][0][currMB->b8mode[k]]
              : img->bipred_mv2[j][i][LIST_1][0][currMB->b8mode[k]];

            enc_picture->ref_idx    [LIST_1][block_y][block_x] = 0;
            enc_picture->ref_pic_id [LIST_1][block_y][block_x] = enc_picture->ref_pic_num[LIST_1 + currMB->list_offset][0];
            enc_picture->mv         [LIST_1][block_y][block_x][0] = cur_mv[0];
            enc_picture->mv         [LIST_1][block_y][block_x][1] = cur_mv[1];
            if(img->MbaffFrameFlag)
              rdopt->refar[LIST_1][j][i] = 0;
          }
          else
          {
            cur_mv = img->all_mv[j][i][LIST_1][(short)l1_refframe[j][i]][currMB->b8mode[k]];

            enc_picture->ref_idx    [LIST_1][block_y][block_x] = l1_refframe[j][i];
            enc_picture->ref_pic_id [LIST_1][block_y][block_x] = enc_picture->ref_pic_num[LIST_1 + currMB->list_offset][(short)l1_refframe[j][i]];
            enc_picture->mv         [LIST_1][block_y][block_x][0] = cur_mv[0];
            enc_picture->mv         [LIST_1][block_y][block_x][1] = cur_mv[1];
            if(img->MbaffFrameFlag)
              rdopt->refar[LIST_1][j][i] = l1_refframe[j][i];
          }
        }
      }
    }
  }

  //==== intra prediction modes ====
  currMB->c_ipred_mode = best_c_imode;
  img->i16offset = best_i16offset;

  if(currMB->mb_type == I8MB)
  {
    memcpy(currMB->intra_pred_modes8x8,b8_intra_pred_modes8x8, MB_BLOCK_PARTITIONS * sizeof(char));
    memcpy(currMB->intra_pred_modes,b8_intra_pred_modes8x8, MB_BLOCK_PARTITIONS * sizeof(char));
    for(j = 0; j < BLOCK_MULTIPLE; j++)
    {
      memcpy(&img->ipredmode[img->block_y+j][img->block_x],b8_ipredmode8x8[j], BLOCK_MULTIPLE * sizeof(char));
      memcpy(&img->ipredmode8x8[img->block_y+j][img->block_x], b8_ipredmode8x8[j], BLOCK_MULTIPLE * sizeof(char));
    }
  }
  else if (mode!=I4MB && mode!=I8MB)
  {
    memset(currMB->intra_pred_modes,DC_PRED, MB_BLOCK_PARTITIONS * sizeof(char));
    for(j = img->block_y; j < img->block_y + BLOCK_MULTIPLE; j++)
      memset(&img->ipredmode[j][img->block_x], DC_PRED, BLOCK_MULTIPLE * sizeof(char));
  }
  // Residue Color Transform
  else if (mode == I4MB)
  {
    memcpy(currMB->intra_pred_modes,b4_intra_pred_modes, MB_BLOCK_PARTITIONS * sizeof(char));
    for(j = 0; j < BLOCK_MULTIPLE; j++)
      memcpy(&img->ipredmode[img->block_y + j][img->block_x],&b4_ipredmode[BLOCK_MULTIPLE * j], BLOCK_MULTIPLE * sizeof(char));
  }

  if(img->MbaffFrameFlag)
  {
    rdopt->c_ipred_mode = currMB->c_ipred_mode;
    rdopt->i16offset = img->i16offset;
    memcpy(rdopt->intra_pred_modes,currMB->intra_pred_modes, MB_BLOCK_PARTITIONS * sizeof(char));
    memcpy(rdopt->intra_pred_modes8x8,currMB->intra_pred_modes8x8, MB_BLOCK_PARTITIONS * sizeof(char));
    for(j = img->block_y; j < img->block_y +BLOCK_MULTIPLE; j++)
      memcpy(&rdopt->ipredmode[j][img->block_x],&ipredmodes[j][img->block_x], BLOCK_MULTIPLE * sizeof(char));
  }

  //==== motion vectors =====
  SetMotionVectorsMB (currMB, bframe);
}



/*!
 *************************************************************************************
 * \brief
 *    Set reference frames and motion vectors
 *************************************************************************************
 */
void SetRefAndMotionVectors (Macroblock *currMB, int block, int mode, int pdir, int fwref, int bwref)
{
  int     i, j=0;
  int     bslice  = (img->type==B_SLICE);
  int     pmode   = (mode==1||mode==2||mode==3?mode:4);
  int     j0      = ((block >> 1)<<1);
  int     i0      = ((block & 0x01)<<1);
  int     j1      = j0 + (input->part_size[pmode][1]);
  int     i1      = i0 + (input->part_size[pmode][0]);
  int     block_x, block_y;
  short   *cur_mv;

  if (pdir<0)
  {
    for (j = img->block_y + j0; j < img->block_y + j1; j++)
    {
      for (i=img->block_x + i0; i<img->block_x +i1; i++)
      {
        enc_picture->ref_pic_id[LIST_0][j][i] = -1;
        enc_picture->ref_pic_id[LIST_1][j][i] = -1;
      }
      memset(&enc_picture->ref_idx[LIST_0][j][img->block_x + i0], -1, (input->part_size[pmode][0]) * sizeof(char));
      memset(&enc_picture->ref_idx[LIST_1][j][img->block_x + i0], -1, (input->part_size[pmode][0]) * sizeof(char));
      memset(enc_picture->mv[LIST_0][j][img->block_x + i0], 0, 2*(input->part_size[pmode][0]) * sizeof(short));
      memset(enc_picture->mv[LIST_1][j][img->block_x + i0], 0, 2*(input->part_size[pmode][0]) * sizeof(short));
    }
    return;
  }

  if (!bslice)
  {
    for (j=j0; j<j1; j++)
    {
      block_y = img->block_y + j;
      memset(&enc_picture->ref_idx   [LIST_0][block_y][img->block_x + i0], fwref, (input->part_size[pmode][0]) * sizeof(char));
      for (i=i0; i<i1; i++)
      {
        block_x = img->block_x + i;
        cur_mv = img->all_mv[j][i][LIST_0][fwref][mode];
        enc_picture->mv        [LIST_0][block_y][block_x][0] = cur_mv[0];
        enc_picture->mv        [LIST_0][block_y][block_x][1] = cur_mv[1];
        enc_picture->ref_pic_id[LIST_0][block_y][block_x] = enc_picture->ref_pic_num[LIST_0+currMB->list_offset][fwref];
      }
    }
    return;
  }
  else
  {
    for (j=j0; j<j1; j++)
    {
      block_y = img->block_y + j;
      for (i=i0; i<i1; i++)
      {
        block_x = img->block_x + i;
        if (mode==0)
        {
          pdir  = direct_pdir[block_y][block_x];
          fwref = direct_ref_idx[LIST_0][block_y][block_x];
          bwref = direct_ref_idx[LIST_1][block_y][block_x];
        }

        if ((pdir==0 || pdir==2))
        {
          if (currMB->bi_pred_me && (pdir == 2) && mode == 1)
          {
            cur_mv = currMB->bi_pred_me == 1
              ? img->bipred_mv1[j][i][LIST_0][0][mode]
              : img->bipred_mv2[j][i][LIST_0][0][mode];

            enc_picture->mv        [LIST_0][block_y][block_x][0] = cur_mv[0];
            enc_picture->mv        [LIST_0][block_y][block_x][1] = cur_mv[1];
            enc_picture->ref_idx   [LIST_0][block_y][block_x]    = 0;
            enc_picture->ref_pic_id[LIST_0][block_y][block_x]    = enc_picture->ref_pic_num[LIST_0+currMB->list_offset][0];
          }
          else
          {
            cur_mv = img->all_mv[j][i][LIST_0][fwref][mode];

            enc_picture->mv        [LIST_0][block_y][block_x][0] = cur_mv[0];
            enc_picture->mv        [LIST_0][block_y][block_x][1] = cur_mv[1];
            enc_picture->ref_idx   [LIST_0][block_y][block_x] = fwref;
            enc_picture->ref_pic_id[LIST_0][block_y][block_x] =
            enc_picture->ref_pic_num[LIST_0+currMB->list_offset][(short)enc_picture->ref_idx[LIST_0][block_y][block_x]];
          }
        }
        else
        {
          enc_picture->mv        [LIST_0][block_y][block_x][0] = 0;
          enc_picture->mv        [LIST_0][block_y][block_x][1] = 0;
          enc_picture->ref_idx   [LIST_0][block_y][block_x]    = -1;
          enc_picture->ref_pic_id[LIST_0][block_y][block_x]    = -1;
        }

        if ((pdir==1 || pdir==2))
        {
          if (currMB->bi_pred_me && (pdir == 2) && mode == 1)
          {
            cur_mv = currMB->bi_pred_me == 1
              ? img->bipred_mv1[j][i][LIST_1][0][mode]
              : img->bipred_mv2[j][i][LIST_1][0][mode];

            enc_picture->mv        [LIST_1][block_y][block_x][0] = cur_mv[0];
            enc_picture->mv        [LIST_1][block_y][block_x][1] = cur_mv[1];
            enc_picture->ref_idx   [LIST_1][block_y][block_x]    = 0;
            enc_picture->ref_pic_id[LIST_1][block_y][block_x]    = enc_picture->ref_pic_num[LIST_1+currMB->list_offset][0];
          }
          else
          {
            cur_mv = img->all_mv[j][i][LIST_1][bwref][mode];

            enc_picture->mv        [LIST_1][block_y][block_x][0] = cur_mv[0];
            enc_picture->mv        [LIST_1][block_y][block_x][1] = cur_mv[1];
            enc_picture->ref_idx   [LIST_1][block_y][block_x] = bwref;
            enc_picture->ref_pic_id[LIST_1][block_y][block_x] =
            enc_picture->ref_pic_num[LIST_1+currMB->list_offset][(short)enc_picture->ref_idx[LIST_1][block_y][block_x]];
          }
        }
        else
        {
          enc_picture->mv        [LIST_1][block_y][block_x][0] = 0;
          enc_picture->mv        [LIST_1][block_y][block_x][1] = 0;
          enc_picture->ref_idx   [LIST_1][block_y][block_x]    = -1;
          enc_picture->ref_pic_id[LIST_1][block_y][block_x]    = -1;
        }
      }
    }
  }
}

/*!
 *************************************************************************************
 * \brief
 *    skip macroblock field inference
 * \return
 *    inferred field flag
 *************************************************************************************
 */
int field_flag_inference(Macroblock *currMB)
{
  int mb_field;

  if (currMB->mbAvailA)
  {
    mb_field = img->mb_data[currMB->mbAddrA].mb_field;
  }
  else
  {
    // check top macroblock pair
    if (currMB->mbAvailB)
      mb_field = img->mb_data[currMB->mbAddrB].mb_field;
    else
      mb_field = 0;
  }

  return mb_field;
}

/*!
 *************************************************************************************
 * \brief
 *    Store motion vectors for 8x8 partition
 *************************************************************************************
 */

void StoreMVBlock8x8(int dir, int block8x8, int mode, int l0_ref, int l1_ref, int pdir8, int bframe)
{
  int i, j, i0, j0, ii, jj;
  short ******all_mv  = img->all_mv;
  short ******pred_mv = img->pred_mv;
  short (*lc_l0_mv8x8)[4][2] = all_mv8x8[dir][LIST_0];
  short (*lc_l1_mv8x8)[4][2] = all_mv8x8[dir][LIST_1];
  short (*lc_pr_mv8x8)[4][2] = NULL;

  i0 = (block8x8 & 0x01) << 1;
  j0 = (block8x8 >> 1) << 1;
  ii = i0+2;
  jj = j0+2;

  if (!bframe)
  {
    if (pdir8>=0) //(mode8!=IBLOCK)&&(mode8!=I16MB))  // && ref != -1)
    {
      lc_pr_mv8x8 = pred_mv8x8[dir][LIST_0];
      for (j=j0; j<jj; j++)
        for (i=i0; i<ii; i++)
        {
          lc_l0_mv8x8[j][i][0] = all_mv [j][i][LIST_0][l0_ref][4][0];
          lc_l0_mv8x8[j][i][1] = all_mv [j][i][LIST_0][l0_ref][4][1];
          lc_pr_mv8x8[j][i][0] = pred_mv[j][i][LIST_0][l0_ref][4][0];
          lc_pr_mv8x8[j][i][1] = pred_mv[j][i][LIST_0][l0_ref][4][1];
        }
    }
  }
  else
  {
    if (pdir8 == 0) // list0
    {
      lc_pr_mv8x8 = pred_mv8x8[dir][LIST_0];
      for (j=j0; j<jj; j++)
        for (i=i0; i<ii; i++)
        {
          lc_l0_mv8x8[j][i][0] = all_mv [j][i][LIST_0][l0_ref][mode][0];
          lc_l0_mv8x8[j][i][1] = all_mv [j][i][LIST_0][l0_ref][mode][1];
          lc_pr_mv8x8[j][i][0] = pred_mv[j][i][LIST_0][l0_ref][mode][0];
          lc_pr_mv8x8[j][i][1] = pred_mv[j][i][LIST_0][l0_ref][mode][1];
        }
    }
    else if (pdir8 == 1) // list1
    {
      lc_pr_mv8x8 = pred_mv8x8[dir][LIST_1];
      for (j=j0; j<jj; j++)
        for (i=i0; i<ii; i++)
        {
          lc_l1_mv8x8[j][i][0] = all_mv [j][i][LIST_1][l1_ref][mode][0];
          lc_l1_mv8x8[j][i][1] = all_mv [j][i][LIST_1][l1_ref][mode][1];
          lc_pr_mv8x8[j][i][0] = pred_mv[j][i][LIST_1][l1_ref][mode][0];
          lc_pr_mv8x8[j][i][1] = pred_mv[j][i][LIST_1][l1_ref][mode][1];
        }
    }
    else if (pdir8==2) // bipred
    {
      lc_pr_mv8x8 = pred_mv8x8[dir][LIST_0];
      for (j=j0; j<jj; j++)
      {
        for (i=i0; i<ii; i++)
        {
          lc_l0_mv8x8[j][i][0] = all_mv [j][i][LIST_0][l0_ref][mode][0];
          lc_l0_mv8x8[j][i][1] = all_mv [j][i][LIST_0][l0_ref][mode][1];
          lc_pr_mv8x8[j][i][0] = pred_mv[j][i][LIST_0][l0_ref][mode][0];
          lc_pr_mv8x8[j][i][1] = pred_mv[j][i][LIST_0][l0_ref][mode][1];
        }
      }
      lc_pr_mv8x8 = pred_mv8x8[dir][LIST_1];
      for (j=j0; j<jj; j++)
      {
        for (i=i0; i<ii; i++)
        {
          lc_l1_mv8x8[j][i][0] = all_mv [j][i][LIST_1][l1_ref][mode][0];
          lc_l1_mv8x8[j][i][1] = all_mv [j][i][LIST_1][l1_ref][mode][1];
          lc_pr_mv8x8[j][i][0] = pred_mv[j][i][LIST_1][l1_ref][mode][0];
          lc_pr_mv8x8[j][i][1] = pred_mv[j][i][LIST_1][l1_ref][mode][1];
        }
      }
    }
    else
    {
      error("invalid direction mode", 255);
    }
  }
}



/*!
 *************************************************************************************
 * \brief
 *    Store motion vectors of 8x8 partitions of one macroblock
 *************************************************************************************
 */
void StoreMV8x8(int dir)
{
  int block8x8;

  int bframe = (img->type == B_SLICE);

  for (block8x8=0; block8x8<4; block8x8++)
    StoreMVBlock8x8(dir, block8x8, tr8x8.part8x8mode[block8x8], tr8x8.part8x8l0ref[block8x8],
    tr8x8.part8x8l1ref[block8x8], tr8x8.part8x8pdir[block8x8], bframe);
}

/*!
*************************************************************************************
* \brief
*    Restore motion vectors for 8x8 partition
*************************************************************************************
*/
void RestoreMVBlock8x8(int dir, int block8x8, RD_8x8DATA tr, int bframe)
{
  int i, j, i0, j0, ii, jj;
  short ******all_mv  = img->all_mv;
  short ******pred_mv = img->pred_mv;
  short (*lc_l0_mv8x8)[4][2] = all_mv8x8[dir][LIST_0];
  short (*lc_l1_mv8x8)[4][2] = all_mv8x8[dir][LIST_1];
  short (*lc_pr_mv8x8)[4][2] = NULL;

  short pdir8  = tr.part8x8pdir [block8x8];
  short mode   = tr.part8x8mode [block8x8];
  short l0_ref = tr.part8x8l0ref[block8x8];
  short l1_ref = tr.part8x8l1ref[block8x8];

  i0 = (block8x8 & 0x01) << 1;
  j0 = (block8x8 >> 1) << 1;
  ii = i0+2;
  jj = j0+2;

  if (!bframe)
  {
    if (pdir8>=0) //(mode8!=IBLOCK)&&(mode8!=I16MB))  // && ref != -1)
    {
      lc_pr_mv8x8 = pred_mv8x8[dir][LIST_0];
      for (j=j0; j<jj; j++)
        for (i=i0; i<ii; i++)
        {
          all_mv [j][i][LIST_0][l0_ref][4][0] = lc_l0_mv8x8[j][i][0] ;
          all_mv [j][i][LIST_0][l0_ref][4][1] = lc_l0_mv8x8[j][i][1] ;
          pred_mv[j][i][LIST_0][l0_ref][4][0] = lc_pr_mv8x8[j][i][0];
          pred_mv[j][i][LIST_0][l0_ref][4][1] = lc_pr_mv8x8[j][i][1];
        }
    }
  }
  else
  {
    if (pdir8==0) // forward
    {
      lc_pr_mv8x8 = pred_mv8x8[dir][LIST_0];
      for (j=j0; j<jj; j++)
      {
        for (i=i0; i<ii; i++)
        {
          all_mv [j][i][LIST_0][l0_ref][mode][0] = lc_l0_mv8x8[j][i][0] ;
          all_mv [j][i][LIST_0][l0_ref][mode][1] = lc_l0_mv8x8[j][i][1] ;
          pred_mv[j][i][LIST_0][l0_ref][mode][0] = lc_pr_mv8x8[j][i][0];
          pred_mv[j][i][LIST_0][l0_ref][mode][1] = lc_pr_mv8x8[j][i][1];
        }
      }
    }
    else if (pdir8==1) // backward
    {
      lc_pr_mv8x8 = pred_mv8x8[dir][LIST_1];
      for (j=j0; j<jj; j++)
      {
        for (i=i0; i<ii; i++)
        {
          all_mv [j][i][LIST_1][l1_ref][mode][0] = lc_l1_mv8x8[j][i][0] ;
          all_mv [j][i][LIST_1][l1_ref][mode][1] = lc_l1_mv8x8[j][i][1] ;
          pred_mv[j][i][LIST_1][l1_ref][mode][0] = lc_pr_mv8x8[j][i][0];
          pred_mv[j][i][LIST_1][l1_ref][mode][1] = lc_pr_mv8x8[j][i][1];
        }
      }
    }
    else if (pdir8==2) // bidir
    {
      lc_pr_mv8x8 = pred_mv8x8[dir][LIST_0];
      for (j=j0; j<jj; j++)
      {
        for (i=i0; i<ii; i++)
        {
          all_mv [j][i][LIST_0][l0_ref][mode][0] = lc_l0_mv8x8[j][i][0] ;
          all_mv [j][i][LIST_0][l0_ref][mode][1] = lc_l0_mv8x8[j][i][1] ;
          pred_mv[j][i][LIST_0][l0_ref][mode][0] = lc_pr_mv8x8[j][i][0];
          pred_mv[j][i][LIST_0][l0_ref][mode][1] = lc_pr_mv8x8[j][i][1];
        }
      }
      lc_pr_mv8x8 = pred_mv8x8[dir][LIST_1];
      for (j=j0; j<jj; j++)
      {
        for (i=i0; i<ii; i++)
        {
          all_mv [j][i][LIST_1][l1_ref][mode][0] = lc_l1_mv8x8[j][i][0] ;
          all_mv [j][i][LIST_1][l1_ref][mode][1] = lc_l1_mv8x8[j][i][1] ;
          pred_mv[j][i][LIST_1][l1_ref][mode][0] = lc_pr_mv8x8[j][i][0];
          pred_mv[j][i][LIST_1][l1_ref][mode][1] = lc_pr_mv8x8[j][i][1];
        }
      }
    }
    else
    {
      error("invalid direction mode", 255);
    }
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Restore motion vectors of 8x8 partitions of one macroblock
 *************************************************************************************
 */
void RestoreMV8x8(int dir)
{
  int block8x8;

  int bframe = (img->type == B_SLICE);

  for (block8x8=0; block8x8<4; block8x8++)
    RestoreMVBlock8x8(dir, block8x8, tr8x8, bframe);
}


/*!
 *************************************************************************************
 * \brief
 *    Store predictors for 8x8 partition
 *************************************************************************************
 */

void StoreNewMotionVectorsBlock8x8(int dir, int block8x8, int mode, int l0_ref, int l1_ref, int pdir8, int bframe)
{
  int i, j, i0, j0, ii, jj;
  short ******all_mv  = img->all_mv;
  short ******pred_mv = img->pred_mv;
  short (*lc_l0_mv8x8)[4][2] = all_mv8x8[dir][LIST_0];
  short (*lc_l1_mv8x8)[4][2] = all_mv8x8[dir][LIST_1];
  short (*lc_pr_mv8x8)[4][2] = NULL;

  i0 = (block8x8 & 0x01) << 1;
  j0 = (block8x8 >> 1) << 1;
  ii = i0+2;
  jj = j0+2;

  if (pdir8<0)
  {
    for (j=j0; j<jj; j++)
    {
      memset(&lc_l0_mv8x8[j][i0], 0, 4 * sizeof(short));
      memset(&lc_l1_mv8x8[j][i0], 0, 4 * sizeof(short));
    }
    return;
  }

  if (!bframe)
  {

    lc_pr_mv8x8 = pred_mv8x8[dir][LIST_0];
    for (j=j0; j<jj; j++)
    {
      for (i=i0; i<ii; i++)
      {
        lc_l0_mv8x8[j][i][0] = all_mv [j][i][LIST_0][l0_ref][4][0];
        lc_l0_mv8x8[j][i][1] = all_mv [j][i][LIST_0][l0_ref][4][1];
        lc_pr_mv8x8[j][i][0] = pred_mv[j][i][LIST_0][l0_ref][4][0];
        lc_pr_mv8x8[j][i][1] = pred_mv[j][i][LIST_0][l0_ref][4][1];
      }
      memset(&lc_l1_mv8x8[j][i0], 0, 4 * sizeof(short));
    }
    return;
  }
  else
  {
    if ((pdir8==0 || pdir8==2))
    {
      lc_pr_mv8x8 = pred_mv8x8[dir][LIST_0];

      for (j=j0; j<jj; j++)
        for (i=i0; i<ii; i++)
        {
          lc_l0_mv8x8[j][i][0] = all_mv [j][i][LIST_0][l0_ref][mode][0];
          lc_l0_mv8x8[j][i][1] = all_mv [j][i][LIST_0][l0_ref][mode][1];
          lc_pr_mv8x8[j][i][0] = pred_mv[j][i][LIST_0][l0_ref][mode][0];
          lc_pr_mv8x8[j][i][1] = pred_mv[j][i][LIST_0][l0_ref][mode][1];
        }
    }
    else
    {
      for (j=j0; j<jj; j++)
        memset(&lc_l0_mv8x8[j][i0], 0, 4 * sizeof(short));
    }

    if ((pdir8==1 || pdir8==2))
    {
      lc_pr_mv8x8 = pred_mv8x8[dir][LIST_1];

      for (j=j0; j<jj; j++)
        for (i=i0; i<ii; i++)
        {
          lc_l1_mv8x8[j][i][0] = all_mv [j][i][LIST_1][l1_ref][mode][0];
          lc_l1_mv8x8[j][i][1] = all_mv [j][i][LIST_1][l1_ref][mode][1];
          lc_pr_mv8x8[j][i][0] = pred_mv[j][i][LIST_1][l1_ref][mode][0];
          lc_pr_mv8x8[j][i][1] = pred_mv[j][i][LIST_1][l1_ref][mode][1];
        }
    }
    else
    {
      for (j=j0; j<jj; j++)
        memset(&lc_l1_mv8x8[j][i0], 0, 4 * sizeof(short));
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    Makes the decision if 8x8 tranform will be used (for RD-off)
 ************************************************************************
 */
int GetBestTransformP8x8()
{
  int    block_y, block_x, pic_pix_y, pic_pix_x, i, j, k;
  int    mb_y, mb_x, block8x8;
  int    cost8x8=0, cost4x4=0;
  int    *diff_ptr;

  if(input->Transform8x8Mode==2) //always allow 8x8 transform
    return 1;

  for (block8x8=0; block8x8<4; block8x8++)
  {
    mb_y = (block8x8 >>   1) << 3;
    mb_x = (block8x8 & 0x01) << 3;
    //===== loop over 4x4 blocks =====
    k=0;
    for (block_y=mb_y; block_y<mb_y+8; block_y+=4)
    {
      pic_pix_y = img->opix_y + block_y;

      //get cost for transform size 4x4
      for (block_x=mb_x; block_x<mb_x+8; block_x+=4)
      {
        pic_pix_x = img->opix_x + block_x;

        //===== get displaced frame difference ======
        diff_ptr=&diff4x4[k];
        for (j=0; j<4; j++)
        {
          for (i=0; i<4; i++, k++)
          {
            //4x4 transform size
            diff4x4[k] = pCurImg[pic_pix_y+j][pic_pix_x+i] - tr4x4.mpr8x8[j+block_y][i+block_x];
            //8x8 transform size
            diff8x8[k] = pCurImg[pic_pix_y+j][pic_pix_x+i] - tr8x8.mpr8x8[j+block_y][i+block_x];
          }
        }

        cost4x4 += distortion4x4 (diff_ptr);
      }
    }
    cost8x8 += distortion8x8 (diff8x8);
  }
  return (cost8x8 < cost4x4);
}

/*!
************************************************************************
* \brief
*    Sets MBAFF RD parameters
************************************************************************
*/
void set_mbaff_parameters(Macroblock  *currMB)
{
  int  i, j, k;
  int  mode         = best_mode;
  int  bframe       = (img->type==B_SLICE);
  char **ipredmodes = img->ipredmode;


  //===== reconstruction values =====
  for (j=0; j < MB_BLOCK_SIZE; j++)
    memcpy(rdopt->rec_mbY[j],&enc_picture->imgY[img->pix_y + j][img->pix_x], MB_BLOCK_SIZE * sizeof(imgpel));

  if (img->yuv_format != YUV400)
  {
    for (j=0; j<img->mb_cr_size_y; j++)
    {
      memcpy(rdopt->rec_mbU[j],&enc_picture->imgUV[0][img->pix_c_y + j][img->pix_c_x], img->mb_cr_size_x * sizeof(imgpel));
      memcpy(rdopt->rec_mbV[j],&enc_picture->imgUV[1][img->pix_c_y + j][img->pix_c_x], img->mb_cr_size_x * sizeof(imgpel));
    }
  }

  //===== coefficients and cbp =====
  rdopt->mode      = mode;
  rdopt->i16offset = img->i16offset;
  rdopt->cbp       = currMB->cbp;
  rdopt->cbp_blk   = currMB->cbp_blk;
  rdopt->mb_type   = currMB->mb_type;

  rdopt->luma_transform_size_8x8_flag = currMB->luma_transform_size_8x8_flag;

  if(rdopt->mb_type == 0 && mode != 0)
  {
    mode=0;
    rdopt->mode=0;
  }

  for(i=0;i<4+img->num_blk8x8_uv;i++)
  {
    for(j=0;j<4;j++)
      for(k=0;k<2;k++)
        memcpy(rdopt->cofAC[i][j][k], img->cofAC[i][j][k], 65 * sizeof(int));
  }

  for(i=0;i<3;i++)
  {
    for(k=0;k<2;k++)
      memcpy(rdopt->cofDC[i][k], img->cofDC[i][k], 18 * sizeof(int));
  }

  memcpy(rdopt->b8mode,currMB->b8mode, BLOCK_MULTIPLE * sizeof(int));
  memcpy(rdopt->b8pdir,currMB->b8pdir, BLOCK_MULTIPLE * sizeof(int));

  //==== reference frames =====
  if (bframe)
  {
    for (j = 0; j < BLOCK_MULTIPLE; j++)
    {
      memcpy(rdopt->refar[LIST_0][j],&enc_picture->ref_idx[LIST_0][img->block_y + j][img->block_x] , BLOCK_MULTIPLE * sizeof(char));
      memcpy(rdopt->refar[LIST_1][j],&enc_picture->ref_idx[LIST_1][img->block_y + j][img->block_x] , BLOCK_MULTIPLE * sizeof(char));
    }
    rdopt->bi_pred_me = currMB->bi_pred_me;
  }
  else
  {
    for (j = 0; j < BLOCK_MULTIPLE; j++)
      memcpy(rdopt->refar[LIST_0][j],&enc_picture->ref_idx[LIST_0][img->block_y + j][img->block_x] , BLOCK_MULTIPLE * sizeof(char));
  }

  memcpy(rdopt->intra_pred_modes,currMB->intra_pred_modes, MB_BLOCK_PARTITIONS * sizeof(char));
  memcpy(rdopt->intra_pred_modes8x8,currMB->intra_pred_modes8x8, MB_BLOCK_PARTITIONS * sizeof(char));
  for (j = img->block_y; j < img->block_y + 4; j++)
  {
    memcpy(&rdopt->ipredmode[j][img->block_x],&ipredmodes[j][img->block_x], BLOCK_MULTIPLE * sizeof(char));
  }
}

/*!
************************************************************************
* \brief
*    store coding state (for rd-optimized mode decision), used for 8x8 transformation
************************************************************************
*/
void store_coding_state_cs_cm(Macroblock *currMB)
{
  store_coding_state(currMB, cs_cm);
}

/*!
************************************************************************
* \brief
*    restore coding state (for rd-optimized mode decision), used for 8x8 transformation
************************************************************************
*/
void reset_coding_state_cs_cm(Macroblock *currMB)
{
  reset_coding_state(currMB, cs_cm);
}

/*!
************************************************************************
* \brief
*    update rounding offsets based on JVT-N011
************************************************************************
*/
void update_offset_params(Macroblock *currMB, int mode, int luma_transform_size_8x8_flag)
{
  int is_inter = (mode != I4MB)&&(mode != I16MB) && (mode != I8MB);
  int luma_pos = AdaptRndPos[(is_inter<<1) + luma_transform_size_8x8_flag][img->type];
  int i,j;
  int qp = currMB->qp + img->bitdepth_luma_qp_scale - MIN_QP;
  int cur_qp = input->AdaptRoundingFixed ? 0 : qp;
  int temp = 0;
  int offsetRange = 1 << (OffsetBits - 1);
  int blk_mask = 0x03 + (luma_transform_size_8x8_flag<<2);
  int blk_shift = 2 + luma_transform_size_8x8_flag;
  short **offsetList = luma_transform_size_8x8_flag ? OffsetList8x8[cur_qp] : OffsetList4x4[cur_qp];

  int **fAdjust = is_inter
    ? (luma_transform_size_8x8_flag ? bestInterFAdjust8x8 : bestInterFAdjust4x4)
    : (luma_transform_size_8x8_flag ? bestIntraFAdjust8x8 : bestIntraFAdjust4x4);

  if( (active_sps->chroma_format_idc == YUV444)&&IS_INDEPENDENT(input) )
  {
    if( luma_transform_size_8x8_flag )  // 8x8
      luma_pos += 5 * img->colour_plane_id;
    else  // 4x4
      luma_pos += img->colour_plane_id;
  }

  for (j=0; j < MB_BLOCK_SIZE; j++)
  {
    int j_pos = ((j & blk_mask)<<blk_shift);
    for (i=0; i < MB_BLOCK_SIZE; i++)
    {
      temp = j_pos + (i & blk_mask);
      offsetList[luma_pos][temp] += fAdjust[j][i];
      offsetList[luma_pos][temp]  = iClip3(0, offsetRange, offsetList[luma_pos][temp]);
    }
  }
  
  if((active_sps->chroma_format_idc == YUV444)&& (IS_INDEPENDENT(input)==0 ))
  { 
    int ***fAdjustCbCr = (int ***) (is_inter
      ? (luma_transform_size_8x8_flag ? bestInterFAdjust8x8Cr : bestInterFAdjust4x4Cr)
      : (luma_transform_size_8x8_flag ? bestIntraFAdjust8x8Cr : bestIntraFAdjust4x4Cr));
    int uv;
    
    for(uv=0; uv<2; uv++)
    {
      luma_pos = AdaptRndPos[(is_inter<<1) + luma_transform_size_8x8_flag][img->type];
      if(luma_transform_size_8x8_flag )  // 8x8
        luma_pos += 5 * (uv+1);
      else  // 4x4
        luma_pos += (uv+1);
      for (j=0; j < MB_BLOCK_SIZE; j++)
      {
        int j_pos = ((j & blk_mask)<<blk_shift);
        for (i=0; i < MB_BLOCK_SIZE; i++)
        {
          temp = j_pos + (i & blk_mask);
          offsetList[luma_pos][temp] += fAdjustCbCr[uv][j][i];
          offsetList[luma_pos][temp] = iClip3(0,offsetRange,offsetList[luma_pos][temp]);
        }
      }
    }
  }
  
  if ((input->yuv_format == YUV420 || input->yuv_format == YUV422 )&&(input->AdaptRndChroma))
  {
    int u_pos = AdaptRndCrPos[is_inter][img->type];
    int v_pos = u_pos + 1;
    int jpos;

    int ***fAdjustCr = is_inter ? bestInterFAdjust4x4Cr : bestIntraFAdjust4x4Cr;
    
    for (j = 0; j < img->mb_cr_size_y; j++)
    {
      jpos = ((j & 0x03)<<2);
      for (i = 0; i < img->mb_cr_size_x; i++)
      {
        temp = jpos + (i & 0x03);
        OffsetList4x4[cur_qp][u_pos][temp] += fAdjustCr[0][j][i];
        OffsetList4x4[cur_qp][u_pos][temp] = iClip3(0,offsetRange,OffsetList4x4[cur_qp][u_pos][temp]);
        OffsetList4x4[cur_qp][v_pos][temp] += fAdjustCr[1][j][i];
        OffsetList4x4[cur_qp][v_pos][temp] = iClip3(0,offsetRange,OffsetList4x4[cur_qp][v_pos][temp]);
      }
    }
  }
}

void assign_enc_picture_params(int mode, char best_pdir, int block, int list_offset, int best_l0_ref, int best_l1_ref, int bframe)
{
  int i,j;
  int block_x, block_y;
  short *cur_mv;

  if (mode==1)
  {
    if (best_pdir==1)
    {
      for (j=img->block_y+(block&2); j<img->block_y+(block&2) + BLOCK_MULTIPLE; j++)
      {
        block_x = img->block_x+(block&1)*2;

        memset(&enc_picture->ref_idx[LIST_0][j][block_x], -1 ,     BLOCK_MULTIPLE * sizeof(char));
        memset(enc_picture->mv      [LIST_0][j][block_x],  0 , 2 * BLOCK_MULTIPLE * sizeof(short));
        for (i=block_x; i<block_x + BLOCK_MULTIPLE; i++)
        {
          enc_picture->ref_pic_id [LIST_0][j][i]    = -1;
        }
      }
    }
    else if (img->bi_pred_me[mode])
    {
      for (j=0; j<BLOCK_MULTIPLE; j++)
      {
        block_y = img->block_y+(block&2)+j;
        block_x = img->block_x+(block&1)*2;
        memset(&enc_picture->ref_idx[LIST_0][block_y][block_x], 0, BLOCK_MULTIPLE * sizeof(char));
        for (i=0; i<BLOCK_MULTIPLE; i++)
        {
          cur_mv = img->bi_pred_me[mode] == 1
            ? img->bipred_mv1[i][j][LIST_0][0][mode]
            : img->bipred_mv2[i][j][LIST_0][0][mode];

          enc_picture->ref_pic_id [LIST_0][block_y][block_x + i]    = enc_picture->ref_pic_num[LIST_0 + list_offset][0];
          enc_picture->mv         [LIST_0][block_y][block_x + i][0] = cur_mv[0];
          enc_picture->mv         [LIST_0][block_y][block_x + i][1] = cur_mv[1];
        }
      }
    }
    else
    {
      for (j=0; j<BLOCK_MULTIPLE; j++)
      {
        block_y = img->block_y+(block&2)+j;
        block_x = img->block_x+(block&1)*2;
        memset(&enc_picture->ref_idx[LIST_0][block_y][block_x], best_l0_ref , BLOCK_MULTIPLE * sizeof(char));
        for (i=0; i<BLOCK_MULTIPLE; i++)
        {
          cur_mv = img->all_mv[j][i][LIST_0][best_l0_ref][mode];

          enc_picture->ref_pic_id [LIST_0][block_y][block_x + i]    = enc_picture->ref_pic_num[LIST_0 + list_offset][best_l0_ref];
          enc_picture->mv         [LIST_0][block_y][block_x + i][0] = cur_mv[0];
          enc_picture->mv         [LIST_0][block_y][block_x + i][1] = cur_mv[1];
        }
      }
    }

    if (bframe)
    {
      if (best_pdir==0)
      {
        for (j=img->block_y+(block&2); j<img->block_y+(block&2) + BLOCK_MULTIPLE; j++)
        {
          block_x = img->block_x+(block&1)*2;
          memset(&enc_picture->ref_idx[LIST_1][j][block_x], -1 , BLOCK_MULTIPLE * sizeof(char));
          memset(enc_picture->mv[LIST_1][j][block_x], 0 , 2 * BLOCK_MULTIPLE * sizeof(short));
          for (i=block_x; i<block_x + BLOCK_MULTIPLE; i++)
          {
            enc_picture->ref_pic_id [LIST_1][j][i] = -1;
          }
        }
      }
      else
      {
        if (img->bi_pred_me[mode])
        {
          for (j=0; j<BLOCK_MULTIPLE; j++)
          {
            block_y = img->block_y+(block&2)+j;
            block_x = img->block_x+(block&1)*2;
            memset(&enc_picture->ref_idx[LIST_1][block_y][block_x], 0, BLOCK_MULTIPLE * sizeof(char));
            for (i=0; i<BLOCK_MULTIPLE; i++)
            {
              cur_mv = img->bi_pred_me[mode] == 1
                ? img->bipred_mv1[i][j][LIST_1][0][mode]
                : img->bipred_mv2[i][j][LIST_1][0][mode];

              enc_picture->ref_pic_id [LIST_1][block_y][block_x + i] =
                enc_picture->ref_pic_num[LIST_1 + list_offset][0];
              enc_picture->mv         [LIST_1][block_y][block_x + i][0] = cur_mv[0];
              enc_picture->mv         [LIST_1][block_y][block_x + i][1] = cur_mv[1];
            }
          }
        }
        else
        {
          for (j=0; j<BLOCK_MULTIPLE; j++)
          {
            block_y = img->block_y+(block&2)+j;
            block_x = img->block_x+(block&1)*2;
            memset(&enc_picture->ref_idx[LIST_1][block_y][block_x], best_l1_ref, BLOCK_MULTIPLE * sizeof(char));
            for (i=0; i<BLOCK_MULTIPLE; i++)
            {

              enc_picture->ref_pic_id [LIST_1][block_y][block_x + i] =
                enc_picture->ref_pic_num[LIST_1 + list_offset][best_l1_ref];
              if(best_l1_ref>=0)
              {
                cur_mv = img->all_mv[j][i][LIST_1][best_l1_ref][mode];
                enc_picture->mv[LIST_1][block_y][block_x + i][0] = cur_mv[0];
                enc_picture->mv[LIST_1][block_y][block_x + i][1] = cur_mv[1];
              }
            }
          }
        }
      }
    }
  }
  else if (mode==2)
  {
    for (j=0; j<2; j++)
    {
      block_y = img->block_y + block * 2 + j;
      for (i=0; i<BLOCK_MULTIPLE; i++)
      {
        block_x = img->block_x + i;
        if (best_pdir==1)
        {
          enc_picture->ref_idx    [LIST_0][block_y][block_x]    = -1;
          enc_picture->ref_pic_id [LIST_0][block_y][block_x]    = -1;
          enc_picture->mv         [LIST_0][block_y][block_x][0] = 0;
          enc_picture->mv         [LIST_0][block_y][block_x][1] = 0;
        }
        else
        {
          cur_mv = img->all_mv[j+block*2][i][LIST_0][best_l0_ref][mode];

          enc_picture->ref_idx    [LIST_0][block_y][block_x]    = best_l0_ref;
          enc_picture->ref_pic_id [LIST_0][block_y][block_x]    =
            enc_picture->ref_pic_num[LIST_0 + list_offset][best_l0_ref];
          enc_picture->mv         [LIST_0][block_y][block_x][0] = cur_mv[0];
          enc_picture->mv         [LIST_0][block_y][block_x][1] = cur_mv[1];
        }

        if (bframe)
        {
          if (best_pdir==0)
          {
            enc_picture->ref_idx    [LIST_1][block_y][block_x]    = -1;
            enc_picture->ref_pic_id [LIST_1][block_y][block_x]    = -1;
            enc_picture->mv         [LIST_1][block_y][block_x][0] = 0;
            enc_picture->mv         [LIST_1][block_y][block_x][1] = 0;
          }
          else
          {
            enc_picture->ref_idx[LIST_1][block_y][block_x] = best_l1_ref;
            if(best_l1_ref>=0)
            {
              cur_mv = img->all_mv[j+ block*2][i][LIST_1][best_l1_ref][mode];

              enc_picture->ref_pic_id [LIST_1][block_y][block_x] =
                enc_picture->ref_pic_num[LIST_1 + list_offset][best_l1_ref];
              enc_picture->mv[LIST_1][block_y][block_x][0] = cur_mv[0];
              enc_picture->mv[LIST_1][block_y][block_x][1] = cur_mv[1];
            }
          }
        }
      }
    }
  }
  else
  {
    for (j=0; j<BLOCK_MULTIPLE; j++)
    {
      block_y = img->block_y+j;
      for (i=0; i<2; i++)
      {
        block_x = img->block_x + block*2 + i;
        if (best_pdir==1)
        {
          enc_picture->ref_idx    [LIST_0][block_y][block_x]    = -1;
          enc_picture->ref_pic_id [LIST_0][block_y][block_x]    = -1;
          enc_picture->mv         [LIST_0][block_y][block_x][0] = 0;
          enc_picture->mv         [LIST_0][block_y][block_x][1] = 0;
        }
        else
        {
          cur_mv = img->all_mv[j][block*2+i][LIST_0][best_l0_ref][mode];

          enc_picture->ref_idx    [LIST_0][block_y][block_x] = best_l0_ref;
          enc_picture->ref_pic_id [LIST_0][block_y][block_x] =
            enc_picture->ref_pic_num[LIST_0 + list_offset][best_l0_ref];
          enc_picture->mv[LIST_0][block_y][block_x][0] = cur_mv[0];
          enc_picture->mv[LIST_0][block_y][block_x][1] = cur_mv[1];
        }

        if (bframe)
        {
          if (best_pdir==0)
          {
            enc_picture->ref_idx    [LIST_1][block_y][block_x]    = -1;
            enc_picture->ref_pic_id [LIST_1][block_y][block_x]    = -1;
            enc_picture->mv         [LIST_1][block_y][block_x][0] = 0;
            enc_picture->mv         [LIST_1][block_y][block_x][1] = 0;
          }
          else
          {
            enc_picture->ref_idx[LIST_1][block_y][block_x] = best_l1_ref;
            if(best_l1_ref>=0)
            {
              cur_mv = img->all_mv[j][block*2+i][LIST_1][best_l1_ref][mode];
              enc_picture->ref_pic_id [LIST_1][block_y][block_x] =
                enc_picture->ref_pic_num[LIST_1 + list_offset][best_l1_ref];

              enc_picture->mv[LIST_1][block_y][block_x][0] = cur_mv[0];
              enc_picture->mv[LIST_1][block_y][block_x][1] = cur_mv[1];
            }
          }
        }
      }
    }
  }
}

void update_refresh_map(int intra, int intra1, Macroblock *currMB)
{
  if (input->RestrictRef==1)
  {
    // Modified for Fast Mode Decision. Inchoon Choi, SungKyunKwan Univ.
    if (input->rdopt<2)
    {
      refresh_map[2*img->mb_y  ][2*img->mb_x  ] = (intra ? 1 : 0);
      refresh_map[2*img->mb_y  ][2*img->mb_x+1] = (intra ? 1 : 0);
      refresh_map[2*img->mb_y+1][2*img->mb_x  ] = (intra ? 1 : 0);
      refresh_map[2*img->mb_y+1][2*img->mb_x+1] = (intra ? 1 : 0);
    }
    else if (input->rdopt==3)
    {
      refresh_map[2*img->mb_y  ][2*img->mb_x  ] = (intra1==0 && (currMB->mb_type==I16MB || currMB->mb_type==I4MB) ? 1 : 0);
      refresh_map[2*img->mb_y  ][2*img->mb_x+1] = (intra1==0 && (currMB->mb_type==I16MB || currMB->mb_type==I4MB) ? 1 : 0);
      refresh_map[2*img->mb_y+1][2*img->mb_x  ] = (intra1==0 && (currMB->mb_type==I16MB || currMB->mb_type==I4MB) ? 1 : 0);
      refresh_map[2*img->mb_y+1][2*img->mb_x+1] = (intra1==0 && (currMB->mb_type==I16MB || currMB->mb_type==I4MB) ? 1 : 0);
    }
  }
  else if (input->RestrictRef==2)
  {
    refresh_map[2*img->mb_y  ][2*img->mb_x  ] = (currMB->mb_type==I16MB || currMB->mb_type==I4MB ? 1 : 0);
    refresh_map[2*img->mb_y  ][2*img->mb_x+1] = (currMB->mb_type==I16MB || currMB->mb_type==I4MB ? 1 : 0);
    refresh_map[2*img->mb_y+1][2*img->mb_x  ] = (currMB->mb_type==I16MB || currMB->mb_type==I4MB ? 1 : 0);
    refresh_map[2*img->mb_y+1][2*img->mb_x+1] = (currMB->mb_type==I16MB || currMB->mb_type==I4MB ? 1 : 0);
  }
}

