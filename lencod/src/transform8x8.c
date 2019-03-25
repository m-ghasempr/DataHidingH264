/*!
 ***************************************************************************
 * \file transform8x8.c
 *
 * \brief
 *    8x8 transform functions
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Yuri Vatis
 *    - Jan Muenster
 *    - Lowell Winger                   <lwinger@lsil.com>
 * \date
 *    12. October 2003
 **************************************************************************
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "global.h"

#include "image.h"
#include "mb_access.h"
#include "elements.h"
#include "cabac.h"
#include "vlc.h"
#include "transform8x8.h"
#include "transform.h"
#include "macroblock.h"
#include "symbol.h"

int   cofAC8x8_chroma[2][4][2][18];
static int diff64[64];

const int quant_coef8[6][8][8] =
{
  {
    {13107, 12222,  16777,  12222,  13107,  12222,  16777,  12222},
    {12222, 11428,  15481,  11428,  12222,  11428,  15481,  11428},
    {16777, 15481,  20972,  15481,  16777,  15481,  20972,  15481},
    {12222, 11428,  15481,  11428,  12222,  11428,  15481,  11428},
    {13107, 12222,  16777,  12222,  13107,  12222,  16777,  12222},
    {12222, 11428,  15481,  11428,  12222,  11428,  15481,  11428},
    {16777, 15481,  20972,  15481,  16777,  15481,  20972,  15481},
    {12222, 11428,  15481,  11428,  12222,  11428,  15481,  11428}
  },
  {
    {11916, 11058,  14980,  11058,  11916,  11058,  14980,  11058},
    {11058, 10826,  14290,  10826,  11058,  10826,  14290,  10826},
    {14980, 14290,  19174,  14290,  14980,  14290,  19174,  14290},
    {11058, 10826,  14290,  10826,  11058,  10826,  14290,  10826},
    {11916, 11058,  14980,  11058,  11916,  11058,  14980,  11058},
    {11058, 10826,  14290,  10826,  11058,  10826,  14290,  10826},
    {14980, 14290,  19174,  14290,  14980,  14290,  19174,  14290},
    {11058, 10826,  14290,  10826,  11058,  10826,  14290,  10826}
  },
  {
    {10082, 9675,   12710,  9675,   10082,  9675, 12710,  9675},
    {9675,  8943,   11985,  8943,   9675,   8943, 11985,  8943},
    {12710, 11985,  15978,  11985,  12710,  11985,  15978,  11985},
    {9675,  8943,   11985,  8943,   9675,   8943, 11985,  8943},
    {10082, 9675,   12710,  9675,   10082,  9675, 12710,  9675},
    {9675,  8943,   11985,  8943,   9675, 8943, 11985,  8943},
    {12710, 11985,  15978,  11985,  12710,  11985,  15978,  11985},
    {9675,  8943,   11985,  8943,   9675, 8943, 11985,  8943}
  },
  {
    {9362,  8931, 11984,  8931, 9362, 8931, 11984,  8931},
    {8931,  8228, 11259,  8228, 8931, 8228, 11259,  8228},
    {11984, 11259,  14913,  11259,  11984,  11259,  14913,  11259},
    {8931,  8228, 11259,  8228, 8931, 8228, 11259,  8228},
    {9362,  8931, 11984,  8931, 9362, 8931, 11984,  8931},
    {8931,  8228, 11259,  8228, 8931, 8228, 11259,  8228},
    {11984, 11259,  14913,  11259,  11984,  11259,  14913,  11259},
    {8931,  8228, 11259,  8228, 8931, 8228, 11259,  8228}
  },
  {
    {8192,  7740, 10486,  7740, 8192, 7740, 10486,  7740},
    {7740,  7346, 9777, 7346, 7740, 7346, 9777, 7346},
    {10486, 9777, 13159,  9777, 10486,  9777, 13159,  9777},
    {7740,  7346, 9777, 7346, 7740, 7346, 9777, 7346},
    {8192,  7740, 10486,  7740, 8192, 7740, 10486,  7740},
    {7740,  7346, 9777, 7346, 7740, 7346, 9777, 7346},
    {10486, 9777, 13159,  9777, 10486,  9777, 13159,  9777},
    {7740,  7346, 9777, 7346, 7740, 7346, 9777, 7346}
  },
  {
    {7282,  6830, 9118, 6830, 7282, 6830, 9118, 6830},
    {6830,  6428, 8640, 6428, 6830, 6428, 8640, 6428},
    {9118,  8640, 11570,  8640, 9118, 8640, 11570,  8640},
    {6830,  6428, 8640, 6428, 6830, 6428, 8640, 6428},
    {7282,  6830, 9118, 6830, 7282, 6830, 9118, 6830},
    {6830,  6428, 8640, 6428, 6830, 6428, 8640, 6428},
    {9118,  8640, 11570,  8640, 9118, 8640, 11570,  8640},
    {6830,  6428, 8640, 6428, 6830, 6428, 8640, 6428}
  }
};


const int dequant_coef8[6][8][8] =
{
  {
    {20,  19, 25, 19, 20, 19, 25, 19},
    {19,  18, 24, 18, 19, 18, 24, 18},
    {25,  24, 32, 24, 25, 24, 32, 24},
    {19,  18, 24, 18, 19, 18, 24, 18},
    {20,  19, 25, 19, 20, 19, 25, 19},
    {19,  18, 24, 18, 19, 18, 24, 18},
    {25,  24, 32, 24, 25, 24, 32, 24},
    {19,  18, 24, 18, 19, 18, 24, 18}
  },
  {
    {22,  21, 28, 21, 22, 21, 28, 21},
    {21,  19, 26, 19, 21, 19, 26, 19},
    {28,  26, 35, 26, 28, 26, 35, 26},
    {21,  19, 26, 19, 21, 19, 26, 19},
    {22,  21, 28, 21, 22, 21, 28, 21},
    {21,  19, 26, 19, 21, 19, 26, 19},
    {28,  26, 35, 26, 28, 26, 35, 26},
    {21,  19, 26, 19, 21, 19, 26, 19}
  },
  {
    {26,  24, 33, 24, 26, 24, 33, 24},
    {24,  23, 31, 23, 24, 23, 31, 23},
    {33,  31, 42, 31, 33, 31, 42, 31},
    {24,  23, 31, 23, 24, 23, 31, 23},
    {26,  24, 33, 24, 26, 24, 33, 24},
    {24,  23, 31, 23, 24, 23, 31, 23},
    {33,  31, 42, 31, 33, 31, 42, 31},
    {24,  23, 31, 23, 24, 23, 31, 23}
  },
  {
    {28,  26, 35, 26, 28, 26, 35, 26},
    {26,  25, 33, 25, 26, 25, 33, 25},
    {35,  33, 45, 33, 35, 33, 45, 33},
    {26,  25, 33, 25, 26, 25, 33, 25},
    {28,  26, 35, 26, 28, 26, 35, 26},
    {26,  25, 33, 25, 26, 25, 33, 25},
    {35,  33, 45, 33, 35, 33, 45, 33},
    {26,  25, 33, 25, 26, 25, 33, 25}
  },
  {
    {32,  30, 40, 30, 32, 30, 40, 30},
    {30,  28, 38, 28, 30, 28, 38, 28},
    {40,  38, 51, 38, 40, 38, 51, 38},
    {30,  28, 38, 28, 30, 28, 38, 28},
    {32,  30, 40, 30, 32, 30, 40, 30},
    {30,  28, 38, 28, 30, 28, 38, 28},
    {40,  38, 51, 38, 40, 38, 51, 38},
    {30,  28, 38, 28, 30, 28, 38, 28}
  },
  {
    {36,  34, 46, 34, 36, 34, 46, 34},
    {34,  32, 43, 32, 34, 32, 43, 32},
    {46,  43, 58, 43, 46, 43, 58, 43},
    {34,  32, 43, 32, 34, 32, 43, 32},
    {36,  34, 46, 34, 36, 34, 46, 34},
    {34,  32, 43, 32, 34, 32, 43, 32},
    {46,  43, 58, 43, 46, 43, 58, 43},
    {34,  32, 43, 32, 34, 32, 43, 32}
  }

};


//! single scan pattern
static const byte SNGL_SCAN8x8[64][2] = {
  {0,0}, {1,0}, {0,1}, {0,2}, {1,1}, {2,0}, {3,0}, {2,1},
  {1,2}, {0,3}, {0,4}, {1,3}, {2,2}, {3,1}, {4,0}, {5,0},
  {4,1}, {3,2}, {2,3}, {1,4}, {0,5}, {0,6}, {1,5}, {2,4},
  {3,3}, {4,2}, {5,1}, {6,0}, {7,0}, {6,1}, {5,2}, {4,3},
  {3,4}, {2,5}, {1,6}, {0,7}, {1,7}, {2,6}, {3,5}, {4,4},
  {5,3}, {6,2}, {7,1}, {7,2}, {6,3}, {5,4}, {4,5}, {3,6},
  {2,7}, {3,7}, {4,6}, {5,5}, {6,4}, {7,3}, {7,4}, {6,5},
  {5,6}, {4,7}, {5,7}, {6,6}, {7,5}, {7,6}, {6,7}, {7,7}
};


//! field scan pattern
static const byte FIELD_SCAN8x8[64][2] = {   // 8x8
  {0,0}, {0,1}, {0,2}, {1,0}, {1,1}, {0,3}, {0,4}, {1,2},
  {2,0}, {1,3}, {0,5}, {0,6}, {0,7}, {1,4}, {2,1}, {3,0},
  {2,2}, {1,5}, {1,6}, {1,7}, {2,3}, {3,1}, {4,0}, {3,2},
  {2,4}, {2,5}, {2,6}, {2,7}, {3,3}, {4,1}, {5,0}, {4,2},
  {3,4}, {3,5}, {3,6}, {3,7}, {4,3}, {5,1}, {6,0}, {5,2},
  {4,4}, {4,5}, {4,6}, {4,7}, {5,3}, {6,1}, {6,2}, {5,4},
  {5,5}, {5,6}, {5,7}, {6,3}, {7,0}, {7,1}, {6,4}, {6,5},
  {6,6}, {6,7}, {7,2}, {7,3}, {7,4}, {7,5}, {7,6}, {7,7}
};


//! array used to find expensive coefficients
static const byte COEFF_COST8x8[2][64] =
{
  {3,3,3,3,2,2,2,2,2,2,2,2,1,1,1,1,
  1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
  {9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
   9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
   9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
   9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9}
};

/*!
 *************************************************************************************
 * \brief
 *    8x8 Intra mode decision for a macroblock
 *************************************************************************************
 */

int Mode_Decision_for_new_Intra8x8Macroblock (Macroblock *currMB, double lambda, double *min_cost)
{
  int cur_cbp = 0, b8;
  double cost8x8;
  int cr_cbp[3] = { 0, 0, 0}; 

  *min_cost = (int)floor(6.0 * lambda + 0.4999);

  cmp_cbp[1] = cmp_cbp[2] = 0;
  for (b8=0; b8<4; b8++)
  {
    if (Mode_Decision_for_new_8x8IntraBlocks (currMB, b8, lambda, &cost8x8, cr_cbp))
    {
      cur_cbp |= (1<<b8);
    }
    *min_cost += cost8x8;

    if(active_sps->chroma_format_idc==YUV444 && !IS_INDEPENDENT(input))
    {
      int k;
      for (k = 1; k < 3; k++)
      {
        if (cr_cbp[k])
        {
          cmp_cbp[k] |= (1<<b8);
          cur_cbp |= cmp_cbp[k];
          cmp_cbp[k] = cur_cbp;
        }
      }
    }
  }

  return cur_cbp;
}

/*!
 *************************************************************************************
 * \brief
 *    8x8 Intra mode decision for a macroblock
 *************************************************************************************
 */

int Mode_Decision_for_new_8x8IntraBlocks (Macroblock *currMB, int b8, double lambda, double *min_cost, int cr_cbp[3])
{
  int     ipmode, best_ipmode = 0, i, j, k, y, dummy;
  double  cost;
  int     c_nz, nonzero = 0;
  static  imgpel  rec8x8[3][8][8];
  double  rdcost = 0.0;
  int     block_x     = (b8 & 0x01) << 3;
  int     block_y     = (b8 >> 1) << 3;
  int     pic_pix_x   = img->pix_x + block_x;
  int     pic_pix_y   = img->pix_y + block_y;
  int     pic_opix_x   = img->opix_x + block_x;
  int     pic_opix_y   = img->opix_y + block_y;
  int     pic_block_x = pic_pix_x >> 2;
  int     pic_block_y = pic_pix_y >> 2;
  double  min_rdcost  = 1e30;
  imgpel    *img_org, *img_prd;
  int       *residual;
  extern  int ****cofAC8x8;
  static int fadjust8x8[2][16][16];
  static int fadjust8x8Cr[2][2][16][16];
  extern  int ****cofAC8x8CbCr[2];
  int uv, c_nzCbCr[3];
  int left_available, up_available, all_available;
  int    (*curr_res)[16] = img->m7[0]; 
  imgpel (*curr_mpr)[16] = img->mpr[0];

  char   upMode;
  char   leftMode;
  int    mostProbableMode;

  PixelPos left_block;
  PixelPos top_block;

  //For residual DPCM
  Boolean lossless_qpprime = ((currMB->qp + img->bitdepth_luma_qp_scale)==0 && img->lossless_qpprime_flag==1);  

  getLuma4x4Neighbour(currMB, block_x - 1, block_y,     &left_block);
  getLuma4x4Neighbour(currMB, block_x,     block_y - 1, &top_block);

  if (input->UseConstrainedIntraPred)
  {
    top_block.available  = top_block.available ? img->intra_block [top_block.mb_addr] : 0;
    left_block.available = left_block.available ? img->intra_block [left_block.mb_addr] : 0;
  }

  if(b8 >> 1)
    upMode    =  top_block.available ? img->ipredmode8x8[top_block.pos_y ][top_block.pos_x ] : -1;
  else
    upMode    =  top_block.available ? img->ipredmode   [top_block.pos_y ][top_block.pos_x ] : -1;
  if(b8 & 0x01)
    leftMode  = left_block.available ? img->ipredmode8x8[left_block.pos_y][left_block.pos_x] : -1;
  else
    leftMode  = left_block.available ? img->ipredmode[left_block.pos_y][left_block.pos_x] : -1;

  mostProbableMode  = (upMode < 0 || leftMode < 0) ? DC_PRED : upMode < leftMode ? upMode : leftMode;

  *min_cost = INT_MAX;

  ipmode_DPCM = NO_INTRA_PMODE; //For residual DPCM

  //===== INTRA PREDICTION FOR 8x8 BLOCK =====
  intrapred_8x8 (currMB, PLANE_Y, pic_pix_x, pic_pix_y, &left_available, &up_available, &all_available);

  if( (active_sps->chroma_format_idc ==YUV444) && (IS_INDEPENDENT(input)==0) )
  { 
    select_plane(PLANE_U);
    intrapred_8x8(currMB, PLANE_U, pic_pix_x, pic_pix_y, &left_available, &up_available, &all_available);
    select_plane(PLANE_V);
    intrapred_8x8(currMB, PLANE_V, pic_pix_x, pic_pix_y, &left_available, &up_available, &all_available);
    select_plane(PLANE_Y);
  }
  
  
  //===== LOOP OVER ALL 8x8 INTRA PREDICTION MODES =====
  for (ipmode = 0; ipmode < NO_INTRA_PMODE; ipmode++)
  {
    if( (ipmode==DC_PRED) ||
        ((ipmode==VERT_PRED||ipmode==VERT_LEFT_PRED||ipmode==DIAG_DOWN_LEFT_PRED) && up_available ) ||
        ((ipmode==HOR_PRED||ipmode==HOR_UP_PRED) && left_available ) ||
        (all_available) )
    {
      if (!input->rdopt)
      {
        for (k=j=0; j<8; j++)
        {
          for (i=0; i<8; i++, k++)
          {
            diff64[k] = pCurImg[pic_opix_y+j][pic_opix_x+i] - img->mpr_8x8[0][ipmode][j][i];
          }
        }
        cost  = (ipmode == mostProbableMode) ? 0 : (int)floor(4 * lambda );
        cost += distortion8x8 (diff64);

        if( (active_sps->chroma_format_idc == YUV444) && (IS_INDEPENDENT(input)==0) )
        {
          int m;
          for (m = PLANE_U; m <= PLANE_V; m++)
          {
            for (k=j=0; j<8; j++)
              for (i=0; i<8; i++, k++)
              {
                diff64[k] = pImgOrg[m][pic_opix_y+j][pic_opix_x+i] - img->mpr_8x8[m][ipmode][j][i]; 
              }
              cost += distortion8x8 (diff64);
          }
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
        for (j=block_y; j < block_y  + 8; j++)
        {
          memcpy(&curr_mpr[j][block_x],img->mpr_8x8[0][ipmode][j - block_y], 8 * sizeof(imgpel));
          img_org  = &pCurImg[img->opix_y+j][pic_opix_x];
          img_prd  = &curr_mpr[j][block_x];
          residual = &curr_res[j][block_x];
          for (i=0; i<8; i++)
          {
            *residual++ = *img_org++ - *img_prd++;
          }
        }
        if( (img->yuv_format == YUV444) && (!IS_INDEPENDENT(input))) 
        {
          for (k = PLANE_U; k <= PLANE_V; k++)
          {
            for (j=0; j<8; j++)
            {
              memcpy(&img->mpr[k][block_y+j][block_x],img->mpr_8x8[k][ipmode][j], 8 * sizeof(imgpel));
              for (i=0; i<8; i++)
              {
                img->m7[k][block_y+j][block_x+i] = (int) (pImgOrg[k][pic_opix_y+j][pic_opix_x+i] - img->mpr_8x8[k][ipmode][j][i]);
              }
            }
          }
        }

        if( (img->yuv_format == YUV444) && (!IS_INDEPENDENT(input))) 
        {
          if((lossless_qpprime)&&(ipmode<2))  //For residual DPCM
          {
            Residual_DPCM_8x8(ipmode, 0, block_y, block_x);
            Residual_DPCM_8x8(ipmode, 1, block_y, block_x);
            Residual_DPCM_8x8(ipmode, 2, block_y, block_x);
            ipmode_DPCM=ipmode;
          }
        }
        else if ( (img->yuv_format == YUV444) && (IS_INDEPENDENT(input)))
        {
          if((lossless_qpprime)&&(ipmode<2))  //For residual DPCM
          {
            Residual_DPCM_8x8(ipmode, 0, block_y, block_x);
            ipmode_DPCM=ipmode;
          }
        }
        //===== store the coding state =====
        // store_coding_state_cs_cm(currMB);
        // get and check rate-distortion cost

        if ((rdcost = RDCost_for_8x8IntraBlocks (currMB, &c_nz, b8, ipmode, lambda, min_rdcost, mostProbableMode, c_nzCbCr)) < min_rdcost)
        {
          //--- set coefficients ---
          for(k=0; k<4; k++) // do 4x now
          {
            for (j=0; j<2; j++)
              memcpy(cofAC8x8[b8][k][j],img->cofAC[b8][k][j], 65 * sizeof(int));
          }

          //--- set reconstruction ---
          for (y=0; y<8; y++)
          {
            memcpy(rec8x8[0][y],&enc_picture->imgY[pic_pix_y+y][pic_pix_x], 8 * sizeof(imgpel));
          }

          if (img->AdaptiveRounding)
          {
            for (j=block_y; j<block_y + 8; j++)
              memcpy(&fadjust8x8[1][j][block_x],&img->fadjust8x8[1][j][block_x], 8 * sizeof(int));

            if(img->yuv_format == YUV444 && !IS_INDEPENDENT(input))
            {
              for (j=block_y; j<block_y + 8; j++)
              {
                memcpy(&fadjust8x8Cr[0][1][j][block_x],&img->fadjust8x8Cr[0][1][j][block_x], 8 * sizeof(int));
                memcpy(&fadjust8x8Cr[1][1][j][block_x],&img->fadjust8x8Cr[1][1][j][block_x], 8 * sizeof(int));
              }
            }            
          }
          if( img->yuv_format == YUV444 && !IS_INDEPENDENT(input) ) 
          { 
            //--- set coefficients ---
            for (uv=0; uv < 2; uv++)
            {
              for(k=0; k<4; k++) // do 4x now
              {
                for (j=0; j<2; j++)
                {
                  memcpy(cofAC8x8CbCr[uv][b8][k][j],img->cofAC[4+b8+4*uv][k][j], 65 * sizeof(int));
                }
              }
              cr_cbp[uv + 1] = c_nzCbCr[uv + 1];
              //--- set reconstruction ---
              for (y=0; y<8; y++)
              {
                memcpy(rec8x8[uv + 1][y],&enc_picture->imgUV[uv][pic_pix_y+y][pic_pix_x], 8 * sizeof(imgpel));
              }
            }
          }

          //--- flag if dct-coefficients must be coded ---
          nonzero = c_nz;

          //--- set best mode update minimum cost ---
          *min_cost   = rdcost;
          min_rdcost  = rdcost;
          best_ipmode = ipmode;
        }
        reset_coding_state_cs_cm(currMB);
      }
    }
  }

  //===== set intra mode prediction =====
  img->ipredmode8x8[pic_block_y][pic_block_x] = (char) best_ipmode;
  ipmode_DPCM = best_ipmode; //For residual DPCM

  if( (active_sps->chroma_format_idc==YUV444) && (IS_INDEPENDENT(input)==0) )
  {
    ColorPlane k;
    CbCr_predmode_8x8[b8] = best_ipmode; 
    for (k = PLANE_U; k <= PLANE_V; k++)
    {
      cr_cbp[k] = 0; 
      select_plane(k);
      for (j=0; j<8; j++)
      {
        for (i=0; i<8; i++)
        {
          img->mpr[k][block_y+j][block_x+i]  = img->mpr_8x8[k][best_ipmode][j][i]; 
          img->m7[k][block_y+j][block_x+i]   = pImgOrg[k][img->pix_y+block_y+j][img->pix_x+block_x+i] - img->mpr_8x8[k][best_ipmode][j][i];
        }
      }
      if(lossless_qpprime)   //For residual DPCM
      {
        ipmode_DPCM=best_ipmode; 
        if((best_ipmode<2))
        {
          Residual_DPCM_8x8(best_ipmode, k, block_y, block_x);
        }
      }

      if (dct_8x8(currMB, k, b8, &dummy, 1))
        cr_cbp[k] = 1;
    }
    select_plane(PLANE_Y);
  }
  
  currMB->intra_pred_modes8x8[4*b8] = (mostProbableMode == best_ipmode)
    ? -1
    : (best_ipmode < mostProbableMode ? best_ipmode : best_ipmode-1);

  for(j = img->mb_y*4+(b8 >> 1)*2; j < img->mb_y*4+(b8 >> 1)*2 + 2; j++)   //loop 4x4s in the subblock for 8x8 prediction setting
   memset(&img->ipredmode8x8[j][img->mb_x*4+(b8 & 0x01)*2], best_ipmode, 2 * sizeof(char));

  if (!input->rdopt)
  {
    // get prediction and prediction error
    for (j = block_y; j < block_y + 8; j++)
    {
      memcpy(&curr_mpr[j][block_x],img->mpr_8x8[0][best_ipmode][j - block_y], 8 * sizeof(imgpel));
      img_org  = &pCurImg[img->opix_y+j][pic_opix_x];
      img_prd  = &curr_mpr[j][block_x];
      residual = &curr_res[j][block_x];
      for (i=0; i<8; i++)
      {
        *residual++ = *img_org++ - *img_prd++;
      }
    }

    if(lossless_qpprime)   //For residual DPCM
    {
      ipmode_DPCM=best_ipmode;

      if((best_ipmode<2))
      {
        Residual_DPCM_8x8(best_ipmode, 0, block_y, block_x);
      }
    }
    nonzero = dct_8x8 (currMB, PLANE_Y, b8, &dummy, 1);    
  }
  else
  {
    //===== restore coefficients =====
    for(k=0; k<4; k++) // do 4x now
    {
      for (j=0; j<2; j++)
        memcpy(img->cofAC[b8][k][j],cofAC8x8[b8][k][j], 65 * sizeof(int));
    }

    if (img->AdaptiveRounding)
    {
      for (j=block_y; j< block_y + 8; j++)
        memcpy(&img->fadjust8x8[1][j][block_x], &fadjust8x8[1][j][block_x], 8 * sizeof(int));
      if(img->yuv_format == YUV444 && !IS_INDEPENDENT (input))
      {
        for (j=0; j<8; j++)
        {
          memcpy(&img->fadjust8x8Cr[0][1][block_y+j][block_x], &fadjust8x8Cr[0][1][block_y+j][block_x], 8 * sizeof(int));
          memcpy(&img->fadjust8x8Cr[1][1][block_y+j][block_x], &fadjust8x8Cr[1][1][block_y+j][block_x], 8 * sizeof(int));
        }
      }
    }

    //===== restore reconstruction and prediction (needed if single coeffs are removed) =====
    for (y=0; y<8; y++)
    {
      memcpy(&enc_picture->imgY[pic_pix_y + y][pic_pix_x], rec8x8[0][y], 8 * sizeof(imgpel));
      memcpy(&curr_mpr[block_y + y][block_x], img->mpr_8x8[0][best_ipmode][y], 8 * sizeof(imgpel));
    }
    if (img->yuv_format==YUV444 && !IS_INDEPENDENT(input) )
    {
      //===== restore coefficients =====
      for(k=0; k<4; k++) // do 4x now    
      {
        for (j=0; j<2; j++)
        {
          memcpy(img->cofAC[4+b8+4*0][k][j], cofAC8x8CbCr[0][b8][k][j], 65 * sizeof(int));
          memcpy(img->cofAC[4+b8+4*1][k][j], cofAC8x8CbCr[1][b8][k][j], 65 * sizeof(int));
        }
      }
      //===== restore reconstruction and prediction (needed if single coeffs are removed) =====
      for (y=0; y<8; y++) 
      {
        memcpy(&enc_picture->imgUV[0][pic_pix_y+y][pic_pix_x], rec8x8[1][y], 8 * sizeof(imgpel));
        memcpy(&enc_picture->imgUV[1][pic_pix_y+y][pic_pix_x], rec8x8[2][y], 8 * sizeof(imgpel));
        memcpy(&img->mpr[1][block_y+y][block_x], img->mpr_8x8[1][best_ipmode][y], 8 * sizeof(imgpel));
        memcpy(&img->mpr[2][block_y+y][block_x], img->mpr_8x8[2][best_ipmode][y], 8 * sizeof(imgpel));
      }
    }
  }

  return nonzero;
}



// Notation for comments regarding prediction and predictors.
// The pels of the 4x4 block are labelled a..p. The predictor pels above
// are labelled A..H, from the left I..P, and from above left X, as follows:
//
//  Z  A  B  C  D  E  F  G  H  I  J  K  L  M   N  O  P
//  Q  a1 b1 c1 d1 e1 f1 g1 h1
//  R  a2 b2 c2 d2 e2 f2 g2 h2
//  S  a3 b3 c3 d3 e3 f3 g3 h3
//  T  a4 b4 c4 d4 e4 f4 g4 h4
//  U  a5 b5 c5 d5 e5 f5 g5 h5
//  V  a6 b6 c6 d6 e6 f6 g6 h6
//  W  a7 b7 c7 d7 e7 f7 g7 h7
//  X  a8 b8 c8 d8 e8 f8 g8 h8


// Predictor array index definitions
#define P_Z (PredPel[0])
#define P_A (PredPel[1])
#define P_B (PredPel[2])
#define P_C (PredPel[3])
#define P_D (PredPel[4])
#define P_E (PredPel[5])
#define P_F (PredPel[6])
#define P_G (PredPel[7])
#define P_H (PredPel[8])
#define P_I (PredPel[9])
#define P_J (PredPel[10])
#define P_K (PredPel[11])
#define P_L (PredPel[12])
#define P_M (PredPel[13])
#define P_N (PredPel[14])
#define P_O (PredPel[15])
#define P_P (PredPel[16])
#define P_Q (PredPel[17])
#define P_R (PredPel[18])
#define P_S (PredPel[19])
#define P_T (PredPel[20])
#define P_U (PredPel[21])
#define P_V (PredPel[22])
#define P_W (PredPel[23])
#define P_X (PredPel[24])

/*!
 ************************************************************************
 * \brief
 *    Make intra 8x8 prediction according to all 9 prediction modes.
 *    The routine uses left and upper neighbouring points from
 *    previous coded blocks to do this (if available). Notice that
 *    inaccessible neighbouring points are signalled with a negative
 *    value in the predmode array .
 *
 *  \par Input:
 *     Starting point of current 8x8 block image posision
 *
 *  \par Output:
 *      none
 ************************************************************************
 */
void intrapred_8x8(Macroblock *currMB, ColorPlane pl, int img_x,int img_y, int *left_available, int *up_available, int *all_available)
{
  int i,j;
  int s0;
  static imgpel PredPel[25];  // array of predictor pels
  imgpel **img_enc = enc_picture->p_curr_img;
  imgpel *img_pel;
  static imgpel (*cur_pred)[16];
  imgpel (*curr_mpr_8x8)[16][16]  = img->mpr_8x8[pl];
  unsigned int dc_pred_value = img->dc_pred_value;

  int ioff = (img_x & 15);
  int joff = (img_y & 15);

  PixelPos pix_a[8];
  PixelPos pix_b, pix_c, pix_d;

  int block_available_up;
  int block_available_left;
  int block_available_up_left;
  int block_available_up_right;

  for (i=0;i<8;i++)
  {
    getNeighbour(currMB, ioff - 1, joff + i , IS_LUMA, &pix_a[i]);
  }

  getNeighbour(currMB, ioff    , joff - 1, IS_LUMA, &pix_b);
  getNeighbour(currMB, ioff + 8, joff - 1, IS_LUMA, &pix_c);
  getNeighbour(currMB, ioff - 1, joff - 1, IS_LUMA, &pix_d);

  pix_c.available = pix_c.available &&!(ioff == 8 && joff == 8);

  if (input->UseConstrainedIntraPred)
  {
    for (i=0, block_available_left=1; i<8;i++)
      block_available_left  &= pix_a[i].available ? img->intra_block[pix_a[i].mb_addr]: 0;

    block_available_up       = pix_b.available ? img->intra_block [pix_b.mb_addr] : 0;
    block_available_up_right = pix_c.available ? img->intra_block [pix_c.mb_addr] : 0;
    block_available_up_left  = pix_d.available ? img->intra_block [pix_d.mb_addr] : 0;
  }
  else
  {
    block_available_left     = pix_a[0].available;
    block_available_up       = pix_b.available;
    block_available_up_right = pix_c.available;
    block_available_up_left  = pix_d.available;
  }

  *left_available = block_available_left;
  *up_available   = block_available_up;
  *all_available  = block_available_up && block_available_left && block_available_up_left;

  i = (img_x & 15);
  j = (img_y & 15);

  // form predictor pels
  // form predictor pels
  if (block_available_up)
  {
    img_pel = &img_enc[pix_b.pos_y][pix_b.pos_x];
    P_A = *(img_pel++);
    P_B = *(img_pel++);
    P_C = *(img_pel++);
    P_D = *(img_pel++);
    P_E = *(img_pel++);
    P_F = *(img_pel++);
    P_G = *(img_pel++);
    P_H = *(img_pel);
  }
  else
  {
    P_A = P_B = P_C = P_D = P_E = P_F = P_G = P_H = dc_pred_value;
  }

  if (block_available_up_right)
  {
    img_pel = &img_enc[pix_c.pos_y][pix_c.pos_x];
    P_I = *(img_pel++);
    P_J = *(img_pel++);
    P_K = *(img_pel++);
    P_L = *(img_pel++);
    P_M = *(img_pel++);
    P_N = *(img_pel++);
    P_O = *(img_pel++);
    P_P = *(img_pel);

  }
  else
  {
    P_I = P_J = P_K = P_L = P_M = P_N = P_O = P_P = P_H;
  }

  if (block_available_left)
  {
    P_Q = img_enc[pix_a[0].pos_y][pix_a[0].pos_x];
    P_R = img_enc[pix_a[1].pos_y][pix_a[1].pos_x];
    P_S = img_enc[pix_a[2].pos_y][pix_a[2].pos_x];
    P_T = img_enc[pix_a[3].pos_y][pix_a[3].pos_x];
    P_U = img_enc[pix_a[4].pos_y][pix_a[4].pos_x];
    P_V = img_enc[pix_a[5].pos_y][pix_a[5].pos_x];
    P_W = img_enc[pix_a[6].pos_y][pix_a[6].pos_x];
    P_X = img_enc[pix_a[7].pos_y][pix_a[7].pos_x];
  }
  else
  {
    P_Q = P_R = P_S = P_T = P_U = P_V = P_W = P_X = dc_pred_value;
  }

  if (block_available_up_left)
  {
    P_Z = img_enc[pix_d.pos_y][pix_d.pos_x];
  }
  else
  {
    P_Z = dc_pred_value;
  }

  for(i=0;i<9;i++)
    curr_mpr_8x8[i][0][0]=-1;

  LowPassForIntra8x8Pred(&(P_Z), block_available_up_left, block_available_up, block_available_left);

  ///////////////////////////////
  // make DC prediction
  ///////////////////////////////
  s0 = 0;
  if (block_available_up && block_available_left)
  {
    // no edge
    s0 = rshift_rnd_sf((P_A + P_B + P_C + P_D + P_E + P_F + P_G + P_H + P_Q + P_R + P_S + P_T + P_U + P_V + P_W + P_X), 4);
  }
  else if (!block_available_up && block_available_left)
  {
    // upper edge
    s0 = rshift_rnd_sf((P_Q + P_R + P_S + P_T + P_U + P_V + P_W + P_X), 3);
  }
  else if (block_available_up && !block_available_left)
  {
    // left edge
    s0 = rshift_rnd_sf((P_A + P_B + P_C + P_D + P_E + P_F + P_G + P_H), 3);
  }
  else //if (!block_available_up && !block_available_left)
  {
    // top left corner, nothing to predict from
    s0 = dc_pred_value;
  }

  // store DC prediction
  cur_pred = curr_mpr_8x8[DC_PRED];
  for (j=0; j < BLOCK_SIZE_8x8; j++)
  {
    for (i=0; i < BLOCK_SIZE_8x8; i++)
    {
      cur_pred[j][i] = (imgpel) s0;
    }
  }


  ///////////////////////////////
  // make horiz and vert prediction
  ///////////////////////////////
  cur_pred = curr_mpr_8x8[VERT_PRED];
  for (i=0; i < BLOCK_SIZE_8x8; i++)
  {
    cur_pred[0][i] =
    cur_pred[1][i] =
    cur_pred[2][i] =
    cur_pred[3][i] =
    cur_pred[4][i] =
    cur_pred[5][i] =
    cur_pred[6][i] =
    cur_pred[7][i] = (imgpel)(&P_A)[i];
  }
  if(!block_available_up)
    cur_pred[0][0]=-1;

  cur_pred = curr_mpr_8x8[HOR_PRED];
  for (i=0; i < BLOCK_SIZE_8x8; i++)
  {
    cur_pred[i][0]  =
    cur_pred[i][1]  =
    cur_pred[i][2]  =
    cur_pred[i][3]  =
    cur_pred[i][4]  =
    cur_pred[i][5]  =
    cur_pred[i][6]  =
    cur_pred[i][7]  = (imgpel) (&P_Q)[i];
  }
  if(!block_available_left)
    cur_pred[0][0]=-1;

  ///////////////////////////////////
  // make diagonal down left prediction
  ///////////////////////////////////
  if (block_available_up)
  {
    // Mode DIAG_DOWN_LEFT_PRED
    cur_pred = curr_mpr_8x8[DIAG_DOWN_LEFT_PRED];
    cur_pred[0][0] = (imgpel) ((P_A + P_C + 2*(P_B) + 2) >> 2);
    cur_pred[0][1] =
    cur_pred[1][0] = (imgpel) ((P_B + P_D + 2*(P_C) + 2) >> 2);
    cur_pred[0][2] =
    cur_pred[1][1] =
    cur_pred[2][0] = (imgpel) ((P_C + P_E + 2*(P_D) + 2) >> 2);
    cur_pred[0][3] =
    cur_pred[1][2] =
    cur_pred[2][1] =
    cur_pred[3][0] = (imgpel) ((P_D + P_F + 2*(P_E) + 2) >> 2);
    cur_pred[0][4] =
    cur_pred[1][3] =
    cur_pred[2][2] =
    cur_pred[3][1] =
    cur_pred[4][0] = (imgpel) ((P_E + P_G + 2*(P_F) + 2) >> 2);
    cur_pred[0][5] =
    cur_pred[1][4] =
    cur_pred[2][3] =
    cur_pred[3][2] =
    cur_pred[4][1] =
    cur_pred[5][0] = (imgpel) ((P_F + P_H + 2*(P_G) + 2) >> 2);
    cur_pred[0][6] =
    cur_pred[1][5] =
    cur_pred[2][4] =
    cur_pred[3][3] =
    cur_pred[4][2] =
    cur_pred[5][1] =
    cur_pred[6][0] = (imgpel) ((P_G + P_I + 2*(P_H) + 2) >> 2);
    cur_pred[0][7] =
    cur_pred[1][6] =
    cur_pred[2][5] =
    cur_pred[3][4] =
    cur_pred[4][3] =
    cur_pred[5][2] =
    cur_pred[6][1] =
    cur_pred[7][0] = (imgpel) ((P_H + P_J + 2*(P_I) + 2) >> 2);
    cur_pred[1][7] =
    cur_pred[2][6] =
    cur_pred[3][5] =
    cur_pred[4][4] =
    cur_pred[5][3] =
    cur_pred[6][2] =
    cur_pred[7][1] = (imgpel) ((P_I + P_K + 2*(P_J) + 2) >> 2);
    cur_pred[2][7] =
    cur_pred[3][6] =
    cur_pred[4][5] =
    cur_pred[5][4] =
    cur_pred[6][3] =
    cur_pred[7][2] = (imgpel) ((P_J + P_L + 2*(P_K) + 2) >> 2);
    cur_pred[3][7] =
    cur_pred[4][6] =
    cur_pred[5][5] =
    cur_pred[6][4] =
    cur_pred[7][3] = (imgpel) ((P_K + P_M + 2*(P_L) + 2) >> 2);
    cur_pred[4][7] =
    cur_pred[5][6] =
    cur_pred[6][5] =
    cur_pred[7][4] = (imgpel) ((P_L + P_N + 2*(P_M) + 2) >> 2);
    cur_pred[5][7] =
    cur_pred[6][6] =
    cur_pred[7][5] = (imgpel) ((P_M + P_O + 2*(P_N) + 2) >> 2);
    cur_pred[6][7] =
    cur_pred[7][6] = (imgpel) ((P_N + P_P + 2*(P_O) + 2) >> 2);
    cur_pred[7][7] = (imgpel) ((P_O + 3*(P_P) + 2) >> 2);

    ///////////////////////////////////
    // make vertical left prediction
    ///////////////////////////////////
    cur_pred = curr_mpr_8x8[VERT_LEFT_PRED];
    cur_pred[0][0] = (imgpel) ((P_A + P_B + 1) >> 1);
    cur_pred[0][1] =
    cur_pred[2][0] = (imgpel) ((P_B + P_C + 1) >> 1);
    cur_pred[0][2] =
    cur_pred[2][1] =
    cur_pred[4][0] = (imgpel) ((P_C + P_D + 1) >> 1);
    cur_pred[0][3] =
    cur_pred[2][2] =
    cur_pred[4][1] =
    cur_pred[6][0] = (imgpel) ((P_D + P_E + 1) >> 1);
    cur_pred[0][4] =
    cur_pred[2][3] =
    cur_pred[4][2] =
    cur_pred[6][1] = (imgpel) ((P_E + P_F + 1) >> 1);
    cur_pred[0][5] =
    cur_pred[2][4] =
    cur_pred[4][3] =
    cur_pred[6][2] = (imgpel) ((P_F + P_G + 1) >> 1);
    cur_pred[0][6] =
    cur_pred[2][5] =
    cur_pred[4][4] =
    cur_pred[6][3] = (imgpel) ((P_G + P_H + 1) >> 1);
    cur_pred[0][7] =
    cur_pred[2][6] =
    cur_pred[4][5] =
    cur_pred[6][4] = (imgpel) ((P_H + P_I + 1) >> 1);
    cur_pred[2][7] =
    cur_pred[4][6] =
    cur_pred[6][5] = (imgpel) ((P_I + P_J + 1) >> 1);
    cur_pred[4][7] =
    cur_pred[6][6] = (imgpel) ((P_J + P_K + 1) >> 1);
    cur_pred[6][7] = (imgpel) ((P_K + P_L + 1) >> 1);
    cur_pred[1][0] = (imgpel) ((P_A + P_C + 2*P_B + 2) >> 2);
    cur_pred[1][1] =
    cur_pred[3][0] = (imgpel) ((P_B + P_D + 2*P_C + 2) >> 2);
    cur_pred[1][2] =
    cur_pred[3][1] =
    cur_pred[5][0] = (imgpel) ((P_C + P_E + 2*P_D + 2) >> 2);
    cur_pred[1][3] =
    cur_pred[3][2] =
    cur_pred[5][1] =
    cur_pred[7][0] = (imgpel) ((P_D + P_F + 2*P_E + 2) >> 2);
    cur_pred[1][4] =
    cur_pred[3][3] =
    cur_pred[5][2] =
    cur_pred[7][1] = (imgpel) ((P_E + P_G + 2*P_F + 2) >> 2);
    cur_pred[1][5] =
    cur_pred[3][4] =
    cur_pred[5][3] =
    cur_pred[7][2] = (imgpel) ((P_F + P_H + 2*P_G + 2) >> 2);
    cur_pred[1][6] =
    cur_pred[3][5] =
    cur_pred[5][4] =
    cur_pred[7][3] = (imgpel) ((P_G + P_I + 2*P_H + 2) >> 2);
    cur_pred[1][7] =
    cur_pred[3][6] =
    cur_pred[5][5] =
    cur_pred[7][4] = (imgpel) ((P_H + P_J + 2*P_I + 2) >> 2);
    cur_pred[3][7] =
    cur_pred[5][6] =
    cur_pred[7][5] = (imgpel) ((P_I + P_K + 2*P_J + 2) >> 2);
    cur_pred[5][7] =
    cur_pred[7][6] = (imgpel) ((P_J + P_L + 2*P_K + 2) >> 2);
    cur_pred[7][7] = (imgpel) ((P_K + P_M + 2*P_L + 2) >> 2);
  }

  ///////////////////////////////////
  // make diagonal down right prediction
  ///////////////////////////////////
  if (block_available_up && block_available_left && block_available_up_left)
  {
    // Mode DIAG_DOWN_RIGHT_PRED
    cur_pred = curr_mpr_8x8[DIAG_DOWN_RIGHT_PRED];
    cur_pred[7][0] = (imgpel) ((P_X + P_V + 2*(P_W) + 2) >> 2);
    cur_pred[6][0] =
    cur_pred[7][1] = (imgpel) ((P_W + P_U + 2*(P_V) + 2) >> 2);
    cur_pred[5][0] =
    cur_pred[6][1] =
    cur_pred[7][2] = (imgpel) ((P_V + P_T + 2*(P_U) + 2) >> 2);
    cur_pred[4][0] =
    cur_pred[5][1] =
    cur_pred[6][2] =
    cur_pred[7][3] = (imgpel) ((P_U + P_S + 2*(P_T) + 2) >> 2);
    cur_pred[3][0] =
    cur_pred[4][1] =
    cur_pred[5][2] =
    cur_pred[6][3] =
    cur_pred[7][4] = (imgpel) ((P_T + P_R + 2*(P_S) + 2) >> 2);
    cur_pred[2][0] =
    cur_pred[3][1] =
    cur_pred[4][2] =
    cur_pred[5][3] =
    cur_pred[6][4] =
    cur_pred[7][5] = (imgpel) ((P_S + P_Q + 2*(P_R) + 2) >> 2);
    cur_pred[1][0] =
    cur_pred[2][1] =
    cur_pred[3][2] =
    cur_pred[4][3] =
    cur_pred[5][4] =
    cur_pred[6][5] =
    cur_pred[7][6] = (imgpel) ((P_R + P_Z + 2*(P_Q) + 2) >> 2);
    cur_pred[0][0] =
    cur_pred[1][1] =
    cur_pred[2][2] =
    cur_pred[3][3] =
    cur_pred[4][4] =
    cur_pred[5][5] =
    cur_pred[6][6] =
    cur_pred[7][7] = (imgpel) ((P_Q + P_A + 2*(P_Z) + 2) >> 2);
    cur_pred[0][1] =
    cur_pred[1][2] =
    cur_pred[2][3] =
    cur_pred[3][4] =
    cur_pred[4][5] =
    cur_pred[5][6] =
    cur_pred[6][7] = (imgpel) ((P_Z + P_B + 2*(P_A) + 2) >> 2);
    cur_pred[0][2] =
    cur_pred[1][3] =
    cur_pred[2][4] =
    cur_pred[3][5] =
    cur_pred[4][6] =
    cur_pred[5][7] = (imgpel) ((P_A + P_C + 2*(P_B) + 2) >> 2);
    cur_pred[0][3] =
    cur_pred[1][4] =
    cur_pred[2][5] =
    cur_pred[3][6] =
    cur_pred[4][7] = (imgpel) ((P_B + P_D + 2*(P_C) + 2) >> 2);
    cur_pred[0][4] =
    cur_pred[1][5] =
    cur_pred[2][6] =
    cur_pred[3][7] = (imgpel) ((P_C + P_E + 2*(P_D) + 2) >> 2);
    cur_pred[0][5] =
    cur_pred[1][6] =
    cur_pred[2][7] = (imgpel) ((P_D + P_F + 2*(P_E) + 2) >> 2);
    cur_pred[0][6] =
    cur_pred[1][7] = (imgpel) ((P_E + P_G + 2*(P_F) + 2) >> 2);
    cur_pred[0][7] = (imgpel) ((P_F + P_H + 2*(P_G) + 2) >> 2);

  ///////////////////////////////////
  // make vertical right prediction
  ///////////////////////////////////
    cur_pred = curr_mpr_8x8[VERT_RIGHT_PRED];
    cur_pred[0][0] =
    cur_pred[2][1] =
    cur_pred[4][2] =
    cur_pred[6][3] = (imgpel) ((P_Z + P_A + 1) >> 1);
    cur_pred[0][1] =
    cur_pred[2][2] =
    cur_pred[4][3] =
    cur_pred[6][4] = (imgpel) ((P_A + P_B + 1) >> 1);
    cur_pred[0][2] =
    cur_pred[2][3] =
    cur_pred[4][4] =
    cur_pred[6][5] = (imgpel) ((P_B + P_C + 1) >> 1);
    cur_pred[0][3] =
    cur_pred[2][4] =
    cur_pred[4][5] =
    cur_pred[6][6] = (imgpel) ((P_C + P_D + 1) >> 1);
    cur_pred[0][4] =
    cur_pred[2][5] =
    cur_pred[4][6] =
    cur_pred[6][7] = (imgpel) ((P_D + P_E + 1) >> 1);
    cur_pred[0][5] =
    cur_pred[2][6] =
    cur_pred[4][7] = (imgpel) ((P_E + P_F + 1) >> 1);
    cur_pred[0][6] =
    cur_pred[2][7] = (imgpel) ((P_F + P_G + 1) >> 1);
    cur_pred[0][7] = (imgpel) ((P_G + P_H + 1) >> 1);
    cur_pred[1][0] =
    cur_pred[3][1] =
    cur_pred[5][2] =
    cur_pred[7][3] = (imgpel) ((P_Q + P_A + 2*P_Z + 2) >> 2);
    cur_pred[1][1] =
    cur_pred[3][2] =
    cur_pred[5][3] =
    cur_pred[7][4] = (imgpel) ((P_Z + P_B + 2*P_A + 2) >> 2);
    cur_pred[1][2] =
    cur_pred[3][3] =
    cur_pred[5][4] =
    cur_pred[7][5] = (imgpel) ((P_A + P_C + 2*P_B + 2) >> 2);
    cur_pred[1][3] =
    cur_pred[3][4] =
    cur_pred[5][5] =
    cur_pred[7][6] = (imgpel) ((P_B + P_D + 2*P_C + 2) >> 2);
    cur_pred[1][4] =
    cur_pred[3][5] =
    cur_pred[5][6] =
    cur_pred[7][7] = (imgpel) ((P_C + P_E + 2*P_D + 2) >> 2);
    cur_pred[1][5] =
    cur_pred[3][6] =
    cur_pred[5][7] = (imgpel) ((P_D + P_F + 2*P_E + 2) >> 2);
    cur_pred[1][6] =
    cur_pred[3][7] = (imgpel) ((P_E + P_G + 2*P_F + 2) >> 2);
    cur_pred[1][7] = (imgpel) ((P_F + P_H + 2*P_G + 2) >> 2);
    cur_pred[2][0] =
    cur_pred[4][1] =
    cur_pred[6][2] = (imgpel) ((P_R + P_Z + 2*P_Q + 2) >> 2);
    cur_pred[3][0] =
    cur_pred[5][1] =
    cur_pred[7][2] = (imgpel) ((P_S + P_Q + 2*P_R + 2) >> 2);
    cur_pred[4][0] =
    cur_pred[6][1] = (imgpel) ((P_T + P_R + 2*P_S + 2) >> 2);
    cur_pred[5][0] =
    cur_pred[7][1] = (imgpel) ((P_U + P_S + 2*P_T + 2) >> 2);
    cur_pred[6][0] = (imgpel) ((P_V + P_T + 2*P_U + 2) >> 2);
    cur_pred[7][0] = (imgpel) ((P_W + P_U + 2*P_V + 2) >> 2);

  ///////////////////////////////////
  // make horizontal down prediction
  ///////////////////////////////////
    cur_pred = img->mpr_8x8[0][HOR_DOWN_PRED];
    cur_pred[0][0] =
    cur_pred[1][2] =
    cur_pred[2][4] =
    cur_pred[3][6] = (imgpel) ((P_Q + P_Z + 1) >> 1);
    cur_pred[1][0] =
    cur_pred[2][2] =
    cur_pred[3][4] =
    cur_pred[4][6] = (imgpel) ((P_R + P_Q + 1) >> 1);
    cur_pred[2][0] =
    cur_pred[3][2] =
    cur_pred[4][4] =
    cur_pred[5][6] = (imgpel) ((P_S + P_R + 1) >> 1);
    cur_pred[3][0] =
    cur_pred[4][2] =
    cur_pred[5][4] =
    cur_pred[6][6] = (imgpel) ((P_T + P_S + 1) >> 1);
    cur_pred[4][0] =
    cur_pred[5][2] =
    cur_pred[6][4] =
    cur_pred[7][6] = (imgpel) ((P_U + P_T + 1) >> 1);
    cur_pred[5][0] =
    cur_pred[6][2] =
    cur_pred[7][4] = (imgpel) ((P_V + P_U + 1) >> 1);
    cur_pred[6][0] =
    cur_pred[7][2] = (imgpel) ((P_W + P_V + 1) >> 1);
    cur_pred[7][0] = (imgpel) ((P_X + P_W + 1) >> 1);
    cur_pred[0][1] =
    cur_pred[1][3] =
    cur_pred[2][5] =
    cur_pred[3][7] = (imgpel) ((P_Q + P_A + 2*P_Z + 2) >> 2);
    cur_pred[1][1] =
    cur_pred[2][3] =
    cur_pred[3][5] =
    cur_pred[4][7] = (imgpel) ((P_Z + P_R + 2*P_Q + 2) >> 2);
    cur_pred[2][1] =
    cur_pred[3][3] =
    cur_pred[4][5] =
    cur_pred[5][7] = (imgpel) ((P_Q + P_S + 2*P_R + 2) >> 2);
    cur_pred[3][1] =
    cur_pred[4][3] =
    cur_pred[5][5] =
    cur_pred[6][7] = (imgpel) ((P_R + P_T + 2*P_S + 2) >> 2);
    cur_pred[4][1] =
    cur_pred[5][3] =
    cur_pred[6][5] =
    cur_pred[7][7] = (imgpel) ((P_S + P_U + 2*P_T + 2) >> 2);
    cur_pred[5][1] =
    cur_pred[6][3] =
    cur_pred[7][5] = (imgpel) ((P_T + P_V + 2*P_U + 2) >> 2);
    cur_pred[6][1] =
    cur_pred[7][3] = (imgpel) ((P_U + P_W + 2*P_V + 2) >> 2);
    cur_pred[7][1] = (imgpel) ((P_V + P_X + 2*P_W + 2) >> 2);
    cur_pred[0][2] =
    cur_pred[1][4] =
    cur_pred[2][6] = (imgpel) ((P_Z + P_B + 2*P_A + 2) >> 2);
    cur_pred[0][3] =
    cur_pred[1][5] =
    cur_pred[2][7] = (imgpel) ((P_A + P_C + 2*P_B + 2) >> 2);
    cur_pred[0][4] =
    cur_pred[1][6] = (imgpel) ((P_B + P_D + 2*P_C + 2) >> 2);
    cur_pred[0][5] =
    cur_pred[1][7] = (imgpel) ((P_C + P_E + 2*P_D + 2) >> 2);
    cur_pred[0][6] = (imgpel) ((P_D + P_F + 2*P_E + 2) >> 2);
    cur_pred[0][7] = (imgpel) ((P_E + P_G + 2*P_F + 2) >> 2);
  }

  ///////////////////////////////////
  // make horizontal up prediction
  ///////////////////////////////////
  if (block_available_left)
  {
    cur_pred = curr_mpr_8x8[HOR_UP_PRED];
    cur_pred[0][0] = (imgpel) ((P_Q + P_R + 1) >> 1);
    cur_pred[1][0] =
    cur_pred[0][2] = (imgpel) ((P_R + P_S + 1) >> 1);
    cur_pred[2][0] =
    cur_pred[1][2] =
    cur_pred[0][4] = (imgpel) ((P_S + P_T + 1) >> 1);
    cur_pred[3][0] =
    cur_pred[2][2] =
    cur_pred[1][4] =
    cur_pred[0][6] = (imgpel) ((P_T + P_U + 1) >> 1);
    cur_pred[4][0] =
    cur_pred[3][2] =
    cur_pred[2][4] =
    cur_pred[1][6] = (imgpel) ((P_U + P_V + 1) >> 1);
    cur_pred[5][0] =
    cur_pred[4][2] =
    cur_pred[3][4] =
    cur_pred[2][6] = (imgpel) ((P_V + P_W + 1) >> 1);
    cur_pred[6][0] =
    cur_pred[5][2] =
    cur_pred[4][4] =
    cur_pred[3][6] = (imgpel) ((P_W + P_X + 1) >> 1);
    cur_pred[4][6] =
    cur_pred[4][7] =
    cur_pred[5][4] =
    cur_pred[5][5] =
    cur_pred[5][6] =
    cur_pred[5][7] =
    cur_pred[6][2] =
    cur_pred[6][3] =
    cur_pred[6][4] =
    cur_pred[6][5] =
    cur_pred[6][6] =
    cur_pred[6][7] =
    cur_pred[7][0] =
    cur_pred[7][1] =
    cur_pred[7][2] =
    cur_pred[7][3] =
    cur_pred[7][4] =
    cur_pred[7][5] =
    cur_pred[7][6] =
    cur_pred[7][7] = (imgpel) P_X;
    cur_pred[6][1] =
    cur_pred[5][3] =
    cur_pred[4][5] =
    cur_pred[3][7] = (imgpel) ((P_W + 3*P_X + 2) >> 2);
    cur_pred[5][1] =
    cur_pred[4][3] =
    cur_pred[3][5] =
    cur_pred[2][7] = (imgpel) ((P_X + P_V + 2*P_W + 2) >> 2);
    cur_pred[4][1] =
    cur_pred[3][3] =
    cur_pred[2][5] =
    cur_pred[1][7] = (imgpel) ((P_W + P_U + 2*P_V + 2) >> 2);
    cur_pred[3][1] =
    cur_pred[2][3] =
    cur_pred[1][5] =
    cur_pred[0][7] = (imgpel) ((P_V + P_T + 2*P_U + 2) >> 2);
    cur_pred[2][1] =
    cur_pred[1][3] =
    cur_pred[0][5] = (imgpel) ((P_U + P_S + 2*P_T + 2) >> 2);
    cur_pred[1][1] =
    cur_pred[0][3] = (imgpel) ((P_T + P_R + 2*P_S + 2) >> 2);
    cur_pred[0][1] = (imgpel) ((P_S + P_Q + 2*P_R + 2) >> 2);
  }
}

/*!
 *************************************************************************************
 * \brief
 *    Prefiltering for Intra8x8 prediction
 *************************************************************************************
 */
void LowPassForIntra8x8Pred(imgpel *PredPel, int block_up_left, int block_up, int block_left)
{
  int i;
  static imgpel LoopArray[25];

  memcpy(LoopArray,PredPel, 25 * sizeof(imgpel));

  if(block_up)
  {
    if(block_up_left)
    {
      LoopArray[1] = (((&P_Z)[0] + ((&P_Z)[1]<<1) + (&P_Z)[2] + 2)>>2);
    }
    else
      LoopArray[1] = (((&P_Z)[1] + ((&P_Z)[1]<<1) + (&P_Z)[2] + 2)>>2);


    for(i = 2; i <16; i++)
    {
      LoopArray[i] = (((&P_Z)[i-1] + ((&P_Z)[i]<<1) + (&P_Z)[i+1] + 2)>>2);
    }
    LoopArray[16] = ((P_P + (P_P<<1) + P_O + 2)>>2);
  }

  if(block_up_left)
  {
    if(block_up && block_left)
    {
      LoopArray[0] = ((P_Q + (P_Z<<1) + P_A +2)>>2);
    }
    else
    {
      if(block_up)
        LoopArray[0] = ((P_Z + (P_Z<<1) + P_A +2)>>2);
      else
        if(block_left)
          LoopArray[0] = ((P_Z + (P_Z<<1) + P_Q +2)>>2);
    }
  }

  if(block_left)
  {
    if(block_up_left)
      LoopArray[17] = ((P_Z + (P_Q<<1) + P_R + 2)>>2);
    else
      LoopArray[17] = ((P_Q + (P_Q<<1) + P_R + 2)>>2);

    for(i = 18; i <24; i++)
    {
      LoopArray[i] = (((&P_Z)[i-1] + ((&P_Z)[i]<<1) + (&P_Z)[i+1] + 2)>>2);
    }
    LoopArray[24] = ((P_W + (P_X<<1) + P_X + 2) >> 2);
  }

  memcpy(PredPel, LoopArray, 25 * sizeof(imgpel));
}





/*!
 *************************************************************************************
 * \brief
 *    R-D Cost for an 8x8 Intra block
 *************************************************************************************
 */

double RDCost_for_8x8IntraBlocks(Macroblock *currMB, int *nonzero, int b8, int ipmode, double lambda, double min_rdcost, int mostProbableMode, int c_nzCbCr[3])
{
  double  rdcost = 0.0;
  int     dummy;
  int     x, y, rate;
  int64   distortion  = 0;
  int     block_x     = (b8 & 0x01) << 3;
  int     block_y     = (b8 >> 1) << 3;
  int     pic_pix_x   = img->pix_x + block_x;
  int     pic_pix_y   = img->pix_y + block_y;
  int     pic_opix_y  = img->opix_y + block_y;
  imgpel  *img_org, *img_enc;

  Slice          *currSlice =  img->currentSlice;
  SyntaxElement  se;
  const int      *partMap   = assignSE2partition[input->partition_mode];
  DataPartition  *dataPart;

  //===== perform DCT, Q, IQ, IDCT, Reconstruction =====
  dummy = 0;

  *nonzero = dct_8x8 (currMB, PLANE_Y, b8, &dummy, 1);

  //===== get distortion (SSD) of 8x8 block =====
  for (y=0; y<8; y++)
  {
    img_org = pCurImg[pic_opix_y+y];
    img_enc = enc_picture->imgY[pic_pix_y+y];
    for (x = pic_pix_x; x < pic_pix_x + 8; x++)
      distortion += iabs2( img_org[x] - img_enc[x]);
  }

    if( img->yuv_format==YUV444 && !IS_INDEPENDENT(input) ) // INDEPENDENT_MODE 2006.07.26
    {
      ColorPlane k;

      for (k = PLANE_U; k <= PLANE_V; k++)
      {
        select_plane(k);
    /*    for (j=0; j<8; j++)   //KHHan, I think these line are not necessary
        {
          for (i=0; i<8; i++)
          {         
            img->mpr[k][block_y+j][block_x+i]  = img->mpr_8x8[k][ipmode][j][i];
            img->m7[k][j][i] = pImgOrg[k][img->pix_y+block_y+j][img->pix_x+block_x+i] - img->mpr_8x8[k][ipmode][j][i];
          }
        }*/
        c_nzCbCr[k]=dct_8x8(currMB, k, b8, &dummy,1);

        for( y = 0; y < 8; y++ )
          for (x=pic_pix_x; x<pic_pix_x+8; x++)
            distortion +=iabs2(pImgOrg[k][pic_opix_y+y][x] - enc_picture->p_curr_img[pic_pix_y+y][x]);
      }
    ipmode_DPCM=NO_INTRA_PMODE;  //For residual DPCM
      select_plane(PLANE_Y);
    }
  else if( img->yuv_format==YUV444 && IS_INDEPENDENT(input) )  //For residual DPCM
  {
    ipmode_DPCM=NO_INTRA_PMODE;  
  }

    
  //===== RATE for INTRA PREDICTION MODE  (SYMBOL MODE MUST BE SET TO CAVLC) =====
  se.value1 = (mostProbableMode == ipmode) ? -1 : ipmode < mostProbableMode ? ipmode : ipmode-1;

  //--- set position and type ---
  se.context = b8;
  se.type    = SE_INTRAPREDMODE;

  //--- choose data partition ---
  if (img->type!=B_SLICE)
    dataPart = &(currSlice->partArr[partMap[SE_INTRAPREDMODE]]);
  else
    dataPart = &(currSlice->partArr[partMap[SE_BFRAME]]);

  //--- encode and update rate ---
  writeIntraPredMode (&se, dataPart);

  rate = se.len;

  //===== RATE for LUMINANCE COEFFICIENTS =====

  if (input->symbol_mode == CAVLC)
  {      
      if ( img->yuv_format==YUV444 && !IS_INDEPENDENT(input) )
      {
        int b4;
        for(b4=0; b4<4; b4++)
        {
          rate  += writeCoeff4x4_CAVLC (currMB, LUMA, b8, b4, 0);
          rate  += writeCoeff4x4_CAVLC (currMB, CB, b8, b4, 0);
          rate  += writeCoeff4x4_CAVLC (currMB, CR, b8, b4, 0);
        }
      }
      else
      {
        rate  += writeCoeff4x4_CAVLC (currMB, LUMA, b8, 0, 0);
        rate  += writeCoeff4x4_CAVLC (currMB, LUMA, b8, 1, 0);
        rate  += writeCoeff4x4_CAVLC (currMB, LUMA, b8, 2, 0);
        rate  += writeCoeff4x4_CAVLC (currMB, LUMA, b8, 3, 0);
      }
  }
  else
  {
    rate  += writeCoeff8x8_CABAC (currMB, PLANE_Y, b8, 1);
    if( img->yuv_format==YUV444 && !IS_INDEPENDENT(input) )
    {
      rate  += writeCoeff8x8_CABAC (currMB, PLANE_U, b8, 1);
      rate  += writeCoeff8x8_CABAC (currMB, PLANE_V, b8, 1);
    }
  }

  rdcost = (double)distortion + lambda*(double)rate;

  return rdcost;
}

/*!
 ************************************************************************
 * \brief
 *    The routine performs transform,quantization,inverse transform, adds the diff.
 *    to the prediction and writes the result to the decoded luma frame. Includes the
 *    RD constrained quantization also.
 *
 * \par Input:
 *    b8: Block position inside a macro block (0,1,2,3).
 *
 * \par Output:
 *    nonzero: 0 if no levels are nonzero.  1 if there are nonzero levels.
 *    coeff_cost: Counter for nonzero coefficients, used to discard expensive levels.
 ************************************************************************
 */

#define MC(coeff) ((coeff)&3)

int dct_8x8(Macroblock *currMB, ColorPlane pl, int b8,int *coeff_cost, int intra)
{
  int i,j,ilev,coeff_ctr;
  int level,scan_pos = 0,run = -1;
  int nonzero = FALSE;  

  int block_x = 8*(b8 & 0x01);
  int block_y = 8*(b8 >> 1);
  int pl_off = pl<<2;
  int*  ACLevel = img->cofAC[b8+pl_off][0][0];
  int*  ACRun   = img->cofAC[b8+pl_off][0][1];  
  const byte *c_cost     = COEFF_COST8x8[input->disthres];  
  imgpel **img_enc       = enc_picture->p_curr_img;
  imgpel (*curr_mpr)[16] = img->mpr[pl];
  int    (*curr_res)[16] = img->m7[pl];   
  
  int max_imgpel_value   = img->max_imgpel_value;
  int scan_poss[4] = { 0 }, runs[4] = { -1, -1, -1, -1 };
  int MCcoeff = 0;
  static imgpel *img_Y, *predY;
  int *m7;
  int qp = currMB->qp_scaled[pl];

  Boolean lossless_qpprime = (Boolean) ((qp == 0) && img->lossless_qpprime_flag==1);
  const byte (*pos_scan)[2] = currMB->is_field_mode ? FIELD_SCAN8x8 : SNGL_SCAN8x8;

  int qp_per = qp_per_matrix[qp];
  int qp_rem = qp_rem_matrix[qp];
  int q_bits = Q_BITS_8 + qp_per;
  int **fadjust8x8 = img->AdaptiveRounding ? (pl ? &img->fadjust8x8Cr[pl-1][intra][block_y] : &img->fadjust8x8[intra][block_y]) :NULL;
  int **levelscale    = LevelScale8x8Comp   [pl][intra][qp_rem];
  int **invlevelscale = InvLevelScale8x8Comp[pl][intra][qp_rem];
  int **leveloffset   = LevelOffset8x8Comp  [pl][intra][qp];

  if (!lossless_qpprime)
  {
    // forward 8x8 transform
    forward8x8(curr_res, curr_res, block_y, block_x);

    // Quant

    for (coeff_ctr = 0; coeff_ctr < 64; coeff_ctr++)
    {
      i=pos_scan[coeff_ctr][0];
      j=pos_scan[coeff_ctr][1];

      run++;
      ilev=0;

      if (currMB->luma_transform_size_8x8_flag && input->symbol_mode == CAVLC)
      {
        MCcoeff = MC(coeff_ctr);
        runs[MCcoeff]++;
      }

      m7 = &curr_res[block_y + j][block_x];
      level = (iabs (m7[i]) * levelscale[j][i] + leveloffset[j][i]) >> q_bits;

      if (level != 0)
      {
        if (img->AdaptiveRounding)
        {
          fadjust8x8[j][block_x + i] = rshift_rnd_sf((AdaptRndWeight * (iabs (m7[i]) * levelscale[j][i] - (level << q_bits))), (q_bits + 1));
        }

        nonzero=TRUE;

        if (currMB->luma_transform_size_8x8_flag && input->symbol_mode == CAVLC)
        {
          *coeff_cost += (level > 1) ? MAX_VALUE : c_cost[runs[MCcoeff]];

          img->cofAC[b8+pl_off][MCcoeff][0][scan_poss[MCcoeff]  ] = isignab(level,m7[i]);
          img->cofAC[b8+pl_off][MCcoeff][1][scan_poss[MCcoeff]++] = runs[MCcoeff];          
          runs[MCcoeff]=-1;
        }
        else
        {
          *coeff_cost += (level > 1) ? MAX_VALUE : c_cost[run];
          ACLevel[scan_pos  ] = isignab(level,m7[i]);
          ACRun  [scan_pos++] = run;
          run=-1;                     // reset zero level counter
        }

        level = isignab(level, m7[i]);

        m7[i] = rshift_rnd_sf(level*invlevelscale[j][i]<<qp_per, 6); // dequantization
      }
      else
      {
        if (img->AdaptiveRounding)
          fadjust8x8[j][block_x + i] = 0;
        m7[i] = 0;
      }
    }

    if (!currMB->luma_transform_size_8x8_flag || input->symbol_mode != CAVLC)
      ACLevel[scan_pos] = 0;
    else
    {
      for(i=0; i<4; i++)
        img->cofAC[b8+pl_off][i][0][scan_poss[i]] = 0;
    }

    if (nonzero)
    {
      //    Inverse Transform
      inverse8x8(curr_res, curr_res, block_y, block_x);

      for( j=block_y; j<block_y + BLOCK_SIZE_8x8; j++)
      {
        img_Y = &img_enc[img->pix_y + j][img->pix_x + block_x];
        predY = &curr_mpr[j][block_x];
        m7 = &curr_res[j][block_x];
        for( i=0; i<BLOCK_SIZE_8x8; i++)
        {
          img_Y[i] = iClip1( max_imgpel_value, rshift_rnd_sf((m7[i]),DQ_BITS_8) + predY[i]);          
        }
      }
    }
    else // no transform coefficients
    {      
      for( j=block_y; j< block_y + BLOCK_SIZE_8x8; j++)
      {
        memcpy(&(img_enc[img->pix_y + j][img->pix_x + block_x]),&(curr_mpr[j][block_x]), BLOCK_SIZE_8x8 * sizeof(imgpel));
      }
    }
  }
  else
  {
    // Quant

    runs[0]=runs[1]=runs[2]=runs[3]=-1;
    scan_poss[0]=scan_poss[1]=scan_poss[2]=scan_poss[3]=0;

    for (coeff_ctr=0; coeff_ctr < 64; coeff_ctr++)
    {
      i=pos_scan[coeff_ctr][0];
      j=pos_scan[coeff_ctr][1];
      
      run++;
      ilev=0;

      if (currMB->luma_transform_size_8x8_flag && input->symbol_mode == CAVLC)
      {
        MCcoeff = MC(coeff_ctr);
        runs[MCcoeff]++;
      }

      m7 = &curr_res[block_y + j][block_x];
      level = iabs (m7[i]);

      if (img->AdaptiveRounding)
      {
        fadjust8x8[j][block_x+i] = 0;
      }

      if (level != 0)
      {
        nonzero = TRUE;

        if (currMB->luma_transform_size_8x8_flag && input->symbol_mode == CAVLC)
        {
          *coeff_cost += MAX_VALUE;

          img->cofAC[b8+pl_off][MCcoeff][0][scan_poss[MCcoeff]  ] = isignab(level,m7[i]);
          img->cofAC[b8+pl_off][MCcoeff][1][scan_poss[MCcoeff]++] = runs[MCcoeff];
          ++scan_pos;
          runs[MCcoeff]=-1;
        }
        else
        {
          *coeff_cost += MAX_VALUE;
          ACLevel[scan_pos  ] = isignab(level,m7[i]);
          ACRun  [scan_pos++] = run;
          run=-1;                     // reset zero level counter
        }

        level = isignab(level, m7[i]);
        ilev = level;
      }
    }
    if (!currMB->luma_transform_size_8x8_flag || input->symbol_mode != CAVLC)
      ACLevel[scan_pos] = 0;
    else
    {
      for(i=0; i<4; i++)
        img->cofAC[b8+pl_off][i][0][scan_poss[i]] = 0;
    }

   if(ipmode_DPCM<2)  //For residual DPCM
   {
   Inv_Residual_DPCM_8x8(pl, block_y, block_x);
   }

    for( j=block_y; j<block_y + BLOCK_SIZE_8x8; j++)
    {            
      for( i=block_x; i< block_x + BLOCK_SIZE_8x8; i++)
      {
        curr_res[j][i] = curr_res[j][i] + curr_mpr[j][i];
        img_enc[img->pix_y + j][img->pix_x + i]= (imgpel) curr_res[j][i];
      }
    }
  }

  //  Decoded block moved to frame memory
  return nonzero;
}

