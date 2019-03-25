
/*!
 *************************************************************************************
 * \file quant8x8_trellis.c
 *
 * \brief
 *    Quantization process for a 4x4 block using trellis based quantization
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *
 *************************************************************************************
 */

#include "contributors.h"

#include <math.h>

#include "global.h"

#include "image.h"
#include "mb_access.h"
#include "vlc.h"
#include "transform.h"
#include "mc_prediction.h"
#include "q_offsets.h"
#include "q_matrix.h"
#include "quant8x8.h"
#include "rdo_quant.h"

const int estErr8x8[6][8][8]={
  {
    {6553600, 6677056, 6400000, 6677056, 6553600, 6677056, 6400000, 6677056}, 
    {6677056, 6765201, 6658560, 6765201, 6677056, 6765201, 6658560, 6765201}, 
    {6400000, 6658560, 6553600, 6658560, 6400000, 6658560, 6553600, 6658560}, 
    {6677056, 6765201, 6658560, 6765201, 6677056, 6765201, 6658560, 6765201}, 
    {6553600, 6677056, 6400000, 6677056, 6553600, 6677056, 6400000, 6677056}, 
    {6677056, 6765201, 6658560, 6765201, 6677056, 6765201, 6658560, 6765201}, 
    {6400000, 6658560, 6553600, 6658560, 6400000, 6658560, 6553600, 6658560}, 
    {6677056, 6765201, 6658560, 6765201, 6677056, 6765201, 6658560, 6765201} 
  },
  {
    {7929856, 8156736, 8028160, 8156736, 7929856, 8156736, 8028160, 8156736}, 
    {8156736, 7537770, 7814560, 7537770, 8156736, 7537770, 7814560, 7537770}, 
    {8028160, 7814560, 7840000, 7814560, 8028160, 7814560, 7840000, 7814560}, 
    {8156736, 7537770, 7814560, 7537770, 8156736, 7537770, 7814560, 7537770}, 
    {7929856, 8156736, 8028160, 8156736, 7929856, 8156736, 8028160, 8156736}, 
    {8156736, 7537770, 7814560, 7537770, 8156736, 7537770, 7814560, 7537770}, 
    {8028160, 7814560, 7840000, 7814560, 8028160, 7814560, 7840000, 7814560}, 
    {8156736, 7537770, 7814560, 7537770, 8156736, 7537770, 7814560, 7537770} 
  },
  {
    {11075584, 10653696, 11151360, 10653696, 11075584, 10653696, 11151360, 10653696}, 
    {10653696, 11045652, 11109160, 11045652, 10653696, 11045652, 11109160, 11045652}, 
    {11151360, 11109160, 11289600, 11109160, 11151360, 11109160, 11289600, 11109160}, 
    {10653696, 11045652, 11109160, 11045652, 10653696, 11045652, 11109160, 11045652}, 
    {11075584, 10653696, 11151360, 10653696, 11075584, 10653696, 11151360, 10653696}, 
    {10653696, 11045652, 11109160, 11045652, 10653696, 11045652, 11109160, 11045652}, 
    {11151360, 11109160, 11289600, 11109160, 11151360, 11109160, 11289600, 11109160}, 
    {10653696, 11045652, 11109160, 11045652, 10653696, 11045652, 11109160, 11045652} 
  },
  {
    {12845056, 12503296, 12544000, 12503296, 12845056, 12503296, 12544000, 12503296}, 
    {12503296, 13050156, 12588840, 13050156, 12503296, 13050156, 12588840, 13050156}, 
    {12544000, 12588840, 12960000, 12588840, 12544000, 12588840, 12960000, 12588840}, 
    {12503296, 13050156, 12588840, 13050156, 12503296, 13050156, 12588840, 13050156}, 
    {12845056, 12503296, 12544000, 12503296, 12845056, 12503296, 12544000, 12503296}, 
    {12503296, 13050156, 12588840, 13050156, 12503296, 13050156, 12588840, 13050156}, 
    {12544000, 12588840, 12960000, 12588840, 12544000, 12588840, 12960000, 12588840}, 
    {12503296, 13050156, 12588840, 13050156, 12503296, 13050156, 12588840, 13050156} 
  },
  {
    {16777216, 16646400, 16384000, 16646400, 16777216, 16646400, 16384000, 16646400}, 
    {16646400, 16370116, 16692640, 16370116, 16646400, 16370116, 16692640, 16370116}, 
    {16384000, 16692640, 16646400, 16692640, 16384000, 16692640, 16646400, 16692640}, 
    {16646400, 16370116, 16692640, 16370116, 16646400, 16370116, 16692640, 16370116}, 
    {16777216, 16646400, 16384000, 16646400, 16777216, 16646400, 16384000, 16646400}, 
    {16646400, 16370116, 16692640, 16370116, 16646400, 16370116, 16692640, 16370116}, 
    {16384000, 16692640, 16646400, 16692640, 16384000, 16692640, 16646400, 16692640}, 
    {16646400, 16370116, 16692640, 16370116, 16646400, 16370116, 16692640, 16370116} 
  },
  {
    {21233664, 21381376, 21667840, 21381376, 21233664, 21381376, 21667840, 21381376}, 
    {21381376, 21381376, 21374440, 21381376, 21381376, 21381376, 21374440, 21381376}, 
    {21667840, 21374440, 21529600, 21374440, 21667840, 21374440, 21529600, 21374440}, 
    {21381376, 21381376, 21374440, 21381376, 21381376, 21381376, 21374440, 21381376}, 
    {21233664, 21381376, 21667840, 21381376, 21233664, 21381376, 21667840, 21381376}, 
    {21381376, 21381376, 21374440, 21381376, 21381376, 21381376, 21374440, 21381376}, 
    {21667840, 21374440, 21529600, 21374440, 21667840, 21374440, 21529600, 21374440}, 
    {21381376, 21381376, 21374440, 21381376, 21381376, 21381376, 21374440, 21381376} 
  }
};

void rdoq_8x8(int (*tblock)[16], int block_y, int block_x,int qp_per, int qp_rem, 
              int **levelscale, int **leveloffset, const byte (*pos_scan)[2], int levelTrellis[64]);
/*!
 ************************************************************************
 * \brief
 *    Quantization process for All coefficients for a 8x8 block
 *
 * \par Input:
 *
 * \par Output:
 *
 ************************************************************************
 */
int quant_8x8_trellis(int (*tblock)[16], int block_y, int block_x, int  qp,                
                      int*  ACLevel, int*  ACRun, 
                      int **fadjust8x8, int **levelscale, int **invlevelscale, int **leveloffset,
                      int *coeff_cost, const byte (*pos_scan)[2], const byte *c_cost)
{
  static int i,j, coeff_ctr;

  static int *m7;

  int   level, run = 0;
  int   nonzero = FALSE;
  int   qp_per = qp_per_matrix[qp];
  int   qp_rem = qp_rem_matrix[qp];
  const byte *p_scan = &pos_scan[0][0];
  int*  ACL = &ACLevel[0];
  int*  ACR = &ACRun[0];

  int levelTrellis[64];

  rdoq_8x8(tblock,block_y,block_x,qp_per,qp_rem,levelscale,leveloffset,pos_scan, levelTrellis);

  // Quantization
  for (coeff_ctr = 0; coeff_ctr < 64; coeff_ctr++)
  {
    i = *p_scan++;  // horizontal position
    j = *p_scan++;  // vertical position

    m7 = &tblock[j][block_x + i];
    /*
    scaled_coeff = iabs (*m7) * levelscale[j][i];
    level = (scaled_coeff + leveloffset[j][i]) >> q_bits;
    */
    level = levelTrellis[coeff_ctr];

    if (level != 0)
    {

      if (params->symbol_mode == CAVLC && img->qp < 10)
        level = imin(level, CAVLC_LEVEL_LIMIT);

      nonzero = TRUE;

      *coeff_cost += (level > 1) ? MAX_VALUE : c_cost[run];

      level  = isignab(level, *m7);
      *m7    = rshift_rnd_sf(((level * invlevelscale[j][i]) << qp_per), 6);
      *ACL++ = level;
      *ACR++ = run; 
      // reset zero level counter
      run    = 0;
    }
    else
    {
      run++;
      *m7 = 0;
    }      
  }

  *ACL = 0;

  return nonzero;
}

/*!
 ************************************************************************
 * \brief
 *    Quantization process for All coefficients for a 8x8 block 
 *    CAVLC version
 *
 * \par Input:
 *
 * \par Output:
 *
 ************************************************************************
 */
int quant_8x8cavlc_trellis(int (*tblock)[16], int block_y, int block_x, int  qp,                 
                   int***  cofAC, 
                   int **fadjust8x8, int **levelscale, int **invlevelscale, int **leveloffset,
                   int *coeff_cost, const byte (*pos_scan)[2], const byte *c_cost)
{
  static int i,j, k, coeff_ctr;

  static int *m7;
  static int scaled_coeff;  

  int level, runs[4] = { 0 };
  int nonzero = FALSE; 
  int qp_per = qp_per_matrix[qp];  
  int q_bits = Q_BITS_8 + qp_per;
  const byte *p_scan = &pos_scan[0][0];
  int*  ACL[4];  
  int*  ACR[4];

  for (k = 0; k < 4; k++)
  {
    ACL[k] = &cofAC[k][0][0];
    ACR[k] = &cofAC[k][1][0];
  }

  // Quantization
  for (coeff_ctr = 0; coeff_ctr < 16; coeff_ctr++)
  {
    for (k = 0; k < 4; k++)
    {
      i = *p_scan++;  // horizontal position
      j = *p_scan++;  // vertical position

      m7 = &tblock[j][block_x + i];
      scaled_coeff = iabs (*m7) * levelscale[j][i];
      level = (scaled_coeff + leveloffset[j][i]) >> q_bits;

      if (level != 0)
      {
        nonzero=TRUE;

        *coeff_cost += (level > 1) ? MAX_VALUE : c_cost[runs[k]];

        level  = isignab(level, *m7);
        *m7    = rshift_rnd_sf(((level * invlevelscale[j][i]) << qp_per), 6);

        *(ACL[k])++ = level;
        *(ACR[k])++ = runs[k];
        // reset zero level counter
        runs[k] = 0;
      }
      else
      {        
        runs[k]++;
        *m7 = 0;      
      }
    }
  }

  for(k = 0; k < 4; k++)
    *(ACL[k]) = 0;

  return nonzero;
}

/*!
************************************************************************
* \brief
*    Rate distortion optimized Quantization process for 
*    all coefficients in a 8x8 block
*
************************************************************************
*/
void rdoq_8x8(int (*tblock)[16], int block_y, int block_x,int qp_per, int qp_rem, 
              int **levelscale, int **leveloffset, const byte (*pos_scan)[2], int levelTrellis[64])
{
  static int i,j, coeff_ctr;

  static int *m7;

  int level;

  int q_bits = Q_BITS_8 + qp_per;

  int m4[8][8], ii;
  levelDataStruct levelData[64];
  double  lambda_md=0, normFact=pow(2,(2 * Q_BITS_8 + 9));
  double err;
  int lowerInt, kStart=0, kStop=0, noCoeff = 0;

  if ((img->type==B_SLICE) && img->nal_reference_idc)
  {
    lambda_md = img->lambda_md[5][img->masterQP];  
  }
  else
  {
    lambda_md = img->lambda_md[img->type][img->masterQP]; 
  }

  for (coeff_ctr = 0; coeff_ctr < 64; coeff_ctr++)
  {
    i = pos_scan[coeff_ctr][0];
    j = pos_scan[coeff_ctr][1];

    m7 = &tblock[j][block_x];
    m4[j][i] = m7[i];
  } 

  for (coeff_ctr=0; coeff_ctr < 64; coeff_ctr++)
  {
    i=pos_scan[coeff_ctr][0];
    j=pos_scan[coeff_ctr][1];

    levelData[coeff_ctr].levelDouble = iabs(m4[j][i] * levelscale[i][j]);
    level=(levelData[coeff_ctr].levelDouble>>q_bits);

    lowerInt=(((int)levelData[coeff_ctr].levelDouble-(level<<q_bits))<(1<<(q_bits-1)))? 1 : 0;

    levelData[coeff_ctr].level[0]=0;
    if (level==0 && lowerInt==1)
    {
      levelData[coeff_ctr].noLevels=1;
    }
    else if (level==0 && lowerInt==0)
    {
      levelData[coeff_ctr].level[1] = level+1;
      levelData[coeff_ctr].noLevels=2;
      kStop=coeff_ctr;
      noCoeff++;
    }
    else if (level>0 && lowerInt==1)
    {
      levelData[coeff_ctr].level[1] = level;
      levelData[coeff_ctr].noLevels=2;
      kStop=coeff_ctr;
      noCoeff++;
    }
    else
    {
      levelData[coeff_ctr].level[1] = level;
      levelData[coeff_ctr].level[2] = level+1;
      levelData[coeff_ctr].noLevels=3;
      kStop=coeff_ctr;
      kStart=coeff_ctr;
      noCoeff++;
    }

    for (ii=0; ii<levelData[coeff_ctr].noLevels; ii++)
    {
      err=(double)(levelData[coeff_ctr].level[ii]<<q_bits)-(double)levelData[coeff_ctr].levelDouble;
      levelData[coeff_ctr].errLevel[ii] = (err * err * (double) estErr8x8[qp_rem][i][j]) / normFact; 
    }
  }

  est_writeRunLevel_CABAC(levelData, levelTrellis, LUMA_8x8, lambda_md, kStart, kStop, noCoeff, 0);
}


