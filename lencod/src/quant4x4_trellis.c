
/*!
 *************************************************************************************
 * \file quant4x4_trellis.c
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
#include "quant4x4.h"
#include "rdo_quant.h"

const int estErr4x4[6][4][4]=
{
  {
    {25600, 27040, 25600, 27040}, 
    {27040, 25600, 27040, 25600}, 
    {25600, 27040, 25600, 27040}, 
    {27040, 25600, 27040, 25600} 
  },
  {
    {30976, 31360, 30976, 31360}, 
    {31360, 32400, 31360, 32400}, 
    {30976, 31360, 30976, 31360}, 
    {31360, 32400, 31360, 32400} 
  },
  {
    {43264, 40960, 43264, 40960}, 
    {40960, 40000, 40960, 40000}, 
    {43264, 40960, 43264, 40960}, 
    {40960, 40000, 40960, 40000} 
  },
  {
    {50176, 51840, 50176, 51840}, 
    {51840, 52900, 51840, 52900}, 
    {50176, 51840, 50176, 51840}, 
    {51840, 52900, 51840, 52900} 
  },
  {
    {65536, 64000, 65536, 64000}, 
    {64000, 62500, 64000, 62500}, 
    {65536, 64000, 65536, 64000}, 
    {64000, 62500, 64000, 62500} 
  },
  {
    {82944, 84640, 82944, 84640}, 
    {84640, 84100, 84640, 84100}, 
    {82944, 84640, 82944, 84640}, 
    {84640, 84100, 84640, 84100} 
  }
};


void rdoq_4x4(int (*tblock)[16], int block_y, int block_x,int qp_per, int qp_rem, 
              int **levelscale, int **leveloffset, const byte (*pos_scan)[2], int levelTrellis[16]);
void rdoq_ac4x4(int (*tblock)[16] , int block_x, int qp_per, int qp_rem, 
                int **levelscale, int **leveloffset, const byte (*pos_scan)[2], int levelTrellis[16]);

/*!
 ************************************************************************
 * \brief
 *    Quantization process for All coefficients for a 4x4 block
 *
 * \par Input:
 *
 * \par Output:
 *
 ************************************************************************
 */
int quant_4x4_trellis(int (*tblock)[16], int block_y, int block_x, int  qp,                
                      int*  ACLevel, int*  ACRun, 
                      int **fadjust4x4, int **levelscale, int **invlevelscale, int **leveloffset,
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
  
  int levelTrellis[16];

  rdoq_4x4(tblock,block_y,block_x,qp_per,qp_rem,levelscale,leveloffset,pos_scan, levelTrellis);

  // Quantization
  for (coeff_ctr = 0; coeff_ctr < 16; coeff_ctr++)
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
      *m7    = rshift_rnd_sf(((level * invlevelscale[j][i]) << qp_per), 4);
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
*    Rate distortion optimized Quantization process for 
*    all coefficients in a 4x4 block
*
************************************************************************
*/
void rdoq_4x4(int (*tblock)[16], int block_y, int block_x,int qp_per, int qp_rem, 
              int **levelscale, int **leveloffset, const byte (*pos_scan)[2], int levelTrellis[16])
{
  static int i,j, coeff_ctr;

  static int *m7;

  int level;

  int q_bits = Q_BITS + qp_per;

  int m4[4][4], ii;
  levelDataStruct levelData[16];
  double  lambda_md=0, normFact=pow(2,(2*DQ_BITS+19));
  double err;
  int lowerInt, kStart=0, kStop=0, noCoeff = 0, estBits;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  if ((img->type==B_SLICE) && img->nal_reference_idc)
  {
    lambda_md = img->lambda_md[5][img->masterQP];  
  }
  else
  {
    lambda_md = img->lambda_md[img->type][img->masterQP]; 
  }

  for (coeff_ctr = 0; coeff_ctr < 16; coeff_ctr++)
  {
    i = pos_scan[coeff_ctr][0];
    j = pos_scan[coeff_ctr][1];

    m7 = &tblock[j][block_x];
    m4[j][i] = m7[i];
  } 
  
  for (coeff_ctr=0;coeff_ctr < 16; coeff_ctr++)
  {
    i=pos_scan[coeff_ctr][0];
    j=pos_scan[coeff_ctr][1];

    levelData[coeff_ctr].levelDouble = iabs(m4[j][i] * levelscale[i][j]);
    level=(levelData[coeff_ctr].levelDouble>>q_bits);

    lowerInt=(((int)levelData[coeff_ctr].levelDouble - (level<<q_bits))<(1<<(q_bits-1)))? 1 : 0;

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
      levelData[coeff_ctr].errLevel[ii] = (err * err * (double) estErr4x4[qp_rem][i][j]) / normFact; 
    }
  }

  estBits=est_write_and_store_CBP_block_bit(currMB, LUMA_4x4, block_y, block_x);
  est_writeRunLevel_CABAC(levelData, levelTrellis, LUMA_4x4, lambda_md, kStart, kStop, noCoeff, estBits);
}

int quant_ac4x4_trellis(int (*tblock)[16], int block_y, int block_x, int qp,                
                        int*  ACLevel, int*  ACRun, 
                        int **fadjust4x4, int **levelscale, int **invlevelscale, int **leveloffset,
                        int *coeff_cost, const byte (*pos_scan)[2], const byte *c_cost)
{
  static int i,j, coeff_ctr;

  static int *m7;
  int   level, run = 0;
  int   nonzero = FALSE;  
  int   qp_per = qp_per_matrix[qp];
  int   qp_rem = qp_rem_matrix[qp];
  const byte *p_scan = &pos_scan[1][0];
  int*  ACL = &ACLevel[0];
  int*  ACR = &ACRun[0];

  int levelTrellis[16]; 

  rdoq_ac4x4(tblock, block_x, qp_per, qp_rem, levelscale, leveloffset, pos_scan, levelTrellis);

  // Quantization
  for (coeff_ctr = 1; coeff_ctr < 16; coeff_ctr++)
  {
    i = *p_scan++;  // horizontal position
    j = *p_scan++;  // vertical position

    m7 = &tblock[j][block_x + i];
    /*
    scaled_coeff = iabs (*m7) * levelscale[j][i];
    level = (scaled_coeff + leveloffset[j][i]) >> q_bits;
    */
    level=levelTrellis[coeff_ctr-1];

    if (level != 0)
    {
      if (params->symbol_mode == CAVLC && img->qp < 10)
        level = imin(level, CAVLC_LEVEL_LIMIT);

      nonzero = TRUE;

      *coeff_cost += (level > 1) ? MAX_VALUE : c_cost[run];

      level  = isignab(level, *m7);
      *m7    = rshift_rnd_sf(((level * invlevelscale[j][i]) << qp_per), 4);
      // inverse scale can be alternative performed as follows to ensure 16bit
      // arithmetic is satisfied.
      // *m7 = (qp_per<4) ? rshift_rnd_sf((level*invlevelscale[j][i]),4-qp_per) : (level*invlevelscale[j][i])<<(qp_per-4);
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

void rdoq_ac4x4(int (*tblock)[16] , int block_x, int qp_per, int qp_rem, 
                int **levelscale, int **leveloffset, const byte (*pos_scan)[2], int levelTrellis[16])
{
  static int i,j, coeff_ctr;

  static int *m7;

  int level;

  int q_bits = Q_BITS + qp_per;
  const byte *p_scan = &pos_scan[1][0];

  int k;
  levelDataStruct levelData[16];
  double  lambda_md=0, normFact=pow(2,(2*DQ_BITS+19));
  double err;
  int lowerInt, kStart=0, kStop=0, noCoeff=0, estBits;
  Macroblock *currMB = &img->mb_data[img->current_mb_nr];

  if ((img->type==B_SLICE) && img->nal_reference_idc)
  {
    lambda_md = img->lambda_md[5][img->masterQP];  
  }
  else
  {
    lambda_md = img->lambda_md[img->type][img->masterQP]; 
  }

  for (coeff_ctr=0;coeff_ctr < 15;coeff_ctr++)
  {
   // i=pos_scan[coeff_ctr+1][0]; // scan is shifted due to DC
    //j=pos_scan[coeff_ctr+1][1]; // scan is shifted due to DC
    i = *p_scan++;  // horizontal position
    j = *p_scan++;  // vertical position
    m7 = &tblock[j][block_x + i];


    levelData[coeff_ctr].levelDouble=iabs(*m7 * levelscale[i][j]);
    level=(levelData[coeff_ctr].levelDouble>>q_bits);

    lowerInt=(((int)levelData[coeff_ctr].levelDouble - (level<<q_bits))<(1<<(q_bits-1)))? 1 : 0;

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

    for (k=0; k<levelData[coeff_ctr].noLevels; k++)
    {
      err=(double)(levelData[coeff_ctr].level[k]<<q_bits)-(double)levelData[coeff_ctr].levelDouble;
      levelData[coeff_ctr].errLevel[k]=(err*err*(double)estErr4x4[qp_rem][i][j])/normFact; 
    }
  }
  estBits=est_write_and_store_CBP_block_bit(currMB, LUMA_16AC, 0, 0);

  est_writeRunLevel_CABAC(levelData, levelTrellis, LUMA_16AC, lambda_md, kStart, kStop, noCoeff, estBits);
}



