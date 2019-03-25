
/*!
 *************************************************************************************
 * \file quant4x4.c
 *
 * \brief
 *    Quantization process for a 4x4 block
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Alexis Michael Tourapis                  <alexismt@ieee.org>
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
int quant_4x4_normal(int (*tblock)[16], int block_y, int block_x, int  qp,
                     int*  ACLevel, int*  ACRun, 
                     int **fadjust4x4, int **levelscale, int **invlevelscale, int **leveloffset,
                     int *coeff_cost, const byte (*pos_scan)[2], const byte *c_cost)
{
  static int i,j, coeff_ctr;

  static int *m7;
  static int scaled_coeff;

  int   level, run = 0;
  int   nonzero = FALSE;  
  int   qp_per = qp_per_matrix[qp];
  int   q_bits = Q_BITS + qp_per;
  const byte *p_scan = &pos_scan[0][0];
  int*  ACL = &ACLevel[0];
  int*  ACR = &ACRun[0];

  // Quantization
  for (coeff_ctr = 0; coeff_ctr < 16; coeff_ctr++)
  {
    i = *p_scan++;  // horizontal position
    j = *p_scan++;  // vertical position

    m7 = &tblock[j][block_x + i];
    scaled_coeff = iabs (*m7) * levelscale[j][i];
    level = (scaled_coeff + leveloffset[j][i]) >> q_bits;

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

int quant_ac4x4_normal(int (*tblock)[16], int block_y, int block_x, int qp,                 
                       int*  ACLevel, int*  ACRun, 
                       int **fadjust4x4, int **levelscale, int **invlevelscale, int **leveloffset,
                       int *coeff_cost, const byte (*pos_scan)[2], const byte *c_cost)
{
  static int i,j, coeff_ctr;

  static int *m7;
  static int scaled_coeff;

  int   level, run = 0;
  int   nonzero = FALSE;  
  int   qp_per = qp_per_matrix[qp];
  int   q_bits = Q_BITS + qp_per;
  const byte *p_scan = &pos_scan[1][0];
  int*  ACL = &ACLevel[0];
  int*  ACR = &ACRun[0];

  // Quantization
  for (coeff_ctr = 1; coeff_ctr < 16; coeff_ctr++)
  {
    i = *p_scan++;  // horizontal position
    j = *p_scan++;  // vertical position

    m7 = &tblock[j][block_x + i];
    scaled_coeff = iabs (*m7) * levelscale[j][i];
    level = (scaled_coeff + leveloffset[j][i]) >> q_bits;

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
 
