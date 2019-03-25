
/*!
 *************************************************************************************
 * \file quantChroma_trellis.c
 *
 * \brief
 *    Quantization process for a Chroma block (trellis based)
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
#include "q_matrix.h"
#include "quant4x4.h"
#include "quantChroma.h"
#include "rdo_quant.h"

/*!
 ************************************************************************
 * \brief
 *    Quantization process for All coefficients for a 2x2 DC block
 *
 * \par Input:
 *
 * \par Output:
 *
 ************************************************************************
 */
int quant_dc2x2_trellis(int (*tblock)[4], int qp, int* DCLevel, int* DCRun, 
                       int **fadjust, int levelscale, int invlevelscale, int **leveloffset,
                       const byte (*pos_scan)[2])
{
  static int i,j, coeff_ctr;

  static int *m7;

  int   level, run = 0;
  int   nonzero = FALSE;  
  int   qp_per = qp_per_matrix[qp];
  int   qp_rem = qp_rem_matrix[qp]; 
  const byte *p_scan = &pos_scan[0][0];
  int*  DCL = &DCLevel[0];
  int*  DCR = &DCRun[0];

  static int levelTrellis[16];

  rdoq_dc(tblock,qp_per,qp_rem, levelscale, leveloffset,pos_scan, levelTrellis, CHROMA_DC);

  for (coeff_ctr=0; coeff_ctr < 4; coeff_ctr++)
  {
    j = *p_scan++;  // note that in this part, coefficients were previously transposed from 2x4 to 4x2.
    i = *p_scan++;  

    m7 = &tblock[j][i];

    if (*m7 != 0)
    {
      level = levelTrellis[coeff_ctr];

      if (level  != 0)
      {
        if (params->symbol_mode == CAVLC)
          level = imin(level, CAVLC_LEVEL_LIMIT);

        level = isignab(level, *m7);

        *m7 = ((level * invlevelscale) << qp_per) >> 5;

        *DCL++  = level;
        *DCR++  = run;
        // reset zero level counter
        run     = 0;
        nonzero = TRUE;
      }
      else
      {
        run++;
        *m7 = 0;
      }
    }
    else
    {
      run++;
    }
  }

  *DCL = 0;

  return nonzero;
}

/*!
 ************************************************************************
 * \brief
 *    Quantization process for All coefficients for a 2x2 DC block
 *
 * \par Input:
 *
 * \par Output:
 *
 ************************************************************************
 */
int quant_dc4x2_trellis(int (*tblock)[4], int qp, int* DCLevel, int* DCRun, 
                       int **fadjust, int levelscale, int invlevelscale, int **leveloffset,
                       const byte (*pos_scan)[2])
{
  static int i,j, coeff_ctr;

  static int *m7;

  int   level, run = 0;
  int   nonzero = FALSE;  
  int   qp_per = qp_per_matrix[qp];
  int   qp_rem = qp_rem_matrix[qp]; 
  const byte *p_scan = &pos_scan[0][0];
  int*  DCL = &DCLevel[0];
  int*  DCR = &DCRun[0];

  static int levelTrellis[16];

  rdoq_dc(tblock,qp_per,qp_rem, levelscale, leveloffset,pos_scan, levelTrellis, CHROMA_DC_2x4);

  for (coeff_ctr=0; coeff_ctr < 4; coeff_ctr++)
  {
    j = *p_scan++;  // note that in this part, somehow coefficients were transposed from 2x4 to 4x2.
    i = *p_scan++;  

    m7 = &tblock[j][i];

    if (*m7 != 0)
    {
      level = levelTrellis[coeff_ctr];

      if (level  != 0)
      {
        if (params->symbol_mode == CAVLC)
          level = imin(level, CAVLC_LEVEL_LIMIT);
        level = isignab(level, *m7);

        *m7 = rshift_rnd_sf(((level * invlevelscale) << qp_per), 6);

        *DCL++  = level;
        *DCR++  = run;
        // reset zero level counter
        run     = 0;
        nonzero = TRUE;
      }
      else
      {
        run++;
        *m7 = 0;
      }
    }
    else
    {
      run++;
    }
  }

  *DCL = 0;

  return nonzero;
}


