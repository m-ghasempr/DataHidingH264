
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
 *
 * \date
 *    12. October 2003
 **************************************************************************
 */

#include "global.h"

#include "image.h"
#include "mb_access.h"
#include "elements.h"
#include "transform8x8.h"
#include "transform.h"
#include "quant.h"



/*!
 ***********************************************************************
 * \brief
 *    Inverse 8x8 transformation
 ***********************************************************************
 */ 
void itrans8x8(Macroblock *currMB,   //!< current macroblock
               ColorPlane pl,        //!< used color plane       
               int ioff,             //!< index to 4x4 block
               int joff)             //!< index to 4x4 block
{
  VideoParameters *p_Vid = currMB->p_Vid;
  Slice *currSlice = currMB->p_Slice;
  int i,j;

  imgpel **mpr    = currSlice->mb_pred[pl];
  imgpel **mb_rec = currSlice->mb_rec[pl];
  int    **m7     = currSlice->mb_rres[pl];
  int     max_imgpel_value = p_Vid->max_pel_value_comp[pl];

  if (currMB->is_lossless == TRUE)
  {
    for( j = joff; j < joff + 8; j++)
    {
      for( i = ioff; i < ioff + 8; i++)
        mb_rec[j][i] = (imgpel) iClip1(max_imgpel_value, (m7[j][i] + (long)mpr[j][i])); 
    }
  }
  else
  {
    inverse8x8(m7, m7, joff, ioff);
    for( j = joff; j < joff + 8; j++)
    {
      for( i = ioff; i < ioff + 8; i++)
        //mb_rec[j][i] = (imgpel) iClip1(max_imgpel_value, rshift_rnd_sf((m7[j][i] + ((long)mpr[j][i] << DQ_BITS_8)), DQ_BITS_8)); 
        mb_rec[j][i] = (imgpel) iClip1(max_imgpel_value, mpr[j][i] + rshift_rnd_sf(m7[j][i], DQ_BITS_8)); 
    }
  }
}
