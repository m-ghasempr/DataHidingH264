
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
void itrans8x8(ImageParameters *img, //!< image parameters
               Macroblock *currMB,   //!< current macroblock
               ColorPlane pl,        //!< used color plane       
               int ioff,             //!< index to 4x4 block
               int joff)             //!< index to 4x4 block
{
  int i,j;

  imgpel **mpr    = img->mb_pred[pl];
  imgpel **mb_rec = img->mb_rec[pl];
  int    **m7     = img->mb_rres[pl];
  int     max_imgpel_value = img->max_imgpel_value_comp[pl];

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
        mb_rec[j][i] = (imgpel) iClip1(max_imgpel_value, rshift_rnd_sf((m7[j][i] + ((long)mpr[j][i] << DQ_BITS_8)), DQ_BITS_8)); 
    }
  }
}
