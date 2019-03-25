
/*!
 ***************************************************************************
 * \file md_distortion.c
 *
 * \brief
 *    Main macroblock mode decision functions and helpers
 *
 **************************************************************************
 */

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <float.h>
#include <memory.h>
#include <string.h>

#include "global.h"
#include "rdopt_coding_state.h"
#include "mb_access.h"
#include "intrarefresh.h"
#include "image.h"
#include "transform8x8.h"
#include "ratectl.h"
#include "mode_decision.h"
#include "fmo.h"
#include "me_umhex.h"
#include "me_umhexsmp.h"
#include "macroblock.h"

/*!
*************************************************************************************
* \brief
*    SSE distortion calculation for a macroblock
*************************************************************************************
*/
int64 distortionSSE(Macroblock *currMB) 
{
  int i, j, k;
  int64 distortionY  = 0;
  int64 distortionCr[2] = {0};
  imgpel *imgOrg, *imgEnc;

  // LUMA
  for (j = 0; j < MB_BLOCK_SIZE; j++)
  {
    imgOrg = &pCurImg[j + img->opix_y][img->opix_x];    
    imgEnc = &enc_picture->p_curr_img[j + img->pix_y][img->pix_x];

    for (i = 0; i < MB_BLOCK_SIZE; i++)
      distortionY += iabs2( *imgOrg++ - *imgEnc++ );
  }

  if (img->yuv_format != YUV400)
  {
    // CHROMA
    for (k = 0; k < 2; k++)
    {
      for (j = 0; j < img->mb_cr_size_y; j++)
      {
        imgOrg = &imgUV_org[k][j + img->opix_c_y][img->opix_c_x];
        imgEnc = &enc_picture->imgUV[k][j + img->pix_c_y][img->pix_c_x];        

        for (i=0; i<img->mb_cr_size_x; i++)
          distortionCr[k] += iabs2( *imgOrg++ - *imgEnc++ );
      }      
    }
  }

  return (int64)( distortionY * input->WeightY + distortionCr[0] * input->WeightCb + distortionCr[1] * input->WeightCr );
}

