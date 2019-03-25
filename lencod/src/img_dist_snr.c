
/*!
 *************************************************************************************
 * \file img_dist_snr.c
 *
 * \brief
 *    Compute signal to noise ratio (SNR) between the encoded image and the reference image
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Woo-Shik Kim                    <wooshik.kim@usc.edu>
 *     - Alexis Michael Tourapis         <alexismt@ieee.org>
 *************************************************************************************
 */
#include "contributors.h"

#include <math.h>

#include "enc_statistics.h"
#include "global.h"
#include "img_distortion.h"
#include "md_distortion.h"

/*!
 ************************************************************************
 * \brief
 *    Find SNR for all three components
 ************************************************************************
 */
void find_snr(VideoParameters *p_Vid, ImageStructure *imgREF, ImageStructure *imgSRC, DistMetric *metricSSE, DistMetric *metricPSNR)
{
  DistortionParams *p_Dist = p_Vid->p_Dist;
  FrameFormat *format = &imgREF->format;
  // Luma.
  metricSSE ->value[0] = (float) compute_SSE(imgREF->data[0], imgSRC->data[0], 0, 0, format->height, format->width);
  metricPSNR->value[0] = psnr(format->max_value_sq[0], format->size_cmp[0], metricSSE->value[0]);

  // Chroma.
  if (format->yuv_format != YUV400)
  {   
    metricSSE ->value[1] = (float) compute_SSE(imgREF->data[1], imgSRC->data[1], 0, 0, format->height_cr, format->width_cr);
    metricPSNR->value[1] = psnr(format->max_value_sq[1], format->size_cmp[1], metricSSE->value[1]);
    metricSSE ->value[2] = (float) compute_SSE(imgREF->data[2], imgSRC->data[2], 0, 0, format->height_cr, format->width_cr);
    metricPSNR->value[2] = psnr(format->max_value_sq[2], format->size_cmp[2], metricSSE->value[2]);
  }
   
  {
    accumulate_average(metricSSE,  p_Dist->frame_ctr);
    accumulate_average(metricPSNR, p_Dist->frame_ctr);

    accumulate_avslice(metricSSE,  p_Vid->type, p_Vid->p_Stats->frame_ctr[p_Vid->type]);
    accumulate_avslice(metricPSNR, p_Vid->type, p_Vid->p_Stats->frame_ctr[p_Vid->type]);
  }
}
