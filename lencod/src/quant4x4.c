
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
#include "global.h"
#include "quant4x4.h"

/*!
************************************************************************
* \brief
*    Quantization initialization function
*
************************************************************************
*/
void init_quant_4x4(InputParameters *params, ImageParameters *img)
{
  if (params->UseRDOQuant == 1)
  {
    quant_4x4     = quant_4x4_trellis;
    quant_dc4x4   = quant_dc4x4_trellis;
    quant_ac4x4   = quant_ac4x4_trellis;
  }
  else if (img->AdaptiveRounding)
  {
    quant_4x4     = quant_4x4_around;
    quant_dc4x4   = quant_dc4x4_normal;
    quant_ac4x4   = quant_ac4x4_around;
  }
  else
  {
    quant_4x4   = quant_4x4_normal;
    quant_dc4x4 = quant_dc4x4_normal;
    quant_ac4x4 = quant_ac4x4_normal;
  }
}

