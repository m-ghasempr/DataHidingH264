
/*!
 ***************************************************************************
 * \file
 *    cabac.h
 *
 * \brief
 *    Headerfile for entropy coding routines
 *
 * \author
 *    Detlev Marpe                                                         \n
 *    Copyright (C) 2000 HEINRICH HERTZ INSTITUTE All Rights Reserved.
 *
 * \date
 *    21. Oct 2000 (Changes by Tobias Oelbaum 28.08.2001)
 ***************************************************************************
 */


#ifndef _CABAC_H_
#define _CABAC_H_

#include "global.h"




/***********************************************************************
 * L O C A L L Y   D E F I N E D   F U N C T I O N   P R O T O T Y P E S
 ***********************************************************************
 */



void unary_bin_encode(EncodingEnvironmentPtr eep_frame,
                      unsigned int symbol,
                      BiContextTypePtr ctx,
                      int ctx_offset);

void unary_bin_max_encode(EncodingEnvironmentPtr eep_frame,
                          unsigned int symbol,
                          BiContextTypePtr ctx,
                          int ctx_offset,
                          unsigned int max_symbol);

void unary_exp_golomb_level_encode( EncodingEnvironmentPtr eep_dp,
                                    unsigned int symbol,
                                    BiContextTypePtr ctx);

void unary_exp_golomb_mv_encode(EncodingEnvironmentPtr eep_dp,
                                unsigned int symbol,
                                BiContextTypePtr ctx,
                                unsigned int max_bin);

#endif  // CABAC_H

