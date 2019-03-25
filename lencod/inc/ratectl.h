
/*!
 ***************************************************************************
 * \file
 *    ratectl.h
 *
 * \author
 *    Zhengguo LI
 *
 * \date
 *    14 Jan 2003
 *
 * \brief
 *    Headerfile for rate control
 **************************************************************************
 */

#ifndef _RATE_CTL_H_
#define _RATE_CTL_H_

#include "global.h"
#include "rc_quadratic.h"


// generic functions
extern int    Qstep2QP          ( double Qstep, int qp_offset );
extern double QP2Qstep          ( int QP );
extern int    ComputeMBMAD      ( int diff[16][16] );
extern double ComputeFrameMAD   ( ImageParameters *p_Img );
extern void   rc_store_mad      ( Macroblock *currMB );

// rate control functions
// init/copy
extern void  rc_alloc_generic           ( ImageParameters *p_Img, RCGeneric **p_quad );
extern void  rc_free_generic            ( RCGeneric **p_quad );
extern void  rc_copy_generic            ( ImageParameters *p_Img, RCGeneric *dst, RCGeneric *src );
extern void  rc_init_gop_params         ( ImageParameters *p_Img, InputParameters *p_Inp );
extern void  rc_init_frame              ( ImageParameters *p_Img, InputParameters *p_Inp);
extern void  rc_init_sequence           ( ImageParameters *p_Img, InputParameters *p_Inp);
extern void  rc_store_slice_header_bits ( ImageParameters *p_Img, InputParameters *p_Inp, int len);

#endif

