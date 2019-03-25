
/*!
 ***************************************************************************
 * \file
 *    rc_quadratic.h
 *
 * \author
 *    Zhengguo LI
 *    Athanasios Leontaris
 *
 * \date
 *    14 Jan 2003
 *
 * \brief
 *    Headerfile for rate control
 **************************************************************************
 */

#ifndef _RC_QUADRATIC_H_
#define _RC_QUADRATIC_H_

#include "rc_types.h"

// rate control functions
// init/copy
extern void rc_alloc_quadratic( ImageParameters *p_Img, InputParameters *p_Inp, RCQuadratic **p_quad );
extern void rc_free_quadratic ( RCQuadratic **p_quad );
extern void rc_copy_quadratic ( ImageParameters *p_Img, InputParameters *p_Inp, RCQuadratic *dst, RCQuadratic *src );

// rate control (externally visible)
extern void rc_init_seq          (ImageParameters *p_Img, InputParameters *p_Inp, RCQuadratic *p_quad, RCGeneric *p_gen);
extern void rc_init_GOP          (ImageParameters *p_Img, InputParameters *p_Inp, RCQuadratic *p_quad, RCGeneric *p_gen, int np, int nb);
extern void rc_update_pict_frame (ImageParameters *p_Img, InputParameters *p_Inp, RCQuadratic *p_quad, RCGeneric *p_gen, int nbits);
extern void rc_init_pict         (ImageParameters *p_Img, InputParameters *p_Inp, 
                           RCQuadratic *p_quad, RCGeneric *p_gen, int fieldpic, int topfield, int targetcomputation, float mult);
extern void rc_update_pict       (ImageParameters *p_Img, InputParameters *p_Inp, RCQuadratic *p_quad, RCGeneric *p_gen, int nbits);
extern void rc_update_picture    (ImageParameters *p_Img, InputParameters *p_Inp, int bits);

extern int  updateQPRC0(ImageParameters *p_Img, InputParameters *p_Inp, RCQuadratic *p_quad, RCGeneric *p_gen, int topfield);
extern int  updateQPRC1(ImageParameters *p_Img, InputParameters *p_Inp, RCQuadratic *p_quad, RCGeneric *p_gen, int topfield);
extern int  updateQPRC2(ImageParameters *p_Img, InputParameters *p_Inp, RCQuadratic *p_quad, RCGeneric *p_gen, int topfield);
extern int  updateQPRC3(ImageParameters *p_Img, InputParameters *p_Inp, RCQuadratic *p_quad, RCGeneric *p_gen, int topfield);

// internal functions
extern void updateQPInterlace   ( RCQuadratic *p_quad, RCGeneric *p_gen );
extern void updateQPNonPicAFF   ( seq_parameter_set_rbsp_t *active_sps, RCQuadratic *p_quad );
extern void updateBottomField   ( InputParameters *p_Inp, RCQuadratic *p_quad );
extern int  updateFirstP        ( ImageParameters *p_Img, InputParameters *p_Inp, RCQuadratic *p_quad, RCGeneric *p_gen, int topfield );
extern int  updateNegativeTarget( ImageParameters *p_Img, InputParameters *p_Inp, RCQuadratic *p_quad, RCGeneric *p_gen, int topfield, int m_Qp );
extern int  updateFirstBU       ( ImageParameters *p_Img, InputParameters *p_Inp, RCQuadratic *p_quad, RCGeneric *p_gen, int topfield );
extern void updateLastBU        ( ImageParameters *p_Img, InputParameters *p_Inp, RCQuadratic *p_quad, RCGeneric *p_gen, int topfield );
extern void predictCurrPicMAD   ( InputParameters *p_Inp, RCQuadratic *p_quad, RCGeneric *p_gen );
extern void updateModelQPBU     ( ImageParameters *p_Img, InputParameters *p_Inp, RCQuadratic *p_quad, int m_Qp );
extern void updateQPInterlaceBU ( RCQuadratic *p_quad, RCGeneric *p_gen );
extern void updateModelQPFrame  ( RCQuadratic *p_quad, int m_Bits );

extern void updateRCModel    (ImageParameters *p_Img, InputParameters *p_Inp, RCQuadratic *p_quad, RCGeneric *p_gen);
extern void updateMADModel   (ImageParameters *p_Img, InputParameters *p_Inp, RCQuadratic *p_quad, RCGeneric *p_gen);
extern void RCModelEstimator (ImageParameters *p_Img, InputParameters *p_Inp, RCQuadratic *p_quad, int n_windowSize, Boolean *m_rgRejected);
extern void MADModelEstimator(ImageParameters *p_Img, InputParameters *p_Inp, RCQuadratic *p_quad, int n_windowSize, Boolean *PictureRejected);
extern int  updateComplexity (ImageParameters *p_Img, RCQuadratic *p_quad, RCGeneric *p_gen, Boolean is_updated, int nbits );
extern void updatePparams    (RCQuadratic *p_quad, RCGeneric *p_gen, int complexity );
extern void updateBparams    (RCQuadratic *p_quad, RCGeneric *p_gen, int complexity );

// external generic functions
extern int  rc_handle_mb         ( Macroblock *currMB, int prev_mb);
extern void rc_init_top_field    ( ImageParameters *p_Img, InputParameters *p_Inp );
extern void rc_init_bottom_field ( ImageParameters *p_Img, InputParameters *p_Inp, int TopFieldBits );
extern void rc_init_frame_rdpic  ( ImageParameters *p_Img, InputParameters *p_Inp, float rateRatio );
extern void rc_allocate_memory   ( ImageParameters *p_Img, InputParameters *p_Inp );
extern void rc_free_memory       ( ImageParameters *p_Img, InputParameters *p_Inp );
extern void rc_update_mb_stats   ( Macroblock *currMB);
extern void rc_save_state        ( ImageParameters *p_Img, InputParameters *p_Inp );
extern void rc_restore_state     ( ImageParameters *p_Img, InputParameters *p_Inp );

#endif
