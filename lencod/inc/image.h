
/*!
 ************************************************************************
 * \file image.h
 *
 * \brief
 *    headers for image processing
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Karsten Sühring                 <suehring@hhi.de> 
 *     - Inge Lille-Langoy               <inge.lille-langoy@telenor.com>
 *     - Alexis Michael Tourapis         <alexismt@ieee.org> 
 *  
 ************************************************************************
 */
#ifndef _IMAGE_H_
#define _IMAGE_H_

#include "mbuffer.h"

extern int     encode_one_frame      ( VideoParameters *p_Vid, InputParameters *p_Inp);
extern Boolean dummy_slice_too_big   ( int bits_slice);
extern void    copy_rdopt_data       ( Macroblock *currMB);       // For MB level field/frame coding tools
extern void    UnifiedOneForthPix    ( VideoParameters *p_Vid, StorablePicture *s);
// For 4:4:4 independent mode
extern void    UnifiedOneForthPix_JV ( VideoParameters *p_Vid, int nplane, StorablePicture *s);
extern void    frame_picture         ( VideoParameters *p_Vid, Picture *frame, ImageData *imgData, int rd_pass);
extern byte    get_idr_flag          ( VideoParameters *p_Vid );
extern void    write_non_vcl_nalu    ( VideoParameters *p_Vid, InputParameters *p_Inp );

#endif

