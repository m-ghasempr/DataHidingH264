
/*!
 *************************************************************************************
 * \file explicit_gop.h
 *
 * \brief
 *    Functions for explicit gop and hierarchy support
 *
 * \author
 *     Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Alexis Michael Tourapis          <alexismt@ieee.org>
 *************************************************************************************
 */

#ifndef _EXPLICIT_GOP_H_
#define _EXPLICIT_GOP_H_

// GOP Hierarchy
extern void init_gop_structure     (ImageParameters *p_Img, InputParameters *p_Inp);
extern void interpret_gop_structure(ImageParameters *p_Img, InputParameters *p_Inp);
extern void create_hierarchy       (ImageParameters *p_Img, InputParameters *p_Inp);
extern void clear_gop_structure    (ImageParameters *p_Img);

#endif
