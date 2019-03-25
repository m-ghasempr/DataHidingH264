/*!
 **************************************************************************************
 * \file
 *    mmco.h
 * \brief
 *    MMCO example operations.
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Alexis Michael Tourapis                     <alexismt@ieee.org>
 *     - Athanasios Leontaris                        <aleon@dolby.com>
 ***************************************************************************************
 */

#ifndef _MMCO_H_
#define _MMCO_H_

extern void mmco_long_term(ImageParameters *p_Img, int current_pic_num);
extern void poc_based_ref_management_frame_pic(ImageParameters *p_Img, int current_pic_num);
extern void poc_based_ref_management_field_pic(ImageParameters *p_Img, int current_pic_num);

#endif 
