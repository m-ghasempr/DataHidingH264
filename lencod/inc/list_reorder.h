/*!
 ***************************************************************************
 * \file
 *    list_reorder.h
 *
 * \date
 *    25 Feb 2009
 *
 * \brief
 *    Headerfile for slice-related functions
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Athanasios Leontaris            <aleon@dolby.com>
 *     - Karsten Sühring                 <suehring@hhi.de> 
 *     - Alexis Michael Tourapis         <alexismt@ieee.org> 
 **************************************************************************
 */

#ifndef _LIST_REORDER_H_
#define _LIST_REORDER_H_

#include "global.h"
#include "mbuffer.h"

extern void init_ref_pic_list_reordering( Slice *currSlice );
extern void reorder_lists               ( Slice *currSlice );
extern void wp_mcprec_reorder_lists     ( Slice *currSlice );


extern void poc_ref_pic_reorder_frame_default( Slice *currSlice, StorablePicture **list, unsigned num_ref_idx_lX_active, 
                                       int *reordering_of_pic_nums_idc, int *abs_diff_pic_num_minus1, int *long_term_pic_idx, int list_no );
extern void poc_ref_pic_reorder_field( Slice *currSlice, StorablePicture **list, unsigned num_ref_idx_lX_active, 
                               int *reordering_of_pic_nums_idc, int *abs_diff_pic_num_minus1, int *long_term_pic_idx, int list_no );

#endif
