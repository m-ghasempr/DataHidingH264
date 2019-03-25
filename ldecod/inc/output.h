
/*!
 **************************************************************************************
 * \file
 *    output.h
 * \brief
 *    Picture writing routine headers
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Karsten Suehring        <suehring@hhi.de>
 ***************************************************************************************
 */

#ifndef _OUTPUT_H_
#define _OUTPUT_H_


extern void write_stored_frame(ImageParameters *p_Img, FrameStore *fs, int p_out);
extern void direct_output     (ImageParameters *p_Img, StorablePicture *p, int p_out);
extern void init_out_buffer   (ImageParameters *p_Img);
extern void uninit_out_buffer (ImageParameters *p_Img);

#if (PAIR_FIELDS_IN_OUTPUT)
extern void flush_pending_output(ImageParameters *p_Img, int p_out);
#endif

#endif //_OUTPUT_H_
