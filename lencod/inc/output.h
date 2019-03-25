
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

extern void flush_direct_output(ImageParameters *p_Img, InputParameters *p_Inp, FrameFormat *output, int p_out);
extern void write_out_picture  ( StorablePicture *p, FrameFormat *output, int p_out);
extern void write_stored_frame (ImageParameters *p_Img, InputParameters *p_Inp, FrameStore *fs, FrameFormat *output, int p_out);
extern void direct_output      (ImageParameters *p_Img, InputParameters *p_Inp, StorablePicture *p, FrameFormat *output, int p_out);
extern void direct_output_paff (ImageParameters *p_Img, InputParameters *p_Inp, StorablePicture *p, FrameFormat *output, int p_out);
extern void init_out_buffer    (ImageParameters *p_Img);
extern void uninit_out_buffer  (ImageParameters *p_Img, InputParameters *p_Inp);


#endif //_OUTPUT_H_
