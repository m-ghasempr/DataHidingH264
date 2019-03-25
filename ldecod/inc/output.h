
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

void write_stored_frame(FrameStore *fs, FILE *p_out);
void direct_output(StorablePicture *p, FILE *p_out);
void init_out_buffer();
void uninit_out_buffer();

#ifdef PAIR_FIELDS_IN_OUTPUT
void flush_pending_output(FILE *p_out);
#endif
