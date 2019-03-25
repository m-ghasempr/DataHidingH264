/*!
 **************************************************************************
 *  \file enc_statistics.h
 *
 *  \brief
 *     statistics reports for the encoding process.
 *
 *  \author
 *      Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Alexis Tourapis                 <alexismt@ieee.org>
 *      - Karsten Sühring                 <suehring@hhi.de>
 *
 **************************************************************************
 */

#ifndef _ENC_STATISTICS_H_
#define _ENC_STATISTICS_H_
#include "global.h"

 typedef struct
 {
   int   quant0;                      //!< quant for the first frame
   int   quant1;                      //!< average quant for the remaining frames
   float bitr;                        //!< bit rate for current frame, used only for output til terminal
   float bitrate;                     //!< average bit rate for the sequence except first frame
   int64 bit_ctr;                     //!< counter for bit usage
   int64 bit_ctr_n;                   //!< bit usage for the current frame
   int   bit_slice;                   //!< number of bits in current slice
   int   stored_bit_slice;            //!< keep number of bits in current slice (to restore status in case of MB re-encoding)
   int   bit_ctr_emulationprevention; //!< stored bits needed to prevent start code emulation
   int   b8_mode_0_use     [NUM_SLICE_TYPES][2];
   int64 mode_use_transform[NUM_SLICE_TYPES][MAXMODE][2];
   int   intra_chroma_mode[4];

   // B pictures
   int     successive_Bframe;
   int     *mode_use_Bframe;
   int     *bit_use_mode_Bframe;

   int     frame_counter;
   int     frame_ctr           [NUM_SLICE_TYPES];
   int64   bit_counter[NUM_SLICE_TYPES];
   float   bitrate_st [NUM_SLICE_TYPES];
   int64   mode_use            [NUM_SLICE_TYPES][MAXMODE]; //!< Macroblock mode usage for Intra frames
   int64   bit_use_mode        [NUM_SLICE_TYPES][MAXMODE]; //!< statistics of bit usage
   int64   bit_use_stuffingBits[NUM_SLICE_TYPES];
   int64   bit_use_mb_type     [NUM_SLICE_TYPES];
   int64   bit_use_header      [NUM_SLICE_TYPES];
   int64   tmp_bit_use_cbp     [NUM_SLICE_TYPES];
   int64   bit_use_coeffC      [NUM_SLICE_TYPES];
   int64   bit_use_coeff    [3][NUM_SLICE_TYPES];  
   int64   bit_use_delta_quant [NUM_SLICE_TYPES];

   int   em_prev_bits_frm;
   int   em_prev_bits_fld;
   int  *em_prev_bits;
   int   bit_ctr_parametersets;
   int   bit_ctr_parametersets_n;
 } StatParameters;

 extern StatParameters *stats;

#endif

