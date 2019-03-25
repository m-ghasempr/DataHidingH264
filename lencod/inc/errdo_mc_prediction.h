/*!
 *************************************************************************************
 * \file errdo_mc_prediction.h
 *
 * \brief
 *    definitions for motion compensated prediction
 *
 * \author
 *      Main contributors (see contributors.h for copyright, 
 *                         address and affiliation details)
 *      - Alexis Michael Tourapis  <alexismt@ieee.org>
 *      - Modified for use in hypothetical decoders at encoder
 *           by Peshala V. Pahalawatta <pesh@ieee.org>
 *
 *************************************************************************************
 */

#ifndef _ERRDO_MC_PREDICTION_H_
#define _ERRDO_MC_PREDICTION_H_

#include "global.h"
#include "mbuffer.h"


extern void get_block_luma  (Macroblock *currMB, int decoder, ColorPlane pl, StorablePicture *dec_picture, StorablePicture *list, int x_pos, int y_pos, int ver_block_size, int hor_block_size, imgpel block[MB_BLOCK_SIZE][MB_BLOCK_SIZE]);
extern void get_block_chroma(Macroblock *currMb, int decoder, int uv, StorablePicture *dec_picture, StorablePicture *list, int x_pos, int y_pos, int hor_block_size, int ver_block_size, imgpel block[MB_BLOCK_SIZE][MB_BLOCK_SIZE]);

extern void perform_mc            (Macroblock* currMB, int decoder, ColorPlane pl, StorablePicture *dec_picture, int pred_dir, int l0_mode, int l1_mode, char*** ref_idx, int i, int j, int block_size_x, int block_size_y, short bipred_me);
extern void perform_mc_concealment(Macroblock* currMB, int decoder, ColorPlane pl, StorablePicture *dec_picture, int pred_dir, int l0_mode, int l1_mode, char*** ref_idx_buf, int i, int j, int block_size_x, int block_size_y);

#endif

