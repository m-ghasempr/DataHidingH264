/*!
 **************************************************************************
 *  \file errdo.h
 *  \brief  Header file for error resilient RDO (name of file should change)
 *
 *  \author 
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Peshala Pahalawatta                     <ppaha@dolby.com>
 *    - Alexis Michael Tourapis                 <atour@ieee.org>
 *
 **************************************************************************
 */

#ifndef _ERRDO_H_
#define _ERRDO_H_

//! Info for the "decoders-in-the-encoder" used for rdoptimization with packet losses
struct decoders
{
  imgpel ***dec_mbY_best;            //!< Best reconstructed macroblock pixel values
  imgpel ****dec_mbY_best8x8;        //!< Best reconstructed 8x8 mode pixel values
  imgpel ****dec_mb_pred_best8x8;    //!< Predicted pixel values for best 8x8 modes
  imgpel ***dec_mb_pred;             //!< Predicted pixel values for macroblock
  int    ***res_img;                 //!< Residual values for macroblock
};

typedef struct decoders Decoders;

//============= rate-distortion opt with packet losses ===========

extern void init_error_conceal      (VideoParameters *p_Vid, int concealment_type);
extern void compute_residue_block   (Macroblock *currMB, imgpel **imgY, int **res_img, imgpel **mb_pred, int b8block, int block_size);
extern void decode_one_b8block      (Macroblock* currMB, StorablePicture *enc_pic, int decoder, int block8x8, short mv_mode, int pred_dir);
extern void errdo_store_best_block  (InputParameters *p_Inp, imgpel*** mbY, imgpel*** dec_img, int block_i, int block_j, int img_i, int img_j, int block_size);
extern void decode_one_mb           (Macroblock* currMB, StorablePicture *enc, int decoder);
extern void UpdateDecoders          (VideoParameters *p_Vid, InputParameters *p_Inp, StorablePicture *enc_pic);
extern void copy_conceal_picture    (VideoParameters *p_Vid, StorablePicture *enc_pic, int decoder);
extern void errdo_get_best_block    (Macroblock *currMB, imgpel*** dec_img, imgpel*** mbY, int j0, int block_size);
extern int  allocate_errdo_mem      (VideoParameters *p_Vid, InputParameters *p_Inp);
extern void free_errdo_mem          (VideoParameters *p_Vid);

#endif

