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
typedef struct
{
  imgpel ***dec_mbY_best;            //!< Best reconstructed macroblock pixel values
  imgpel ****dec_mbY_best8x8;        //!< Best reconstructed 8x8 mode pixel values
  imgpel ****dec_mb_pred_best8x8;    //!< Predicted pixel values for best 8x8 modes
  imgpel ***dec_mb_pred;             //!< Predicted pixel values for macroblock
  int    ***res_img;                 //!< Residual values for macroblock
} Decoders;

extern Decoders *decs;

//============= rate-distortion opt with packet losses ===========

void init_error_conceal(int concealment_type);

void compute_residue_block (ImageParameters *image, imgpel **imgY, int **res_img, imgpel **mb_pred, int b8block, int block_size);

void decode_one_b8block (ImageParameters *image, StorablePicture *enc_pic, Macroblock* currMB, int decoder, int block8x8, short mv_mode, int pred_dir);
void errdo_store_best_block(imgpel*** mbY, imgpel*** dec_img, int block_i, int block_j, int img_i, int img_j, int block_size);

void decode_one_mb  (ImageParameters *image, StorablePicture *enc, int decoder, Macroblock* currMB);
void UpdateDecoders (InputParameters *params, ImageParameters *image, StorablePicture *enc_pic);

void (*error_conceal_picture)(ImageParameters *image, StorablePicture *enc_pic, int decoder);
void copy_conceal_picture(ImageParameters *image, StorablePicture *enc_pic, int decoder);

void errdo_get_best_block(ImageParameters* image, imgpel*** dec_img, imgpel*** mbY, int j0, int block_size);

int  allocate_errdo_mem(InputParameters *pparams);
void free_errdo_mem(void);

#endif

