
/*!
************************************************************************
* \file quantChroma.h
*
* \brief
*    Quantization process for chroma header file
*
* \author
*    Limin Liu                       <lliu@dolby.com>
*    Alexis Michael Tourapis         <alexismt@ieee.org>                
*
************************************************************************
*/

#ifndef _QUANT_CR_H_
#define _QUANT_CR_H_

void init_quant_Chroma(InputParameters *params, ImageParameters *img, Slice *currSlice);

int (*quant_4x4cr)   (int **tblock, int block_y, int block_x, int qp, 
                      int*  ACLevel, int*  ACRun, 
                      int **fadjust4x4, int **levelscale, int **invlevelscale, int **leveloffset,
                      int *coeff_cost, const byte (*pos_scan)[2], const byte *c_cost, int is_cavlc);

int quant_dc2x2_normal (int **tblock, int qp, int* DCLevel, int* DCRun, 
                        int **fadjust2x2, int levelscale, int invlevelscale, int leveloffset,
                        const byte (*pos_scan)[2], int is_cavlc);

int quant_dc2x2_around (int **tblock, int qp, int* DCLevel, int* DCRun, 
                        int **fadjust2x2, int levelscale, int invlevelscale, int leveloffset,
                        const byte (*pos_scan)[2], int is_cavlc);

int quant_dc2x2_trellis(int **tblock, int qp, int* DCLevel, int* DCRun, 
                        int **fadjust, int levelscale, int invlevelscale, int leveloffset, 
                        const byte (*pos_scan)[2], int is_cavlc);

int quant_dc4x2_normal (int **tblock, int qp, int* DCLevel, int* DCRun, 
                        int **fadjust, int levelscale, int invlevelscale, int leveloffset,
                        const byte (*pos_scan)[2], int is_cavlc);

int quant_dc4x2_around (int **tblock, int qp, int* DCLevel, int* DCRun, 
                        int **fadjust, int levelscale, int invlevelscale, int leveloffset,
                        const byte (*pos_scan)[2], int is_cavlc);

int quant_dc4x2_trellis(int **tblock, int qp, int* DCLevel, int* DCRun, 
                        int **fadjust, int levelscale, int invlevelscale, int leveloffset,
                        const byte (*pos_scan)[2], int is_cavlc);

int (*quant_dc_cr)     (int **tblock, int qp, int* DCLevel, int* DCRun, 
                        int **fadjust, int levelscale, int invlevelscale, int leveloffset,
                        const byte (*pos_scan)[2], int is_cavlc);

int (*quant_ac4x4cr)   (int **tblock, int block_y, int block_x, int qp,                 
                        int*  ACLevel, int*  ACRun, 
                        int **fadjust4x4, int **levelscale, int **invlevelscale, int **leveloffset,
                        int *coeff_cost, const byte (*pos_scan)[2], const byte *c_cost, int type, int is_cavlc);

void rdoq_dc_cr_CAVLC  (int **tblock, int qp_per, int qp_rem, int levelscale, int leveloffset, 
                        const byte (*pos_scan)[2], int levelTrellis[16], int type);

void rdoq_dc_cr_CABAC  (int **tblock, int qp_per, int qp_rem, int levelscale, int leveloffset, 
                        const byte (*pos_scan)[2], int levelTrellis[16], int type);

void (*rdoq_dc_cr)     (int **tblock, int qp_per, int qp_rem, int levelscale, int leveloffset, 
                        const byte (*pos_scan)[2], int levelTrellis[16], int type);

#endif

