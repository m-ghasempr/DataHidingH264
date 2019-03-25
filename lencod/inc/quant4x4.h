
/*!
 ************************************************************************
 * \file quant4x4.h
 *
 * \brief
 *    Quantization process header file
 *
 * \author
 *    Alexis Michael Tourapis         <alexismt@ieee.org>                
 *
 ************************************************************************
 */

#ifndef _QUANT4x4_H_
#define _QUANT4x4_H_

void init_quant_4x4(ImageParameters *img);

int quant_4x4_normal(int (*tblock)[16], int block_y, int block_x, int qp, 
                     int*  ACLevel, int*  ACRun, 
                     int **fadjust4x4, int **levelscale, int **invlevelscale, int **leveloffset,
                     int *coeff_cost, const byte (*pos_scan)[2], const byte *c_cost);

int quant_4x4_around(int (*tblock)[16], int block_y, int block_x, int qp, 
                     int*  ACLevel, int*  ACRun, 
                     int **fadjust4x4, int **levelscale, int **invlevelscale, int **leveloffset,
                     int *coeff_cost, const byte (*pos_scan)[2], const byte *c_cost);

int quant_4x4_trellis (int (*tblock)[16], int block_y, int block_x, int qp,
                       int*  ACLevel, int*  ACRun, 
                       int **fadjust4x4, int **levelscale, int **invlevelscale, int **leveloffset,
                       int *coeff_cost, const byte (*pos_scan)[2], const byte *c_cost);

int (*quant_4x4)    (int (*tblock)[16], int block_y, int block_x, int qp, 
                     int*  ACLevel, int*  ACRun, 
                     int **fadjust4x4, int **levelscale, int **invlevelscale, int **leveloffset,
                     int *coeff_cost, const byte (*pos_scan)[2], const byte *c_cost);

int (*quant_4x4cr)    (int (*tblock)[16], int block_y, int block_x, int qp, 
                     int*  ACLevel, int*  ACRun, 
                     int **fadjust4x4, int **levelscale, int **invlevelscale, int **leveloffset,
                     int *coeff_cost, const byte (*pos_scan)[2], const byte *c_cost);


int quant_ac4x4_normal(int (*tblock)[16], int block_y, int block_x, int qp,                 
                       int*  ACLevel, int*  ACRun, 
                       int **fadjust4x4, int **levelscale, int **invlevelscale, int **leveloffset,
                       int *coeff_cost, const byte (*pos_scan)[2], const byte *c_cost);

int quant_ac4x4_around(int (*tblock)[16], int block_y, int block_x, int qp,                 
                       int*  ACLevel, int*  ACRun, 
                       int **fadjust4x4, int **levelscale, int **invlevelscale, int **leveloffset,
                       int *coeff_cost, const byte (*pos_scan)[2], const byte *c_cost);

int quant_ac4x4_trellis(int (*tblock)[16], int block_y, int block_x, int qp,                
                        int*  ACLevel, int*  ACRun, 
                        int **fadjust4x4, int **levelscale, int **invlevelscale, int **leveloffset,
                        int *coeff_cost, const byte (*pos_scan)[2], const byte *c_cost);


int (*quant_ac4x4)    (int (*tblock)[16], int block_y, int block_x, int qp,                 
                       int*  ACLevel, int*  ACRun, 
                       int **fadjust4x4, int **levelscale, int **invlevelscale, int **leveloffset,
                       int *coeff_cost, const byte (*pos_scan)[2], const byte *c_cost);

int (*quant_ac4x4cr)  (int (*tblock)[16], int block_y, int block_x, int qp,                 
                       int*  ACLevel, int*  ACRun, 
                       int **fadjust4x4, int **levelscale, int **invlevelscale, int **leveloffset,
                       int *coeff_cost, const byte (*pos_scan)[2], const byte *c_cost);


#endif

