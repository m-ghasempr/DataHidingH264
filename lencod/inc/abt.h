/*
***********************************************************************
* COPYRIGHT AND WARRANTY INFORMATION
*
* Copyright 2001, International Telecommunications Union, Geneva
*
* DISCLAIMER OF WARRANTY
*
* These software programs are available to the user without any
* license fee or royalty on an "as is" basis. The ITU disclaims
* any and all warranties, whether express, implied, or
* statutory, including any implied warranties of merchantability
* or of fitness for a particular purpose.  In no event shall the
* contributor or the ITU be liable for any incidental, punitive, or
* consequential damages of any kind whatsoever arising from the
* use of these programs.
*
* This disclaimer of warranty extends to the user of these programs
* and user's customers, employees, agents, transferees, successors,
* and assigns.
*
* The ITU does not represent or warrant that the programs furnished
* hereunder are free of infringement of any third-party patents.
* Commercial implementations of ITU-T Recommendations, including
* shareware, may be subject to royalty fees to patent holders.
* Information regarding the ITU-T patent policy is available from
* the ITU Web site at http://www.itu.int.
*
* THIS IS NOT A GRANT OF PATENT RIGHTS - SEE THE ITU-T PATENT POLICY.
************************************************************************
*/

/*
 *************************************************************************************
 * \file
 *    abt.h
 *
 * \brief
 *    Description
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *     -  Mathias Wien        <wien@ient.rwth-aachen.de>
 *
 * \date
 *    Fri Mar 8 2002
 *
 *  copyright : (C) 2002 by   Mathias Wien
 *                            Institut und Lehrstuhl für Nachrichtentechnik
 *                            RWTH Aachen University
 *                            52072 Aachen
 *                            Germany
 *************************************************************************************
 */

#ifndef _ABT_H_
#define _ABT_H_

#include "global.h"
//#include "defines.h"

#define NUM_2D_TABLES 4
#define CODE2D_ESCAPE_SYMBOL 59

// ========================================================
// external variables
// ========================================================
extern const int ABT_matrix[2][8][8];        //!< ABT transform matrices
extern const int ABT_TRSIZE[4][2];           //!< ABT transform block sizes by abt_mode
extern const int ABT_NUMTR[4][2];            //!< ABT number of hor./ver. transforms by abt_mode
extern const int ABT_TRIDX[4][2];            //!< ABT transform index  by abt_mode
extern const int ABT_Q[4][QUANT_PERIOD][3];  //!< ABT Quantization table
extern const int ABT_R[4][QUANT_PERIOD][3];  //!< ABT De-Quantization table
extern const int ABT_QMAP[4][2][2];          //!< ABT mapping the proper Q value for coefficient position
extern const int ABT_QF[4][2];               //!< ABT rounding factor
extern const int ABT_N[4][2];                //!< ABT quantization normalization
extern const int ABT_SHIFT0[4][2][2];        //!< ABT bit shift needed after transform step to stay inside 16 bits.
extern const int ABT_SHIFT1[4];              //!< ABT bit shift needed after inverse transform step to stay inside 16 bits.
extern const int ABT_SCAN[2][4][64][2];      //!< ABT scan positions. Positions are stored as (pix,lin).
extern const int ABT_COEFF_COST[4][64];      //!< TML COEFF_COST 'stretched' for ABT
extern const char ABT_2D_VLC[NUM_2D_TABLES][16][8];  //   Inter, Intra0-13, Intra14-21, Intra22-31
extern const unsigned short int cbp_blk_masks[4][4];


// ========================================================
// typedefs
// ========================================================
typedef enum
{
  PIX,
  LIN
} Direction;


// ========================================================
// functions
// ========================================================
void transform_ABT_B8      (int abt_mode, int blk_off_x, int blk_off_y, int curr_blk[B8_SIZE][B8_SIZE]);
void inv_transform_ABT_B8  (int abt_mode, int blk_off_x, int blk_off_y, int curr_blk[B8_SIZE][B8_SIZE]);
void quant_abt_B8          (int qp, int abt_mode, int blk_off_x, int blk_off_y, int curr_blk[B8_SIZE][B8_SIZE]);
int  scanquant_ABT_B8      (int qp, int abt_mode, int b8, int blk_off_x, int blk_off_y, int curr_blk[B8_SIZE][B8_SIZE], int scrFlag, int *cbp, int *cbp_blk);
int  trans_scanquant_ABT_sp(int abt_mode, int block8x8, int blk_off_x,int blk_off_y, int curr_blk[B8_SIZE][B8_SIZE], int scrFlag, int *cbp, int *cbp_blk);
int  find_sad_abt          (int iMode, int iSizeX, int iSizeY, int iOffX, int iOffY, int m7[MB_BLOCK_SIZE][MB_BLOCK_SIZE]);
int  sad_hadamard          (int iSizeX, int iSizeY, int iOffX, int iOffY, int m7[MB_BLOCK_SIZE][MB_BLOCK_SIZE]);
int  writeLumaCoeffABT_B8  (int b8,int intra,int blk_off_x,int blk_off_y);
void setDirectModeABT      (int block8x8);
int  getDirectModeABT      (int block8x8);
void copyblock_SP_ABT      (int abt_mode,int b8,int blk_off_x, int blk_off_y);
void get_quant_consts      (int abt_mode,int qp,int intra,int Q[2][2],int Qrshift[2][2],int qp_const[2][2]);
void get_dequant_consts    (int abt_mode,int qp,int R[2][2]);

int Mode_Decision_for_ABT_IntraBlocks(int b8,int b4,double lambda,int *min_cost,int bs_x,int bs_y);
double RDCost_for_ABTIntraBlocks(int *nonzero,int b8,int b4,int ipmode,double lambda,double  min_rdcost,int bs_x,int bs_y, int mostProbableMode);

#endif // _ABT_H_
