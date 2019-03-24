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

/*!
 ***************************************************************************
 * \file
 *    biariencode.h
 *
 * \brief
 *    Headerfile for binary arithmetic encoding routines
 *
 * \author
 *    Detlev Marpe,
 *    Gabi Blaettermann                                                     \n
 *    Copyright (C) 2000 HEINRICH HERTZ INSTITUTE All Rights Reserved.
 *
 * \date
 *    21. Oct 2000
 **************************************************************************
 */


#ifndef _BIARIENCOD_H_
#define _BIARIENCOD_H_



/************************************************************************
 * D e f i n i t i o n s
 ***********************************************************************
 */

// some definitions to increase the readability of the source code

#define Elow            (eep->Elow)
#define Erange          (eep->Erange)
#define Ebits_to_follow (eep->Ebits_to_follow)
#define Ebuffer         (eep->Ebuffer)
#define Ebits_to_go     (eep->Ebits_to_go)
#define Ecodestrm       (eep->Ecodestrm)
#define Ecodestrm_len   (eep->Ecodestrm_len)
#define Ecodestrm_laststartcode   (eep->Ecodestrm_laststartcode)

/* Only necessary for new AC */
#define B_BITS    16
#define F_BITS    14
#define CACM99_HALF       (1 << (B_BITS-1))
#define CACM99_QUARTER    (1 << (B_BITS-2))

#define MAX_FREQ                    16 


const unsigned short rLPS_table_64x4[64][4]=
{
{9216,  11264,  13312,  15360},
{8832,  10816,  12800,  14720},
{8512,  10368,  12288,  14144},
{8128,  9920,   11712,  13504},
{7680,  9344,   11072,  12736},
{7168,  8768,   10368,  11968},
{6912,  8448,   9984,   11520},
{6336,  7808,   9216,   10624},
{5888,  7232,   8512,   9856},
{5440,  6656,   7872,   9088},
{5120,  6208,   7360,   8512},
{4608,  5632,   6656,   7680},
{4224,  5184,   6144,   7104},
{3968,  4800,   5696,   6592},
{3712,  4480,   5312,   6144},
{3456,  4224,   4992,   5760},
{3072,  3776,   4416,   5120},
{2816,  3456,   4096,   4736},
{2624,  3200,   3776,   4416},
{2432,  3008,   3520,   4096},
{2304,  2816,   3328,   3840},
{2048,  2496,   2944,   3392},
{1856,  2240,   2688,   3072},
{1664,  2048,   2432,   2816},
{1536,  1856,   2240,   2560},
{1408,  1728,   2048,   2368},
{1344,  1600,   1920,   2176},
{1216,  1472,   1792,   2048},
{1152,  1408,   1664,   1920},
{1088,  1344,   1536,   1792},
{1024,  1280,   1472,   1728},
{960,   1216,   1408,   1600},
{896,   1152,   1344,   1536},
{896,   1088,   1280,   1472},
{832,   1024,   1216,   1408},
{832,   960,    1152,   1344},
{768,   960,    1088,   1280},
{768,   896,    1088,   1216},
{704,   896,    1024,   1152},
{704,   832,    960,    1152},
{640,   832,    960,    1088},
{640,   768,    896,    1088},
{640,   768,    896,    1024},
{576,   704,    832,    960},
{576,   704,    832,    960},
{576,   704,    832,    960},
{512,   640,    768,    896},
{512,   640,    768,    896},
{512,   640,    768,    832},
{512,   640,    704,    832},
{512,   576,    704,    832},
{448,   576,    704,    768},
{448,   576,    640,    768},
{448,   512,    640,    704},
{448,   512,    640,    704},
{448,   512,    576,    704},
{384,   512,    576,    704},
{384,   512,    576,    640},
{384,   448,    576,    640},
{384,   448,    576,    640},
{384,   448,    512,    640},
{384,   448,    512,    640},
{384,   448,    512,    576},
{384,   448,    512,    576}

};

unsigned short AC_next_state_MPS_64_INTRA[64] =    
{
  //0 - MPS
        1,2,3,4,5,6,7,8,9,10,
        11,12,13,14,15,16,17,18,19,20,
        21,22,23,23,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,
        0,0,0,63
};  

unsigned short AC_next_state_MPS_64_INTER[64] =    
{
  //0 - MPS
        1,2,3,4,5,6,7,8,9,10,
        11,12,13,14,15,16,17,18,19,20,
        21,22,23,24,25,26,27,28,29,30,
        31,32,33,34,35,36,37,38,39,40,
        41,42,43,44,45,46,47,38,39,50,
        51,52,53,54,55,56,57,58,59,60,
        61,62,63,63
};  

const unsigned short AC_next_state_LPS_64[64] =    
{
  //0 - LPS
        0, 0, 1, 2, 2, 3, 4, 5, 6, 7, 8, 8, 10, 10, 10,
        11, 12, 13, 14, 14, 14, 14, 15, 16, 17, 18, 19, 
        19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 
        25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 30, 
        31, 32, 33, 33, 34, 34, 35, 35, 36, 37, 37, 38,
        38
};


#endif  // BIARIENCOD_H

