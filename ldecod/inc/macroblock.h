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
 ************************************************************************
 * \file macroblock.h
 *
 * \author
 *  Inge Lille-Langøy               <inge.lille-langoy@telenor.com>
 *  Copyright (C) 1999  Telenor Satellite Services, Norway
 ************************************************************************
 */

#ifndef _MACROBLOCK_H_
#define _MACROBLOCK_H_

#define SINGLE_SCAN 0
#define DOUBLE_SCAN 1

/*
Array used to decide type of prediction:
0 = DC prediction
1 = vertical prediction
2 = horisontal prediction
3 = diagonal
4 = diagonal
all not appl. choises is set to 0
*/

const byte PRED_IPRED[7][7][6]=
{
  {
    {0,0,0,0,0,0},
    {0,4,5,0,0,0},
    {0,4,5,0,0,0},
    {0,4,5,0,0,0},
    {0,4,5,0,0,0},
    {4,0,5,0,0,0},
    {5,0,4,0,0,0},
  },
  {
    {0,2,1,0,0,0},
    {0,4,1,3,5,2},
    {0,1,4,3,2,5},
    {0,1,2,3,4,5},
    {3,0,4,1,5,2},
    {4,0,3,5,1,2},
    {5,4,0,3,1,2},
  },
  {
    {1,0,2,0,0,0},
    {1,0,4,3,2,5},
    {1,0,2,4,3,5},
    {1,0,2,3,4,5},
    {3,1,0,4,2,5},
    {4,0,1,5,3,2},
    {0,1,5,4,3,2},
  },
  {
    {2,0,1,0,0,0},
    {2,3,0,4,1,5},
    {2,0,3,1,4,5},
    {2,1,0,3,4,5},
    {2,3,1,0,5,4},
    {2,4,0,3,5,1},
    {2,0,1,4,5,3},
  },
  {
    {0,1,2,0,0,0},
    {3,0,4,2,1,5},
    {0,3,2,1,4,5},
    {3,0,2,1,4,5},
    {3,0,4,2,1,5},
    {4,3,0,5,1,2},
    {5,3,0,4,1,2},
  },
  {
    {0,1,2,0,0,0},
    {0,4,3,1,5,2},
    {0,4,1,3,2,5},
    {0,4,2,1,3,5},
    {4,0,3,5,1,2},
    {4,0,3,5,1,2},
    {4,5,0,3,1,2},
  },
  {
    {0,1,2,0,0,0},
    {0,4,3,5,1,2},
    {0,1,4,3,5,2},
    {0,1,3,2,4,5},
    {3,0,5,4,1,2},
    {4,0,5,3,1,2},
    {5,0,4,1,3,2},
  }
};

const byte IPRED_ORDER[36][2]=
{
  {0,0},{1,0},{0,1},{0,2},{1,1},{2,0},
  {3,0},{2,1},{1,2},{0,3},{0,4},{1,3},
  {2,2},{3,1},{4,0},{5,0},{4,1},{3,2},
  {2,3},{1,4},{0,5},{1,5},{2,4},{3,3},
  {4,2},{5,1},{5,2},{4,3},{3,4},{2,5},
  {3,5},{4,4},{5,3},{5,4},{4,5},{5,5}
};

//! single scan pattern
const byte SNGL_SCAN[16][2] =
{
  {0,0},{1,0},{0,1},{0,2},
  {1,1},{2,0},{3,0},{2,1},
  {1,2},{0,3},{1,3},{2,2},
  {3,1},{3,2},{2,3},{3,3}
};

//! double scan pattern
const byte DBL_SCAN[8][2][2] =
{
  {{0,0},{0,1}},
  {{1,0},{0,2}},
  {{2,1},{0,1}},
  {{2,1},{1,2}},
  {{2,0},{2,3}},
  {{3,1},{0,3}},
  {{3,2},{1,3}},
  {{3,3},{2,3}},
};

//! gives CBP value from codeword number, both for intra and inter
const byte NCBP[48][2]=
{
  {47, 0},{31,16},{15, 1},{ 0, 2},{23, 4},{27, 8},{29,32},{30, 3},{ 7, 5},{11,10},{13,12},{14,15},
  {39,47},{43, 7},{45,11},{46,13},{16,14},{ 3, 6},{ 5, 9},{10,31},{12,35},{19,37},{21,42},{26,44},
  {28,33},{35,34},{37,36},{42,40},{44,39},{ 1,43},{ 2,45},{ 4,46},{ 8,17},{17,18},{18,20},{20,24},
  {24,19},{ 6,21},{ 9,26},{22,28},{25,23},{32,27},{33,29},{34,30},{36,22},{40,25},{38,38},{41,41},
};

//! used to control block sizes : Not used/16x16/16x8/8x16/8x8/8x4/4x8/4x4
const int BLOCK_STEP[8][2]=
{
  {0,0},{4,4},{4,2},{2,4},{2,2},{2,1},{1,2},{1,1}
};

//! QP dependent scale factor for coefficients
const int JQ1[]=
{
  3881,  4351, 4890, 5481,  6154,  6914,  7761,  8718,
  9781, 10987,12339,13828, 15523, 17435, 19561, 21873,
  24552,27656,30847,34870, 38807, 43747, 49103, 54683,
  61694,68745,77615,89113,100253,109366,126635,141533,
};
const int JQ[32] =
{
   620,  553, 492,  439,   391, 348, 310,  276,
   246,  219, 195,  174,   155, 138, 123,  110,
    98,   87,  78,   69,    62,  55,  49,   44,
    39,   35,  31,   27,    24,  22,  19,   17,
};
// gives chroma QP from QP
const byte QP_SCALE_CR[32] =
{
  0 , 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,
  16,17,17,18,19,20,20,21,22,22,23,23,24,24,25,25
};

#endif

