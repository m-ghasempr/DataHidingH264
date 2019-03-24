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
 * \file
 *    macroblock.h
 *
 * \author
 *    Inge Lille-Langøy               <inge.lille-langoy@telenor.com>     \n
 *    Telenor Satellite Services                                          \n
 *    P.O.Box 6914 St.Olavs plass                                         \n
 *    N-0130 Oslo, Norway
 *
 ************************************************************************/
#ifndef _MACROBLOCK_H_
#define _MACROBLOCK_H_


//! just to make new temp intra mode table
const int  MODTAB[3][2]=
{
  { 0, 4},
  {16,12},
  { 8,20}
};

//! gives codeword number from CBP value, both for intra and inter
const int NCBP[48][2]=
{
  { 3, 0},{29, 2},{30, 3},{17, 7},{31, 4},{18, 8},{37,17},{ 8,13},{32, 5},{38,18},{19, 9},{ 9,14},
  {20,10},{10,15},{11,16},{ 2,11},{16, 1},{33,32},{34,33},{21,36},{35,34},{22,37},{39,44},{ 4,40},
  {36,35},{40,45},{23,38},{ 5,41},{24,39},{ 6,42},{ 7,43},{ 1,19},{41, 6},{42,24},{43,25},{25,20},
  {44,26},{26,21},{46,46},{12,28},{45,27},{47,47},{27,22},{13,29},{28,23},{14,30},{15,31},{ 0,12},
};

/*!
Return prob.(0-8) for the input intra prediction mode(0-8),
depending on previous (right/above) pred mode(0-9).

NA values are set to 9 in the array.
The modes for the neighbour blocks are signalled:
0 = DC prediction
1 = vertical prediction
2 = horizontal prediction
3 = diagonal SE
4 = diagonal NE
5 = diagonal SSE
6 = diagonal NNE
7 = diagonal ENE
8 = diagonal ESE

prob order=PRED_IPRED[A(block above)][B(block left)][intra mode] */
const byte PRED_IPRED[10][10][9]=
{
  { // A=outside
    {0,9,9,9,9,9,9,9,9}, // B=outside
    {0,9,1,9,9,9,9,9,9}, // B=mode0   
    {0,9,1,9,9,9,9,9,9}, // B=mode1
    {1,9,0,9,9,9,9,9,9}, // B=mode2
    {1,9,0,9,9,9,9,9,9}, // B=mode3
    {1,9,0,9,9,9,9,9,9}, // B=mode4
    {1,9,0,9,9,9,9,9,9}, // B=mode5
    {1,9,0,9,9,9,9,9,9}, // B=mode6
    {0,9,1,9,9,9,9,9,9}, // B=mode7
    {0,9,1,9,9,9,9,9,9}, // B=mode8
  },
  { // A=mode0
    {0,1,9,9,9,9,9,9,9}, // outside
    {0,2,1,8,4,6,3,7,5},
    {1,0,2,6,5,4,3,8,7},
    {2,3,0,6,5,8,7,4,1},
    {1,2,0,3,6,5,8,7,4},
    {1,2,0,7,3,8,4,5,6},
    {0,1,3,5,7,2,4,8,6},
    {0,1,3,8,4,6,2,5,7},
    {2,3,0,7,4,8,6,1,5},
    {2,3,0,5,6,7,8,4,1},
  },
  { // A=mode1
    {1,0,9,9,9,9,9,9,9}, // outside
    {5,0,1,4,6,2,3,8,7},
    {5,0,2,4,6,3,1,8,7},
    {7,1,0,5,8,6,3,2,4},
    {8,0,1,3,6,2,4,7,5},
    {3,0,2,8,4,5,1,7,6},
    {7,0,2,4,6,1,3,8,5},
    {2,0,3,7,4,5,1,6,8},
    {4,1,0,8,7,6,3,2,5},
    {8,0,1,4,5,6,7,2,3},
  },
  { // A=mode2
    {1,0,9,9,9,9,9,9,9}, // outside
    {0,2,1,8,7,6,5,4,3},
    {2,0,1,8,6,4,3,5,7},
    {4,3,0,7,6,8,5,2,1},
    {1,3,0,4,7,6,8,5,2},
    {1,3,0,7,2,8,6,4,5},
    {1,2,0,8,5,3,6,7,4},
    {1,3,0,8,4,7,2,5,6},
    {4,3,0,8,6,5,7,1,2},
    {4,3,0,6,5,8,7,2,1},
  },
  { // A=mode3
    {1,0,9,9,9,9,9,9,9}, // outside
    {0,2,1,3,7,4,6,8,5},
    {1,0,2,4,6,3,5,8,7},
    {3,2,0,4,8,5,7,6,1},
    {8,4,1,0,5,2,6,7,3},
    {2,4,1,7,0,5,3,8,6},
    {7,2,3,1,6,0,5,8,4},
    {2,0,3,7,4,5,1,8,6},
    {2,3,0,8,5,4,7,1,6},
    {5,4,0,2,8,3,7,6,1},
  },
  { // A=mode4
    {1,0,9,9,9,9,9,9,9}, // outside
    {1,3,0,8,4,6,2,5,7},
    {3,0,2,6,4,5,1,7,8},
    {5,6,0,8,4,7,3,2,1},
    {3,2,1,6,0,7,4,8,5},
    {3,4,1,8,0,5,2,6,7},
    {3,0,1,5,6,2,4,7,8},
    {2,3,4,7,1,6,0,5,8},
    {4,5,0,8,2,7,3,1,6},
    {8,5,0,6,3,7,4,2,1},
  },
  { // A=mode5
    {1,0,9,9,9,9,9,9,9}, // outside
    {6,1,2,3,7,0,4,8,5},
    {5,0,4,3,6,1,2,8,7},
    {8,1,0,3,7,2,4,6,5},
    {8,2,3,1,6,0,4,7,5},
    {6,0,2,7,3,4,1,8,5},
    {5,1,4,2,7,0,3,8,6},
    {4,0,3,6,5,2,1,7,8},
    {5,2,0,7,8,4,3,1,6},
    {7,2,0,3,6,1,4,8,5},
  },
  { // A=mode6
    {1,0,9,9,9,9,9,9,9}, // outside
    {3,0,2,6,5,4,1,7,8},
    {6,0,4,5,3,2,1,7,8},
    {8,1,0,7,4,6,2,3,5},
    {7,1,0,6,5,3,2,8,4},
    {4,2,3,8,1,5,0,6,7},
    {5,0,3,4,6,2,1,8,7},
    {3,1,4,7,2,6,0,5,8},
    {6,4,0,8,2,5,3,1,7},
    {8,1,0,6,4,7,2,5,3},
  },
  { // A=mode7
    {1,0,9,9,9,9,9,9,9}, // outside
    {1,5,0,8,2,7,4,3,6},
    {3,1,2,8,4,6,0,5,7},
    {3,4,0,7,6,8,5,1,2},
    {2,5,0,4,1,8,7,6,3},
    {3,5,1,8,0,7,4,2,6},
    {2,1,0,8,7,4,5,6,3},
    {5,3,0,8,2,6,1,4,7},
    {3,6,0,8,2,7,5,1,4},
    {5,4,0,6,3,8,7,2,1},
  },
  { // A=mode8
    {1,0,9,9,9,9,9,9,9}, // outside
    {1,3,0,4,5,7,6,8,2},
    {2,0,1,7,8,5,3,6,4},
    {4,3,0,5,8,7,6,2,1},
    {5,4,2,1,6,3,8,7,0},
    {1,5,0,8,2,4,7,6,3},
    {2,1,0,5,7,4,6,8,3},
    {3,1,0,8,5,6,2,7,4},
    {4,6,0,8,3,7,5,1,2},
    {3,6,0,2,5,8,7,4,1},
  }
};



/*
  return codeword number from two combined intra prediction blocks 
*/

const int IPRED_ORDER[9][9]= 
{
  { 0, 1, 4, 6,10,14,19,28,36},
  { 2, 3, 8,13,17,23,30,38,45},
  { 5, 9,16,21,25,32,41,47,53},
  { 7,12,22,27,34,43,51,58,64},
  {11,18,26,35,40,49,55,62,66},
  {15,24,33,44,50,57,60,68,73},
  {20,31,42,52,56,61,70,71,75},
  {29,39,48,59,63,69,72,77,78},
  {37,46,54,65,67,74,76,79,80}
};


extern int QP2QUANT[40];

#endif

