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
 * \file image.h
 *
 * \author
 *  Inge Lille-Langøy               <inge.lille-langoy@telenor.com>
 *  Copyright (C) 1999  Telenor Satellite Services, Norway
 ************************************************************************
 */

#ifndef _IMAGE_H_
#define _IMAGE_H_

//! TAPs used in the oneforthpix()routine
const int ONE_FOURTH_TAP[3][2] =
{
  {20,20},
  {-5,-4},
  { 1, 0},
};

//! for new loopfilter
const byte FILTER_STR[32][4] =
{
  {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},
  {0,0,0,0},{0,0,0,1},{0,0,0,1},{0,0,1,1},{0,0,1,1},{0,0,1,1},{0,1,1,1},{0,1,1,1},
  {0,1,1,1},{0,1,1,2},{0,1,1,2},{0,1,1,2},{0,1,1,2},{0,1,2,3},{0,1,2,3},{0,1,2,3},
  {0,1,2,4},{0,2,3,4},{0,2,3,5},{0,2,3,5},{0,2,3,5},{0,2,4,7},{0,3,5,8},{0,3,5,9},
};

const byte LIM[32] =
{
   7,   8,  9, 10, 11, 12, 14, 16,
  18,  20, 22, 25, 28, 31, 35, 39,
  44,  49, 55, 62, 69, 78, 88, 98,
  110,124,139,156,175,197,221,248
};

//! convert from H.263 QP to H.26L quant given by: quant=pow(2,QP/6)
const int QP2QUANT[32]=
{
   1, 1, 1, 1, 2, 2, 2, 2,
   3, 3, 3, 4, 4, 4, 5, 6,
   6, 7, 8, 9,10,11,13,14,
  16,17,20,23,25,29,32,36
};


int two[6]  =  {64,-320,3328,1280,-320,64};
int three[6] = {128,-640,2560,2560,-640,128};

int five[6][6]  = {{  4, -20,  80,  80, -20,  4},
                   {-20, 100,-400,-400, 100,-20},
                   { 80,-400,1600,1600,-400, 80},
                   { 80,-400,1600,1600,-400, 80},
                   {-20, 100,-400,-400, 100,-20},
                   {  4, -20,  80,  80, -20,  4}
                  };
int six[6][6]  =  {{ 1,  -5,  52,  20,  -5,  1},
                   {-5,  25,-260,-100,  25, -5},
                   {52,-260,2704,1040,-260, 52},
                   {20,-100,1040, 400,-100, 20},
                   {-5,  25,-260,-100,  25, -5},
                   { 1,  -5,  52,  20,  -5,  1}
                  };
int seven[6][6] = {{  2, -10,  40,  40, -10,  2},
                   {-10,  50,-200,-200,  50,-10},
                   {104,-520,2080,2080,-520,104},
                   { 40,-200, 800, 800,-200, 40},
                   {-10,  50,-200,-200,  50,-10},
                   {  2, -10,  40,  40, -10,  2}
                  };

#endif

