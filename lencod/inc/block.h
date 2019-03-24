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
 * \file block.h
 *
 * \author
 *  Inge Lille-Langøy               <inge.lille-langoy@telenor.com>    \n
 *  Telenor Satellite Services                                         \n
 *  P.O.Box 6914 St.Olavs plass                                        \n
 *  N-0130 Oslo, Norway
 *
 ************************************************************************
 */

#ifndef _BLOCK_H_
#define _BLOCK_H_

#include "global.h"


//! make chroma QP from quant
const byte QP_SCALE_CR[52]=
{
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,
   12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,
   28,29,29,30,31,32,32,33,34,34,35,35,36,36,37,37,
   37,38,38,38,39,39,39,39
};


//! single scan pattern
const byte SNGL_SCAN[16][2] =
{
  {0,0},{1,0},{0,1},{0,2},
  {1,1},{2,0},{3,0},{2,1},
  {1,2},{0,3},{1,3},{2,2},
  {3,1},{3,2},{2,3},{3,3}
};

//! array used to find expencive coefficients
const byte COEFF_COST[16] =
{
  3,2,2,1,1,1,0,0,0,0,0,0,0,0,0,0
};



//! bit cost for coefficients
const byte COEFF_BIT_COST[3][16][16]=
{
  { // 2x2 scan (corrested per Gisle's Email 11/23/2000 by StW
    { 3, 5, 7, 9, 9,11,11,11,11,13,13,13,13,13,13,13},
    { 5, 7, 9, 9,11,11,11,11,13,13,13,13,13,13,13,13},
    { 7, 9, 9,11,11,11,11,13,13,13,13,13,13,13,13,15},
    { 7, 9, 9,11,11,11,11,13,13,13,13,13,13,13,13,15},
    { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
    { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
    { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
    { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
    { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
    { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
    { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
    { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
    { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
    { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
    { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
    { 7, 7, 9, 9, 9, 9,11,11,11,11,11,11,11,11,13,13},
  },
  {  // double scan
    { 3, 5, 7, 7, 7, 9, 9, 9, 9,11,11,13,13,13,13,15},
    { 5, 9, 9,11,11,13,13,13,13,15,15,15,15,15,15,15},
    { 7,11,11,13,13,13,13,15,15,15,15,15,15,15,15,17},
    { 9,11,11,13,13,13,13,15,15,15,15,15,15,15,15,17},
    { 9,11,11,13,13,13,13,15,15,15,15,15,15,15,15,17},
    {11,11,13,13,13,13,15,15,15,15,15,15,15,15,17,17},
    {11,11,13,13,13,13,15,15,15,15,15,15,15,15,17,17},
    {11,11,13,13,13,13,15,15,15,15,15,15,15,15,17,17},
    {11,11,13,13,13,13,15,15,15,15,15,15,15,15,17,17},
    {11,11,13,13,13,13,15,15,15,15,15,15,15,15,17,17},
    {11,11,13,13,13,13,15,15,15,15,15,15,15,15,17,17},
    {11,11,13,13,13,13,15,15,15,15,15,15,15,15,17,17},
    {11,11,13,13,13,13,15,15,15,15,15,15,15,15,17,17},
    {11,11,13,13,13,13,15,15,15,15,15,15,15,15,17,17},
  },
  {    // single scan
    { 3, 7, 9, 9,11,13,13,15,15,15,15,17,17,17,17,17},
    { 5, 9,11,13,13,15,15,15,15,17,17,17,17,17,17,17},
    { 5, 9,11,13,13,15,15,15,15,17,17,17,17,17,17,17},
    { 7,11,13,13,15,15,15,15,17,17,17,17,17,17,17,17},
    { 7,11,13,13,15,15,15,15,17,17,17,17,17,17,17,17},
    { 7,11,13,13,15,15,15,15,17,17,17,17,17,17,17,17},
    { 9,11,13,13,15,15,15,15,17,17,17,17,17,17,17,17},
    { 9,11,13,13,15,15,15,15,17,17,17,17,17,17,17,17},
    { 9,11,13,13,15,15,15,15,17,17,17,17,17,17,17,17},
    { 9,11,13,13,15,15,15,15,17,17,17,17,17,17,17,17},
    {11,13,13,15,15,15,15,17,17,17,17,17,17,17,17,19},
    {11,13,13,15,15,15,15,17,17,17,17,17,17,17,17,19},
    {11,13,13,15,15,15,15,17,17,17,17,17,17,17,17,19},
    {11,13,13,15,15,15,15,17,17,17,17,17,17,17,17,19},
    {11,13,13,15,15,15,15,17,17,17,17,17,17,17,17,19},
    {11,13,13,15,15,15,15,17,17,17,17,17,17,17,17,19},
  },
};

#endif

