/*
***********************************************************************
*  COPYRIGHT  AND  WARRANTY INFORMATION
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
 *    golomb_dec.h
 *
 * \brief
 *    Description
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *     -  Achim Dahlhoff      <dahlhoff@ient.rwth-aachen.de>
 *
 * \date
 *    Fri Mar 8 2002
 *
 *  copyright : (C) 2002      Institut und Lehrstuhl für Nachrichtentechnik
 *                            RWTH Aachen University
 *                            52072 Aachen
 *                            Germany
 *************************************************************************************
 */

#ifndef GOLOMB_H
#define GOLOMB_H

#include "global.h"

unsigned int decode_golomb_word(const unsigned char **buffer,unsigned int *bitoff,unsigned int grad0,unsigned int max_levels);
unsigned int decode_multilayer_golomb_word(const unsigned char **buffer,unsigned int *bitoff,const unsigned int *grad0,const unsigned int *max_levels);



int  readSyntaxElement_GOLOMB(SyntaxElement *sym, struct img_par *img, struct inp_par *inp, struct datapartition *dp);


#endif
