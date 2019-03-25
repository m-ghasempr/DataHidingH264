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
 * \file  nalu.c
 *
 * \brief
 *    Common NALU support functions
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Stephan Wenger   <stewe@cs.tu-berlin.de>
 ************************************************************************
 */

#include <stdio.h>
#include <assert.h>
#include <malloc.h>
#include "memory.h"

#include "nalu.h"
#include "memalloc.h"


/*! 
 *************************************************************************************
 * \brief
 *    Allocates memory for a NALU
 *
 * \param none
 *
 * \return
 *    pointer to a NALU
 *************************************************************************************
 */
 

NALU_t *AllocNALU()
{
  NALU_t *n;

  if ((n = calloc (1, sizeof (NALU_t))) == NULL) no_mem_exit ("AllocNALU");
  return n;
}


/*! 
 *************************************************************************************
 * \brief
 *    Frees a NALU
 *
 * \param NALU to be freed
 *
 * \return
 *    none
 *************************************************************************************
 */

void FreeNALU(NALU_t *n)
{
  free (n);
}


