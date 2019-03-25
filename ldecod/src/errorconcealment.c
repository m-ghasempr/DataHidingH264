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
 ***********************************************************************
 * \file errorconcealment.c
 *
 * \brief
 *    Implements error concealment scheme for H.26L decoder
 *
 * \date
 *    6.10.2000
 *
 * \version
 *    1.0
 *
 * \note
 *    This simple error concealment implemented in this decoder uses
 *    the existing dependencies of syntax elements.
 *    In case that an element is detected as false this elements and all
 *    dependend elements are marked as elements to conceal in the ec_flag[]
 *    array. If the decoder requests a new element by the function
 *    readSyntaxElement_xxxx() this array is checked first if an error concealment has
 *    to be applied on this element.
 *    In case that an error occured a concealed element is given to the
 *    decoding function in macroblock().
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Sebastian Purreiter   <sebastian.purreiter@mch.siemens.de>
 ***********************************************************************
 */

#include "contributors.h"
#include "global.h"
#include "elements.h"

static int ec_flag[SE_MAX_ELEMENTS];        //!< array to set errorconcealment
/*
static char SEtypes[][25] =
{
    "SE_HEADER",
    "SE_PTYPE",
    "SE_MBTYPE",
    "SE_REFFRAME",
    "SE_INTRAPREDMODE",
    "SE_MVD",
    "SE_CBP_INTRA",
    "SE_LUM_DC_INTRA",
    "SE_CHR_DC_INTRA",
    "SE_LUM_AC_INTRA",
    "SE_CHR_AC_INTRA",
    "SE_CBP_INTER",
    "SE_LUM_DC_INTER",
    "SE_CHR_DC_INTER",
    "SE_LUM_AC_INTER",
    "SE_CHR_AC_INTER",
    "SE_DELTA_QUANT_INTER",
    "SE_DELTA_QUANT_INTRA",
    "SE_BFRAME",
    "SE_EOS"
};
*/

/*!
 ***********************************************************************
 * \brief
 *    set concealment for all elements in same partition
 *    and dependend syntax elements
 * \return
 *    EC_REQ, elements of same type or depending type need error concealment. \n
 *    EX_SYNC   sync on next header
 ***********************************************************************
 */
int set_ec_flag(
  int se)   //!< type of syntax element to conceal
{

  /*
  if (ec_flag[se] == NO_EC)
    printf("Error concealment on element %s\n",SEtypes[se]);
  */
  switch (se)
  {
  case SE_HEADER :
    ec_flag[SE_HEADER] = EC_REQ;
  case SE_PTYPE :
    ec_flag[SE_PTYPE] = EC_REQ;
  case SE_MBTYPE :
    ec_flag[SE_MBTYPE] = EC_REQ;

  case SE_REFFRAME :
    ec_flag[SE_REFFRAME] = EC_REQ;
    ec_flag[SE_MVD] = EC_REQ; // set all motion vectors to zero length
    se = SE_CBP_INTER;      // conceal also Inter texture elements
    break;

  case SE_INTRAPREDMODE :
    ec_flag[SE_INTRAPREDMODE] = EC_REQ;
    se = SE_CBP_INTRA;      // conceal also Intra texture elements
    break;
  case SE_MVD :
    ec_flag[SE_MVD] = EC_REQ;
    se = SE_CBP_INTER;      // conceal also Inter texture elements
    break;

  default:
    break;
  }

  switch (se)
  {
  case SE_CBP_INTRA :
    ec_flag[SE_CBP_INTRA] = EC_REQ;
  case SE_LUM_DC_INTRA :
    ec_flag[SE_LUM_DC_INTRA] = EC_REQ;
  case SE_CHR_DC_INTRA :
    ec_flag[SE_CHR_DC_INTRA] = EC_REQ;
  case SE_LUM_AC_INTRA :
    ec_flag[SE_LUM_AC_INTRA] = EC_REQ;
  case SE_CHR_AC_INTRA :
    ec_flag[SE_CHR_AC_INTRA] = EC_REQ;
    break;

  case SE_CBP_INTER :
    ec_flag[SE_CBP_INTER] = EC_REQ;
  case SE_LUM_DC_INTER :
    ec_flag[SE_LUM_DC_INTER] = EC_REQ;
  case SE_CHR_DC_INTER :
    ec_flag[SE_CHR_DC_INTER] = EC_REQ;
  case SE_LUM_AC_INTER :
    ec_flag[SE_LUM_AC_INTER] = EC_REQ;
  case SE_CHR_AC_INTER :
    ec_flag[SE_CHR_AC_INTER] = EC_REQ;
    break;
  case SE_DELTA_QUANT_INTER :
    ec_flag[SE_DELTA_QUANT_INTER] = EC_REQ;
    break;
  case SE_DELTA_QUANT_INTRA :
    ec_flag[SE_DELTA_QUANT_INTRA] = EC_REQ;
    break;
  default:
    break;

  }
  return EC_REQ;
}

/*!
 ***********************************************************************
 * \brief
 *    resets EC_Flags called at the start of each slice
 *
 ***********************************************************************
 */
void reset_ec_flags()
{
  int i;
  for (i=0; i<SE_MAX_ELEMENTS; i++)
    ec_flag[i] = NO_EC;
}


/*!
 ***********************************************************************
 * \brief
 *    get error concealed element in dependence of syntax
 *    element se.                                                          \n
 *    This function implements the error concealment.
 * \return
 *    NO_EC if no error concealment required                               \n
 *    EC_REQ if element requires error concealment
 ***********************************************************************
 */
int get_concealed_element(SyntaxElement *sym)
{
  if (ec_flag[sym->type] == NO_EC)
    return NO_EC;
/*
#if TRACE
  printf("TRACE: get concealed element for %s!!!\n", SEtypes[sym->type]);
#endif
*/
  switch (sym->type)
  {
  case SE_HEADER :
    sym->len = 31;
    sym->inf = 0; // Picture Header
    break;

  case SE_PTYPE : // inter_img_1
  case SE_MBTYPE : // set COPY_MB
  case SE_REFFRAME :
    sym->len = 1;
    sym->inf = 0;
    break;

  case SE_INTRAPREDMODE :
  case SE_MVD :
    sym->len = 1;
    sym->inf = 0;  // set vector to zero length
    break;

  case SE_CBP_INTRA :
    sym->len = 5;
    sym->inf = 0; // codenumber 3 <=> no CBP information for INTRA images
    break;

  case SE_LUM_DC_INTRA :
  case SE_CHR_DC_INTRA :
  case SE_LUM_AC_INTRA :
  case SE_CHR_AC_INTRA :
    sym->len = 1;
    sym->inf = 0;  // return EOB
    break;

  case SE_CBP_INTER :
    sym->len = 1;
    sym->inf = 0; // codenumber 1 <=> no CBP information for INTER images
    break;

  case SE_LUM_DC_INTER :
  case SE_CHR_DC_INTER :
  case SE_LUM_AC_INTER :
  case SE_CHR_AC_INTER :
    sym->len = 1;
    sym->inf = 0;  // return EOB
    break;

  case SE_DELTA_QUANT_INTER:
    sym->len = 1;
    sym->inf = 0;
    break;
  case SE_DELTA_QUANT_INTRA:
    sym->len = 1;
    sym->inf = 0;
    break;
  default:
    break;
  }

  return EC_REQ;
}

