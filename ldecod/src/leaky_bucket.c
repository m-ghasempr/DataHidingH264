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
 * \file  leaky_bucket.c
 *
 * \brief
 *   Calculate if decoder leaky bucket parameters meets HRD constraints specified by encoder. 
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Shankar Regunathan                   <shanre@microsoft.com>
 ************************************************************************
 */

#include "contributors.h"
#include "global.h"
#include "stdlib.h"

#ifdef _LEAKYBUCKET_
/*!
 ***********************************************************************
 * \brief
 *   Function to get unsigned long word from a file.
 * \param fp
 *    Filepointer
 * \return
 *    unsigned long double word
 * \par SideEffects
 *     None.
 *  \par Notes
 *     File should be opened to read in binary format.
 * \author
 *    Shankar Regunathan                   shanre@microsoft.com
 *  \date 
 *      December 06, 2001.
 ***********************************************************************
 */
/* gets unsigned double stored in Big Endian Order */
unsigned long GetBigDoubleWord(FILE *fp)
{
  register unsigned long dw;
  dw =  (unsigned long) (fgetc(fp) & 0xFF);
  dw = ((unsigned long) (fgetc(fp) & 0xFF)) | (dw << 0x08);
  dw = ((unsigned long) (fgetc(fp) & 0xFF)) | (dw << 0x08);
  dw = ((unsigned long) (fgetc(fp) & 0xFF)) | (dw << 0x08);
  return(dw);
}

/*!
 ***********************************************************************
 * \brief
 *   Calculates if decoder leaky bucket parameters meets HRD constraints specified by encoder.
 * \param inp
 *    Structure which contains decoder leaky bucket parameters.   
 * \return
 *    None
 * \par SideEffects
 *     None.
 * \par Notes
 *     Failure if LeakyBucketParam file is missing or if it does not have
 *     the correct number of entries.
 * \author
 *    Shankar Regunathan                   shanre@microsoft.com
 *  \date 
 *      December 06, 2001.
 ***********************************************************************
 */

/* Main Routine to verify HRD compliance */
void calc_buffer(struct inp_par *inp)
{
  unsigned long NumberLeakyBuckets, *Rmin, *Bmin, *Fmin;
  float B_interp,  F_interp;
  unsigned long iBucket;
  float dnr, frac1, frac2;
  unsigned long R_decoder, B_decoder, F_decoder;
  FILE *outf;
        
  if ((outf=fopen(inp->LeakyBucketParamFile,"rb"))==NULL)
    {
    snprintf(errortext, ET_SIZE, "Error open file %s \n",inp->LeakyBucketParamFile);
    error(errortext,1);
    }

  NumberLeakyBuckets = GetBigDoubleWord(outf);
  printf(" Number Leaky Buckets: %8ld \n\n", NumberLeakyBuckets);
  Rmin = calloc(sizeof(unsigned long), NumberLeakyBuckets);
  Bmin = calloc(sizeof(unsigned long), NumberLeakyBuckets);
  Fmin = calloc(sizeof(unsigned long), NumberLeakyBuckets);

  for(iBucket =0; iBucket < NumberLeakyBuckets; iBucket++) 
  {
    Rmin[iBucket] = GetBigDoubleWord(outf);
    Bmin[iBucket] = GetBigDoubleWord(outf);
    Fmin[iBucket] = GetBigDoubleWord(outf);
    printf(" %8ld %8ld %8ld \n", Rmin[iBucket], Bmin[iBucket], Fmin[iBucket]);
  }
  fclose(outf);

  R_decoder = inp->R_decoder;
  F_decoder = inp->F_decoder;
  B_decoder = inp->B_decoder;

  for( iBucket =0; iBucket < NumberLeakyBuckets; iBucket++) 
  {
    if(R_decoder < Rmin[iBucket])
      break;
  }

  printf("\n");
  if(iBucket > 0 ) {
    if(iBucket < NumberLeakyBuckets) {
      dnr = (float) (Rmin[iBucket] - Rmin[iBucket-1]);
      frac1 = (float) (R_decoder - Rmin[iBucket-1]);
      frac2 = (float) (Rmin[iBucket] - R_decoder);
      B_interp = (float) (Bmin[iBucket] * frac1 + Bmin[iBucket-1] * frac2) /dnr;
      F_interp = (float) (Fmin[iBucket] * frac1 + Fmin[iBucket-1] * frac2) /dnr;
    }
    else {
      B_interp = (float) Bmin[iBucket-1];
      F_interp = (float) Fmin[iBucket-1];
    }
    printf(" Min.buffer %8.2f Decoder buffer size %ld \n Minimum Delay %8.2f DecoderDelay %ld \n", B_interp, B_decoder, F_interp, F_decoder);
    if(B_decoder > B_interp && F_decoder > F_interp)
      printf(" HRD Compliant \n"); 
    else
      printf(" HRD Non Compliant \n");
  }
  else { // (iBucket = 0)
    printf(" Decoder Rate is too small; HRD cannot be verified \n");
  }
  
  free(Rmin);
  free(Bmin);
  free(Fmin);
  return;
}
#endif
