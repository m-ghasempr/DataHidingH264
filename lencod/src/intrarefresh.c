
/*!
 *****************************************************************************
 *
 * \file intrarefresh.c
 *
 * \brief
 *    Encoder support for pseudo-random intra macroblock refresh
 *
 * \date
 *    16 June 2002
 *
 * \author
 *    Stephan Wenger   stewe@cs.tu-berlin.de
 *****************************************************************************/

#include "global.h"


/*!
 ************************************************************************
 * \brief
 *    RandomIntraInit: Initializes Random Intra module.  Should be called
 *    only after initialization (or changes) of the picture size or the
 *    random intra refresh value.  In version jm2.1 it is impossible to
 *    change those values on-the-fly, hence RandomIntraInit should be
 *    called immediately after the parsing of the config file
 *
 * \par Input:
 *    xsize, ysize: size of the picture (in MBs)
 *    refresh     : refresh rate in MBs per picture
 ************************************************************************
 */
void RandomIntraInit(ImageParameters *p_Img, int xsize, int ysize, int refresh)
{
  int i, pos;

  srand (1);      // A fixed random initializer to make things reproducible
  p_Img->NumberOfMBs = xsize * ysize;
  p_Img->NumberIntraPerPicture = refresh;

  if (refresh != 0)
  {
    p_Img->RefreshPattern = malloc (sizeof (int) * p_Img->NumberOfMBs);
    if (p_Img->RefreshPattern == NULL) no_mem_exit("RandomIntraInit: p_Img->RefreshPattern");

    p_Img->IntraMBs = malloc (sizeof (int) * refresh);
    if (p_Img->IntraMBs == NULL) no_mem_exit("RandomIntraInit: p_Img->IntraMBs");

    for (i= 0; i<p_Img->NumberOfMBs; i++)
      p_Img->RefreshPattern[i] = -1;

    for (i=0; i<p_Img->NumberOfMBs; i++)
    {
      do
      {
        pos = rand() % p_Img->NumberOfMBs;
      } while (p_Img->RefreshPattern [pos] != -1);
      p_Img->RefreshPattern [pos] = i;
    }
    /*
    for (i=0; i<p_Img->NumberOfMBs; i++) printf ("%d\t", p_Img->RefreshPattern[i]);
    getchar();
    */
  }
  else
  {
    p_Img->RefreshPattern = NULL;
    p_Img->IntraMBs = NULL;
  }
}

/*!
 ************************************************************************
 * \brief
 *    RandomIntra: Code an MB as Intra?
 *
 * \par Input
 *    MacroblockNumberInScanOrder
 * \par Output
 *    1 if an MB should be forced to Intra, according the the
 *      RefreshPattern
 *    0 otherwise
 *
 ************************************************************************
 */
int RandomIntra (ImageParameters *p_Img, int mb)
{
  int i;

  for (i=0; i<p_Img->NumberIntraPerPicture; i++)
    if (p_Img->IntraMBs[i] == mb)
      return 1;
  return 0;
}


/*!
 ************************************************************************
 * \brief
 *    RandomIntraNewPicture: Selects new set of MBs for forced Intra
 *
 * \par
 *    This function should be called exactly once per picture, and
 *    requires a finished initialization
 *
 ************************************************************************
 */
void RandomIntraNewPicture (ImageParameters *p_Img)
{
  int i, j;

  p_Img->WalkAround += p_Img->NumberIntraPerPicture;
  for (j=0, i = p_Img->WalkAround; j<p_Img->NumberIntraPerPicture; j++, i++)
    p_Img->IntraMBs[j] = p_Img->RefreshPattern [i%p_Img->NumberOfMBs];
}

void RandomIntraUninit(ImageParameters *p_Img)
{
  if (p_Img->NumberIntraPerPicture >0 )
  {
    free(p_Img->RefreshPattern);
    free(p_Img->IntraMBs);
  }
}
