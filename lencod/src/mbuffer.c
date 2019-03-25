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
 *  \file
 *      mbuffer.c
 *
 *  \brief
 *      Frame buffer functions
 *
 *  \author
 *      Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Karsten Sühring          <suehring@hhi.de>
 ***********************************************************************
 */

#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "global.h"
#include "mbuffer.h"
#include "memalloc.h"

/*!
 ************************************************************************
 * \brief
 *    allocate memory for frame buffer
 *
 * \par Input:
 *    Input Parameters struct inp_par *inp, Image Parameters struct img_par *img
 *
 ************************************************************************
 */
void init_frame_buffers(InputParameters *inp, ImageParameters *img)
{
  int i;
  int bufsize=img->buf_cycle+1;    // Tian: PLUS1, June 6, 2002

  // Tian: the following line should be deleted! June 6, 2002
//  if ((fb=(FrameBuffer*)calloc(1,sizeof (FrameBuffer)))==NULL) no_mem_exit("init_frame_buffers: fb");

//frame buffers
  if ((frm=(FrameBuffer*)calloc(1,sizeof (FrameBuffer)))==NULL) no_mem_exit("init_frame_buffers: frm");

  if ((frm->picbuf_short=(Frame**)calloc(bufsize,sizeof (Frame*)))==NULL) no_mem_exit("init_frame_buffers: frm->picbuf_short");

  for (i=0;i<bufsize;i++)
    if ((frm->picbuf_short[i]=(Frame*)calloc(1,sizeof (Frame)))==NULL) no_mem_exit("init_frame_buffers: frm->picbuf_short");

  for (i=0;i<bufsize;i++)
  {
    get_mem2D(&(frm->picbuf_short[i]->mref), (img->height+2*IMG_PAD_SIZE)*4, (img->width+2*IMG_PAD_SIZE)*4);
    get_mem3D(&(frm->picbuf_short[i]->mcef), 2, img->height_cr, img->width_cr);
    if (NULL == (frm->picbuf_short[i]->Refbuf11 = malloc ((img->width * img->height) * sizeof (pel_t))))
      no_mem_exit ("init_frame_buffers: Refbuf11");

    if (input->WeightedPrediction || input->WeightedBiprediction)
    {
      get_mem2D(&(frm->picbuf_short[i]->mref_w), (img->height+2*IMG_PAD_SIZE)*4, (img->width+2*IMG_PAD_SIZE)*4);

      if (NULL == (frm->picbuf_short[i]->Refbuf11_w = malloc ((img->width * img->height) * sizeof (pel_t))))
        no_mem_exit ("init_frame_buffers: Refbuf11_w");
    }
  }


  if (input->WeightedPrediction > 0 || input->WeightedBiprediction > 0)
  {
    get_mem3Dint(&wp_weight, 2, MAX_REFERENCE_PICTURES, 3);
    get_mem3Dint(&wp_offset, 2, MAX_REFERENCE_PICTURES, 3);
    get_mem4Dint(&wbp_weight, 2, MAX_REFERENCE_PICTURES, MAX_REFERENCE_PICTURES, 3);

    get_mem2Dint(&weight,  MAX_REFERENCE_PICTURES, 3);
    get_mem2Dint(&offset, MAX_REFERENCE_PICTURES, 3);
  }

  frm->short_size=bufsize;
  frm->short_used=0;
  frm->long_size=0;
  frm->long_used=0;

//field buffers
  if(input->InterlaceCodingOption != FRAME_CODING)
  {
    bufsize *=2;
    if ((fld=(FrameBuffer*)calloc(1,sizeof (FrameBuffer)))==NULL) no_mem_exit("init_field_buffers: fld");
    
    if ((fld->picbuf_short=(Frame**)calloc(bufsize,sizeof (Frame*)))==NULL) no_mem_exit("init_field_buffers: fld->picbuf_short");
    
    for (i=0;i<bufsize;i++)
      if ((fld->picbuf_short[i]=(Frame*)calloc(1,sizeof (Frame)))==NULL) no_mem_exit("init_field_buffers: fld->picbuf_short");
      
    for (i=0;i<bufsize;i++)
    {
      get_mem2D(&(fld->picbuf_short[i]->mref), (img->height/2+2*IMG_PAD_SIZE)*4, (img->width+2*IMG_PAD_SIZE)*4);
      get_mem3D(&(fld->picbuf_short[i]->mcef), 2, img->height_cr/2, img->width_cr);
      if (NULL == (fld->picbuf_short[i]->Refbuf11 = malloc ((img->width * img->height / 2) * sizeof (pel_t))))
        no_mem_exit ("init_field_buffers: Refbuf11");
      if (input->WeightedPrediction || input->WeightedBiprediction)
      {
        get_mem2D(&(fld->picbuf_short[i]->mref_w), (img->height/2+2*IMG_PAD_SIZE)*4, (img->width+2*IMG_PAD_SIZE)*4);
        if (NULL == (fld->picbuf_short[i]->Refbuf11_w = malloc ((img->width * img->height / 2) * sizeof (pel_t))))
          no_mem_exit ("init_field_buffers: Refbuf11_w");
      }
    }
    
    fld->short_size=bufsize;
    fld->short_used=0;
    fld->long_size=0;
    fld->long_used=0;
    
    if(input->InterlaceCodingOption >= MB_CODING)
    {
      if ((mref_mbfld = (byte***)calloc(fld->short_size+fld->long_size,sizeof(byte**))) == NULL)
        no_mem_exit("alloc_mref: mref_mbfld");
      for(i=0;i<bufsize;i++)
        get_mem2D(&mref_mbfld[i], (img->height/2+2*IMG_PAD_SIZE)*4, (img->width+2*IMG_PAD_SIZE)*4);
      if (input->WeightedPrediction || input->WeightedBiprediction)
      {
        if ((mref_mbfld_w = (byte***)calloc(fld->short_size+fld->long_size,sizeof(byte**))) == NULL)
          no_mem_exit("alloc_mref: mref_mbfld_w");
        for(i=0;i<bufsize;i++)
          get_mem2D(&mref_mbfld_w[i], (img->height/2+2*IMG_PAD_SIZE)*4, (img->width+2*IMG_PAD_SIZE)*4);
      }
    }
    //  bufsize=img->buf_cycle+1;
  }
}

/*!
 ************************************************************************
 * \brief
 *    free frame buffer memory
 *
 * \par Input:
 *    Input Parameters struct inp_par *inp, Image Parameters struct img_par *img
 *
 ************************************************************************
 */
void free_frame_buffers(InputParameters *inp, ImageParameters *img)
{
  int i;

  for (i=0;i<frm->short_size;i++)
  {
    free_mem2D(frm->picbuf_short[i]->mref);
    free_mem3D(frm->picbuf_short[i]->mcef,2);
    free(frm->picbuf_short[i]->Refbuf11);
    if (input->WeightedPrediction || input->WeightedBiprediction)
    {
      free_mem2D(frm->picbuf_short[i]->mref_w);
      free(frm->picbuf_short[i]->Refbuf11_w);
    }
  }


  for (i=0;i<frm->short_size;i++)
    free (frm->picbuf_short[i]);

  free (frm->picbuf_short);

  if (frm->picbuf_long)
  {
    for (i=0;i<frm->long_size;i++)
    {
      free_mem2D(frm->picbuf_long[i]->mref);
      free_mem3D(frm->picbuf_long[i]->mcef,2);
      free(frm->picbuf_long[i]->Refbuf11);
      if (input->WeightedPrediction || input->WeightedBiprediction)
      {
        free_mem2D(frm->picbuf_long[i]->mref_w);
        free(frm->picbuf_long[i]->Refbuf11_w);
      }
    } 
    for (i=0;i<frm->long_size;i++)
      free (frm->picbuf_long[i]);

  }

  free (frm);

  if(input->InterlaceCodingOption != FRAME_CODING)
  {
    for (i=0;i<fld->short_size;i++)
    {
      free_mem2D(fld->picbuf_short[i]->mref);
      free_mem3D(fld->picbuf_short[i]->mcef,2);
      free(fld->picbuf_short[i]->Refbuf11);
      if (input->WeightedPrediction || input->WeightedBiprediction)
      {
        free_mem2D(fld->picbuf_short[i]->mref_w);
        free(fld->picbuf_short[i]->Refbuf11_w);
      }
    }


    for (i=0;i<fld->short_size;i++)
      free (fld->picbuf_short[i]);

    free (fld->picbuf_short);

    if (fld->picbuf_long)
    {
      for (i=0;i<fld->long_size;i++)
      {
        free_mem2D(fld->picbuf_long[i]->mref);
        free_mem3D(fld->picbuf_long[i]->mcef,2);
        free(fld->picbuf_long[i]->Refbuf11);
        if (input->WeightedPrediction || input->WeightedBiprediction)
        {
          free_mem2D(fld->picbuf_long[i]->mref_w);
          free(fld->picbuf_long[i]->Refbuf11_w);
        }
      } 
      for (i=0;i<fld->long_size;i++)
        free (fld->picbuf_long[i]);
    }

    free (fld);
  }
}

/*!
 ************************************************************************
 * \brief
 *    allocate or reallocate long term buffer memory
 *
 * \param size
 *    number of frames
 *
 ************************************************************************
 */
void init_long_term_buffer(int size, ImageParameters *img)
{
  Frame **pb;
  int i;

  /* nothing to do if size unchanged */
  if (fb->long_size==size) return;

  if (fb->long_size>size)
  {
    /* shrink size */

    pb=fb->picbuf_long;
    if ((fb->picbuf_long=(Frame**)calloc(size,sizeof (Frame*)))==NULL) no_mem_exit("init_long_term_buffer: fb->picbuf_long");

    for (i=0;i<size;i++)
      fb->picbuf_long[i]=pb[i];

    for (i=size;i<fb->long_size;i++)
    {
      free_mem2D(pb[i]->mref);
      free_mem3D(pb[i]->mcef, 2);
      free(pb[i]->Refbuf11);
      if (input->WeightedPrediction || input->WeightedBiprediction)
      {
        free_mem2D(pb[i]->mref_w);
        free(pb[i]->Refbuf11_w);
      }
      free (pb[i]);
    }

    free (pb);

  } else
  {
    /* grow size */

    pb=fb->picbuf_long;

    if ((fb->picbuf_long=(Frame**)calloc(size,sizeof (Frame*)))==NULL) no_mem_exit("init_frame_buffers: fb->picbuf_long");

    for (i=0;i<fb->long_size;i++)
      fb->picbuf_long[i]=pb[i];

    for (i=fb->long_size;i<size;i++)
    {
      if ((fb->picbuf_long[i]=(Frame*)calloc(1,sizeof (Frame)))==NULL) no_mem_exit("init_frame_buffers: fb->picbuf_long");
      get_mem2D(&(fb->picbuf_long[i]->mref), (img->height+2*IMG_PAD_SIZE)*4, (img->width+2*IMG_PAD_SIZE)*4);
      if (input->WeightedPrediction || input->WeightedBiprediction)
        get_mem2D(&(fb->picbuf_long[i]->mref_w), (img->height+2*IMG_PAD_SIZE)*4, (img->width+2*IMG_PAD_SIZE)*4);
      get_mem3D(&(fb->picbuf_long[i]->mcef), 2, img->height_cr, img->width_cr);
      if (NULL == (fb->picbuf_long[i]->Refbuf11 = malloc ((img->width * img->height) * sizeof (pel_t))))
        no_mem_exit ("init_long_term_buffer: Refbuf11");
      if (input->WeightedPrediction || input->WeightedBiprediction)
        if (NULL == (fb->picbuf_long[i]->Refbuf11_w = malloc ((img->width * img->height) * sizeof (pel_t))))
          no_mem_exit ("init_long_term_buffer: Refbuf11_w");
    }

    free (pb);

  }

  // finally set the size
  fb->long_size=size;

  // and correct the size of mref/mcef, Refbuf11
  alloc_mref(img);
  alloc_Refbuf(img);
}

/*!
 ************************************************************************
 * \brief
 *    assign a long term index to a short term picture
 *
 * \param shortID
 *    short term ID
 *
 * \param longID
 *    long term ID
 *    
 ************************************************************************
 */
void assign_long_term_id(int shortID, int longID, ImageParameters *img)
{
  int i, uv;
  int idx=-1;
  Frame *f;
  int j,idx1,id;

  // try to find short term picture in buffer
  for (i=0;i<fb->short_used;i++)
  {
    if (fb->picbuf_short[i]->picID==shortID)
      idx=i;
  }

  // short term picture not found: error
  if (idx==-1)
    error ("Trying to assign a long term ID to a not existent short term picture",400);

  // return if reassignment (no error)
  if (fb->picbuf_short[idx]->lt_picID==longID)
    return;

  // copy the frame to the long term buffer
  // this may not be the best solution, but avoids messing up pointers

  printf("Storing long term ID=%d\n",longID);

  /* delete frames with same short term ID */
  remove_long_term(longID);


  /* cycle frame buffer */
  f=fb->picbuf_long[fb->long_size-1];

  for (i=fb->long_size-2;i>=0;i--)
    fb->picbuf_long[i+1]=fb->picbuf_long[i];

  fb->picbuf_long[0]=f;


  fb->picbuf_long[0]->used=1;
  fb->picbuf_long[0]->picID=shortID;
  fb->picbuf_long[0]->lt_picID=longID;

  fb->picbuf_short[idx]->lt_picID=longID;

  (fb->long_used)++;
  if (fb->long_used>fb->long_size)
    fb->long_used=fb->long_size;
  
  memcpy (fb->picbuf_long[0]->mref[0], fb->picbuf_short[idx]->mref[0], ((img->height+2*IMG_PAD_SIZE)*4)*((img->width+2*IMG_PAD_SIZE)*4)); 
  if (input->WeightedPrediction || input->WeightedBiprediction)
     memcpy (fb->picbuf_long[0]->mref_w[0], fb->picbuf_short[idx]->mref_w[0], ((img->height+2*IMG_PAD_SIZE)*4)*((img->width+2*IMG_PAD_SIZE)*4)); 

  for (uv=0; uv < 2; uv++)
    memcpy (fb->picbuf_long[0]->mcef[uv][0], fb->picbuf_short[idx]->mcef[uv][0], img->width_cr*img->height_cr); 

  memcpy (fb->picbuf_long[0]->Refbuf11, fb->picbuf_short[idx]->Refbuf11, img->width*img->height); 
  if (input->WeightedPrediction || input->WeightedBiprediction)
        memcpy (fb->picbuf_long[0]->Refbuf11, fb->picbuf_short[idx]->Refbuf11, img->width*img->height); 

  // sort the long term buffer by lt_IDs
  // Note: Just a  simple bubble sort implementation. The buffer is not very large
  //       so this should not take too much time. Feel free to implement something better.

  for (i=0,idx1=0;i<fb->long_used-1;i++)
  {
    id=fb->picbuf_long[i]->lt_picID;

    for (j=i+1;j<fb->long_used;j++)
    {
      if (fb->picbuf_long[j]->lt_picID<id) 
      {
        id=fb->picbuf_long[j]->lt_picID;
        idx1=j;
      }
    }
    if (id<fb->picbuf_long[i]->lt_picID)
    {
      f=fb->picbuf_long[i];
      fb->picbuf_long[i]=fb->picbuf_long[idx1];
      fb->picbuf_long[idx1]=f;
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    remove long term frame from buffer
 *
 * \param longID
 *    long term ID
 *    
 ************************************************************************
 */
void remove_long_term(int longID)
{
  int i,j;
  Frame *f;

  for (i=0;i<fb->long_used;i++)
  {
    if ((fb->picbuf_long[i]->lt_picID)==longID)
    {
      fb->picbuf_long[i]->used=0;
      fb->picbuf_long[i]->picID=-1;
      fb->picbuf_long[i]->lt_picID=-1;
      
      fb->long_used--;

      if (i<fb->long_used) 
      {
        f=fb->picbuf_long[i];
        for (j=i;j<fb->long_used;j++)
          fb->picbuf_long[j]=fb->picbuf_long[j+1];
        fb->picbuf_long[fb->long_used]=f;
      } 
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    remove short term frame from buffer
 *
 * \param shortID
 *    short term ID
 *    
 ************************************************************************
 */
void remove_short_term(int shortID)
{
  int i,j;
  Frame *f;

  for (i=0;i<fb->short_used;i++)
  {
    if ((fb->picbuf_short[i]->picID)==shortID)
    {
      fb->picbuf_short[i]->used=0;
      fb->picbuf_short[i]->lt_picID=-1;
      fb->picbuf_short[i]->picID=-1;
      
      fb->short_used--;

      // Tian Dong: June 15, 2002
      // Use short_size to replace short_used.
      if (i<fb->short_size) 
      {
        f=fb->picbuf_short[i];
        for (j=i;j<fb->short_size-1;j++)
          fb->picbuf_short[j]=fb->picbuf_short[j+1];
        fb->picbuf_short[fb->short_size-1]=f;
      } 
      // old lines:
/*      if (i<fb->short_size) 
      {
        f=fb->picbuf_short[i];
        for (j=i;j<fb->short_used-1;j++);
          fb->picbuf_short[j]=fb->picbuf_short[j+1];
        fb->picbuf_short[fb->short_used-1]=f;
      } */
    }
  }
}

/*!
 ************************************************************************
 * \brief
 *    (re)-allocate memory for mref pointer structure
 ************************************************************************
 */
void alloc_mref(ImageParameters *img)
{

  if (mref_frm) 
    free (mref_frm);

  if((mref_frm = (byte***)calloc(frm->short_size+frm->long_size,sizeof(byte**))) == NULL)
    no_mem_exit("alloc_mref: mref_frm");


  if (mcef_frm) 
    free (mcef_frm);

  if((mcef_frm = (byte****)calloc(frm->short_size+frm->long_size,sizeof(byte***))) == NULL)
    no_mem_exit("alloc_mref: mcef_frm");

  if(input->InterlaceCodingOption != FRAME_CODING)
  {
    if (mref_fld) 
      free (mref_fld);

    if ((mref_fld = (byte***)calloc(fld->short_size+fld->long_size,sizeof(byte**))) == NULL)
      no_mem_exit("alloc_mref: mref_fld");

    if (mcef_fld) 
      free (mcef_fld);

    if ((mcef_fld = (byte****)calloc(fld->short_size+fld->long_size,sizeof(byte***))) == NULL)
      no_mem_exit("alloc_mref: mcef_fld");
  }


  if (input->WeightedPrediction || input->WeightedBiprediction)
  {
    if (mref_frm_w) 
      free (mref_frm_w);

    if((mref_frm_w = (byte***)calloc(frm->short_size+frm->long_size,sizeof(byte**))) == NULL)
      no_mem_exit("alloc_mref: mref_frm_w");

    if(input->InterlaceCodingOption != FRAME_CODING)
    {
      if (mref_fld_w) 
        free (mref_fld_w);

      if ((mref_fld_w = (byte***)calloc(fld->short_size+fld->long_size,sizeof(byte**))) == NULL)
        no_mem_exit("alloc_mref: mref_fld_w");
    }
  }

}

/*!
 ************************************************************************
 * \brief
 *    init mref pointer structure
 ************************************************************************
 */
void init_mref()
{
  int i,j;

  for (i=0,j=0;i<fb->short_used;i++,j++)
  {
    mref[j]=fb->picbuf_short[i]->mref;
    if (input->WeightedPrediction || input->WeightedBiprediction)
      mref_w[j]=fb->picbuf_short[i]->mref_w;
    mcef[j]=fb->picbuf_short[i]->mcef;
  }
  for (i=0;i<fb->long_used;i++,j++)
  {
    mref[j]=fb->picbuf_long[i]->mref;
    if (input->WeightedPrediction || input->WeightedBiprediction)
      mref_w[j]=fb->picbuf_long[i]->mref_w;
    mcef[j]=fb->picbuf_long[i]->mcef;
  }

  // set all other mref pointers to NULL !KS!
  for (;j<fb->long_size+fb->short_size;j++)
  {
    mref[j]=NULL;
    if (input->WeightedPrediction || input->WeightedBiprediction)
      mref_w[j]=NULL;
    mcef[j]=NULL;
  }
  
}

/*!
 ************************************************************************
 * \brief
 *    process RMPNI command list
 ************************************************************************
 */
void reorder_mref(ImageParameters *img)
{
  RMPNIbuffer_t *r;

  int pnp = img->pn;
  int pnq = 0;
  int index = 0,i;
  int size;
  int islong=0;
  int found=0;
  
  Frame **fr;
  Frame *f;

  // if nothing to do update mref and return
  if (img->currentSlice->rmpni_buffer== NULL) 
  {
    init_mref();
    init_Refbuf(img);
    return;
  }

  // We have to do the reordering with complete frame structures
  // so we need a temporary frame buffer of ordered references

  size=fb->short_used+fb->long_used;
  if ((fr=(Frame**)calloc(size,sizeof (Frame*)))==NULL) no_mem_exit("init_frame_buffers: fb->picbuf_long");

  for (i=0;i<fb->short_used;i++)
  {
    fr[i]=fb->picbuf_short[i];
    fr[i]->islong=0;
  }

  for (i=0;i<fb->long_used;i++)
  {
    fr[i+fb->short_used]=fb->picbuf_long[i];
    fr[i+fb->short_used]->islong=1;
  }

  r=img->currentSlice->rmpni_buffer;
  while (r)
  {
    switch (r->RMPNI)
    {

    case 0:
      // ADPN correspondig to negative difference
      if ((pnp-r->Data)<0)
        pnq=pnp-r->Data+fb->short_size;
      else
        pnq=pnp-r->Data;

      islong=0;

      break;

    case 1:
      // ADPN correspondig to positive difference
      if ((pnp+r->Data)>fb->short_size-1)
        pnq=pnp+r->Data-fb->short_size;
      else
        pnq=pnp+r->Data;

      islong=0;

      break;

    case 2:
      // LPIR specifying long term index
      pnq=r->Data;
      islong=1;

      break;

    case 3:
//      return; // Tian: I don't think it should return here. June 6, 2002
      assert(r->Next == NULL);  // Tian: for debug purpose, try to assert ... June 6, 2002
      break;

    }

    // now scan the frame buffer for the needed frame
    if ( r->RMPNI == 0 || r->RMPNI == 1 || r->RMPNI == 2 )  // the IF is added by Tian
    {
      found=0;
      i=0;

      while ((!found)&&(i<size))
      {
        if (((!islong)&&(!fr[i]->islong) && (fr[i]->picID==pnq)) ||
            ((islong)&&(fr[i]->islong) && (fr[i]->lt_picID==pnq)))
        {
            found=1;
            index=i;
        }
        i++;
      }

      if (!found) error ("tried to remap non-existent picture",400);

      // now do the actual reordering
      /* cycle frame buffer */
      f=fr[index];
      for (i=index-1;i>=0;i--)
      {
        fr[i+1]=fr[i];
      }
      fr[0]=f;

      // set the picture number prediction correctly
      pnp=pnq;
    }

    img->currentSlice->rmpni_buffer = r->Next;
    free( r );
    r = img->currentSlice->rmpni_buffer;
  }

  // at last init mref, mcef and Refbuf from temporary structure
  for (i=0;i<size;i++)
  {
    mref[i]=fr[i]->mref;
    mcef[i]=fr[i]->mcef;
    Refbuf11[i]=fr[i]->Refbuf11;
    if (input->WeightedPrediction || input->WeightedBiprediction)
    {
      mref_w[i]=fr[i]->mref_w;
      Refbuf11_w[i]=fr[i]->Refbuf11_w;
    }
  }

  // free temporary memory
  free (fr);

}

/*!
 ************************************************************************
 * \brief
 *    (re)-allocate memory for Refbuf11 pointer structure
 ************************************************************************
 */
void alloc_Refbuf (ImageParameters *img)
{
 // int num_frames = fb->long_size+fb->short_size; 
 //   if (Refbuf11) free (Refbuf11);

  int num_frames = frm->long_size+frm->short_size;
 
  if (Refbuf11_frm) free (Refbuf11_frm);
  if (input->WeightedBiprediction || input->WeightedPrediction)
        if (Refbuf11_frm_w) free (Refbuf11_frm_w);

  if (NULL == (Refbuf11_frm = malloc (num_frames * sizeof (pel_t *))))
    no_mem_exit ("alloc_Refbuf: Refbuf11_frm");
  
  if (input->WeightedBiprediction || input->WeightedPrediction)
          if (NULL == (Refbuf11_frm_w = malloc (num_frames * sizeof (pel_t *))))
                no_mem_exit ("alloc_Refbuf: Refbuf11_frm_w");
  
  if(input->InterlaceCodingOption != FRAME_CODING)
  {
    num_frames = fld->long_size+fld->short_size;
    if (Refbuf11_fld) free (Refbuf11_fld);

    if (NULL == (Refbuf11_fld = malloc (num_frames * sizeof (pel_t *))))
      no_mem_exit ("alloc_Refbuf: Refbuf11_fld");
        if (input->WeightedBiprediction || input->WeightedPrediction)
            if (NULL == (Refbuf11_fld_w = malloc (num_frames * sizeof (pel_t *))))
                        no_mem_exit ("alloc_Refbuf: Refbuf11_fld_w");
  }
}

/*!
 ************************************************************************
 * \brief
 *    init Refbuf pointer structure
 ************************************************************************
 */
void init_Refbuf(ImageParameters *img)
{
  int i,j;

  for (i=0,j=0;i<fb->short_used;i++)
  {
    Refbuf11[j]=fb->picbuf_short[i]->Refbuf11;
        if (input->WeightedBiprediction || input->WeightedPrediction)
                Refbuf11_w[j]=fb->picbuf_short[i]->Refbuf11_w;
    j++;
  }
  for (i=0;i<fb->long_size;i++)
  {
    Refbuf11[j]=fb->picbuf_long[i]->Refbuf11;
        if (input->WeightedBiprediction || input->WeightedPrediction)
            Refbuf11_w[j]=fb->picbuf_long[i]->Refbuf11_w;
    j++;
  }
  for (;j<fb->long_size+fb->short_size;j++)
  {
    Refbuf11[j]=NULL;
        if (input->WeightedBiprediction || input->WeightedPrediction)
                Refbuf11_w[j]=NULL;
  }
}


/*!
 ************************************************************************
 * \brief
 *    mark all frames except last decoded unused
 ************************************************************************
 */
void reset_buffers()
{
  int i;

  /* mark short term pictures unused */
  for (i=1;i<fb->short_used;i++)
    fb->picbuf_short[i]->used=0;
  fb->short_used=1;

  /* mark short term pictures unused */
  for (i=0;i<fb->long_used;i++)
    fb->picbuf_short[i]->used=0;
  fb->long_used=0;

}


/*!
 ************************************************************************
 * \brief
 *    store reconstructed frame in multiple reference buffer
 ************************************************************************
 */
void add_frame(ImageParameters *img)
{
  int i;
  Frame *f;
  Boolean remove_last = FALSE;

  /* delete frames with same short term ID */
  if(mref==mref_frm) remove_short_term(img->pn);

  if (fb->short_used == fb->short_size && fb==frm) remove_last = TRUE;  // Tian: PLUS1

  /* cycle frame buffer */
  f=fb->picbuf_short[fb->short_size-1];

  for (i=fb->short_size-2;i>=0;i--)
    fb->picbuf_short[i+1]=fb->picbuf_short[i];

  fb->picbuf_short[0]=f;


  fb->picbuf_short[0]->used=1;
  fb->picbuf_short[0]->picID=img->pn;
  fb->picbuf_short[0]->frame_num_256=img->number % 256;   // Tian Dong (Sept 2002)
  fb->picbuf_short[0]->lt_picID=-1;

  // indicate which layer and sub-seq current ref frame comes from.
//  fb->picbuf_short[0]->layer_no=currPictureInfo.refFromLayerNumber;
//  fb->picbuf_short[0]->sub_seq_no=currPictureInfo.refFromSubSequenceIdentifier;

  (fb->short_used)++;
  if (fb->short_used>fb->short_size)
    fb->short_used=fb->short_size;
  if (remove_last)
  {
    fb->picbuf_short[fb->short_size-1]->used=0;
    fb->picbuf_short[fb->short_size-1]->picID=-1;
    fb->picbuf_short[fb->short_size-1]->frame_num_256=-1;   // Tian Dong (Sept 2002)
    fb->picbuf_short[fb->short_size-1]->lt_picID=-1;
    fb->short_used--;
    printf("remove last........................\n");
  }
}

void copy_mref()
{
  int i,j,k;

  // For MB level field/frame coding -- 
  // temporary solution till we remove mref==mref_fld from pic level 
  // then we can use mref as the reference field buffer                                     
  if(input->InterlaceCodingOption >= MB_CODING)
  {
    for (i=0;i<fld->short_used;i++)
      for (j=0;j<4*(input->img_height/2 + 2*IMG_PAD_SIZE);j++)
        for(k=0;k<4*(img->width + 2*IMG_PAD_SIZE);k++)
        {
          mref_mbfld[i][j][k] = fld->picbuf_short[i]->mref[j][k];
          if (input->WeightedPrediction || input->WeightedBiprediction)
            mref_mbfld_w[i][j][k] = fld->picbuf_short[i]->mref_w[j][k];
        }

    // still need to do this regardless of mref==mref_fld is used or not
    if(mref != mref_fld)
    {
      for (i=0,j=0;i<fld->short_used;i++)
      {
        mref_fld[j]=fld->picbuf_short[i]->mref;
        if (input->WeightedPrediction || input->WeightedBiprediction)
          mref_fld_w[j]=fld->picbuf_short[i]->mref_w;
        mcef_fld[j]=fld->picbuf_short[i]->mcef;
        j++;
      }
      for (i=0;i<fld->long_used;i++)
      {
        mref_fld[j]=fld->picbuf_long[i]->mref;
        if (input->WeightedPrediction || input->WeightedBiprediction)
          mref_fld_w[j]=fld->picbuf_long[i]->mref_w;
        mcef_fld[j]=fld->picbuf_long[i]->mcef;
        j++;
      }
      
      // set all other mref pointers to NULL !KS!
      for (;j<fld->long_size+fld->short_size;j++)
      {
        mref_fld[j]=NULL;
        if (input->WeightedPrediction || input->WeightedBiprediction)
          mref_fld_w[j]=NULL;
        mcef_fld[j]=NULL;
      }
      
      // need to do this, because MB level frame/field overwrites the current Refbuf11 
      for (i=0,j=0;i<fld->short_used;i++)
      {
        Refbuf11_fld[j]=fld->picbuf_short[i]->Refbuf11;
        if (input->WeightedBiprediction || input->WeightedPrediction)
                   Refbuf11_fld_w[j]=fld->picbuf_short[i]->Refbuf11_w;
        j++;
      }
      for (i=0;i<fld->long_size;i++)
      {
        Refbuf11_fld[j]=fld->picbuf_long[i]->Refbuf11;
        if (input->WeightedBiprediction || input->WeightedPrediction)
              Refbuf11_fld_w[j]=fld->picbuf_long[i]->Refbuf11_w;
        j++;
      }
      for (;j<fld->long_size+fld->short_size;j++)
      {
        Refbuf11_fld[j]=NULL;
        if (input->WeightedBiprediction || input->WeightedPrediction)
                Refbuf11_fld_w[j]=NULL;
      }
    }
  }
}


void alloc_ref_pic_list_reordering_buffer(Slice *currSlice)
{
  int size = img->num_ref_pic_active_fwd_minus1+1;

  if (img->type!=INTRA_IMG /* && img->type!=SI_IMG */)
  {
    if ((currSlice->remapping_of_pic_nums_idc_l0 = calloc(size,sizeof(int)))==NULL) no_mem_exit("alloc_ref_pic_list_reordering_buffer: remapping_of_pic_nums_idc_l0");
    if ((currSlice->abs_diff_pic_num_minus1_l0 = calloc(size,sizeof(int)))==NULL) no_mem_exit("alloc_ref_pic_list_reordering_buffer: abs_diff_pic_num_minus1_l0");
    if ((currSlice->long_term_pic_idx_l0 = calloc(size,sizeof(int)))==NULL) no_mem_exit("alloc_ref_pic_list_reordering_buffer: long_term_pic_idx_l0");
  }
  else
  {
    currSlice->remapping_of_pic_nums_idc_l0 = NULL;
    currSlice->abs_diff_pic_num_minus1_l0 = NULL;
    currSlice->long_term_pic_idx_l0 = NULL;
  }
  
  size = img->num_ref_pic_active_bwd_minus1+1;

  if (img->type!=B_IMG)
  {
    if ((currSlice->remapping_of_pic_nums_idc_l1 = calloc(size,sizeof(int)))==NULL) no_mem_exit("alloc_ref_pic_list_reordering_buffer: remapping_of_pic_nums_idc_l1");
    if ((currSlice->abs_diff_pic_num_minus1_l1 = calloc(size,sizeof(int)))==NULL) no_mem_exit("alloc_ref_pic_list_reordering_buffer: abs_diff_pic_num_minus1_l1");
    if ((currSlice->long_term_pic_idx_l1 = calloc(size,sizeof(int)))==NULL) no_mem_exit("alloc_ref_pic_list_reordering_buffer: long_term_pic_idx_l1");
  }
  else
  {
    currSlice->remapping_of_pic_nums_idc_l1 = NULL;
    currSlice->abs_diff_pic_num_minus1_l1 = NULL;
    currSlice->long_term_pic_idx_l1 = NULL;
  }
}


void free_ref_pic_list_reordering_buffer(Slice *currSlice)
{

  if (currSlice->remapping_of_pic_nums_idc_l0) 
    free(currSlice->remapping_of_pic_nums_idc_l0);
  if (currSlice->abs_diff_pic_num_minus1_l0)
    free(currSlice->abs_diff_pic_num_minus1_l0);
  if (currSlice->long_term_pic_idx_l0)
    free(currSlice->long_term_pic_idx_l0);

  currSlice->remapping_of_pic_nums_idc_l0 = NULL;
  currSlice->abs_diff_pic_num_minus1_l0 = NULL;
  currSlice->long_term_pic_idx_l0 = NULL;
  
  if (currSlice->remapping_of_pic_nums_idc_l1)
    free(currSlice->remapping_of_pic_nums_idc_l1);
  if (currSlice->abs_diff_pic_num_minus1_l1)
    free(currSlice->abs_diff_pic_num_minus1_l1);
  if (currSlice->long_term_pic_idx_l1)
    free(currSlice->long_term_pic_idx_l1);
  
  currSlice->remapping_of_pic_nums_idc_l1 = NULL;
  currSlice->abs_diff_pic_num_minus1_l1 = NULL;
  currSlice->long_term_pic_idx_l1 = NULL;
}

