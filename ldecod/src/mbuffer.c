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
void init_frame_buffers(struct inp_par *inp, ImageParameters *img)
{
  int i;

  img->buf_cycle = inp->buf_cycle+1;

  if (img->structure != FRAME)
  {
    img->height *= 2;
    img->height_cr *= 2;
  }

  //frame buffers
  if ((frm=(FrameBuffer*)calloc(1,sizeof (FrameBuffer)))==NULL) no_mem_exit("init_frame_buffers: frm");

  if ((frm->picbuf_short=(Frame**)calloc(img->buf_cycle,sizeof (Frame*)))==NULL) no_mem_exit("init_frame_buffers: frm->picbuf_short");

  for (i=0;i<img->buf_cycle;i++)
    if ((frm->picbuf_short[i]=(Frame*)calloc(1,sizeof (Frame)))==NULL) no_mem_exit("init_frame_buffers: frm->picbuf_short");

  for (i=0;i<img->buf_cycle;i++)
  {
    get_mem2D(&(frm->picbuf_short[i]->mref), img->height, img->width);
    get_mem3D(&(frm->picbuf_short[i]->mcef), 2, img->height_cr, img->width_cr);

    frm->picbuf_short[i]->picID = -1;
    frm->picbuf_short[i]->lt_picID = -1;
    frm->picbuf_short[i]->used = 0;
  }


  frm->short_size=img->buf_cycle;
  frm->short_used=0;
  frm->long_size=0;
  frm->long_used=0;

//field buffers
  img->buf_cycle=2*img->buf_cycle;
  if ((fld=(FrameBuffer*)calloc(1,sizeof (FrameBuffer)))==NULL) no_mem_exit("init_field_buffers: fld");

  if ((fld->picbuf_short=(Frame**)calloc(img->buf_cycle,sizeof (Frame*)))==NULL) no_mem_exit("init_field_buffers: fld->picbuf_short");

  for (i=0;i<img->buf_cycle;i++)
    if ((fld->picbuf_short[i]=(Frame*)calloc(1,sizeof (Frame)))==NULL) no_mem_exit("init_field_buffers: fld->picbuf_short");

  for (i=0;i<img->buf_cycle;i++)
  {
    get_mem2D(&(fld->picbuf_short[i]->mref), img->height/2, img->width);
    get_mem3D(&(fld->picbuf_short[i]->mcef), 2, img->height_cr/2, img->width_cr);
  }

//  bufsize=img->buf_cycle+1;

  fld->short_size=img->buf_cycle;
  fld->short_used=0;
  fld->long_size=0;
  fld->long_used=0;

  if (img->structure != FRAME)
  {
    img->height /= 2;
    img->height_cr /= 2;
  }
    img->buf_cycle = inp->buf_cycle+1;
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
void free_frame_buffers(struct inp_par *inp, ImageParameters *img)
{
  int i;

  img->buf_cycle = inp->buf_cycle+1;

  for (i=0;i<img->buf_cycle;i++)
  {
    free_mem2D(frm->picbuf_short[i]->mref);
    free_mem3D(frm->picbuf_short[i]->mcef,2);
//    free(frm->picbuf_short[i]->Refbuf11);
  }


  for (i=0;i<img->buf_cycle;i++)
    free (frm->picbuf_short[i]);

  free (frm->picbuf_short);		//Access violation here may mean input file is empty

  if (frm->picbuf_long)
  {
    for (i=0;i<frm->long_size;i++)
    {
      free_mem2D(frm->picbuf_long[i]->mref);
      free_mem3D(frm->picbuf_long[i]->mcef,2);
//      free(frm->picbuf_long[i]->Refbuf11);
    } 
    for (i=0;i<frm->long_size;i++)
      free (frm->picbuf_long[i]);

  }

  free (frm);

  for (i=0;i<2 * img->buf_cycle;i++)
  {
    free_mem2D(fld->picbuf_short[i]->mref);
    free_mem3D(fld->picbuf_short[i]->mcef,2);
//    free(fld->picbuf_short[i]->Refbuf11);
  }


  for (i=0;i<2*img->buf_cycle;i++)
    free (fld->picbuf_short[i]);

  free (fld->picbuf_short);

  if (fld->picbuf_long)
  {
    for (i=0;i<fld->long_size;i++)
    {
      free_mem2D(fld->picbuf_long[i]->mref);
      free_mem3D(fld->picbuf_long[i]->mcef,2);
//      free(fld->picbuf_long[i]->Refbuf11);
    } 
    for (i=0;i<fld->long_size;i++)
      free (fld->picbuf_long[i]);
  }

  free (fld);  

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

  printf ("init long term buffer size = %d\n",size);

  /* nothing to do if size unchanged */
  if (fb->long_size==size) return;

  if (fb->long_size>size)
  {
    /* shrink size */

    pb=fb->picbuf_long;
    if ((fb->picbuf_long=(Frame**)calloc(size,sizeof (Frame*)))==NULL) no_mem_exit("init_frame_buffers: fb->picbuf_long");

    for (i=0;i<size;i++)
      fb->picbuf_long[i]=pb[i];

    for (i=size;i<fb->long_size;i++)
    {
      free_mem2D(pb[i]->mref);
      free_mem3D(pb[i]->mcef, 2);
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
      get_mem2D(&(fb->picbuf_long[i]->mref), img->height, img->width);
      get_mem3D(&(fb->picbuf_long[i]->mcef), 2, img->height_cr, img->width_cr);
    }

    free (pb);

  }

  // finally set the size
  fb->long_size=size;

  // and correct the size of mref/mcef
  alloc_mref(img);
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
  int j,id, idx1;
  Frame *f;

  printf ("assign long term long term id %d to short term id %d \n",longID, shortID);

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

  /* delete frames with same long term ID */
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
  
  memcpy (fb->picbuf_long[0]->mref[0], fb->picbuf_short[idx]->mref[0], img->width*img->height); 
  
  for (uv=0; uv < 2; uv++)
    memcpy (fb->picbuf_long[0]->mcef[uv][0], fb->picbuf_short[idx]->mcef[uv][0],img->width_cr*img->height_cr); 

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
      printf("removing long term ID %d\n",longID);

      fb->picbuf_long[i]->used=0;
      fb->picbuf_long[i]->lt_picID=-1;
      fb->picbuf_long[i]->picID=-1;
      
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
/*      if (i<fb->short_used) 
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

  if (mref_fld) 
    free (mref_fld);

  if((mref_fld = (byte***)calloc(fld->short_size+fld->long_size,sizeof(byte**))) == NULL)
    no_mem_exit("alloc_mref: mref_fld");

  if (mcef_fld) 
    free (mcef_fld);

  if((mcef_fld = (byte****)calloc(fld->short_size+fld->long_size,sizeof(byte***))) == NULL)
    no_mem_exit("alloc_mref: mcef_fld");

}

/*!
 ************************************************************************
 * \brief
 *    init mref pointer structure
 ************************************************************************
 */
void init_mref(ImageParameters *img)
{
  int i,j;

  for (i=0,j=0;i<frm->short_used;i++,j++)
  {
    mref_frm[j]=frm->picbuf_short[i]->mref;
    mcef_frm[j]=frm->picbuf_short[i]->mcef;
  }
  for (i=0;i<frm->long_size;i++,j++)
  {
    mref_frm[j]=frm->picbuf_long[i]->mref;
    mcef_frm[j]=frm->picbuf_long[i]->mcef;
  }

  // Tian Dong. June 15, 2002. Let the others be NULL
  for (i=frm->short_used;i<frm->short_size;i++,j++)
  {
    mref_frm[j]=NULL;
    mcef_frm[j]=NULL;
  }
  for (i=frm->long_used;i<frm->long_size;i++,j++)
  {
    mref_frm[j]=NULL;
    mcef_frm[j]=NULL;
  }

  for (i=0,j=0;i<fld->short_used;i++,j++)
  {
    mref_fld[j]=fld->picbuf_short[i]->mref;
    mcef_fld[j]=fld->picbuf_short[i]->mcef;
  }
  for (i=0;i<fld->long_size;i++,j++)
  {
    mref_fld[j]=fld->picbuf_long[i]->mref;
    mcef_fld[j]=fld->picbuf_long[i]->mcef;
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

  //!KS: do nothing but freeing the buffers, needs to be rewritten for list0 and list1
  init_mref(img);
  free_ref_pic_list_reordering_buffer(img->currentSlice);
  return;

#if (0)
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
  if (img->currentSlice->rmpni_buffer==NULL) 
  {
    init_mref(img);
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
//      return;   // Tian: I don't think it is a return point, so I add the judgment (*) below.
      assert( r->Next == NULL );  // Tian assert ...
      break;

    }

    // now scan the frame buffer for the needed frame

    if ( r->RMPNI == 0 || r->RMPNI == 1 || r->RMPNI == 2 )  // Tian: judgment (*)
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
      if ( fr[index]->valid == 0 )  // Tian Dong. June 10, 2002
        printf("Warning: an invalid reference frame is remapped.\n");

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

    img->currentSlice->rmpni_buffer=r->Next;
    free (r);
    r=img->currentSlice->rmpni_buffer;
  }

  // at last init mref and mcef from temporary structure
  for (i=0;i<size;i++)
  {
    mref[i]=fr[i]->mref;
    mcef[i]=fr[i]->mcef;
  }

  // free temporary memory
  free (fr);
#endif
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
 *    calculate the short term picture ID from DPN
 ************************************************************************
 */
int get_short_ID(int DPN, ImageParameters *img)
{
  if ((img->pn-DPN)<0)
    return (img->pn-DPN+img->buf_cycle);
  else
    return (img->pn-DPN);
}

/*!
 ************************************************************************
 * \brief
 *    store reconstructed frame in multiple reference buffer
 ************************************************************************
 */
void copy2fb(ImageParameters *img)
{
  int i,j,uv;
  Frame *f;
  MMCObuffer_t *tmp_mmco;

  /* delete frames with same short term ID */
  if(imgY == imgY_frm||imgY == imgY_top) remove_short_term(img->pn);

  /* cycle frame buffer */
  f=fb->picbuf_short[fb->short_size-1];

  for (i=fb->short_size-2;i>=0;i--)
    fb->picbuf_short[i+1]=fb->picbuf_short[i];

  fb->picbuf_short[0]=f;


  fb->picbuf_short[0]->used=1;
  fb->picbuf_short[0]->picID=img->pn;
  fb->picbuf_short[0]->lt_picID=-1;
  fb->picbuf_short[0]->valid=1;     // Tian Dong: This is a normal reference frame.

  (fb->short_used)++;
  if (fb->short_used>fb->short_size)
    fb->short_used=fb->short_size;
  
  for (j=0; j < img->height; j++)
    for (i=0; i < img->width; i++)
      fb->picbuf_short[0]->mref[j][i]=imgY[j][i]; 
  
  for (uv=0; uv < 2; uv++)
    for (i=0; i < img->width_cr; i++)
      for (j=0; j < img->height_cr; j++)
        fb->picbuf_short[0]->mcef[uv][j][i]=imgUV[uv][j][i];// just copy 1/1 pix, interpolate "online"  

//  fb->picbuf_short[0]->layerNumber = currPictureInfo.layerNumber;
//  fb->picbuf_short[0]->subSequenceIdentifier = currPictureInfo.subSequenceIdentifier;

  // MMCO will be done after storing the decoded frame
  while (img->mmco_buffer)
  {
    tmp_mmco=img->mmco_buffer;
    switch (tmp_mmco->MMCO)
    {
    case 1: // Mark Short-Term Picture as "Unused"
      remove_short_term(get_short_ID(tmp_mmco->DPN,img));
      break;
    case 2: // Mark Long-Term Picture as "Unused"
      remove_long_term(tmp_mmco->LPIN);
      break;
    case 3: // Assign a Long-Term Index to a Picture
      assign_long_term_id(get_short_ID(tmp_mmco->DPN,img),tmp_mmco->LPIN,img);
      break;
    case 4: // Specify the Maximum Long-Term Picture Index
      init_long_term_buffer(tmp_mmco->MLIP1-1,img);
      break;
    case 5: // Reset
      reset_buffers();
      break;
    }

    img->mmco_buffer=tmp_mmco->Next;
    free (tmp_mmco);
  }
}

/*!
 ************************************************************************
 * \brief
 *      Tian Dong, June 13, 2002
 *      If a PN gap is found, try to fill the gap
 * \param pn_expected
 *      which picture will be added to the buffer
 * \param valid
 *      1: a normal reference frame is added to the buffer
 *      0: a blank reference frame is inserted to the buffer for filling
 *          the PN gap
 ************************************************************************
 */
void add_frame(int pn_expected, int valid)
{
  Frame *f;
  int i;

  f=fb->picbuf_short[fb->short_size-1];

  for (i=fb->short_size-2;i>=0;i--)
    fb->picbuf_short[i+1]=fb->picbuf_short[i];

  fb->picbuf_short[0]=f;


  fb->picbuf_short[0]->used=1;
  fb->picbuf_short[0]->picID=pn_expected;
  fb->picbuf_short[0]->lt_picID=-1;
  fb->picbuf_short[0]->valid=valid;

  (fb->short_used)++;
  if (fb->short_used>fb->short_size)
    fb->short_used=fb->short_size;
}

/*!
 ************************************************************************
 * \brief
 *      Tian Dong, June 13, 2002
 *      If a PN gap is found, try to fill the gap
 * \param img
 *      
 ************************************************************************
 */
void fill_PN_gap(ImageParameters *img)
{
  int pn_expected, pn_new;

  if (fb==0 || fb!=frm) return;

  if ( img->type == B_IMG_MULT || img->type ==  B_IMG_1 ) return;

  pn_expected = (fb->picbuf_short[0]->picID+1) % (img->buf_cycle);
  pn_new = img->pn;

  while ( pn_new != pn_expected)
  {
    fb = frm;
    add_frame(pn_expected, 0);

    fb = fld;
    add_frame(pn_expected, 0);
    add_frame(pn_expected, 0);

    pn_expected = (fb->picbuf_short[0]->picID + 1) % (img->buf_cycle);
    img->number++;
  }
  fb = frm;
}


void alloc_ref_pic_list_reordering_buffer(Slice *currSlice)
{
  int size = img->num_ref_pic_active_fwd;

  if (img->type!=INTRA_IMG && img->type!=SI_IMG)
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
  
  size = img->num_ref_pic_active_bwd;

  if (img->type!=B_IMG_1 || img->type!=B_IMG_MULT)
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

