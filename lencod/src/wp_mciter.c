
/*!
*************************************************************************************
* \file wp_mciter.c
*
* \brief
*    Estimate weights for WP using iterative MC method
*
* \author
*    Main contributors (see contributors.h for copyright, address and affiliation details)
*************************************************************************************
*/
#include "contributors.h"

#include <math.h>
#include <float.h>

#include "global.h"
#include "image.h"
#include "wp.h"

/*!
************************************************************************
* \brief
*    Estimates reference picture weighting factors for P slices
************************************************************************
*/

void EstimateWPPSliceAlg2(ImageParameters *img, InputParameters *params, int select_offset)
{
  double dc_org = 0.0;
  double dc_org_UV[2] = {0.0};
  double dc_ref[MAX_REFERENCE_PICTURES] = { 0.0 };
  double dc_ref_UV[MAX_REFERENCE_PICTURES][2] = { {0.0}};

  int i, n, k;
  int default_weight[3];
  int list_offset   = ((img->MbaffFrameFlag)&&(img->mb_data[img->current_mb_nr].mb_field))? (img->current_mb_nr & 0x01) ? 4 : 2 : 0;
  int weight[2][MAX_REFERENCE_PICTURES][3];
  int offset[2][MAX_REFERENCE_PICTURES][3];
  int clist;

  imgpel **tmpPtr;

  luma_log_weight_denom   = 5;
  chroma_log_weight_denom = 5;

  wp_luma_round           = 1 << (luma_log_weight_denom - 1);
  wp_chroma_round         = 1 << (chroma_log_weight_denom - 1);
  default_weight[0]       = 1 << luma_log_weight_denom;
  default_weight[1]       = default_weight[2] = 1 << chroma_log_weight_denom;
  
  dc_org = ComputeImgSum(pCurImg, img->height, img->width);

  if (params->ChromaWeightSupport == 1)
  {
    for (k = 0; k < 2; k++)
    {
      dc_org_UV[k] = ComputeImgSum(pImgOrg[k + 1], img->height_cr, img->width_cr);
    } 
  }

  for (clist = 0; clist < 2 + list_offset; clist++)
  {
    for (n = 0; n < listXsize[clist]; n++)
    {
      if ( wpxDetermineWP( params, img, clist, n ) )
      {
        /* set all values to defaults */
        for (i = 0; i < 3; i++)
        {
          weight[clist][n][i]    = default_weight[i];
          wp_weight[clist][n][i] = default_weight[i];
          wp_offset[clist][n][i] = 0;
          offset[clist][n][i]    = 0;
        }

        // Y
        tmpPtr = listX[clist][n]->p_curr_img;      
        dc_ref[n] = ComputeImgSum(tmpPtr, img->height, img->width);

        if (params->ChromaWeightSupport == 1)
        {
          for (k = 0; k < 2; k++)
          {
            // UV
            tmpPtr = listX[clist][n]->imgUV[k];
            dc_ref_UV[n][k] = ComputeImgSum(tmpPtr, img->height_cr, img->width_cr);
          }        
        }

        if (select_offset == 0)
        {
          if (dc_ref[n] != 0.0)
            weight[clist][n][0] = (int) (default_weight[0] * dc_org / dc_ref[n] + 0.5);
          else
            weight[clist][n][0] = default_weight[0];  // only used when reference picture is black
          weight[clist][n][0] = iClip3(-128, 127, weight[clist][n][0]);
          if (params->ChromaWeightSupport == 1)
          {
            if (dc_ref_UV[n][0] != 0)
              weight[clist][n][1] = (int) (default_weight[1] * dc_org_UV[0] / dc_ref_UV[n][0] + 0.5);
            else
              weight[clist][n][1] = default_weight[1];  // only used when reference picture is black
            weight[clist][n][1] = iClip3(-128, 127, weight[clist][n][1]);

            if (dc_ref_UV[n][1] != 0)
              weight[clist][n][2] = (int) (default_weight[2] * dc_org_UV[1] / dc_ref_UV[n][1] + 0.5);
            else
              weight[clist][n][2] = default_weight[2];  // only used when reference picture is black
            weight[clist][n][2] = iClip3(-64, 128, weight[clist][n][2]);
          }
        }
        else
        {
          offset[clist][n][0] = img->frameOffset[clist][n];

          offset[clist][n][0] = (offset[clist][n][0]+((img->bitdepth_luma-8)>>1))>>(img->bitdepth_luma-8);
          offset[clist][n][0] = iClip3( -128, 127, offset[clist][n][0]);
          offset[clist][n][0] = offset[clist][n][0]<<(img->bitdepth_luma-8);
          weight[clist][n][0] = default_weight[0];

          if (params->ChromaWeightSupport == 1)
          {
            offset[clist][n][1] = (int) ((dc_org_UV[0] - dc_ref_UV[n][0])/(img->size_cr)+0.5);
            offset[clist][n][1] = (offset[clist][n][1] + ((img->bitdepth_chroma - 8)>>1))>>(img->bitdepth_chroma-8);
            offset[clist][n][1] = iClip3( -128, 127, offset[clist][n][1]);
            offset[clist][n][1] = offset[clist][n][1]<<(img->bitdepth_chroma - 8);
            
            weight[clist][n][1] = default_weight[1];

            offset[clist][n][2] = (int) ((dc_org_UV[1] - dc_ref_UV[n][1])/(img->size_cr)+0.5);
            offset[clist][n][2] = (offset[clist][n][2] + ((img->bitdepth_chroma - 8)>>1))>>(img->bitdepth_chroma-8);
            offset[clist][n][2] = iClip3( -128, 127, offset[clist][n][2]);
            offset[clist][n][2] = offset[clist][n][2]<<(img->bitdepth_chroma - 8);

            weight[clist][n][2] = default_weight[2];
          }
        }

        for (i=0; i < 3; i ++)
        {
          wp_weight[clist][n][i] = weight[clist][n][i];
          wp_offset[clist][n][i] = offset[clist][n][i];
#if DEBUG_WP
          printf("index %d component %d weight %d offset %d\n",n,i,weight[0][n][i],offset[0][n][i]);
#endif
        }
      }
    }
  }
}

/*!
************************************************************************
* \brief
*    Estimates reference picture weighting factors for B slices
************************************************************************
*/
void EstimateWPBSliceAlg2(ImageParameters *img, InputParameters *params)
{
  int i, j, k, n;

  int tx,DistScaleFactor;

  int index;
  int comp;
  double dc_org = 0.0;
  double dc_org_UV[2] = { 0.0 };
  double dc_ref_UV[6][MAX_REFERENCE_PICTURES][2] = { {{0.0}} };

  int default_weight[3];
  int list_offset   = ((img->MbaffFrameFlag)&&(img->mb_data[img->current_mb_nr].mb_field))? (img->current_mb_nr & 0x01) ? 4 : 2 : 0;
  int weight[6][MAX_REFERENCE_PICTURES][3];
  int offset[6][MAX_REFERENCE_PICTURES][3];
  int im_weight[6][MAX_REFERENCE_PICTURES][MAX_REFERENCE_PICTURES][3];
  int clist;
  int wf_weight, wf_offset;
  imgpel **tmpPtr;

  if (active_pps->weighted_bipred_idc == 2) //! implicit mode. Values are fixed and it is important to show it here
  {
    luma_log_weight_denom = 5;
    chroma_log_weight_denom = 5;
  }
  else                                     //! explicit mode. Values can be changed for higher precision.
  {
    luma_log_weight_denom = 5;
    chroma_log_weight_denom = 5;
  }

  wp_luma_round     = 1 << (luma_log_weight_denom - 1);
  wp_chroma_round   = 1 << (chroma_log_weight_denom - 1);
  default_weight[0] = 1 << luma_log_weight_denom;
  default_weight[1] = 1 << chroma_log_weight_denom;
  default_weight[2] = 1 << chroma_log_weight_denom;

  if (active_pps->weighted_bipred_idc == 2) //! implicit mode
  {
    for (i = 0; i < listXsize[LIST_0]; i++)
    {
      for (j = 0; j < listXsize[LIST_1]; j++)
      {
        int td, tb;
        td = iClip3(-128, 127,(listX[LIST_1][j]->poc - listX[LIST_0][i]->poc));
        tb = iClip3(-128, 127,(enc_picture->poc - listX[LIST_0][i]->poc));
        for (comp = 0; comp < 3; comp++)
        {
          // implicit weights
          if (td == 0)
          {
            im_weight[0][i][j][comp] = default_weight[comp];
            im_weight[1][i][j][comp] = default_weight[comp];
          }
          else
          {
            tx = (16384 + iabs(td/2))/td;
            DistScaleFactor = iClip3(-1024, 1023, (tx*tb + 32 )>>6);
            im_weight[1][i][j][comp] = DistScaleFactor>>2;
            if (im_weight[1][i][j][comp] < -64 || im_weight[1][i][j][comp] >128)
              im_weight[1][i][j][comp] = default_weight[comp];
            im_weight[0][i][j][comp] = 64 - im_weight[1][i][j][comp];
          }
        }
#if DEBUG_WP
        printf ("%d imp weight[%d][%d] = %d  , %d (%d %d %d) (%d %d) (%d %d)\n",enc_picture->poc, i, j,  im_weight[0][i][j][0], im_weight[1][i][j][0],
          enc_picture->poc,listX[LIST_0][i]->poc, listX[LIST_1][j]->poc,
          DistScaleFactor ,tx,td,tb);
#endif
      }
    }

    for (k = 0; k < 2; k++)
    {
      for (i = 0; i < listXsize[LIST_0]; i++)
      {
        for (j = 0; j < listXsize[LIST_1]; j++)
        {
          for (comp = 0; comp < 3; comp++)
          {
            wbp_weight[k][i][j][comp] = im_weight[k][i][j][comp];
          }
        }
      }
    }

    for (clist=0; clist<2 + list_offset; clist++)
    {
      for (index = 0; index < listXsize[clist]; index++)
      {
        for (comp = 0; comp < 3; comp++)
        {
          wp_weight[clist][index][comp] = default_weight[comp];
          wp_offset[clist][index][comp] = 0;
        }
      }
    }
  }
  else
  {
    dc_org = ComputeImgSum(pCurImg, img->height, img->width);

    if (params->ChromaWeightSupport == 1)
    {
      for (k = 0; k < 2; k++)
      {
        dc_org_UV[k] = ComputeImgSum(pImgOrg[k + 1], img->height_cr, img->width_cr);
      } 
    }

    for (clist=0; clist<2 + list_offset; clist++)
    {
      for (n = 0; n < listXsize[clist]; n++)
      {
        if ( wpxDetermineWP( params, img, clist, n ) )
        {
          /* set all values to defaults */
          for (i = 0; i < 3; i++)
          {
            wp_weight[clist][n][i] = default_weight[i];
            wp_offset[clist][n][i] = 0;
            offset   [clist][n][i] = 0;
            weight   [clist][n][i] = default_weight[i];
          }

          offset[clist][n][0] = wf_offset = img->frameOffset[clist][n];
          weight[clist][n][0] = wf_weight = default_weight[0];


          // UV
          if (params->ChromaWeightSupport == 1)
          {          
            for (k = 0; k < 2; k++)
            {        	
              tmpPtr = listX[clist][n]->imgUV[k];
              dc_ref_UV[clist][n][k] = ComputeImgSum(tmpPtr, img->height_cr, img->width_cr);

              if (dc_ref_UV[clist][n][k] != 0.0)
                wf_weight = (int) (default_weight[k + 1] * dc_org_UV[k] / dc_ref_UV[clist][n][k] + 0.5);
              else
                wf_weight = default_weight[k + 1];  // only used when reference picture is black
              wf_weight = iClip3(-128, 127, wf_weight);
              wf_offset = 0;

              weight[clist][n][k + 1] = wf_weight;
              offset[clist][n][k + 1] = wf_offset;
            }
          }
          else
          {
            weight[clist][n][1] = default_weight[1];
            weight[clist][n][2] = default_weight[2];        
            offset[clist][n][1] = 0;
            offset[clist][n][2] = 0;
          }

          for (i = 0; i < 3; i++)
          {
            wp_weight[clist][n][i] = weight[clist][n][i];
            wp_offset[clist][n][i] = offset[clist][n][i];
#if DEBUG_WP
            printf("%d %d\n",wp_weight[clist][index][comp],wp_offset[clist][index][comp]);
#endif
          }
        }
      }
    }

    if (active_pps->weighted_bipred_idc != 1)
    {
      for (clist=0; clist<2 + list_offset; clist++)
      {
        for (index = 0; index < listXsize[clist]; index++)
        {
          memcpy(wp_weight[clist][index], default_weight, 3 * sizeof(int));
          memset(wp_offset[clist][index], 0, 3 * sizeof(int));
        }
      }
    }


    for (i = 0; i < listXsize[LIST_0]; i++)
    {
      for (j = 0; j < listXsize[LIST_1]; j++)
      {
        for (comp = 0; comp < 3; comp++)
        {
          wbp_weight[0][i][j][comp] = wp_weight[0][i][comp];
          wbp_weight[1][i][j][comp] = wp_weight[1][j][comp];
        }
#if DEBUG_WP
        printf ("bpw weight[%d][%d] = %d  , %d (%d %d %d) (%d %d) (%d %d)\n", i, j, wbp_weight[0][i][j][0], wbp_weight[1][i][j][0],
          enc_picture->poc,listX[LIST_0][i]->poc, listX[LIST_1][j]->poc,
          DistScaleFactor ,tx,tx,tx);
#endif
      }
    }
  }
}


/*!
************************************************************************
* \brief
*    Tests P slice weighting factors to perform or not WP RD decision
************************************************************************
*/

int TestWPPSliceAlg2(ImageParameters *img, InputParameters *params, int select_offset)
{
  int i, j, k, n;

  int index;
  int comp;
  double dc_org = 0.0;
  double dc_org_UV[2] = {0.0};  
  double dc_ref[MAX_REFERENCE_PICTURES] = { 0.0 };
  double dc_ref_UV[MAX_REFERENCE_PICTURES][2] = { {0.0}};

  int default_weight[3];

  int list_offset   = ((img->MbaffFrameFlag)&&(img->mb_data[img->current_mb_nr].mb_field))? (img->current_mb_nr & 0x01) ? 4 : 2 : 0;
  int weight[2][MAX_REFERENCE_PICTURES][3];
  int offset[2][MAX_REFERENCE_PICTURES][3];
  int clist;
  int perform_wp = 0;
  imgpel **tmpPtr;

  luma_log_weight_denom = 5;
  chroma_log_weight_denom = 5;
  wp_luma_round = 1 << (luma_log_weight_denom - 1);
  wp_chroma_round = 1 << (chroma_log_weight_denom - 1);
  default_weight[0] = 1 << luma_log_weight_denom;
  default_weight[1] = default_weight[2] = 1 << chroma_log_weight_denom;

  /* set all values to defaults */
  for (i = 0; i < 2 + list_offset; i++)
  {
    for (j = 0; j < listXsize[i]; j++)
    {
      for (n = 0; n < 3; n++)
      {
        weight[i][j][n] = default_weight[n];
        wp_weight[i][j][n] = default_weight[n];
        wp_offset[i][j][n] = 0;
        offset[i][j][n] = 0;
      }
    }
  }

  dc_org = ComputeImgSum(pCurImg, img->height, img->width);

  if (params->ChromaWeightSupport == 1)
  {
    for (k = 0; k < 2; k++)
    {
      dc_org_UV[k] = ComputeImgSum(pImgOrg[k + 1], img->height_cr, img->width_cr);
    } 
  }

  for (clist = 0; clist < 2 + list_offset; clist++)
  {
    for (n = 0; n < listXsize[clist]; n++)
    {
      tmpPtr = listX[clist][n]->p_curr_img;
      dc_ref[n] = ComputeImgSum(tmpPtr, img->height, img->width);

      if (params->ChromaWeightSupport == 1)
      {
        for (k = 0; k < 2; k++)
        {
          tmpPtr = listX[clist][n]->imgUV[k];
          dc_ref_UV[n][k] = ComputeImgSum(tmpPtr, img->height_cr, img->width_cr);
        }        
      }

      if (select_offset == 0)
      {
        if (dc_ref[n] != 0.0)
          weight[clist][n][0] = (int) (default_weight[0] * dc_org / dc_ref[n] + 0.5);
        else
          weight[clist][n][0] = default_weight[0];  // only used when reference picture is black
        weight[clist][n][0] = iClip3(-128, 127, weight[clist][n][0]);
        if (params->ChromaWeightSupport == 1)
        {
          if (dc_ref_UV[n][0] != 0)
            weight[clist][n][1] = (int) (default_weight[1] * dc_org_UV[0] / dc_ref_UV[n][0] + 0.5);
          else
            weight[clist][n][1] = default_weight[1];  // only used when reference picture is black
          weight[clist][n][1] = iClip3(-128, 127, weight[clist][n][1]);

          if (dc_ref_UV[n][1] != 0)
            weight[clist][n][2] = (int) (default_weight[2] * dc_org_UV[1] / dc_ref_UV[n][1] + 0.5);
          else
            weight[clist][n][2] = default_weight[2];  // only used when reference picture is black
          weight[clist][n][2] = iClip3(-64, 128, weight[clist][n][2]);
        }
      }
      else
      {
        offset[clist][n][0] = img->frameOffset[clist][n]; 

        offset[clist][n][0] = (offset[clist][n][0]+((img->bitdepth_luma-8)>>1))>>(img->bitdepth_luma-8);
        offset[clist][n][0] = iClip3( -128, 127, offset[clist][n][0]);
        offset[clist][n][0] = offset[clist][n][0]<<(img->bitdepth_luma-8);
        weight[clist][n][0] = default_weight[0];

        if (params->ChromaWeightSupport == 1)
        {
            offset[clist][n][1] = (int) ((dc_org_UV[0] - dc_ref_UV[n][0])/(img->size_cr)+0.5);
            offset[clist][n][1] = (offset[clist][n][1] + ((img->bitdepth_chroma - 8)>>1))>>(img->bitdepth_chroma-8);
            offset[clist][n][1] = iClip3( -128, 127, offset[clist][n][1]);
            offset[clist][n][1] = offset[clist][n][1]<<(img->bitdepth_chroma - 8);
            
            weight[clist][n][1] = default_weight[1];

            offset[clist][n][2] = (int) ((dc_org_UV[1] - dc_ref_UV[n][1])/(img->size_cr)+0.5);
            offset[clist][n][2] = (offset[clist][n][2] + ((img->bitdepth_chroma - 8)>>1))>>(img->bitdepth_chroma-8);
            offset[clist][n][2] = iClip3( -128, 127, offset[clist][n][2]);
            offset[clist][n][2] = offset[clist][n][2]<<(img->bitdepth_chroma - 8);

            weight[clist][n][2] = default_weight[2];
        }
      }
    }
  }

  for (clist=0; clist<2 + list_offset; clist++)
  {
    for (index = 0; index < listXsize[clist]; index++)
    {
      for (comp=0; comp < 3; comp ++)
      {
        int offset_test = params->RDPSliceBTest && active_sps->profile_idc != 66
          ? iabs(offset[clist][index][comp]) > 2
          : offset[clist][index][comp] != 0;

        if (weight[clist][index][comp] != default_weight[0] || offset_test)
        {
          perform_wp = 1;
          break;
        }
      }
      if (perform_wp == 1) break;
    }
    if (perform_wp == 1) break;
  }

  return perform_wp;
}

/*!
************************************************************************
* \brief
*    TestWPBSliceAlg2:
*    Tests B slice weighting prediction
************************************************************************
*/
int TestWPBSliceAlg2(ImageParameters *img, InputParameters *params, int select_method)
{
  int i, j, k, n;

  int tx,DistScaleFactor;

  int index;
  int comp;
  double dc_org = 0.0;
  double dc_org_UV[2] = { 0.0 };    
  double dc_ref_UV[6][MAX_REFERENCE_PICTURES][2] = { {{0.0}} };

  int default_weight[3];
  // this needs to be fixed.
  int list_offset   = ((img->MbaffFrameFlag)&&(img->mb_data[img->current_mb_nr].mb_field))? (img->current_mb_nr & 0x01) ? 4 : 2 : 0;
  int weight[6][MAX_REFERENCE_PICTURES][3];
  int offset[6][MAX_REFERENCE_PICTURES][3];
  int im_weight[6][MAX_REFERENCE_PICTURES][MAX_REFERENCE_PICTURES][3];
  int clist;
  int wf_weight, wf_offset;
  int perform_wp = 0;
  imgpel **tmpPtr;

  if (select_method == 1) //! implicit mode
  {
    luma_log_weight_denom = 5;
    chroma_log_weight_denom = 5;
  }
  else
  {
    luma_log_weight_denom = 5;
    chroma_log_weight_denom = 5;
  }

  wp_luma_round     = 1 << (luma_log_weight_denom - 1);
  wp_chroma_round   = 1 << (chroma_log_weight_denom - 1);
  default_weight[0] = 1 << luma_log_weight_denom;
  default_weight[1] = 1 << chroma_log_weight_denom;
  default_weight[2] = 1 << chroma_log_weight_denom;

  /* set all values to defaults */
  for (i = 0; i < 2 + list_offset; i++)
  {
    for (j = 0; j < listXsize[i]; j++)
    {
      for (n = 0; n < 3; n++)
      {
        wp_weight[i][j][n] = default_weight[n];
        wp_offset[i][j][n] = 0;
        offset   [i][j][n] = 0;
        weight   [i][j][n] = default_weight[n];
      }
    }
  }

  for (i = 0; i < listXsize[LIST_0]; i++)
  {
    for (j = 0; j < listXsize[LIST_1]; j++)
    {
      int td, tb;
      td = iClip3(-128, 127,(listX[LIST_1][j]->poc - listX[LIST_0][i]->poc));
      tb = iClip3(-128, 127,(enc_picture->poc - listX[LIST_0][i]->poc));
      for (comp = 0; comp < 3; comp++)
      {
        // implicit weights
        if (td == 0)
        {
          im_weight[0][i][j][comp] = default_weight[comp];
          im_weight[1][i][j][comp] = default_weight[comp];
        }
        else
        {
          tx = (16384 + iabs(td/2))/td;
          DistScaleFactor = iClip3(-1024, 1023, (tx*tb + 32 )>>6);
          im_weight[1][i][j][comp] = DistScaleFactor >> 2;
          if (im_weight[1][i][j][comp] < -64 || im_weight[1][i][j][comp] >128)
            im_weight[1][i][j][comp] = default_weight[comp];
          im_weight[0][i][j][comp] = 64 - im_weight[1][i][j][comp];
        }
      }
    }
  }


  if (select_method == 1) //! implicit mode
  {
    for (i = 0; i < listXsize[LIST_0]; i++)
    {
      for (j = 0; j < listXsize[LIST_1]; j++)
      {
        for (comp = 0; comp < 3; comp++)
        {
          wbp_weight[1][i][j][comp] = im_weight[1][i][j][comp] ;
          wbp_weight[0][i][j][comp] = im_weight[0][i][j][comp];
        }
      }
    }

    for (clist=0; clist<2 + list_offset; clist++)
    {
      for (index = 0; index < listXsize[clist]; index++)
      {
        for (comp = 0; comp < 3; comp++)
        {
          wp_weight[clist][index][comp] = default_weight[comp];
          wp_offset[clist][index][comp] = 0;
        }
      }
    }
  }
  else
  {
    dc_org = ComputeImgSum(pCurImg, img->height, img->width);

    if (params->ChromaWeightSupport == 1)
    {
      for (k = 0; k < 2; k++)
      {
        dc_org_UV[k] = ComputeImgSum(pImgOrg[k + 1], img->height_cr, img->width_cr);
      } 
    }

    for (clist=0; clist<2 + list_offset; clist++)
    {
      for (n = 0; n < listXsize[clist]; n++)
      {
        offset[clist][n][0] = wf_offset = img->frameOffset[clist][n];          
        weight[clist][n][0] = wf_weight = default_weight[0];         

        // UV
        if (params->ChromaWeightSupport == 1)
        {          
          for (k = 0; k < 2; k++)
          {
            tmpPtr = listX[clist][n]->imgUV[k];
            dc_ref_UV[clist][n][k] = ComputeImgSum(tmpPtr, img->height_cr, img->width_cr);

            if (dc_ref_UV[clist][n][k] != 0.0)
              wf_weight = (int) (default_weight[k + 1] * dc_org_UV[k] / dc_ref_UV[clist][n][k] + 0.5);
            else
              wf_weight = default_weight[k + 1];  // only used when reference picture is black
            wf_weight = iClip3(-128, 127, wf_weight);
            wf_offset = 0;

            weight[clist][n][k + 1] = wf_weight;
            offset[clist][n][k + 1] = wf_offset;
          }
        }
        else
        {
          weight[clist][n][1] = default_weight[1];
          weight[clist][n][2] = default_weight[2];        
          offset[clist][n][1] = 0;
          offset[clist][n][2] = 0;
        }
      }
    }

    if (select_method == 0) //! explicit mode
    {
      for (clist=0; clist<2 + list_offset; clist++)
      {
        for (index = 0; index < listXsize[clist]; index++)
        {
          memcpy(wp_weight[clist][index], weight[clist][index], 3 * sizeof(int));
          memcpy(wp_offset[clist][index], offset[clist][index], 3 * sizeof(int));
        }
      }
    }
    else
    {
      for (clist=0; clist<2 + list_offset; clist++)
      {
        for (index = 0; index < listXsize[clist]; index++)
        {
          memcpy(wp_weight[clist][index], default_weight, 3 * sizeof(int));
          memset(wp_offset[clist][index], 0, 3 * sizeof(int));
        }
      }
    }

    for (i = 0; i < listXsize[LIST_0]; i++)
    {
      for (j = 0; j < listXsize[LIST_1]; j++)
      {
        for (comp = 0; comp < 3; comp++)
        {
          wbp_weight[0][i][j][comp] = wp_weight[0][i][comp];
          wbp_weight[1][i][j][comp] = wp_weight[1][j][comp];
        }
#if DEBUG_WP
        printf ("bpw weight[%d][%d] = %d  , %d (%d %d %d) (%d %d) (%d %d)\n", i, j, wbp_weight[0][i][j][0], wbp_weight[1][i][j][0],
        enc_picture->poc,listX[LIST_0][i]->poc, listX[LIST_1][j]->poc,
        DistScaleFactor ,tx,tx,tx);
#endif
      }
    }
  }

  if (select_method == 0) //! implicit mode
  {
    int active_refs[2];

    active_refs[0] = (params->B_List0_refs == 0 ? listXsize[0] : imin(params->B_List0_refs, listXsize[0]));
    active_refs[1] = (params->B_List1_refs == 0 ? listXsize[1] : imin(params->B_List1_refs, listXsize[1]));

    for (clist=0; clist<2 + list_offset; clist++)
    {
      for (index = 0; index < active_refs[clist]; index++)
      {
        for (comp=0; comp < 3; comp ++)
        {
          if (wp_weight[clist][index][comp] != default_weight[comp] || (params->WPIterMC && img->nal_reference_idc && wp_offset[clist][index][comp] != 0) )
          {
            perform_wp = 1;
            break;
          }
        }
        if (perform_wp == 1) break;
      }
      if (perform_wp == 1) break;
    }
  }
  return perform_wp;
}

void compute_offset()
{
  Macroblock *currMB;
  int i, j, x, y, xj, yi, temp, valOrg;
  int mvx=0,  mvy=0;
  int ref_frame=0;
  int subblock=0;
  int x_orig, y_orig;       
  int x_pos, y_pos;
  int out4Y_width  = (img->width  + 2 * IMG_PAD_SIZE) * 4 - 1;  
  int out4Y_height = (img->height + 2 * IMG_PAD_SIZE) * 4 - 1;

  int frame, list, offset; 
  double dtemp;
  int numlists  = (img->type == B_SLICE) ? 2 : 1;


  for(list = 0; list < 2; list++)
  {
    for(frame = 0; frame < MAX_REFERENCE_PICTURES; frame++)
    {
      img->frameOffsetTotal[list][frame] = 0;
      img->frameOffsetCount[list][frame] = 0;
    }
  }


  for(i=0; i<img->height/16; i++) //y
  {
    for(j=0; j<img->width/16; j++)  //x
    {
      { 
        currMB = &img->mb_data[i*img->width/16+j];
        if(IS_INTRA(currMB)) //intra macroblocks are not used for calculation of the filter coeffs.
          continue;

        x_orig = MB_BLOCK_SIZE*j;
        y_orig = MB_BLOCK_SIZE*i;

        for(subblock = 0; subblock < 16; subblock++)
        {
          //List 0
          x = x_orig+4*(subblock%4);
          y = y_orig+4*(subblock/4);
          mvx = enc_picture->motion.mv[LIST_0][y/4][x/4][0];
          mvy = enc_picture->motion.mv[LIST_0][y/4][x/4][1];
          ref_frame = enc_picture->motion.ref_idx[LIST_0][y/4][x/4];

          if(ref_frame != -1)
          {
            for(yi = 0; yi < 4; yi++)
            {    //y
              for(xj = 0; xj < 4; xj++)
              {  //x
                valOrg   = pCurImg[y+yi][x+xj];

                y_pos = imax(0,imin(out4Y_height,4*(y+yi)+4*IMG_PAD_SIZE+mvy));
                x_pos = imax(0,imin(out4Y_width, 4*(x+xj)+4*IMG_PAD_SIZE+mvx));

                temp=listX[LIST_0][ref_frame]->p_curr_img_sub[(y_pos & 0x03)][(x_pos & 0x03)][y_pos >> 2][x_pos >> 2];
                img->frameOffsetTotal[LIST_0][ref_frame]+=(valOrg-temp);
                img->frameOffsetCount[LIST_0][ref_frame]++;          
              }
            }
          } 


          //List 1
          mvx = enc_picture->motion.mv[LIST_1][y/4][x/4][0];
          mvy = enc_picture->motion.mv[LIST_1][y/4][x/4][1];
          ref_frame = enc_picture->motion.ref_idx[LIST_1][y/4][x/4];
          if(ref_frame != -1)
          {
            for(yi = 0; yi < 4; yi++)
            {    //y
              for(xj = 0; xj < 4; xj++)
              {  //x
                valOrg   = pCurImg[y+yi][x+xj];

                y_pos = imax(0,imin(out4Y_height,4*(y+yi)+4*IMG_PAD_SIZE+mvy));
                x_pos = imax(0,imin(out4Y_width, 4*(x+xj)+4*IMG_PAD_SIZE+mvx));

                temp=listX[LIST_0][ref_frame]->p_curr_img_sub[(y_pos & 0x03)][(x_pos & 0x03)][y_pos >> 2][x_pos >> 2];
                img->frameOffsetTotal[LIST_1][ref_frame]+=(valOrg-temp);
                img->frameOffsetCount[LIST_1][ref_frame]++;          
              }
            }
          }
        }//  for(subblock = 0; subblock < 16; subblock++)
      }
    }

    for(list = 0; list < numlists; list++)
    {
      for(frame = 0; frame < listXsize[list]; frame++)
      {
        dtemp=(double)img->frameOffsetTotal[list][frame];

        if (img->frameOffsetCount[list][frame]>0)
        {
          offset=(int)(fabs(dtemp)/(double)img->frameOffsetCount[list][frame]+0.5);
          if (img->frameOffsetTotal[list][frame]>=0)
          {
            img->frameOffset[list][frame]=offset;
          }
          else
          {
            img->frameOffset[list][frame]=-offset;
          }
          //printf("list %d frame %d offset %d frameOffsetCount %d, frameOffsetTotal %d\n", list, frame, img->frameOffset[list][frame], img->frameOffsetCount[list][frame], img->frameOffsetTotal[list][frame]);
        }
        //else
        //{
        //printf("list %d frame %d, frameOffsetCount=0\n", list, frame);
        //}
      }
    }
  }
}
