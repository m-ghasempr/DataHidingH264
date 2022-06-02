
/*!
 ***********************************************************************
 *  \file
 *     decoder_test.c
 *  \brief
 *     H.264/AVC decoder test 
 *  \author
 *     Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Yuwen He       <yhe@dolby.com>
 ***********************************************************************
 */

#include "contributors.h"

#include <sys/stat.h>

#include "win32.h"
#include "h264decoder.h"
#include "configfile.h"
#include "Enc_Entropy.h"
#include "memalloc.h"

#define DECOUTPUT_TEST      0

#define PRINT_OUTPUT_POC    0
#define BITSTREAM_FILENAME  "test.264"
#define DECRECON_FILENAME   "test_dec.yuv"
#define ENCRECON_FILENAME   "test_rec.yuv"
#define FCFR_DEBUG_FILENAME "fcfr_dec_rpu_stats.txt"
#define DECOUTPUT_VIEW0_FILENAME  "H264_Decoder_Output_View0.yuv"
#define DECOUTPUT_VIEW1_FILENAME  "H264_Decoder_Output_View1.yuv"


static void Configure(InputParameters *p_Inp, int ac, char *av[])
{
  //char *config_filename=NULL;
  //char errortext[ET_SIZE];
  memset(p_Inp, 0, sizeof(InputParameters));
  strcpy(p_Inp->infile, BITSTREAM_FILENAME); //! set default bitstream name
  strcpy(p_Inp->outfile, DECRECON_FILENAME); //! set default output file name
  strcpy(p_Inp->reffile, ENCRECON_FILENAME); //! set default reference file name
  
#ifdef _LEAKYBUCKET_
  strcpy(p_Inp->LeakyBucketParamFile,"leakybucketparam.cfg");    // file where Leaky Bucket parameters (computed by encoder) are stored
#endif

  ParseCommand(p_Inp, ac, av);

  fprintf(stdout,"----------------------------- JM %s %s -----------------------------\n", VERSION, EXT_VERSION);
  //fprintf(stdout," Decoder config file                    : %s \n",config_filename);
  if(!p_Inp->bDisplayDecParams)
  {
    fprintf(stdout,"--------------------------------------------------------------------------\n");
    fprintf(stdout," Input H.264 bitstream                  : %s \n",p_Inp->infile);
    fprintf(stdout," Output decoded YUV                     : %s \n",p_Inp->outfile);
    //fprintf(stdout," Output status file                     : %s \n",LOGFILE);
    fprintf(stdout," Input reference file                   : %s \n",p_Inp->reffile);

    fprintf(stdout,"--------------------------------------------------------------------------\n");
  #ifdef _LEAKYBUCKET_
    fprintf(stdout," Rate_decoder        : %8ld \n",p_Inp->R_decoder);
    fprintf(stdout," B_decoder           : %8ld \n",p_Inp->B_decoder);
    fprintf(stdout," F_decoder           : %8ld \n",p_Inp->F_decoder);
    fprintf(stdout," LeakyBucketParamFile: %s \n",p_Inp->LeakyBucketParamFile); // Leaky Bucket Param file
    calc_buffer(p_Inp);
    fprintf(stdout,"--------------------------------------------------------------------------\n");
  #endif
  }
  
}

/*********************************************************
if bOutputAllFrames is 1, then output all valid frames to file onetime; 
else output the first valid frame and move the buffer to the end of list;
*********************************************************/
static int WriteOneFrame(DecodedPicList *pDecPic, int hFileOutput0, int hFileOutput1, int bOutputAllFrames)
{
  int iOutputFrame=0;
  DecodedPicList *pPic = pDecPic;

  if(pPic && (((pPic->iYUVStorageFormat==2) && pPic->bValid==3) || ((pPic->iYUVStorageFormat!=2) && pPic->bValid==1)) )
  {
    int i, iWidth, iHeight, iStride, iWidthUV, iHeightUV, iStrideUV;
    byte *pbBuf;    
    int hFileOutput;
    int res;

    iWidth = pPic->iWidth*((pPic->iBitDepth+7)>>3);
    iHeight = pPic->iHeight;
    iStride = pPic->iYBufStride;
    if(pPic->iYUVFormat != YUV444)
      iWidthUV = pPic->iWidth>>1;
    else
      iWidthUV = pPic->iWidth;
    if(pPic->iYUVFormat == YUV420)
      iHeightUV = pPic->iHeight>>1;
    else
      iHeightUV = pPic->iHeight;
    iWidthUV *= ((pPic->iBitDepth+7)>>3);
    iStrideUV = pPic->iUVBufStride;
    
    do
    {
      if(pPic->iYUVStorageFormat==2)
        hFileOutput = (pPic->iViewId&0xffff)? hFileOutput1 : hFileOutput0;
      else
        hFileOutput = hFileOutput0;
      if(hFileOutput >=0)
      {
        //Y;
        pbBuf = pPic->pY;
        for(i=0; i<iHeight; i++)
        {
          res = write(hFileOutput, pbBuf+i*iStride, iWidth);
          if (-1==res)
          {
            error ("error writing to output file.", 600);
          }
        }

        if(pPic->iYUVFormat != YUV400)
        {
         //U;
         pbBuf = pPic->pU;
         for(i=0; i<iHeightUV; i++)
         {
           res = write(hFileOutput, pbBuf+i*iStrideUV, iWidthUV);
           if (-1==res)
           {
             error ("error writing to output file.", 600);
           }
}
         //V;
         pbBuf = pPic->pV;
         for(i=0; i<iHeightUV; i++)
         {
           res = write(hFileOutput, pbBuf+i*iStrideUV, iWidthUV);
           if (-1==res)
           {
             error ("error writing to output file.", 600);
           }
         }
        }

        iOutputFrame++;
      }

      if (pPic->iYUVStorageFormat == 2)
      {
        hFileOutput = ((pPic->iViewId>>16)&0xffff)? hFileOutput1 : hFileOutput0;
        if(hFileOutput>=0)
        {
          int iPicSize =iHeight*iStride;
          //Y;
          pbBuf = pPic->pY+iPicSize;
          for(i=0; i<iHeight; i++)
          {
            res = write(hFileOutput, pbBuf+i*iStride, iWidth);
            if (-1==res)
            {
              error ("error writing to output file.", 600);
            }
          }

          if(pPic->iYUVFormat != YUV400)
          {
           iPicSize = iHeightUV*iStrideUV;
           //U;
           pbBuf = pPic->pU+iPicSize;
           for(i=0; i<iHeightUV; i++)
           {
             res = write(hFileOutput, pbBuf+i*iStrideUV, iWidthUV);
             if (-1==res)
             {
               error ("error writing to output file.", 600);
             }
           }
           //V;
           pbBuf = pPic->pV+iPicSize;
           for(i=0; i<iHeightUV; i++)
           {
             res = write(hFileOutput, pbBuf+i*iStrideUV, iWidthUV);
             if (-1==res)
             {
               error ("error writing to output file.", 600);
             }
           }
          }

          iOutputFrame++;
        }
      }

#if PRINT_OUTPUT_POC
      fprintf(stdout, "\nOutput frame: %d/%d\n", pPic->iPOC, pPic->iViewId);
#endif
      pPic->bValid = 0;
      pPic = pPic->pNext;
    }while(pPic != NULL && pPic->bValid && bOutputAllFrames);
  }
#if PRINT_OUTPUT_POC
  else
    fprintf(stdout, "\nNone frame output\n");
#endif

  return iOutputFrame;
}

/*!
 ***********************************************************************
 * \brief
 *    main function for JM decoder
 ***********************************************************************
 */
int main(int argc, char **argv)
{
  int iRet;
  DecodedPicList *pDecPicList;
  int hFileDecOutput0=-1, hFileDecOutput1=-1;
  int iFramesOutput=0, iFramesDecoded=0;
  InputParameters InputParams;
  
  unsigned char fileN[256];
  unsigned char sizeH = 0;
  int OffsetEnd = 5, BitOffsetEnd = 8;
  unsigned char threeEnd = 0x3;

#if DECOUTPUT_TEST
  hFileDecOutput0 = open(DECOUTPUT_VIEW0_FILENAME, OPENFLAGS_WRITE, OPEN_PERMISSIONS);
  fprintf(stdout, "Decoder output view0: %s\n", DECOUTPUT_VIEW0_FILENAME);
  hFileDecOutput1 = open(DECOUTPUT_VIEW1_FILENAME, OPENFLAGS_WRITE, OPEN_PERMISSIONS);
  fprintf(stdout, "Decoder output view1: %s\n", DECOUTPUT_VIEW1_FILENAME);
#endif

  init_time();

  //get input parameters;
  Configure(&InputParams, argc, argv);
  Enc_MB = (Enc_Macroblock *)calloc(1, sizeof(Enc_Macroblock));
  Enc_currSlice = (Enc_Slice *)calloc(1, sizeof(Enc_Slice));
  Enc_dp = (Enc_DataPartition *)calloc(1, sizeof(Enc_DataPartition));
  Enc_Vid = (Enc_VideoParameters *)calloc(1, sizeof(Enc_VideoParameters));
  Enc_Stream = (EncodingEnvironment *)calloc(1, sizeof(EncodingEnvironment));
  Enc_VStream = (Enc_Bitstream *)calloc(1, sizeof(Enc_Bitstream));
  Enc_Vid->active_pps = (pic_parameter_set_rbsp_t*)calloc(1, sizeof(pic_parameter_set_rbsp_t));
  Enc_Vid->PicPos = (BlockPos *)calloc(1, sizeof(BlockPos));
  Enc_currSlice->tex_ctx = (TextureInfoContexts *)calloc(1, sizeof(TextureInfoContexts));

  Enc_currSlice->p_Vid = Enc_Vid;
  Enc_currSlice->partArr = Enc_dp;
  Enc_currSlice->partArr->ee_cabac = *Enc_Stream;
  Enc_MB->p_Slice = Enc_currSlice;
  Enc_MB->p_Vid = Enc_Vid;

  strcpy(fileN, InputParams.infile);
  strcat(fileN, ".264");
  Input_File = fopen(InputParams.infile, "rb");
  Output_File = fopen(fileN, "wb");

  InsertingSlice = Find_Slice_type(Input_File, fileN);
  //read Size
  fseek(Input_File, -1, SEEK_END);
  fread(&sizeH, 1, 1, Input_File);

  fseek(Input_File, (sizeH * -1), SEEK_END);
  fread(fileN, 1, 5, Input_File);

  if ((fileN[1] == 0x0) && (fileN[2] == 0x0) && (fileN[3] == 0x1) && (fileN[4] == 28))				// fileN[0] bekhater ffmpeg 3 bayti check nashode
  {
	  Inserted = 1;
	  for (int j = 5; j <= sizeH; j++)
	  {
		  fread(&fileN[j], 1, 1, Input_File);
		  if ((fileN[j - 2] == 0) && (fileN[j - 1] == 0) && (fileN[j] == 3))
		  {
			  fread(&fileN[j], 1, 1, Input_File);
			  sizeH--;
		  }
	  }
	  endInfo.FrameNum = ReadExpGlomb(fileN, &OffsetEnd, &BitOffsetEnd);
	  endInfo.SliceMbNum = ReadExpGlomb(fileN, &OffsetEnd, &BitOffsetEnd);
	  endInfo.Threshold = ReadExpGlomb(fileN, &OffsetEnd, &BitOffsetEnd);
	  endInfo.MetaDataNum = ReadExpGlomb(fileN, &OffsetEnd, &BitOffsetEnd);
	  endInfo.FrameType = ReadExpGlomb(fileN, &OffsetEnd, &BitOffsetEnd);
	  if (THR != 1)
	  {
		  if (endInfo.Threshold != THR)
		  {
			  fprintf(stderr, "THN\n");
			  return -1;
		  }
	  }
	  if (InsFrame != 4)
	  {
		  if (endInfo.FrameType != InsFrame)
		  {
			  fprintf(stderr, "FTN\n");
			  return -1;
		  }
	  }
  }
  else
  {
	  endInfo.FrameNum = 0;
	  endInfo.SliceMbNum = -1;
	  endInfo.Threshold = THR;
	  endInfo.MetaDataNum = 0;
	  endInfo.FrameType = InsFrame;
	  if (InsFrame == 4)
	  {
		  fprintf(stderr, "BFN\n");
		  return -1;
	  }
	  if (THR == 0)
	  {
		  fprintf(stderr, "ZTH\n");
		  return -1;
	  }
	  AllowInsertion = 1;
  }
  THR = endInfo.Threshold;
  if ((endInfo.Threshold < 1) || (endInfo.Threshold > 14))
  {
	  fprintf(stderr, "HMH\n");
	  return -1;
  }


  if (InsFrame == 0)
	  InsertingSlice = B_SLICE;
  /*else if ((InsFrame != 4) && (InsFrame != 1))
	  InsertingSlice = InsFrame;
  else if (InsFrame == 1)
	  InsertingSlice = InsertingSlice;*/
  else if (InsFrame == 1)
	  InsertingSlice = I_SLICE;
  else if (InsFrame == 2)
	  InsertingSlice = P_SLICE;
  else if (InsFrame == 3)
	  InsertingSlice = SP_SLICE; //limited P
  else if (InsFrame == 4)
	  InsertingSlice = (endInfo.FrameType == B_SLICE) ? InsertingSlice : endInfo.FrameType;

  //open decoder;
  iRet = OpenDecoder(&InputParams);
  if(iRet != DEC_OPEN_NOERR)
  {
    fprintf(stderr, "Open encoder failed: 0x%x!\n", iRet);
    return -1; //failed;
  }

  Input_MD = Init_Output_Buffer(G_File_MDIn, NULL, &Input_Len, &large_File);
  Input_MD[1] = endInfo.MetaDataNum >> 8; Input_MD[2] = endInfo.MetaDataNum & 0xFF;
  pInput_MD = 0;
  Start_Slice_Input += 4;
  get_mem3Dint(&Enc_Vid->nz_coeff, 512 * 256, 4, 4);				//frame MB no.



  //decoding;
  do
  {
    iRet = DecodeOneFrame(&pDecPicList);
	if (I_Finish)
	{
		fprintf(stderr, "Frame : %d\n", fna);
		break;
	}
	fna++;
	if ((fna % 10) == 0)
	{
		fprintf(stderr, "Frame : %d\n", fna);
		fprintf(stderr, "Percent : %5.2f\n", (float)((pInput_MD + Pmd) * 99) / (float)(total_len << 3));
	}


    if(iRet==DEC_EOS || iRet==DEC_SUCCEED)
    {
      //process the decoded picture, output or display;
      iFramesOutput += WriteOneFrame(pDecPicList, hFileDecOutput0, hFileDecOutput1, 0);
      iFramesDecoded++;
    }
    else
    {
      //error handling;
      fprintf(stderr, "Error in decoding process: 0x%x\n", iRet);
    }
  }while((iRet == DEC_SUCCEED) && ((p_Dec->p_Inp->iDecFrmNum==0) || (iFramesDecoded<p_Dec->p_Inp->iDecFrmNum)));
  if (Enc_Vid->nz_coeff)
	  free_mem3Dint(Enc_Vid->nz_coeff);
  Copy2End(Input_File_PTR, Input_File, Output_File);
  fprintf(stderr, "Percent : 100.0\n");
  fflush(stderr);
  if (Inserted)
  {
	  fseek(Input_File, (sizeH * -1) + 1, SEEK_END);					//vase ffmpeg
	  fseek(Output_File, (sizeH * -1) + 1, SEEK_END);					//vase ffmpeg
  }
  else
	  fseek(Input_File, 0, SEEK_END);

  if (I_Finish)
  {
	  memset(fileN, 0, 256);
	  OffsetEnd = 5; BitOffsetEnd = 8;
	  fileN[0] = fileN[1] = fileN[2] = 0x0;
	  fileN[3] = 0x1;
	  fileN[4] = 28;
	  WriteExpGlomb(fileN, &OffsetEnd, &BitOffsetEnd, endInfo.FrameNum);
	  WriteExpGlomb(fileN, &OffsetEnd, &BitOffsetEnd, endInfo.SliceMbNum);
	  WriteExpGlomb(fileN, &OffsetEnd, &BitOffsetEnd, endInfo.Threshold);
	  WriteExpGlomb(fileN, &OffsetEnd, &BitOffsetEnd, endInfo.MetaDataNum);
	  WriteExpGlomb(fileN, &OffsetEnd, &BitOffsetEnd, endInfo.FrameType);
	  fprintf(stderr, "Id : %d\n", endInfo.MetaDataNum);
	  sizeH = OffsetEnd + 1;
	  if (Inserted)
		  fwrite(fileN + 1, 1, 4, Output_File);
	  else
		  fwrite(fileN, 1, 5, Output_File);
	  for (int j = 5; j <= OffsetEnd; j++)
	  {
		  fwrite(&fileN[j], 1, 1, Output_File);
		  if ((fileN[j] == 0) && (fileN[j - 1] == 0) && (fileN[j + 1] <= 3))
		  {
			  fwrite(&threeEnd, 1, 1, Output_File);
			  sizeH++;
		  }
	  }
	  sizeH++;
	  fwrite(&sizeH, 1, 1, Output_File);
	  fprintf(stderr, "Last : %ld\n", endInfo.FrameNum);
	  fprintf(stderr, "Completed.");
	  fflush(stderr);
  }
  else
  {
	  fprintf(stderr, "VoL\n");
	  fflush(stderr);
  }

  iRet = FinitDecoder(&pDecPicList);
  iFramesOutput += WriteOneFrame(pDecPicList, hFileDecOutput0, hFileDecOutput1 , 1);
  iRet = CloseDecoder();

  //quit;
  if(hFileDecOutput0>=0)
  {
    close(hFileDecOutput0);
  }
  if(hFileDecOutput1>=0)
  {
    close(hFileDecOutput1);
  }

  printf("%d frames are decoded.\n", iFramesDecoded);
  Enc_MB = NULL;
  Enc_currSlice = NULL;
  Enc_dp = NULL;
  Enc_Vid = NULL;
  free(Enc_MB);
  free(Enc_currSlice);
  free(Enc_dp);
  free(Enc_Vid);

  return 0;
}


