
/*!
 *************************************************************************************
 * \file q_matrix.c
 *
 * \brief
 *    read q_matrix parameters from input file: q_matrix.cfg
 *
 *************************************************************************************
 */
#include <stdlib.h>
#include <string.h>

#include "global.h"

extern char *GetConfigFileContent (char *Filename, int error_type);

#define MAX_ITEMS_TO_PARSE  1000

int matrix4x4_check[6] = {0, 0, 0, 0, 0, 0};
int matrix8x8_check[2] = {0, 0};

/*!
 ***********************************************************************
 * \brief
 *    Check the parameter name.
 * \param s
 *    parameter name string
 * \param type
 *    4x4 or 8x8 matrix type
 * \return
 *    the index number if the string is a valid parameter name,         \n
 *    -1 for error
 ***********************************************************************
 */
int CheckParameterName (char *s, int *type)
{
  int i = 0;

  *type = 0;
  while ((MatrixType4x4[i] != NULL) && (i<6))
  {
    if (0==strcmp (MatrixType4x4[i], s))
      return i;
    else
      i++;
  }

  i = 0;
  *type = 1;
  while ((MatrixType8x8[i] != NULL) && (i<2))
  {
    if (0==strcmp (MatrixType8x8[i], s))
      return i;
    else
      i++;
  }

  return -1;
};

/*!
 ***********************************************************************
 * \brief
 *    Parse the Q matrix values read from cfg file.
 * \param buf
 *    buffer to be parsed
 * \param bufsize
 *    buffer size of buffer
 ***********************************************************************
 */
void ParseMatrix (char *buf, int bufsize)
{
  char *items[MAX_ITEMS_TO_PARSE];
  int MapIdx;
  int item = 0;
  int InString = 0, InItem = 0;
  char *p = buf;
  char *bufend = &buf[bufsize];
  int IntContent;
  int i, j, range, type, cnt;
  short *ScalingList;

  while (p < bufend)
  {
    switch (*p)
    {
      case 13:
        p++;
        break;
      case '#':                 // Found comment
        *p = '\0';              // Replace '#' with '\0' in case of comment immediately following integer or string
        while (*p != '\n' && p < bufend)  // Skip till EOL or EOF, whichever comes first
          p++;
        InString = 0;
        InItem = 0;
        break;
      case '\n':
        InItem = 0;
        InString = 0;
        *p++='\0';
        break;
      case ' ':
      case '\t':              // Skip whitespace, leave state unchanged
        if (InString)
          p++;
        else
        {                     // Terminate non-strings once whitespace is found
          *p++ = '\0';
          InItem = 0;
        }
        break;

      case '"':               // Begin/End of String
        *p++ = '\0';
        if (!InString)
        {
          items[item++] = p;
          InItem = ~InItem;
        }
        else
          InItem = 0;
        InString = ~InString; // Toggle
        break;

      case ',':
        p++;
        InItem = 0;
        break;

      default:
        if (!InItem)
        {
          items[item++] = p;
          InItem = ~InItem;
        }
        p++;
    }
  }

  item--;

  for (i=0; i<item; i+=cnt)
  {
    cnt=0;
    if (0 > (MapIdx = CheckParameterName (items[i+cnt], &type)))
    {
      snprintf (errortext, ET_SIZE, " Parsing error in config file: Parameter Name '%s' not recognized.", items[i+cnt]);
      error (errortext, 300);
    }
    cnt++;
    if (strcmp ("=", items[i+cnt]))
    {
      snprintf (errortext, ET_SIZE, " Parsing error in config file: '=' expected as the second token in each item.");
      error (errortext, 300);
    }
    cnt++;

    if (!type) //4x4 Matrix
    {
      range = 16;
      ScalingList = ScalingList4x4input[MapIdx];
      matrix4x4_check[MapIdx] = 1; //to indicate matrix found in cfg file
    }
    else //8x8 matrix
    {
      range = 64;
      ScalingList = ScalingList8x8input[MapIdx];
      matrix8x8_check[MapIdx] = 1; //to indicate matrix found in cfg file
    }

    for(j=0; j<range; j++)
    {
      if (1 != sscanf (items[i+cnt+j], "%d", &IntContent))
      {
        snprintf (errortext, ET_SIZE, " Parsing error: Expected numerical value for Parameter of %s, found '%s'.", items[i], items[i+cnt+j]);
        error (errortext, 300);
      }

      ScalingList[j] = (short)IntContent; //save value in matrix
    }
    cnt+=j;
    printf (".");
  }
}

/*!
 ***********************************************************************
 * \brief
 *    Check Q Matrix values. If invalid values found in matrix,
 *    whole matrix will be patch with default value 16.
 ***********************************************************************
 */
void PatchMatrix(void)
{
  short *ScalingList;
  int i, cnt, fail;

  for(i=0; i<6; i++)
  {
    if(input->ScalingListPresentFlag[i])
    {
      ScalingList=ScalingList4x4input[i];
      if(matrix4x4_check[i])
      {
        fail=0;
        for(cnt=0; cnt<16; cnt++)
        {
          if(ScalingList[cnt]<0 || ScalingList[cnt]>255) // ScalingList[0]=0 to indicate use default matrix
          {
            fail=1;
            break;
          }
        }

        if(fail) //value of matrix exceed range
        {
          printf("\n%s value exceed range. (Value must be 1 to 255)\n", MatrixType4x4[i]);
          printf("Padding with default values for this matrix.\n\n");
          if(i>2)
            memcpy(ScalingList, Quant_inter_default, sizeof(short)*16);
          else
            memcpy(ScalingList, Quant_intra_default, sizeof(short)*16);
        }
      }
      else //matrix not found, pad with default value
      {
        printf("\n%s matrix not found... padding with default values.\n\n", MatrixType4x4[i]);
        if(i>2)
          memcpy(ScalingList, Quant_inter_default, sizeof(short)*16);
        else
          memcpy(ScalingList, Quant_intra_default, sizeof(short)*16);
      }
    }

    if((i<2) && input->ScalingListPresentFlag[i+6])
    {
      ScalingList=ScalingList8x8input[i];
      if(matrix8x8_check[i])
      {
        fail=0;
        for(cnt=0; cnt<64; cnt++)
        {
          if(ScalingList[cnt]<0 || ScalingList[cnt]>255) // ScalingList[0]=0 to indicate use default matrix
          {
            fail=1;
            break;
          }
        }

        if(fail) //value of matrix exceed range
        {
          printf("\n%s value exceed range. (Value must be 1 to 255)\n", MatrixType8x8[i]);
          printf("Padding with default values for this matrix.\n\n");
          if(i==7)
            memcpy(ScalingList, Quant8_inter_default, sizeof(short)*64);
          else
            memcpy(ScalingList, Quant8_intra_default, sizeof(short)*64);
        }
      }
      else //matrix not found, pad with default value
      {
        printf("\n%s matrix not found... padding with default values.\n\n", MatrixType8x8[i]);
        if(i==7)
          memcpy(ScalingList, Quant8_inter_default, sizeof(short)*64);
        else
          memcpy(ScalingList, Quant8_intra_default, sizeof(short)*64);
      }
    }
  }
}

/*!
 ***********************************************************************
 * \brief
 *    Initialise Q matrix values.
 ***********************************************************************
 */
void Init_QMatrix (void)
{
  char *content;

  if(input->ScalingMatrixPresentFlag)
  {
    printf ("Parsing Configfile %s ", input->QmatrixFile);
    content = GetConfigFileContent(input->QmatrixFile, 0);
    if(content!='\0')
      ParseMatrix(content, strlen (content));
    else
      printf("\nError...file not found, proceeding with default values for all the matrix.");
    PatchMatrix();
    printf("\n");

    memset(UseDefaultScalingMatrix4x4Flag, 0, sizeof(short)*6);
    UseDefaultScalingMatrix8x8Flag[0]=UseDefaultScalingMatrix8x8Flag[1]=0;
  
    free(content);
  }
}
