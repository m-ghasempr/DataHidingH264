#include "Enc_Entropy.h"

void Alloc_Init_CABAC()
{
	FirstMB = 0;
	Enc_Stream->C = 0;
	Enc_Stream->E = 0;
	Enc_Stream->Elow = 0;
	Enc_Stream->Epbuf = -1;
	Enc_Stream->Ebuffer = 0;
	Enc_Stream->Erange = 510;
	Enc_Stream->Ebits_to_go = 17;
	Enc_Stream->Ecodestrm = Enc_Buffer;
	Enc_Stream->Echunks_outstanding = 0;
	Enc_Stream->Ecodestrm_len = &Enc_Buffer_len;
}
void Alloc_Init_CAVLC()
{
	outputBuffer.byteBuffFrame = 0;
	outputBuffer.bits2Go = 8;
	outputBuffer.BuffFrame[0] = '\0';
	outputBuffer.pBuffFrame = 0;
	writepoint = 0;
}
SliceType Find_Slice_type(FILE *Input_File, unsigned char *fileN)
{
	SliceType slice_type = 5;
	int startcode = 0, IFrame = 0, PFrame = 0;
	int frame_count = 0;
	while (slice_type == 5)
	{
		if (!fread(fileN, 1, 256, Input_File))
			return (PFrame ? P_SLICE : I_SLICE);
		for (int i = 0; i < 256; i++)
		{
			if (startcode == 5)
			{
				if (((fileN[i] & 0x7c) == 0x18) || ((fileN[i] & 0x40) == 0x40))					//Pframe
				{
					frame_count++;
					PFrame = 1;
				}
				else if (((fileN[i] & 0x7c) == 0x1c) || ((fileN[i] & 0x70) == 0x20))			//Bframe
				{
					frame_count++;
					slice_type = B_SLICE;
					break;
				}
				else if (((fileN[i] & 0x7F) == 0x8) || ((fileN[i] & 0x70) == 0x30))			//Iframe
				{
					frame_count++;
					if (!IFrame)
						IFrame = 1;
					else
					{
						slice_type = PFrame ? P_SLICE : I_SLICE;
						break;
					}
				}
				startcode = 0;
			}
			if ((startcode == 4) && (((fileN[i] & 0x1F) == 0x5) || ((fileN[i] & 0x1F) == 0x1)))		//IDR or non-IDR
				startcode++;
			if ((fileN[i] == 0) && (startcode < 3))
				startcode++;
			else if ((fileN[i] != 0) && (fileN[i] != 1) && (startcode != 5))
				startcode = 0;
			else if ((fileN[i] == 1) && ((startcode == 3) || (startcode == 2)))
				startcode = 4;
		}
	}
	GoP_Size = frame_count - 1;
	return slice_type;
}
unsigned char* Init_Output_Buffer(char *File_name, char *Input_string, int *size, unsigned char *large)
{
	FILE *Input_file;
	unsigned char *Output_string;
	if (File_name)
	{
		Input_file = fopen(File_name, "rb");			//Open Metadata file (binary)

		fseek(Input_file, 0L, SEEK_END);				//Find Size of Metadata
		*size = ftell(Input_file);
		total_len = *size;
		rewind(Input_file);
		if (*size > LARGE_FILE)
		{
			*size = LARGE_FILE + 8;
			*large = 1;
			Output_string = (unsigned char *)calloc(1, *size);
			fread(Output_string + 8, 1, LARGE_FILE, Input_file);
		}
		else
		{
			*size += 8;
			Output_string = (unsigned char *)calloc(1, *size);
			fread(Output_string + 8, 1, (*size) - 8, Input_file);
		}

		fclose(Input_file);
	}
	else
	{
		Output_string = (unsigned char *)calloc(1, *size + 8);
		for (int i = 0; i < *size; i++)
		{
			Output_string[i + 4] = Input_string[i];
		}
		*size += 8;
	}
	Output_string[0] = Output_string[3] = '$';
	Output_string[1] = Output_string[2] = '0';
	Output_string[4] = (total_len >> 24) & 0xFF;
	Output_string[5] = (total_len >> 16) & 0xFF;
	Output_string[6] = (total_len >> 8) & 0xFF;
	Output_string[7] = total_len & 0xFF;
	return Output_string;
}
void Handle_Large_File(char *File_name, int *size, char *large, char *Output_string, unsigned int num)
{
	FILE *Input_file;
	Input_file = fopen(File_name, "rb");			//Open Metadata file (binary)

	fseek(Input_file, 0L, SEEK_END);				//Find Size of Metadata
	*size = ftell(Input_file) - (*large * LARGE_FILE);
	if (*size > LARGE_FILE)
	{
		*size = LARGE_FILE;
		fseek(Input_file, (*large * LARGE_FILE), SEEK_SET);
		fread(Output_string, 1, LARGE_FILE, Input_file);
		(*large)++;
	}
	else
	{
		fseek(Input_file, (*large * LARGE_FILE), SEEK_SET);
		fread(Output_string, 1, LARGE_FILE, Input_file);
		*large = 0;
	}
	fclose(Input_file);
}
int ReadExpGlomb(char *Bit_stream, int *Offset, int *Bit_offset_to_go)
{
	int M = 0, i = *Offset;
	int Output = 0;
	while (!ReadFlag(Bit_stream, &i, Bit_offset_to_go))
	{
		M++;
	}
	for (int j = 0; j < M; j++)
	{
		Output |= ReadFlag(Bit_stream, &i, Bit_offset_to_go) << (M - 1 - j);
	}
	*Offset = i;
	return Output + (int)pow(2, M) - 1;
}
boolean ReadFlag(char *Bit_stream, int *Offset, int *Bit_offset_to_go)
{
	boolean Output;
	int i = *Offset;
	if (*Bit_offset_to_go == 0)
	{
		*Bit_offset_to_go = 8;
		i++;
	}
	Output = ((Bit_stream[i] & (1 << ((*Bit_offset_to_go) - 1))) != 0);
	(*Bit_offset_to_go)--;
	*Offset = i;
	return Output;
}
void WriteExpGlomb(char *Bit_stream, int *Offset, int *Bit_offset_to_go, int data)
{
	int M = 0, i = *Offset;
	int D = data + 1;
	while (D > 1)
	{
		M++;
		D = D >> 1;
		WriteFlag(Bit_stream, &i, Bit_offset_to_go, 0);
	}
	WriteFlag(Bit_stream, &i, Bit_offset_to_go, 1);
	D = data + 1 - (int)pow(2, M);
	for (int j = 0; j < M; j++)
	{
		WriteFlag(Bit_stream, &i, Bit_offset_to_go, D & (1 << (M - j - 1)));
	}
	*Offset = i;
}
void WriteFlag(char *Bit_stream, int *Offset, int *Bit_offset_to_go, Boolean data)
{
	int i = *Offset;
	if (*Bit_offset_to_go == 0)
	{
		*Bit_offset_to_go = 8;
		i++;
	}
	if (data)
		Bit_stream[i] |= 1 << ((*Bit_offset_to_go) - 1);

	(*Bit_offset_to_go)--;
	*Offset = i;
}

void Copy2End(unsigned long int S_point, FILE *INP, FILE *OUTP)
{
	size_t size = 0;
	double End = 0, Cur = 0;
	fseek(INP, 0L, SEEK_END);
	End = ftell(INP);

	fseek(INP, S_point, SEEK_SET);
	Cur = S_point;
	while ((size = fread(Enc_Buffer, 1, 4 * 1024 * 1024, INP)) > 0)
	{
		fwrite(Enc_Buffer, 1, size, OUTP);
		Cur += size;
		fprintf(stderr, "Percent : 99.%d\n", (int)((Cur * 99) / End));
	}
}