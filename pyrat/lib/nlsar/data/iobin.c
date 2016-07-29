/*
** iobin.c: I/O of SAR images in PolSARPro format
**
** This file is part of NL-SAR Toolbox version 0.6.
**
** Copyright Charles-Alban Deledalle (2013)
** Email charles-alban.deledalle@math.u-bordeaux1.fr
**
** This software is a computer program whose purpose is to provide a
** suite of tools to manipulate SAR images.
**
** This software is governed by the CeCILL license under French law and
** abiding by the rules of distribution of free software. You can use,
** modify and/ or redistribute the software under the terms of the CeCILL
** license as circulated by CEA, CNRS and INRIA at the following URL
** "http://www.cecill.info".
**
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided only
** with a limited warranty and the software's author, the holder of the
** economic rights, and the successive licensors have only limited
** liability.
**
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards their
** requirements in conditions enabling the security of their systems and/or
** data to be ensured and, more generally, to use and operate it in the
** same conditions as regards security.
**
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL license and that you accept its terms.
**
**
** Started on  Wed Jul 24 15:54:05 2013 Charles-Alban Deledalle
** Last update Tue Oct  8 09:37:48 2013 Charles-Alban Deledalle
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <math.h>
#include "tools/sarerror.h"
#include "tools/mathtools.h"
#include "sardata.h"

#define PolSARPro_ID 100

#define SQRT2	1.414213562373095145474621858739

typedef enum {
  PSP_Unknown,
  PSP_Sinclair,
  PSP_Coherency,
  PSP_Covariance,
  PSP_Incoherent
} POLSARPro_format;

typedef struct
{
  char id;
  char polartype[32];
  char polarcase[32];
  char prefix[32];
  POLSARPro_format format;
} PolSARPro_info;


static void swapbyte(char* bytes, int N, unsigned int S)
{
  unsigned long int i;
  int k;
  char tmp;

  for (i = 0; i < S; ++i)
    for (k = 0; k < N/2; ++k)
      {
	tmp = bytes[i * N + k];
	bytes[i * N + k] = bytes[i * N + N - 1 - k];
	bytes[i * N + N - 1 - k] = tmp;
      }
}

sardata* binread_header(const char* dirname, sardata* sar, char swap)
{
    char filename[1024];
    PolSARPro_info info;
    FILE *f;

    sprintf(filename, "%s/config.txt", dirname);
    if (!(f = fopen(filename, "r")))
      {
	sarerror_msg_perror("Cannot open config.txt file");
	return NULL;
      }
    if (fscanf(f, "Nrow\n") != 0 ||
	fscanf(f, "%d\n", &(sar->M)) != 1 ||
	fscanf(f, "---------\n") != 0 ||
	fscanf(f, "Ncol\n") != 0 ||
	fscanf(f, "%d\n", &(sar->N)) != 1 ||
	fscanf(f, "---------\n") != 0 ||
	fscanf(f, "PolarCase\n") != 0 ||
	fscanf(f, "%s\n", info.polarcase) != 1 ||
	fscanf(f, "---------\n") != 0 ||
	fscanf(f, "PolarType\n") != 0 ||
	fscanf(f, "%s\n", info.polartype) != 1)
      {
	sarerror_msg("Bad format in config.txt");
	return NULL;
      }
    fclose(f);

    info.id = PolSARPro_ID;
    info.format = PSP_Unknown;

    // Sinclair?
    if (info.format == PSP_Unknown)
      {
	sprintf(filename, "%s/s11.bin", dirname);
	if ((f = fopen(filename, "r")))
	  {
	    fclose(f);
	    sprintf(info.prefix, "s");
	    info.format = PSP_Sinclair;
	  }
      }
    if (info.format == PSP_Unknown)
      {
	sprintf(filename, "%s/S11.bin", dirname);
	if ((f = fopen(filename, "r")))
	  {
	    fclose(f);
	    sprintf(info.prefix, "S");
	    info.format = PSP_Sinclair;
	  }
      }

    // Coherency?
    if (info.format == PSP_Unknown)
      {
	sprintf(filename, "%s/T11.bin", dirname);
	if ((f = fopen(filename, "r")))
	  {
	    fclose(f);
	    sprintf(info.prefix, "T");
	    info.format = PSP_Coherency;
	  }
      }
    if (info.format == PSP_Unknown)
      {
	sprintf(filename, "%s/T%d/T11.bin", dirname, sar->D);
	if ((f = fopen(filename, "r")))
	  {
	    fclose(f);
	    sprintf(info.prefix, "T%d/T", sar->D);
	    info.format = PSP_Coherency;
	  }
      }

    // Covariance?
    if (info.format == PSP_Unknown)
      {
	sprintf(filename, "%s/C11.bin", dirname);
	if ((f = fopen(filename, "r")))
	  {
	    fclose(f);
	    sprintf(info.prefix, "C");
	    info.format = PSP_Covariance;
	  }
      }
    if (info.format == PSP_Unknown)
      {
	sprintf(filename, "%s/C%d/C11.bin", dirname, sar->D);
	if ((f = fopen(filename, "r")))
	  {
	    fclose(f);
	    sprintf(info.prefix, "C%d/C", sar->D);
	    info.format = PSP_Covariance;
	  }
      }

    // Incoherent?
    if (info.format == PSP_Unknown)
      {
	sprintf(filename, "%s/I11.bin", dirname);
	if ((f = fopen(filename, "r")))
	  {
	    fclose(f);
	    sprintf(info.prefix, "I");
	    info.format = PSP_Incoherent;
	  }
      }

    sar->D = 0;
    if (!strcmp(info.polarcase, "monostatic") && !strcmp(info.polartype, "full"))
      sar->D = 3;
    if (!strcmp(info.polarcase, "bistatic") && !strcmp(info.polartype, "full"))
      sar->D = 4;
    if (!strcmp(info.polarcase, "monostatic") &&
	(!strcmp(info.polartype, "pp1") || !strcmp(info.polartype, "pp2") || !strcmp(info.polartype, "pp3")))
      sar->D = 2;
    if (!strcmp(info.polartype, "pp0"))
      sar->D = 1;
    if (sar->D == 0)
      {
	sarerror_msg("Case %s, %s not implemented", info.polarcase, info.polartype);
	return NULL;
      }
    memcpy(&(sar->extra), &info, sizeof(PolSARPro_info));

    swap = swap;

    return sar;
}

sardata* binread_extract(const char* dirname, sardata* sar, char swap,
			 long int xoffset, long int yoffset,
			 long int width, long int height)
{
  char filename[1024];
  unsigned long int old_size, total_size, target_size;
  complex float tmp;
  FILE* f;
  int i, j, k ,l, m, n;

  old_size = sar->M * sar->N * sar->D * sar->D;
  if (!(sar = binread_header(dirname, sar, swap)))
    return NULL;

  if (xoffset + width > sar->M || yoffset + height > sar->N)
    {
      sarerror_msg("Limits out of bounds");
      return NULL;
    }
  if (width < 0)
    width = sar->M - xoffset;
  if (height < 0)
    height = sar->N - yoffset;
  int filewidth = sar->M;
  int fileheight = sar->N;
  sar->M = width;
  sar->N = height;
  total_size = sar->M * sar->N * sar->D * sar->D;
  total_size = total_size;
  target_size = width * height * sar->D * sar->D;
  if (old_size != target_size)
    if (!(sar->array = realloc(sar->array, target_size * sizeof(float complex))))
      {
	sarerror_msg_perror("Cannot allocate array field");
	return NULL;
      }
  PolSARPro_info* info = (PolSARPro_info *) &(sar->extra);
  switch (info->format)
    {
      case PSP_Sinclair:
	for (i = 0; i < width; ++i)
	  for (j = 0; j < height; ++j)
	    for (k = 0; k < sar->D; ++k)
	      for (l = k; l < sar->D; ++l)
		SARDATA_ACCESS(sar, i, j, k, l) = 1;
	k = 0;
	for (m = 0; m < MIN(2, sar->D); ++m)
	  for (n = (sar->D == 3 ? m : 0); n < MIN(2, sar->D); ++n)
	    {
	      sprintf(filename, "%s/%s%d%d.bin", dirname, info->prefix, m + 1, n + 1);
	      if (!(f = fopen(filename, "r")))
		{
		  sarerror_msg_perror("Cannot open %s file", filename);
		  return NULL;
		}
	      fseek(f, (xoffset * fileheight + yoffset) * sizeof(complex float), SEEK_CUR);
	      for (i = 0; i < width; ++i)
		{
		  for (j = 0; j < height; ++j)
		    {
		      if (fread(&tmp, sizeof(complex float), 1, f) != 1)
			{
			  sarerror_msg("Cannot read data in %s file", filename);
			  return NULL;
			}
		      if (swap)
			swapbyte((char*) &tmp, sizeof(float), 2);
		      if (sar->D == 3 && k == 1)
			tmp = SQRT2 * tmp;
		      for (l = 0; l < sar->D; ++l)
			{
			  SARDATA_ACCESS(sar, i, j, k, l) *= tmp;
			  SARDATA_ACCESS(sar, i, j, l, k) *= conj(tmp);
			}
		    }
		  if (i + xoffset < filewidth - 1)
		    fseek(f, (fileheight - height) * sizeof(complex float), SEEK_CUR);
		}
	      fclose(f);
	      ++k;
	    }
	sprintf(info->prefix, "C");
	info->format = PSP_Covariance;
	break;
      case PSP_Coherency:
      case PSP_Covariance:
	for (i = 0; i < width; ++i)
	  for (j = 0; j < height; ++j)
	    for (k = 0; k < sar->D; ++k)
	      for (l = k; l < sar->D; ++l)
		SARDATA_ACCESS(sar, i, j, k, l) = 0;
	for (k = 0; k < sar->D; ++k)
	  for (l = k; l < sar->D; ++l)
	    {
	      if (l == k)
		sprintf(filename, "%s/%s%d%d.bin", dirname, info->prefix, k + 1, l + 1);
	      else
		sprintf(filename, "%s/%s%d%d_real.bin", dirname, info->prefix, k + 1, l + 1);
	      if (!(f = fopen(filename, "r")))
		{
		  sarerror_msg_perror("Cannot open %s file", filename);
		  return NULL;
		}
	      fseek(f, (xoffset * fileheight + yoffset) * sizeof(float), SEEK_CUR);
	      for (i = 0; i < width; ++i)
		{
		  for (j = 0; j < height; ++j)
		    {
		      if (fread(&tmp, sizeof(float), 1, f) != 1)
			{
			  sarerror_msg("Cannot read data in %s file", filename);
			  return NULL;
			}
		      if (swap)
			swapbyte((char*) &tmp, sizeof(float), 1);
		      SARDATA_ACCESS(sar, i, j, k, l) = tmp;
		    }
		  if (i + xoffset < filewidth - 1)
		    fseek(f, (fileheight - height) * sizeof(float), SEEK_CUR);
		}
	      fclose(f);
	      if (l > k)
		{
		  sprintf(filename, "%s/%s%d%d_imag.bin", dirname, info->prefix, k + 1, l + 1);
		  if (!(f = fopen(filename, "r")))
		    {
		      sarerror_msg_perror("Cannot open %s file", filename);
		      return NULL;
		    }
		  fseek(f, (xoffset * fileheight + yoffset) * sizeof(float), SEEK_CUR);
		  for (i = 0; i < width; ++i)
		    {
		      for (j = 0; j < height; ++j)
			{
			  if (fread(&tmp, sizeof(float), 1, f) != 1)
			    {
			      sarerror_msg("Cannot read data in %s file", filename);
			      return NULL;
			    }
			  if (swap)
			    swapbyte((char*) &tmp, sizeof(float), 1);
			  *(((float*) &SARDATA_ACCESS(sar, i, j, k, l)) + 1) = tmp;
			  SARDATA_ACCESS(sar, i, j, l, k) = conjf(SARDATA_ACCESS(sar, i, j, k, l));
			}
		      if (i + xoffset < filewidth - 1)
			fseek(f, (fileheight - height) * sizeof(float), SEEK_CUR);
		    }
		  fclose(f);
		}
	    }
	break;
      case PSP_Incoherent:
	sarerror_msg("Incoherent matrices not implemented yet");
	return NULL;
	break;
      default:
	sarerror_msg("Unknown format");
	return NULL;
	break;
    }

  return sar;
}

sardata* binread(const char* dirname, sardata* sar, char swap)
{
  return binread_extract(dirname, sar, swap, 0, 0, -1, -1);
}

int binwrite(const sardata* sar, const char* dirname, char swap)
{
  char filename[1024];
  FILE* f;
  int i, j, k ,l;
  float tmp;

  if (mkdir(dirname, S_IRWXU) == -1 && errno != EEXIST)
    {
      sarerror_perror();
      return 0;
    }
  sprintf(filename, "%s/config.txt", dirname);
  if ((f = fopen(filename, "w")) == NULL)
    {
      sarerror_msg_perror("Cannot create config.txt file");
      return 0;
    }
  char* prefix = "C";
  char* polarcase = "unknown";
  char* polartype = "unknown";
  if (sar->extra && *((char*) sar->extra) == PolSARPro_ID)
    {
      PolSARPro_info* info = (PolSARPro_info*) &(sar->extra);
      polarcase = info->polarcase;
      polartype = info->polartype;
      prefix = info->prefix;
    }
  else
    {
      if (sar->D == 3)
	{
	  polarcase = "monostatic";
	  polartype = "full";
	}
      if (sar->D == 4)
	{
	  polarcase = "bistatic";
	  polartype = "full";
	}
      if (sar->D == 2)
	{
	  polarcase = "monostatic";
	  polartype = "pp1";
	}
      if (sar->D == 1)
	{
	  polarcase = "monostatic";
	  polartype = "pp0";
	}
    }
  fprintf(f, "Nrow\n");
  fprintf(f, "%d\n", sar->M);
  fprintf(f, "---------\n");
  fprintf(f, "Ncol\n");
  fprintf(f, "%d\n", sar->N);
  fprintf(f, "---------\n");
  fprintf(f, "PolarCase\n");
  fprintf(f, "%s\n", polarcase);
  fprintf(f, "---------\n");
  fprintf(f, "PolarType\n");
  fprintf(f, "%s\n", polartype);
  fclose(f);

  for (k = 0; k < sar->D; ++k)
    for (l = k; l < sar->D; ++l)
      {
	if (l == k)
	  sprintf(filename, "%s/%s%d%d.bin", dirname, prefix, k + 1, l + 1);
	else
	  sprintf(filename, "%s/%s%d%d_real.bin", dirname, prefix, k + 1, l + 1);
	if (!(f = fopen(filename, "w")))
	  {
	    sarerror_msg_perror("Cannot create %s file", filename);
	    return 0;
	  }
	for (i = 0; i < sar->M; ++i)
	  for (j = 0; j < sar->N; ++j)
	    {
	      tmp = SARDATA_ACCESS(sar, i, j, k, l);
	      if (swap)
		swapbyte((char*) &tmp, sizeof(float), 1);
	      fwrite(&tmp, sizeof(float), 1, f);
	    }
	fclose(f);
	if (l > k)
	  {
	    sprintf(filename, "%s/%s%d%d_imag.bin", dirname, prefix, k + 1, l + 1);
	    if (!(f = fopen(filename, "w")))
	      {
		sarerror_msg_perror("Cannot create %s file", filename);
		return 0;
	      }
	    for (i = 0; i < sar->M; ++i)
	      for (j = 0; j < sar->N; ++j)
		{
		  tmp = *(((float*)&SARDATA_ACCESS(sar, i, j, k, l)) + 1);
		  if (swap)
		    swapbyte((char*) &tmp, sizeof(float), 1);
		  fwrite(&tmp, sizeof(float), 1, f);
		}
	    fclose(f);
	  }
      }
  return 1;
}

