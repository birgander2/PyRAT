/*
** iorat.c: I/O of SAR images in RAT formats
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
** Started on  Wed Jul 24 15:53:26 2013 Charles-Alban Deledalle
** Last update Fri Aug 23 16:47:48 2013 Charles-Alban Deledalle
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include "tools/sarerror.h"
#include "sardata.h"

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

static sardata* ratread_header2_(FILE* f, sardata* sar, int* var, int* is_xdr);

static sardata* ratread_header_(FILE* f, sardata* sar, int* var, int* is_xdr)
{
  int dim;
  int* size;
  int type;
  char dummy[19];
  char info[80];

  if (fread(&dummy, sizeof(char), 4, f) == 4)
    {
      if (!strncmp(dummy, "RAT2", 4))
	return ratread_header2_(f, sar, var, is_xdr);
      fseek(f, 0, SEEK_SET);
    }
  else
    {
      sarerror_msg("Invalid read of the first 4 bytes");
      return NULL;
    }
  if (fread(&dim, sizeof(int), 1, f) != 1)
    {
      sarerror_msg("Invalid dim field");
      return NULL;
    }
  *is_xdr = (dim < 0 || dim > 100);
  if (*is_xdr)
    swapbyte((char*) &dim, sizeof(int), 1);
  if (!(size = malloc(dim * sizeof(int))))
    {
      sarerror_msg_perror("Cannot allocate size field");
      return NULL;
    }
  if (fread(size, sizeof(int), dim, f) != (unsigned) dim)
    {
      sarerror_msg("Failed to read RAT header: read of size field");
      return NULL;
    }
  if (*is_xdr)
    swapbyte((char*) size, sizeof(int), dim);

  if (!sar)
    sar = sardata_alloc();
  sar->M = size[dim-1];
  sar->N = size[dim-2];
  if (dim <= 2)
    sar->D = 1;
  else
    sar->D = size[0];
  free(size);

  if (fread(var, sizeof(int), 1, f) != 1)
    {
      sarerror_msg("Invalid var field");
      return NULL;
    }
  if (*is_xdr)
    swapbyte((char*) var, sizeof(int), 1);
  if (fread(&type, sizeof(int), 1, f) != 1)
    {
      sarerror_msg("Invalid type field");
      return NULL;
    }
  if (*is_xdr)
    swapbyte((char*) &type, sizeof(int), 1);
  if (fread(dummy, sizeof(char), 12, f) != 12)
    {
      sarerror_msg("Invalid dummy 12 bytes");
      return NULL;
    }
  if (fread(info, sizeof(char), 80, f) != 80)
    {
      sarerror_msg("Invalid info bytes");
      return NULL;
    }
  if (*is_xdr)
    swapbyte((char*) &info, sizeof(char), 80);
  return sar;
}

static sardata* ratread_header2_(FILE* f, sardata* sar, int* var, int* is_xdr)
{
  int dim;
  int channel;
  int size[8];
  int sub[2];
  int type;
  int dummy[150];
  char info[100];
  int geo[25];
  int stat[25];

  if (fread(&dummy, sizeof(int), 1, f) != 1)
    {
      sarerror_msg("Invalid version field");
      return NULL;
    }
  if (fread(&dim, sizeof(int), 1, f) != 1)
    {
      sarerror_msg("Invalid dim field");
      return NULL;
    }
  *is_xdr = (dim < 0 || dim > 100);
  if (*is_xdr)
    swapbyte((char*) &dim, sizeof(int), 1);
  if (fread(&channel, sizeof(int), 1, f) != 1)
    {
      sarerror_msg("Invalid channel field");
      return NULL;
    }
  if (*is_xdr)
    swapbyte((char*) &channel, sizeof(int), 1);
  if (fread(size, sizeof(int), 8, f) != 8)
    {
      sarerror_msg("Failed to read RAT header: read of size field");
      return NULL;
    }
  if (*is_xdr)
    swapbyte((char*) size, sizeof(int), 8);

  if (dim <= 2 && channel != 1)
    {
      sarerror_msg("Cannot interpret dim=%d and channel=%d", dim, channel);
      return NULL;
    }
  if (dim > 2 && channel != size[0] && channel != size[1])
    {
      sarerror_msg("Cannot interpret dim=%d, size[0]=%d, size[1]=%d and channel=%d",
		   dim, size[0], size[1], channel);
      return NULL;
    }

  if (!sar)
    sar = sardata_alloc();
  sar->M = size[dim-1];
  sar->N = size[dim-2];
  sar->D = channel;
  if (fread(var, sizeof(int), 1, f) != 1)
    {
      sarerror_msg("Invalid var field");
      return NULL;
    }
  if (*is_xdr)
    swapbyte((char*) var, sizeof(int), 1);
  if (fread(sub, sizeof(int), 2, f) != 2)
    {
      sarerror_msg("Invalid sub field");
      return NULL;
    }
  if (fread(&type, sizeof(int), 1, f) != 1)
    {
      sarerror_msg("Invalid type field");
      return NULL;
    }
  if (*is_xdr)
    swapbyte((char*) &type, sizeof(int), 1);
  if (fread(dummy, sizeof(int), 9, f) != 9)
    {
      sarerror_msg("Invalid dummy 9 bytes");
      return NULL;
    }
  if (fread(info, sizeof(char), 100, f) != 100)
    {
      sarerror_msg("Invalid info bytes");
      return NULL;
    }
  if (*is_xdr)
    swapbyte((char*) &info, sizeof(char), 100);
  if (fread(geo, sizeof(int), 25, f) != 25)
    {
      sarerror_msg("Invalid info bytes");
      return NULL;
    }
  if (*is_xdr)
    swapbyte((char*) &geo, sizeof(int), 25);
  if (fread(stat, sizeof(int), 25, f) != 25)
    {
      sarerror_msg("Invalid info bytes");
      return NULL;
    }
  if (*is_xdr)
    swapbyte((char*) &stat, sizeof(int), 25);
  if (fread(dummy, sizeof(int), 150, f) != 150)
    {
      sarerror_msg("Invalid dummy 150*4 bytes");
      return NULL;
    }
  return sar;
}

sardata* ratread_header(const char* filename, sardata* sar)
{
  int var, is_xdr;
  FILE* f;

  if (!(f = fopen(filename, "r")))
    {
      sarerror_perror();
      return NULL;
    }
  sar = ratread_header_(f, sar, &var, &is_xdr);
  fclose(f);
  if (!sar)
    return NULL;
  return sar;
}

sardata* ratread_extract(const char* filename, sardata* sar,
			 long int xoffset, long int yoffset,
			 long int width, long int height)
{
  int var, is_xdr;
  unsigned long int input_size, input_bytes, total_size, total_bytes, target_size, target_bytes, k;
  FILE* f;

  input_size = sar->M * sar->N * sar->D * sar->D;
  input_bytes = input_size * sizeof(float complex);
  input_bytes = input_bytes;
  if (!(f = fopen(filename, "r")))
    {
      sarerror_perror();
      return NULL;
    }
  if (!(sar = ratread_header_(f, sar, &var, &is_xdr)))
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

  total_size = sar->M * sar->N * sar->D * sar->D;
  target_size = width * height * sar->D * sar->D;
  long curpos = ftell(f);
  fseek(f, 0, SEEK_END);
  long endpos = ftell(f);
  fseek(f, curpos, SEEK_SET);
  switch (var)
    {
      case 1:
	total_bytes = total_size * sizeof(char);
	target_bytes = target_size * sizeof(char);
	break;
      case 2:
	total_bytes = total_size * sizeof(short);
	target_bytes = target_size * sizeof(short);
	break;
      case 3:
	total_bytes = total_size * sizeof(int);
	target_bytes = target_size * sizeof(int);
	break;
      case 4:
	total_bytes = total_size * sizeof(float);
	target_bytes = target_size * sizeof(float);
	break;
      case 5:
	total_bytes = total_size * sizeof(double);
	target_bytes = target_size * sizeof(double);
	break;
      case 6:
	total_bytes = total_size * sizeof(float complex);
	target_bytes = target_size * sizeof(float complex);
	break;
      case 9:
	total_bytes = total_size * sizeof(double complex);
	target_bytes = target_size * sizeof(double complex);
	break;
      default:
	sarerror_msg("Case var=%d unknown", var);
        return 0;
    }
  if ((unsigned long int) (endpos - curpos) >= total_bytes + (long) (4 * sizeof(char)))
    {
      fseek(f, 4 * sizeof(char), SEEK_CUR);
      curpos += 4;
    }
  if ((unsigned long int) (endpos - curpos) < total_bytes)
    {
      sarerror_msg("Expected %ld bytes but found %ld bytes",
		   total_bytes,
		   endpos - curpos);
      return NULL;
    }
  if (!(sar->array = realloc(sar->array, target_bytes)))
    {
      sarerror_msg_perror("Cannot allocate array field");
      return NULL;
    }
  int x;
  switch (var)
    {
      case 4:
	if (height == sar->N && width == sar->M)
	  {
	    if (fread(sar->array, sizeof(float), target_size, f) != target_size)
	      {
		sarerror_msg("Invalid array field", filename);
		return NULL;
	      }
	  }
	else
	  {
	    sarerror_msg("Extraction for case var=%d not implemented\n", var);
	    return 0;
	  }
	if (is_xdr)
	  swapbyte((char*) sar->array, sizeof(float), target_size);
	if (!(sar->array = realloc(sar->array, target_size * sizeof(float complex))))
	  {
	    sarerror_msg_perror("Cannot allocate array field");
	    return NULL;
	  }
	for (k = target_size - 1; ; k--)
	  {
	    ((float complex*) sar->array)[k] = ((float*) sar->array)[k];
	    if (k == 0)
	      break;
	  }
	break;
      case 6:
	if (height == sar->N)
	  {
	    fseek(f, (xoffset * sar->N * sar->D * sar->D) * sizeof(float complex), SEEK_CUR);
	    if (fread(sar->array, sizeof(float complex), target_size, f) != target_size)
	      {
		sarerror_msg("Invalid array field", filename);
		return NULL;
	      }
	  }
	else
	  {
	    fseek(f, (xoffset * sar->N + yoffset) * sar->D * sar->D * sizeof(float complex), SEEK_CUR);
	    for (x = 0; x < width; ++x)
	      {
		if (fread(sar->array + x * height * sar->D * sar->D, sizeof(float complex), height * sar->D * sar->D, f)
		    != (unsigned long) height * sar->D * sar->D)
		  {
		    sarerror_msg("Invalid array field", filename);
		    return NULL;
		  }
		if (x + xoffset < sar->M - 1)
		  fseek(f, (sar->N - height) * sar->D * sar->D * sizeof(float complex), SEEK_CUR);
	      }
	  }
	if (is_xdr)
	  swapbyte((char*) sar->array, sizeof(float), 2 * target_size);
	break;
      case 9:
	if (height == sar->N && width == sar->M)
	  {
	    if (fread(sar->array, sizeof(double complex), target_size, f) != target_size)
	      {
		sarerror_msg("Invalid array field", filename);
		return NULL;
	      }
	  }
	else
	  {
	    sarerror_msg("Extraction for case var=%d not implemented\n", var);
	    return 0;
	  }
	if (is_xdr)
	  swapbyte((char*) sar->array, sizeof(double), 2 * target_size);
	for (k = 0; k < target_size; ++k)
	  ((float complex*) sar->array)[k] = ((double complex*) sar->array)[k];
	if (!(sar->array = realloc(sar->array, target_size * sizeof(float complex))))
	  {
	    sarerror_msg_perror("Cannot allocate array field");
	    return NULL;
	  }
	break;
      default:
	sarerror_msg("Case var=%d not implemented\n", var);
        return 0;
    }
  fclose(f);
  sar->M = width;
  sar->N = height;
  return sar;
}

sardata* ratread(const char* filename, sardata* sar)
{
  return ratread_extract(filename, sar, 0, 0, -1, -1);
}

int ratwrite(const sardata* sar, const char* filename)
{
  int dim;
  int size[4];
  unsigned long int total_size = 0;
  int var;
  int type;
  char dummy[12] = { 0, };
  char info[80] = { ' ', };
  FILE* f;
  if (!(f = fopen(filename, "w")))
    {
      sarerror_perror();
      return 0;
    }
  dim = 4;
  fwrite(&dim, sizeof(int), 1, f);
  size[3] = sar->M;
  size[2] = sar->N;
  size[1] = sar->D;
  size[0] = sar->D;
  fwrite(size, sizeof(int), dim, f);
  var = 6;
  fwrite(&var, sizeof(int), 1, f);
  type = 0;
  fwrite(&type, sizeof(int), 1, f);
  fwrite(dummy, sizeof(char), 12, f);
  fwrite(info, sizeof(char), 80, f);
  fwrite(dummy, sizeof(char), 4, f);
  switch (var)
    {
      case 6:
	total_size = sar->M * sar->N * sar->D * sar->D;
	fwrite(sar->array, sizeof(float complex), total_size, f);
	break;
      default:
	sarerror_msg("Case var=%d not implemented\n", var);
        return 0;
    }
  fclose(f);
  return 1;
}
