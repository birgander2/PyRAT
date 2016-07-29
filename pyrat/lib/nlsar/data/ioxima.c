/*
** ioxima.c: I/O of SAR images in XIMA formats
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
** Started on  Wed Jul 24 15:53:00 2013 Charles-Alban Deledalle
** Last update Thu Oct 10 18:06:41 2013 Charles-Alban Deledalle
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

sardata* ximaread_header(const char* filename, sardata* sar, char swap)
{
    char filename_dim[1024];
    FILE *f;
    int len = strlen(filename);

    if (len < 3)
      {
	sarerror_msg("Cannot add .dim extension: Filename too short");
	return NULL;
      }
    strcpy(filename_dim, filename);
    filename_dim[len - 3] = 'd';
    filename_dim[len - 2] = 'i';
    filename_dim[len - 1] = 'm';
    if ((f = fopen(filename_dim, "r")) == NULL)
      {
	sarerror_msg_perror("Cannot open %s file", filename_dim);
	return NULL;
      }
    if (fscanf(f, "%d %d\n", &(sar->N), &(sar->M)) != 2)
      {
	sarerror_msg("Bad format in %s", filename_dim);
	return NULL;
      }
    fclose(f);
    sar->D = 1;
    swap = swap;
    return sar;
}

sardata* imaread_extract(const char* filename, sardata* sar, char swap,
			 long int xoffset, long int yoffset,
			 long int width, long int height)
{
  unsigned long int input_size, total_size, target_size;
  FILE* f;
  int k;

  input_size = sar->M * sar->N * sar->D * sar->D;
  if (!(sar = ximaread_header(filename, sar, swap)))
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
  total_size = total_size;
  target_size = width * height * sar->D * sar->D;
  if (input_size != target_size)
    if (!(sar->array = realloc(sar->array, target_size * sizeof(float complex))))
      {
	sarerror_msg_perror("Cannot allocate array field");
	return NULL;
      }
  if (!(f = fopen(filename, "r")))
    {
      sarerror_msg_perror("Cannot open %s file", filename);
      return NULL;
    }
  if (height == sar->N)
    {
      fseek(f, (xoffset * sar->N * sar->D * sar->D) * sizeof(unsigned char), SEEK_CUR);
      if (fread(sar->array, sizeof(unsigned char), target_size, f) != target_size)
	{
	  sarerror_msg("Cannot read data in %s file", filename);
	  return NULL;
	}
    }
  else
    {
      int x;
      fseek(f, (xoffset * sar->N + yoffset) * sar->D * sar->D * sizeof(unsigned char), SEEK_CUR);
      for (x = 0; x < width; ++x)
	{
	  if (fread((unsigned char *) sar->array + x * height * sar->D * sar->D, sizeof(unsigned char), height * sar->D * sar->D, f)
	      != (unsigned long) height * sar->D * sar->D)
	    {
	      sarerror_msg("Invalid array field", filename);
	      return NULL;
	    }
	  if (x + xoffset < sar->M - 1)
	    fseek(f, (sar->N - height) * sar->D * sar->D * sizeof(unsigned char), SEEK_CUR);
	}
    }
  fclose(f);
  sar->M = width;
  sar->N = height;
  if (swap)
    swapbyte((char*) sar->array, sizeof(unsigned char), target_size);
  for (k = target_size - 1; ; k--)
    {
      ((float complex*) sar->array)[k] = ((unsigned char*) sar->array)[k];
      if (k == 0)
	break;
    }
  return sar;
}

sardata* imwread_extract(const char* filename, sardata* sar, char swap,
			 long int xoffset, long int yoffset,
			 long int width, long int height)
{
  unsigned long int input_size, total_size, target_size;
  FILE* f;
  int k;

  input_size = sar->M * sar->N * sar->D * sar->D;
  if (!(sar = ximaread_header(filename, sar, swap)))
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
  total_size = total_size;
  target_size = width * height * sar->D * sar->D;
  if (input_size != target_size)
    if (!(sar->array = realloc(sar->array, target_size * sizeof(float complex))))
      {
	sarerror_msg_perror("Cannot allocate array field");
	return NULL;
      }
  if (!(f = fopen(filename, "r")))
    {
      sarerror_msg_perror("Cannot open %s file", filename);
      return NULL;
    }
  if (height == sar->N)
    {
      fseek(f, (xoffset * sar->N * sar->D * sar->D) * sizeof(unsigned short), SEEK_CUR);
      if (fread(sar->array, sizeof(unsigned short), target_size, f) != target_size)
	{
	  sarerror_msg("Cannot read data in %s file", filename);
	  return NULL;
	}
    }
  else
    {
      int x;
      fseek(f, (xoffset * sar->N + yoffset) * sar->D * sar->D * sizeof(unsigned short), SEEK_CUR);
      for (x = 0; x < width; ++x)
	{
	  if (fread((unsigned short *) sar->array + x * height * sar->D * sar->D, sizeof(unsigned short), height * sar->D * sar->D, f)
	      != (unsigned long) height * sar->D * sar->D)
	    {
	      sarerror_msg("Invalid array field", filename);
	      return NULL;
	    }
	  if (x + xoffset < sar->M - 1)
	    fseek(f, (sar->N - height) * sar->D * sar->D * sizeof(unsigned short), SEEK_CUR);
	}
    }
  fclose(f);
  sar->M = width;
  sar->N = height;
  if (swap)
    swapbyte((char*) sar->array, sizeof(unsigned short), target_size);
  for (k = target_size - 1; ; k--)
    {
      ((float complex*) sar->array)[k] = ((unsigned short*) sar->array)[k];
      if (k == 0)
	break;
    }
  return sar;
}

sardata* imsread_extract(const char* filename, sardata* sar, char swap,
			 long int xoffset, long int yoffset,
			 long int width, long int height)
{
  unsigned long int input_size, total_size, target_size;
  FILE* f;
  int k;

  input_size = sar->M * sar->N * sar->D * sar->D;
  if (!(sar = ximaread_header(filename, sar, swap)))
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
  total_size = total_size;
  target_size = width * height * sar->D * sar->D;
  if (input_size != target_size)
    if (!(sar->array = realloc(sar->array, target_size * sizeof(float complex))))
      {
	sarerror_msg_perror("Cannot allocate array field");
	return NULL;
      }
  if (!(f = fopen(filename, "r")))
    {
      sarerror_msg_perror("Cannot open %s file", filename);
      return NULL;
    }
  if (height == sar->N)
    {
      fseek(f, (xoffset * sar->N * sar->D * sar->D) * sizeof(short), SEEK_CUR);
      if (fread(sar->array, sizeof(short), target_size, f) != target_size)
	{
	  sarerror_msg("Cannot read data in %s file", filename);
	  return NULL;
	}
    }
  else
    {
      int x;
      fseek(f, (xoffset * sar->N + yoffset) * sar->D * sar->D * sizeof(short), SEEK_CUR);
      for (x = 0; x < width; ++x)
	{
	  if (fread((short *) sar->array + x * height * sar->D * sar->D, sizeof(short), height * sar->D * sar->D, f)
	      != (unsigned long) height * sar->D * sar->D)
	    {
	      sarerror_msg("Invalid array field", filename);
	      return NULL;
	    }
	  if (x + xoffset < sar->M - 1)
	    fseek(f, (sar->N - height) * sar->D * sar->D * sizeof(short), SEEK_CUR);
	}
    }
  fclose(f);
  sar->M = width;
  sar->N = height;
  if (swap)
    swapbyte((char*) sar->array, sizeof(short), target_size);
  for (k = target_size - 1; ; k--)
    {
      ((float complex*) sar->array)[k] = ((short*) sar->array)[k];
      if (k == 0)
	break;
    }
  return sar;
}

sardata* imfread_extract(const char* filename, sardata* sar, char swap,
			 long int xoffset, long int yoffset,
			 long int width, long int height)
{
  unsigned long int input_size, total_size, target_size;
  FILE* f;
  int k;

  input_size = sar->M * sar->N * sar->D * sar->D;
  if (!(sar = ximaread_header(filename, sar, swap)))
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
  total_size = total_size;
  target_size = width * height * sar->D * sar->D;
  if (input_size != target_size)
    if (!(sar->array = realloc(sar->array, target_size * sizeof(float complex))))
      {
	sarerror_msg_perror("Cannot allocate array field");
	return NULL;
      }
  if (!(f = fopen(filename, "r")))
    {
      sarerror_msg_perror("Cannot open %s file", filename);
      return NULL;
    }
  if (height == sar->N)
    {
      fseek(f, (xoffset * sar->N * sar->D * sar->D) * sizeof(float), SEEK_CUR);
      if (fread(sar->array, sizeof(float), target_size, f) != target_size)
	{
	  sarerror_msg("Cannot read data in %s file", filename);
	  return NULL;
	}
    }
  else
    {
      int x;
      fseek(f, (xoffset * sar->N + yoffset) * sar->D * sar->D * sizeof(float), SEEK_CUR);
      for (x = 0; x < width; ++x)
	{
	  if (fread((float *) sar->array + x * height * sar->D * sar->D, sizeof(float), height * sar->D * sar->D, f)
	      != (unsigned long) height * sar->D * sar->D)
	    {
	      sarerror_msg("Invalid array field", filename);
	      return NULL;
	    }
	  if (x + xoffset < sar->M - 1)
	    fseek(f, (sar->N - height) * sar->D * sar->D * sizeof(float), SEEK_CUR);
	}
    }
  fclose(f);
  sar->M = width;
  sar->N = height;
  if (swap)
    swapbyte((char*) sar->array, sizeof(float), 2 * target_size);
  for (k = target_size - 1; ; k--)
    {
      ((float complex*) sar->array)[k] = ((float*) sar->array)[k];
      if (k == 0)
	break;
    }
  return sar;
}


sardata* imdread_extract(const char* filename, sardata* sar, char swap,
			 long int xoffset, long int yoffset,
			 long int width, long int height)
{
  unsigned long int input_size, total_size, target_size;
  FILE* f;
  unsigned int k;

  input_size = sar->M * sar->N * sar->D * sar->D;
  if (!(sar = ximaread_header(filename, sar, swap)))
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
  total_size = total_size;
  target_size = width * height * sar->D * sar->D;
  if (input_size != target_size)
    if (!(sar->array = realloc(sar->array, target_size * sizeof(double complex))))
      {
	sarerror_msg_perror("Cannot allocate array field");
	return NULL;
      }
  if (!(f = fopen(filename, "r")))
    {
      sarerror_msg_perror("Cannot open %s file", filename);
      return NULL;
    }
  if (height == sar->N)
    {
      fseek(f, (xoffset * sar->N * sar->D * sar->D) * sizeof(double), SEEK_CUR);
      if (fread(sar->array, sizeof(double), target_size, f) != target_size)
	{
	  sarerror_msg("Cannot read data in %s file", filename);
	  return NULL;
	}
    }
  else
    {
      int x;
      fseek(f, (xoffset * sar->N + yoffset) * sar->D * sar->D * sizeof(double), SEEK_CUR);
      for (x = 0; x < width; ++x)
	{
	  if (fread((double *) sar->array + x * height * sar->D * sar->D, sizeof(double), height * sar->D * sar->D, f)
	      != (unsigned long) height * sar->D * sar->D)
	    {
	      sarerror_msg("Invalid array field", filename);
	      return NULL;
	    }
	  if (x + xoffset < sar->M - 1)
	    fseek(f, (sar->N - height) * sar->D * sar->D * sizeof(double), SEEK_CUR);
	}
    }
  fclose(f);
  sar->M = width;
  sar->N = height;
  if (swap)
    swapbyte((char*) sar->array, sizeof(double), 2 * target_size);
  for (k = 0; k < target_size; k++)
    ((float complex*) sar->array)[k] = ((double*) sar->array)[k];
  sar->array = realloc(sar->array, target_size * sizeof(float complex));
  return sar;
}

sardata* cxbread_extract(const char* filename, sardata* sar, char swap,
			 long int xoffset, long int yoffset,
			 long int width, long int height)
{
  unsigned long int input_size, total_size, target_size;
  FILE* f;
  int k;

  input_size = sar->M * sar->N * sar->D * sar->D;
  if (!(sar = ximaread_header(filename, sar, swap)))
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
  total_size = total_size;
  target_size = width * height * sar->D * sar->D;
  if (input_size != target_size)
    if (!(sar->array = realloc(sar->array, target_size * sizeof(float complex))))
      {
	sarerror_msg_perror("Cannot allocate array field");
	return NULL;
      }
  if (!(f = fopen(filename, "r")))
    {
      sarerror_msg_perror("Cannot open %s file", filename);
      return NULL;
    }
  if (height == sar->N)
    {
      fseek(f, (xoffset * sar->N * sar->D * sar->D) * sizeof(char complex), SEEK_CUR);
      if (fread(sar->array, sizeof(char complex), target_size, f) != target_size)
	{
	  sarerror_msg("Cannot read data in %s file", filename);
	  return NULL;
	}
    }
  else
    {
      int x;
      fseek(f, (xoffset * sar->N + yoffset) * sar->D * sar->D * sizeof(char complex), SEEK_CUR);
      for (x = 0; x < width; ++x)
	{
	  if (fread(sar->array + x * height * sar->D * sar->D, sizeof(char complex), height * sar->D * sar->D, f)
	      != (unsigned long) height * sar->D * sar->D)
	    {
	      sarerror_msg("Invalid array field", filename);
	      return NULL;
	    }
	  if (x + xoffset < sar->M - 1)
	    fseek(f, (sar->N - height) * sar->D * sar->D * sizeof(char complex), SEEK_CUR);
	}
    }
  fclose(f);
  sar->M = width;
  sar->N = height;
  if (swap)
    swapbyte((char*) sar->array, sizeof(char), 2 * target_size);
  for (k = target_size - 1; ; k--)
    {
      ((float*) sar->array)[2*k+1] = ((char*) sar->array)[2*k+1];
      ((float*) sar->array)[2*k] = ((char*) sar->array)[2*k];
      if (k == 0)
	break;
    }
  return sar;
}

sardata* cxsread_extract(const char* filename, sardata* sar, char swap,
			 long int xoffset, long int yoffset,
			 long int width, long int height)
{
  unsigned long int input_size, total_size, target_size;
  FILE* f;
  int k;

  input_size = sar->M * sar->N * sar->D * sar->D;
  if (!(sar = ximaread_header(filename, sar, swap)))
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
  total_size = total_size;
  target_size = width * height * sar->D * sar->D;
  if (input_size != target_size)
    if (!(sar->array = realloc(sar->array, target_size * sizeof(float complex))))
      {
	sarerror_msg_perror("Cannot allocate array field");
	return NULL;
      }
  if (!(f = fopen(filename, "r")))
    {
      sarerror_msg_perror("Cannot open %s file", filename);
      return NULL;
    }
  if (height == sar->N)
    {
      fseek(f, (xoffset * sar->N * sar->D * sar->D) * sizeof(short complex), SEEK_CUR);
      if (fread(sar->array, sizeof(short complex), target_size, f) != target_size)
	{
	  sarerror_msg("Cannot read data in %s file", filename);
	  return NULL;
	}
    }
  else
    {
      int x;
      fseek(f, (xoffset * sar->N + yoffset) * sar->D * sar->D * sizeof(short complex), SEEK_CUR);
      for (x = 0; x < width; ++x)
	{
	  if (fread(sar->array + x * height * sar->D * sar->D, sizeof(short complex), height * sar->D * sar->D, f)
	      != (unsigned long) height * sar->D * sar->D)
	    {
	      sarerror_msg("Invalid array field", filename);
	      return NULL;
	    }
	  if (x + xoffset < sar->M - 1)
	    fseek(f, (sar->N - height) * sar->D * sar->D * sizeof(short complex), SEEK_CUR);
	}
    }
  fclose(f);
  sar->M = width;
  sar->N = height;
  if (swap)
    swapbyte((char*) sar->array, sizeof(short), 2 * target_size);
  for (k = target_size - 1; ; k--)
    {
      ((float*) sar->array)[2*k+1] = ((short*) sar->array)[2*k+1];
      ((float*) sar->array)[2*k] = ((short*) sar->array)[2*k];
      if (k == 0)
	break;
    }
  return sar;
}

sardata* cxfread_extract(const char* filename, sardata* sar, char swap,
			 long int xoffset, long int yoffset,
			 long int width, long int height)
{
  unsigned long int input_size, total_size, target_size;
  FILE* f;

  input_size = sar->M * sar->N * sar->D * sar->D;
  if (!(sar = ximaread_header(filename, sar, swap)))
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
  total_size = total_size;
  target_size = width * height * sar->D * sar->D;
  if (input_size != target_size)
    if (!(sar->array = realloc(sar->array, target_size * sizeof(float complex))))
      {
	sarerror_msg_perror("Cannot allocate array field");
	return NULL;
      }
  if (!(f = fopen(filename, "r")))
    {
      sarerror_msg_perror("Cannot open %s file", filename);
      return NULL;
    }
  if (height == sar->N)
    {
      fseek(f, (xoffset * sar->N * sar->D * sar->D) * sizeof(float complex), SEEK_CUR);
      if (fread(sar->array, sizeof(float complex), target_size, f) != target_size)
	{
	  sarerror_msg("Cannot read data in %s file", filename);
	  return NULL;
	}
    }
  else
    {
      int x;
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
  fclose(f);
  sar->M = width;
  sar->N = height;
  if (swap)
    swapbyte((char*) sar->array, sizeof(float), 2 * target_size);
  return sar;
}

sardata* cxdread_extract(const char* filename, sardata* sar, char swap,
			 long int xoffset, long int yoffset,
			 long int width, long int height)
{
  unsigned long int input_size, total_size, target_size;
  FILE* f;
  unsigned int k;

  input_size = sar->M * sar->N * sar->D * sar->D;
  if (!(sar = ximaread_header(filename, sar, swap)))
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
  total_size = total_size;
  target_size = width * height * sar->D * sar->D;
  if (input_size != target_size)
    if (!(sar->array = realloc(sar->array, target_size * sizeof(double complex))))
      {
	sarerror_msg_perror("Cannot allocate array field");
	return NULL;
      }
  if (!(f = fopen(filename, "r")))
    {
      sarerror_msg_perror("Cannot open %s file", filename);
      return NULL;
    }
  if (height == sar->N)
    {
      fseek(f, (xoffset * sar->N * sar->D * sar->D) * sizeof(double complex), SEEK_CUR);
      if (fread(sar->array, sizeof(double complex), target_size, f) != target_size)
	{
	  sarerror_msg("Cannot read data in %s file", filename);
	  return NULL;
	}
    }
  else
    {
      int x;
      fseek(f, (xoffset * sar->N + yoffset) * sar->D * sar->D * sizeof(double complex), SEEK_CUR);
      for (x = 0; x < width; ++x)
	{
	  if (fread(sar->array + x * height * sar->D * sar->D, sizeof(double complex), height * sar->D * sar->D, f)
	      != (unsigned long) height * sar->D * sar->D)
	    {
	      sarerror_msg("Invalid array field", filename);
	      return NULL;
	    }
	  if (x + xoffset < sar->M - 1)
	    fseek(f, (sar->N - height) * sar->D * sar->D * sizeof(double complex), SEEK_CUR);
	}
    }
  fclose(f);
  sar->M = width;
  sar->N = height;
  if (swap)
    swapbyte((char*) sar->array, sizeof(double), 2 * target_size);
  for (k = 0; k < target_size; k++)
    {
      ((float*) sar->array)[2*k] = ((double*) sar->array)[2*k];
      ((float*) sar->array)[2*k+1] = ((double*) sar->array)[2*k+1];
    }
  sar->array = realloc(sar->array, target_size * sizeof(float complex));
  return sar;
}

sardata* imaread(const char* filename, sardata* sar, char swap)
{
  return imaread_extract(filename, sar, swap, 0, 0, -1, -1);
}

sardata* imwread(const char* filename, sardata* sar, char swap)
{
  return imwread_extract(filename, sar, swap, 0, 0, -1, -1);
}

sardata* imsread(const char* filename, sardata* sar, char swap)
{
  return imsread_extract(filename, sar, swap, 0, 0, -1, -1);
}

sardata* imfread(const char* filename, sardata* sar, char swap)
{
  return imfread_extract(filename, sar, swap, 0, 0, -1, -1);
}

sardata* imdread(const char* filename, sardata* sar, char swap)
{
  return imfread_extract(filename, sar, swap, 0, 0, -1, -1);
}

sardata* cxbread(const char* filename, sardata* sar, char swap)
{
  return cxbread_extract(filename, sar, swap, 0, 0, -1, -1);
}

sardata* cxsread(const char* filename, sardata* sar, char swap)
{
  return cxsread_extract(filename, sar, swap, 0, 0, -1, -1);
}

sardata* cxfread(const char* filename, sardata* sar, char swap)
{
  return cxfread_extract(filename, sar, swap, 0, 0, -1, -1);
}

sardata* cxdread(const char* filename, sardata* sar, char swap)
{
  return cxfread_extract(filename, sar, swap, 0, 0, -1, -1);
}

int cxfwrite(const sardata* sar, const char* filename)
{
  unsigned long int total_size;
  FILE* f;
  char* fn_dim;
  int len;

  if (sar->D > 1)
    {
      sarerror_msg("Case D>1 not implemented\n");
      return 0;
    }
  len = strlen(filename);
  if (len < 3)
    {
      sarerror_msg("Cannot add .dim extension: Filename too short");
      return 0;
    }
  if (!(f = fopen(filename, "w")))
    {
      sarerror_perror();
      return 0;
    }
  total_size = sar->M * sar->N;
  fwrite(sar->array, sizeof(float complex), total_size, f);
  fclose(f);

  fn_dim = strdup(filename);
  fn_dim[len - 3] = 'd';
  fn_dim[len - 2] = 'i';
  fn_dim[len - 1] = 'm';
  if (!(f = fopen(fn_dim, "w")))
    {
      sarerror_perror();
      return 0;
    }
  fprintf(f, "%d %d\n", sar->N, sar->M);
  fclose(f);

  return 1;
}



