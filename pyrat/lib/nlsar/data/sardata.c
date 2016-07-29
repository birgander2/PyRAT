/*
** sardata.c: manipulation of SAR data
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
** Started on  Wed Jul 24 15:51:36 2013 Charles-Alban Deledalle
** Last update Mon Oct  7 18:59:57 2013 Charles-Alban Deledalle
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools/sarerror.h"
#include "tools/mathtools.h"
#include "sardata.h"

sardata* sardata_alloc()
{
  sardata* new = malloc(sizeof(sardata));
  int k;
  new->M = 0;
  new->N = 0;
  new->D = 0;
  new->array = NULL;
  for (k = 0; k < SAR_EXTRA_SIZE; ++k)
    new->extra[k] = 0;
  return new;
}

sardata* sardata_alloc_size(int M, int N, int D)
{
  sardata* new = malloc(sizeof(sardata));
  int k;
  new->M = M;
  new->N = N;
  new->D = D;
  new->array = malloc(M * N * D * D * sizeof(float complex));
  for (k = 0; k < SAR_EXTRA_SIZE; ++k)
    new->extra[k] = 0;
  return new;
}

sardata* sardata_calloc_size(int M, int N, int D)
{
  sardata* new = malloc(sizeof(sardata));
  int k;
  new->M = M;
  new->N = N;
  new->D = D;
  new->array = calloc(M * N * D * D, sizeof(float complex));
  for (k = 0; k < SAR_EXTRA_SIZE; ++k)
    new->extra[k] = 0;
  return new;
}

sardata* sardata_realloc_size(sardata* sar, int M, int N, int D)
{
  if (!sar)
    return sardata_alloc_size(M, N, D);
  sar->array =
    realloc(sar->array, D * D * M * N * sizeof(float complex));
  sar->D = D;
  sar->M = M;
  sar->N = N;
  return sar;
}

sardata* sardata_free(sardata* old)
{
  if (old)
    {
      if (old->array)
	free(old->array);
      free(old);
    }
  return NULL;
}

sardata* sardata_copy(const sardata* src, sardata* dst)
{
  int k;
  if (dst == src)
    return dst;
  if (dst->M != src->M || dst->N != src->N || dst->D != src->D)
    {
      dst->M = src->M;
      dst->N = src->N;
      dst->D = src->D;
      dst->array = realloc(dst->array,
			   dst->M * dst->N * dst->D * dst->D * sizeof(float complex));
    }
  for (k = 0; k < SAR_EXTRA_SIZE; ++k)
    dst->extra[k] = src->extra[k];
  memcpy(dst->array, src->array, dst->M * dst->N * dst->D * dst->D * sizeof(float complex));
  return dst;
}

sardata* sardata_dup(const sardata* src)
{
  sardata* new = sardata_alloc();
  return sardata_copy(src, new);
}

sardata* sardata_extract(const sardata* src, sardata* dst,
			 long int xoffset, long int yoffset,
			 long int width, long int height,
			 long int step)
{
  long int i, j;
  int k;

  if (xoffset < 0 ||
      yoffset < 0 ||
      xoffset + width > src->M ||
      yoffset + height > src->N)
    {
      sarerror_msg("Limits out of bounds");
      return NULL;
    }
  dst->M = width / step;
  dst->N = height / step;
  dst->D = src->D;
  if (!(dst->array = realloc(dst->array,
			     dst->M * dst->N * dst->D * dst->D * sizeof(float complex))))
    {
      sarerror_perror();
      return NULL;
    }
  if (step <= 1)
    for (i = 0; i < width; ++i)
      {
	memcpy(&SARDATA_ACCESS(dst, i, 0, 0, 0),
	       &SARDATA_ACCESS(src, xoffset + i, yoffset, 0, 0),
	       1 * dst->N * dst->D * dst->D * sizeof(float complex));
      }
  else
    for (i = 0; i < width && i / step < dst->M; i += step)
      for (j = 0; j < height && j / step < dst->N; j += step)
	memcpy(&SARDATA_ACCESS(dst, i / step, j / step, 0, 0),
	       &SARDATA_ACCESS(src, xoffset + i, yoffset + j, 0, 0),
	       1 * dst->D * dst->D * sizeof(float complex));
  for (k = 0; k < SAR_EXTRA_SIZE; ++k)
    dst->extra[k] = src->extra[k];
  return dst;
}

sardata* sardata_cat(sardata* const* src_list, sardata* dst,
		     int nx, int ny)
{
  int i, j, k, x, offset;
  int M, N, D;
  for (i = 0; i < nx; ++i)
    for (j = 0; j < ny; ++j)
      if (src_list[i * ny + j]->M != src_list[i * ny]->M ||
	  src_list[i * ny + j]->N != src_list[j]->N ||
	  src_list[i * ny + j]->D != src_list[0]->D)
	{
	  sarerror_msg("Dimensions are inconsistent");
	  return NULL;
	}
  M = 0;
  N = 0;
  if (nx * ny == 0)
    D = 1;
  else
    D = src_list[0]->D;
  for (i = 0; i < nx; ++i)
    M = M + src_list[i * ny]->M;
  for (j = 0; j < ny; ++j)
    N = N + src_list[j]->N;
  if (dst->M != M || dst->N != N || dst->D != D)
    {
      dst->M = M;
      dst->N = N;
      dst->D = D;
      dst->array = realloc(dst->array,
			   dst->M * dst->N * dst->D * dst->D * sizeof(float complex));
    }
  offset = 0;
  for (i = 0; i < nx; ++i)
    {
      for (j = 0; j < ny; ++j)
	{
	  for (x = 0; x < src_list[i * ny + j]->M; ++x)
	    memcpy(dst->array + (offset + x * dst->N) * D * D,
		   src_list[i * ny + j]->array + x * src_list[i * ny + j]->N * D * D,
		   src_list[i * ny + j]->N * D * D * sizeof(float complex));
	  offset += src_list[i * ny + j]->N;
	}
      offset += dst->N * src_list[i * ny]->M - dst->N;
    }

  for (k = 0; k < SAR_EXTRA_SIZE; ++k)
    dst->extra[k] = src_list[0]->extra[k];
  return dst;
}

sardata* sardata_join(sardata* const* src_list, sardata* dst, int D)
{
  int i, j, k, l;
  unsigned long int old_size, total_size;
  int M, N;

  old_size = dst->M * dst->N * dst->D * dst->D;
  M = src_list[0]->M;
  N = src_list[0]->N;
  for (k = 1; k < D; ++k)
    {
      if (src_list[k]->M != M || src_list[k]->N != N)
	{
	  sarerror_msg("SAR input images dimensions mismatch");
	  return NULL;
	}
      if (src_list[k]->D != 1)
	{
	  sarerror_msg("At least one of the input image is not a mono-dimensional complex SAR images");
	  return NULL;
	}
    }
  dst->M = M;
  dst->N = N;
  dst->D = D;
  total_size = dst->M * dst->N * dst->D * dst->D;
  if (old_size != total_size)
    if (!(dst->array = realloc(dst->array, total_size * sizeof(float complex))))
      {
	sarerror_msg_perror("Cannot allocate array field");
	return NULL;
      }
  for (i = 0; i < M; ++i)
    for (j = 0; j < N; ++j)
      for (k = 0; k < D; ++k)
	for (l = 0; l < D; ++l)
	  SARDATA_ACCESS(dst, i, j, k, l) =
	    SARDATA_ACCESS(src_list[k], i, j, 0, 0) *
	    conjf(SARDATA_ACCESS(src_list[l], i, j, 0, 0));
  for (k = 0; k < SAR_EXTRA_SIZE; ++k)
    dst->extra[k] = 0;
  return dst;
}
