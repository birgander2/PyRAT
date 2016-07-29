/*
** fltdata.c: manipulation of FLT data
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
** Started on  Wed Jul 24 15:54:17 2013 Charles-Alban Deledalle
** Last update Wed Jul 24 15:54:21 2013 Charles-Alban Deledalle
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools/sarerror.h"
#include "tools/mathtools.h"
#include "fltdata.h"

fltdata* fltdata_alloc()
{
  fltdata* new = malloc(sizeof(fltdata));
  new->M = 0;
  new->N = 0;
  new->D = 0;
  new->array = NULL;
  return new;
}

fltdata* fltdata_alloc_size(int M, int N, int D)
{
  fltdata* new = malloc(sizeof(fltdata));
  new->M = M;
  new->N = N;
  new->D = D;
  new->array = malloc(M * N * D * sizeof(float));
  return new;
}

fltdata* fltdata_calloc_size(int M, int N, int D)
{
  fltdata* new = malloc(sizeof(fltdata));
  new->M = M;
  new->N = N;
  new->D = D;
  new->array = calloc(M * N * D, sizeof(float));
  return new;
}

fltdata* fltdata_realloc_size(fltdata* flt, int M, int N, int D)
{
  int i, j;
  if (!flt)
    return fltdata_alloc_size(M, N, D);
  if (flt->N == N)
    {
      flt->M = M;
      return flt;
    }
  if (M * N > flt->M * flt->N)
    flt->array = realloc(flt->array, D * M * N * sizeof(unsigned char));
  if (N > flt->N)
    for (i = MIN(M, flt->M) - 1; i >= 0; --i)
      for (j = MIN(N, flt->N) - 1; j >= 0; --j)
	memcpy(&flt->array[(i * N + j) * D],
	       &FLTDATA_ACCESS(flt, i, j, 0),
	       D * sizeof(unsigned char));
  else
    for (i = 0; i < MIN(M, flt->M); ++i)
      for (j = 0; j < MIN(N, flt->N); ++j)
	memcpy(&flt->array[(i * N + j) * D],
	       &FLTDATA_ACCESS(flt, i, j, 0),
	       D * sizeof(unsigned char));
  if (M * N <= flt->M * flt->N)
    flt->array = realloc(flt->array, D * M * N * sizeof(unsigned char));
  flt->M = M;
  flt->N = N;
  return flt;
}

fltdata* fltdata_free(fltdata* old)
{
  if (old)
    {
      if (old->array)
	free(old->array);
      free(old);
    }
  return NULL;
}

fltdata* fltdata_copy(const fltdata* src, fltdata* dst)
{
  if (dst == src)
    return dst;
  if (dst->M != src->M || dst->N != src->N || dst->D != src->D)
    {
      dst->M = src->M;
      dst->N = src->N;
      dst->D = src->D;
      dst->array = realloc(dst->array,
			   dst->M * dst->N * dst->D * sizeof(float));
    }
  memcpy(dst->array, src->array, dst->M * dst->N * dst->D *  sizeof(float));
  return dst;
}

fltdata* fltdata_extract(const fltdata* src, fltdata* dst,
			 long int xoffset, long int yoffset,
			 long int width, long int height,
			 long int step)
{
  long int i, j;

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
			     dst->M * dst->N * dst->D * sizeof(float))))
    {
      sarerror_perror();
      return NULL;
    }
  if (step <= 1)
    for (i = 0; i < width; ++i)
      {
	memcpy(&FLTDATA_ACCESS(dst, i, 0, 0),
	       &FLTDATA_ACCESS(src, xoffset + i, yoffset, 0),
	       1 * dst->N * dst->D * sizeof(float));
      }
  else
    for (i = 0; i < width && i / step < dst->M; i += step)
      for (j = 0; j < height && j / step < dst->N; j += step)
	memcpy(&FLTDATA_ACCESS(dst, i / step, j / step, 0),
	       &FLTDATA_ACCESS(src, xoffset + i, yoffset + j, 0),
	       1 * dst->D * sizeof(float));
  return dst;
}

fltdata* fltdata_dup(const fltdata* src)
{
  fltdata* new = fltdata_alloc();
  return fltdata_copy(src, new);
}
