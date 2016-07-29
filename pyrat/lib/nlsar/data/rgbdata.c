/*
** rgbdata.c: manipulation of RGB data
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
** Started on  Wed Jul 24 15:52:10 2013 Charles-Alban Deledalle
** Last update Wed Jul 24 15:52:14 2013 Charles-Alban Deledalle
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools/sarerror.h"
#include "tools/mathtools.h"
#include "rgbdata.h"

rgbdata* rgbdata_alloc()
{
  rgbdata* new = malloc(sizeof(rgbdata));
  new->M = 0;
  new->N = 0;
  new->array = NULL;
  return new;
}

rgbdata* rgbdata_alloc_size(int M, int N)
{
  rgbdata* new = malloc(sizeof(rgbdata));
  new->M = M;
  new->N = N;
  new->array = malloc(3 * M * N * sizeof(unsigned char));
  return new;
}


rgbdata* rgbdata_calloc_size(int M, int N)
{
  rgbdata* new = malloc(sizeof(rgbdata));
  new->M = M;
  new->N = N;
  new->array = calloc(3 * M * N, sizeof(unsigned char));
  return new;
}

rgbdata* rgbdata_realloc_size(rgbdata* rgb, int M, int N)
{
  int i, j;
  if (!rgb)
    return rgbdata_alloc_size(M, N);
  if (rgb->N == N)
    {
      rgb->M = M;
      return rgb;
    }
  if (M * N > rgb->M * rgb->N)
    rgb->array = realloc(rgb->array, 3 * M * N * sizeof(unsigned char));
  if (N > rgb->N)
    for (i = MIN(M, rgb->M) - 1; i >= 0; --i)
      for (j = MIN(N, rgb->N) - 1; j >= 0; --j)
	memcpy(&rgb->array[(i * N + j) * 3],
	       &RGBDATA_ACCESS(rgb, i, j, 0),
	       3 * sizeof(unsigned char));
  else
    for (i = 0; i < MIN(M, rgb->M); ++i)
      for (j = 0; j < MIN(N, rgb->N); ++j)
	memcpy(&rgb->array[(i * N + j) * 3],
	       &RGBDATA_ACCESS(rgb, i, j, 0),
	       3 * sizeof(unsigned char));
  if (M * N <= rgb->M * rgb->N)
    rgb->array = realloc(rgb->array, 3 * M * N * sizeof(unsigned char));
  rgb->M = M;
  rgb->N = N;
  return rgb;
}

rgbdata* rgbdata_free(rgbdata* old)
{
  if (old)
    {
      if (old->array)
	free(old->array);
      free(old);
    }
  return NULL;
}

rgbdata* rgbdata_copy(const rgbdata* src, rgbdata* dst)
{
  if (dst == src)
    return dst;
  if (dst->M != src->M || dst->N != src->N)
    {
      dst->M = src->M;
      dst->N = src->N;
      dst->array = realloc(dst->array, dst->M * dst->N * 3 * sizeof(unsigned char));
    }
  memcpy(dst->array, src->array, dst->M * dst->N * 3 * sizeof(unsigned char));
  return dst;
}

rgbdata* rgbdata_extract(const rgbdata* src, rgbdata* dst,
			 int xoffset, int yoffset, int width, int height,
			 int step)
{
  int i, j;

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
  if (!(dst->array = realloc(dst->array,
			     dst->M * dst->N * 3 % sizeof(unsigned char))))
    {
      sarerror_perror();
      return NULL;
    }
  if (step <= 1)
    for (i = 0; i < width; ++i)
      memcpy(&RGBDATA_ACCESS(dst, i, 0, 0),
	     &RGBDATA_ACCESS(src, xoffset + i, yoffset, 0),
	     1 * dst->N * 3 * sizeof(unsigned char));
  else
    for (i = 0; i < width && i / step < dst->M; i += step)
      for (j = 0; j < height && j / step < dst->N; j += step)
	memcpy(&RGBDATA_ACCESS(dst, i / step, j / step, 0),
	       &RGBDATA_ACCESS(src, xoffset + i, yoffset + j, 0),
	       1 * 3 * sizeof(unsigned char));
  return dst;
}

rgbdata* rgbdata_dup(const rgbdata* src)
{
  rgbdata* new = rgbdata_alloc();
  return rgbdata_copy(src, new);
}
