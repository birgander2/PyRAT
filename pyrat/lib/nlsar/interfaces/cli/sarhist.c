/*
** sarhist.c: CLI for histograms of the 1st diagonal element of SAR data
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
** Started on  Wed Jul 24 14:45:55 2013 Charles-Alban Deledalle
** Last update Wed Jul 24 14:46:38 2013 Charles-Alban Deledalle
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools/sarerror.h"
#include "data/sardata.h"
#include "data/iosar.h"

static int usage(const char* argv0)
{
  fprintf(stderr, "usage: %s filein alpha_min alpha_max nb_bins [xoffset yoffset width height]\n", argv0);
  return 1;
}

static int compf(const void *p1, const void *p2)
{
  if ((* (float*) p1) > (* (float*) p2))
    return 1;
  if ((* (float*) p1) < (* (float*) p2))
    return -1;
  return 0;
}


int main(int argc, char* argv[])
{
  if (argc < 5)
    return usage(argv[0]);
  char* fn_in  = argv[1];
  int xoffset, yoffset, width, height, nb_bins;
  float alpha_min, alpha_max;
  if ((sscanf(argv[2], "%f", &alpha_min) != 1) |
      (sscanf(argv[3], "%f", &alpha_max) != 1) |
      (sscanf(argv[4], "%d", &nb_bins) != 1))
    return usage(argv[0]);

  if (!(0 <= alpha_min && alpha_min < alpha_max && alpha_max <= 1))
    {
      sarerror_msg("Quantiles sould verify 0 <= alpha_min < alpha_max <= 1");
      fprintf(stderr, "%s\n", sarerror);
      return 2;
    }
  sardata* sar_in = sardata_alloc();
  sardata* sar_out = sardata_alloc();
  if (!(sar_in = sarread(fn_in, sar_in)))
    {
      sarerror_msg_msg("Cannot open file %s", fn_in);
      fprintf(stderr, "%s\n", sarerror);
      return 2;
    }
  xoffset = 0;
  yoffset = 0;
  width = sar_in->M;
  height = sar_in->N;
  if ((argc > 8  && sscanf(argv[5], "%d", &xoffset) != 1) |
      (argc > 9  && sscanf(argv[6], "%d", &yoffset) != 1) |
      (argc > 10 && sscanf(argv[7], "%d", &width) != 1) |
      (argc > 11 && sscanf(argv[8], "%d", &height) != 1))
    {
      sardata_free(sar_in);
      usage(argv[0]);
      return 2;
    }
  if (!(sardata_extract(sar_in, sar_out, xoffset, yoffset, width, height, 1)))
    {
      sarerror_msg_msg("Cannot extract region");
      fprintf(stderr, "%s\n", sarerror);
      return 2;
    }

  int M = sar_out->M;
  int N = sar_out->N;
  int i, j, k;
  float q_min, q_max, q;
  float* values = malloc(M * N * sizeof(float));
  for (i = 0; i < M; ++i)
    for (j = 0; j < N; ++j)
      values[i * N + j] = sqrtf(SARDATA_ACCESS(sar_in, i, j, 0, 0));
  qsort(values, M * N, sizeof(float), compf);
  q_min = values[(int) ((M * N - 1) * alpha_min)];
  q_max = values[(int) ((M * N - 1) * alpha_max)];
  int* hist = malloc(nb_bins * sizeof(int));
  int total = 0;
  i = M * N * alpha_min;
  for (k = 0; k < nb_bins; ++k)
    {
      q = q_min + k * (q_max - q_min) / (nb_bins - 1);
      for (j = 0; i + j < M * N && values[i + j] <= q; ++j)
	;
      hist[k] = j;
      if (k == 1)
	++hist[k];
      total += j;
      printf("%.2e\t%d\n", q, hist[k]);
      i += j;
    }
  free(values);
  free(hist);
  sardata_free(sar_in);
  sardata_free(sar_out);

  return 0;
}
