/*
** sarsimstats.c: estimation of the statistis of similarity criteria
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
** Started on  Wed Jul 24 15:45:50 2013 Charles-Alban Deledalle
** Last update Wed Jul 24 15:46:00 2013 Charles-Alban Deledalle
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "data/sardata.h"
#include "tools/sarerror.h"
#include "tools/mathtools.h"
#include "sarsim.h"
#include "sarsimstats.h"

static int compf(const void *p1, const void *p2)
{
  if ((* (float*) p1) > (* (float*) p2))
    return 1;
  if ((* (float*) p1) < (* (float*) p2))
    return -1;
  return 0;
}

sarsimstats* sarsimstats_create(const sardata* input,
				const sarsimfuncs* funcs,
				sarsimstats* stats,
				int L,
				int hP,
				float alpha_min,
				float alpha_max,
				int nb)
{
  int M = input->M;
  int N = input->N;
  int D = input->D;
  int Mext, Next;
  int x, y, x_dx, y_dy, k;
  float q_min, q_max, q;

  Mext = M+2*hP+1;
  Next = N+2*hP+1;

  if (!(0 <= alpha_min && alpha_min < alpha_max && alpha_max <= 1))
    {
      sarerror_msg("Quantiles sould verify 0 <= alpha_min < alpha_max <= 1");
      return NULL;
    }
  if (nb <= 1)
    {
      sarerror_msg("Number of bins should be at least 2");
      return NULL;
    }
  if (!stats)
    {
      if (!(stats = malloc(sizeof(sarsimstats))))
	{
	  sarerror_msg("Cannot allocate memory");
	  return NULL;
	}
      stats->alpha = NULL;
      stats->quantile = NULL;
    }
  stats->N = nb;
  float* diff_cum = malloc(Mext * Next * sizeof(float));
  float* diff = malloc(M * N * sizeof(float));
  if (!diff && !diff_cum)
    {
      sarerror_msg("Cannot allocate memory");
      return NULL;
    }
  sarsimdata* inputsarsim = funcs->create(L, input);
  if (!inputsarsim)
    {
      sarerror_perror();
      return NULL;
    }
  int dx = 2 * hP + 1;
  int dy = 2 * hP + 1;
  float mean, m2, min, max;

  // Compute similarities
  for (x = -hP-1; x < M+hP; ++x)
    {
      for (y = -hP-1; y < N+hP; ++y)
	{
	  x_dx = MOD(x + dx, M);
	  y_dy = MOD(y + dy, N);
	  diff_cum[(x+1+hP) * Next + (y+1+hP)] =
	    funcs->lsarsim(D, L,
			   SARSIMDATA_ACCESS(inputsarsim, MOD(x, M), MOD(y, N)),
			   SARSIMDATA_ACCESS(inputsarsim, x_dx, y_dy));
	}
    }
  // Compute commulative sums
  for (y = 1; y < Next; ++y)
    diff_cum[y] += diff_cum[y-1];
  for (x = 1; x < Mext; ++x)
    diff_cum[x * Next] += diff_cum[(x-1) * Next];
  for (x = 1; x < Mext; ++x)
    for (y = 1; y < Next; ++y)
      diff_cum[x * Next + y] +=
	+ diff_cum[(x-1) * Next + y]
	+ diff_cum[x     * Next + (y-1)]
	- diff_cum[(x-1) * Next + (y-1)];

  // Compute sums
  for (x = 0; x < M; ++x)
    for (y = 0; y < N; ++y)
      {
	x_dx = MOD(x + dx, M);
	y_dy = MOD(y + dy, N);

	diff[x * N + y] =
	  + diff_cum[(x+1+hP+hP)   * Next + (y+1+hP+hP)]
	  - diff_cum[(x+1+hP+hP)   * Next + (y+1+hP-hP-1)]
	  - diff_cum[(x+1+hP-hP-1) * Next + (y+1+hP+hP)]
	  + diff_cum[(x+1+hP-hP-1) * Next + (y+1+hP-hP-1)];
	diff[x * N + y] /= (2*hP+1) * (2*hP+1);
      }

  // Compute statistics
  stats->selfsim =
    funcs->lsarsim(D, L,
		   SARSIMDATA_ACCESS(inputsarsim, MOD(dx, M), MOD(dy, N)),
		   SARSIMDATA_ACCESS(inputsarsim, MOD(dx, M), MOD(dy, N)));
  funcs->free(inputsarsim);
  float avNorm = 0.0;
  mean = 0.0 ;
  m2 = 0.0;
  min = INFINITY;
  max = -INFINITY;
  for (k = 0; k < M * N; ++k)
    {
      //if (isnan(diff[k]) || isinf(diff[k]))
      //continue;
      mean += diff[k];
      m2 += diff[k] * diff[k];
      min  = MIN(diff[k], min);
      max  = MAX(diff[k], max);
      avNorm += 1;
    }
  mean /= avNorm;
  m2 /= avNorm;
  stats->mean = mean;
  stats->std = sqrtf(m2 - mean * mean);
  stats->min = min;
  stats->max = max;
  qsort(diff, M * N, sizeof(float), compf);
  q_min = diff[(int) ((M * N - 1) * alpha_min)];
  q_max = diff[(int) ((M * N - 1) * alpha_max)];
  q_min = q_min;
  q_max = q_max;
  if (!stats->alpha)
    if (!(stats->alpha = malloc(nb * sizeof(float))))
      {
	sarerror_perror();
	return NULL;
      }
  if (!stats->quantile)
    if (!(stats->quantile = malloc(nb * sizeof(float))))
      {
	sarerror_perror();
	return NULL;
      }
  float alpha;
  for (k = 0; k < nb; ++k)
    {
      alpha = alpha_min + k * (alpha_max - alpha_min) / (nb - 1);
      q = diff[(int) ((M * N - 1) * alpha)];
      stats->alpha[k] = alpha;
      stats->quantile[k] = q;
    }
  free(diff);
  free(diff_cum);

  return stats;
}


sarsimstats* sarsimstats_free(sarsimstats* stats)
{
  if (stats)
    {
      if (stats->quantile)
	free(stats->quantile);
      if (stats->alpha)
	free(stats->alpha);
      free(stats);
    }
  return NULL;
}
