/*
** noisegen.c: Noise generators
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
** Started on  Mon Aug 19 09:39:49 2013 Charles-Alban Deledalle
** Last update Mon Aug 19 09:46:38 2013 Charles-Alban Deledalle
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "tools/sarerror.h"
#include "tools/mathtools.h"
#include "algos/noisegen/noisegen.h"
#include "data/sardata.h"

#define PI	3.141592653589793115997963468544
#define SQRT2	1.414213562373095145474621858739

static float randn()
{
  return sqrtf(-2.0*logf((float) rand() / RAND_MAX)) *
    cosf(2.0 * PI * (float) rand() / RAND_MAX);
}

sardata* wishartrnd(int M, int N, int D, int L, float rho)
{
  sardata* output;
  int i, j, k, l, n;
  float* chol = NULL;
  float a;

  if (!(output = sardata_alloc_size(M, N, D)))
    {
      sarerror_perror();
      return NULL;
    }

  chol = calloc(D * D, sizeof(float complex));
  chol[0] = 1;
  for (k = 1; k < D; ++k)
    {
      a = 0;
      for (l = 0; l <= k - 2; ++l)
	{
	  chol[l + D * k] = chol[l + D * (k-1)];
	  a += chol[l + D * k] * chol[l + D * k];
	}
      chol[(k-1) + D * k] = (rho - a) / chol[(k-1) + D * (k-1)];
      a += chol[(k-1) + D * k] * chol[(k-1) + D * k];
      chol[k + D * k] = sqrtf(1 - a);
    }

  float complex* z = malloc(D * sizeof(float complex));
  float complex* Lz = malloc(D * sizeof(float complex));
  for (i = 0; i < M; ++i)
    for (j = 0; j < N; ++j)
      {
	for (k = 0; k < D; ++k)
	  for (l = 0; l < D; ++l)
	    SARDATA_ACCESS(output, i, j, k, l) = 0;
	for (n = 0; n < L; ++n)
	  {
	    for (k = 0; k < D; ++k)
	      z[k] = (randn() + I * randn()) / SQRT2;
	    for (k = 0; k < D; ++k)
	      {
		Lz[k] = 0;
		for (l = 0; l < D; ++l)
		  Lz[k] += chol[k + D * l] * z[l];
	      }
	    for (k = 0; k < D; ++k)
	      for (l = 0; l < D; ++l)
		SARDATA_ACCESS(output, i, j, k, l) += Lz[k] * conjf(Lz[l]);
	  }
	for (k = 0; k < D; ++k)
	  for (l = 0; l < D; ++l)
	    SARDATA_ACCESS(output, i, j, k, l) /= L;
      }
  free(z);

  return output;
}
