/*
** sarmire.c: CLI to generate a resolution target with speckle
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
** Started on  Wed Jul 24 14:45:02 2013 Charles-Alban Deledalle
** Last update Fri Aug 23 16:04:34 2013 Charles-Alban Deledalle
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "tools/sarerror.h"
#include "data/sardata.h"
#include "data/iosar.h"

#define PI	3.141592653589793115997963468544
#define SQRT2	1.414213562373095145474621858739

static float randn()
{
  return sqrtf(-2.0*logf((float) rand() / RAND_MAX)) *
    cosf(2.0 * PI * (float) rand() / RAND_MAX);
}

static int usage(const char* argv0)
{
  fprintf(stderr, "usage: %s fileout width height dim look\n", argv0);
  return 1;
}

int main(int argc, char* argv[])
{
  srand(time(0));
  if (argc < 6)
    return usage(argv[0]);

  char* fn_out  = argv[1];
  int M, N, D, L;
  if ((sscanf(argv[2], "%d", &M) != 1) |
      (sscanf(argv[3], "%d", &N) != 1) |
      (sscanf(argv[4], "%d", &D) != 1) |
      (sscanf(argv[5], "%d", &L) != 1))
    return usage(argv[0]);

  sardata* sar_out = sardata_alloc_size(M, N, D);
  int i, j, k, l, n;
  float complex* z = malloc(D * sizeof(float complex));

  for (i = 0; i < M; ++i)
    for (j = 0; j < N; ++j)
      {	for (k = 0; k < D; ++k)
	  for (l = 0; l < D; ++l)
	    SARDATA_ACCESS(sar_out, i, j, k, l) = 0;
	for (n = 0; n < L; ++n)
	  {
	    for (k = 0; k < D; ++k)
	      {
		if ((j*1.0/N < 0.5 - 0.25 * (i*1.0/M - 0.5) * (i*1.0/M - 0.5)  &&
		     ((j*1.0/N - 3./16) * (j*1.0/N - 3./16) +
		      (i*1.0/M - 3./16) * (i*1.0/M - 3./16)
		      > (1./16) * (1./16)) &&
		     ((j*1.0/N - 2.5/16) * (j*1.0/N - 2.5/16) +
		      (i*1.0/M - 6.5/16) * (i*1.0/M - 6.5/16)
		      > (0.5/16) * (0.5/16)) &&
		     ((j*1.0/N - 2.25/16) * (j*1.0/N - 2.25/16) +
		      (i*1.0/M - 10.25/16) * (i*1.0/M - 10.25/16)
		      > (0.25/16) * (0.25/16)) &&
		     ((j*1.0/N - 2.125/16) * (j*1.0/N - 2.125/16) +
		      (i*1.0/M - 14.125/16) * (i*1.0/M - 14.125/16)
		      > (0.125/16) * (0.125/16))))
		  z[k] = (randn() + I * randn()) / SQRT2;
		else
		  z[k] = 4 * (randn() + I * randn());
		if (3.0 / 4 + 0.025*cosf(8*2*PI*i*1.0/M) < j*1.0/N &&
		    j*1.0/N < 3.0 / 4 + 7.0/N + 0.025*cosf(8*2*PI*i*1.0/M))
		  z[k] = 0.5 * (randn() + I * randn());
		if ((i == 7 * M / 8) && (j == 3 * N / 8))
		  z[k] = 160 * (randn() + I * randn());
	      }
	    for (k = 0; k < D; ++k)
	      for (l = 0; l < D; ++l)
		SARDATA_ACCESS(sar_out, i, j, k, l) += z[k] * conjf(z[l]);
	  }
	for (k = 0; k < D; ++k)
	  for (l = 0; l < D; ++l)
	    SARDATA_ACCESS(sar_out, i, j, k, l) /= L;
      }
  free(z);
  if (!(sarwrite(sar_out, fn_out)))
    {
      sarerror_msg_msg("Cannot create file %s", fn_out);
      fprintf(stderr, "%s\n", sarerror);
      return 2;
    }
  sardata_free(sar_out);

  return 0;
}
