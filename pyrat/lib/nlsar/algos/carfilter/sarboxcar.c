/*
** sarboxcar.c: implementation of the boxcar filter
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
** Started on  Wed Jul 24 14:58:54 2013 Charles-Alban Deledalle
** Last update Wed Jul 24 14:58:55 2013 Charles-Alban Deledalle
*/

#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include "data/sardata.h"
#include "tools/mathtools.h"
#include "carfilter.h"

sardata* sarboxcar(const sardata* input, sardata* output, int hW)
{
  int M = input->M;
  int N = input->N;
  int D = input->D;
  int di, dj, t, i_di, j_dj;
  int norm;

  fftwf_plan plan;

  int sizes[4];
  sizes[0] = M;
  sizes[1] = N;
  sizes[2] = D;
  sizes[3] = D;

  output = sardata_realloc_size(output, M, N, D);
  output = sardata_copy(input, output);

  fftwf_complex* kernel = fftwf_malloc(M * N * D * D * sizeof(fftwf_complex));
  for (t = 0; t < M * N * D * D; ++t)
    kernel[t] = 0;
  norm = 0;
  for (di = -hW; di <= hW; ++di)
    for (dj = -hW; dj <= hW; ++dj)
      {
	i_di = MOD(di, M);
	j_dj = MOD(dj, N);
	kernel[((i_di * N + j_dj) * D) * D] = 1;
	norm++;
      }
  plan = fftwf_plan_dft(4, sizes,
			kernel, kernel,
			FFTW_FORWARD, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);

  plan = fftwf_plan_dft(4, sizes,
			output->array, output->array,
			FFTW_FORWARD, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);

  for (t = 0; t < M * N * D * D; ++t)
    output->array[t] *= conj(kernel[t]) / (M * N * D * D) / norm;

  plan = fftwf_plan_dft(4, sizes,
			output->array, output->array,
			FFTW_BACKWARD, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
  fftwf_free(kernel);
  fftwf_cleanup();

  return output;
}
