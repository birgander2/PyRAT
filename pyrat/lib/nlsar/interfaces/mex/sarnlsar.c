/*
** sarnlsar.c: Matlab interface of the non-local SAR filter
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
** Started on  Wed Jul 24 16:03:32 2013 Charles-Alban Deledalle
** Last update Fri Aug 23 15:29:40 2013 Charles-Alban Deledalle
*/

#include <mex.h>
#include <string.h>
#include "data/sardata.h"
#include "data/fltdata.h"
#include "tools/sarerror.h"
#include "algos/nlsar/nlsar.h"

static void usage()
{
  char str[1024];
  sprintf(str, "usage: [sarout lookmap] = sarnlsar(sarin, look, verbose, hW, hP, noise)\n");
  mexErrMsgTxt(str);
}

#define isScalarValue(a) (mxIsNumeric(a) && !mxIsComplex(a) && mxGetM(a) == 1 && mxGetN(a) == 1)
#define isString(a) (mxIsChar(a) && mxGetM(a) == 1)

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  float L;
  int hW = 12, hP = 5, verbose = 1;
  int M, N, D;
  sardata* input;
  sardata* sar_noise;
  sardata* output;
  fltdata* look = NULL;
  const mwSize* sizes;
  float* Pr;
  float* Pi;
  int i, j, k, l;

  if (nrhs < 2 || nrhs > 6 || nlhs > 2)
    {
      usage();
      return;
    }
  if (mxIsSingle(prhs[0]) != 1 || mxGetNumberOfDimensions(prhs[0]) != 4 ||
      !isScalarValue(prhs[1]) ||
      (nrhs > 2 && !isScalarValue(prhs[2])) ||
      (nrhs > 3 && !isScalarValue(prhs[3])) ||
      (nrhs > 4 && !isScalarValue(prhs[4])) ||
      (nrhs > 5 && (mxIsSingle(prhs[5]) != 1 || mxGetNumberOfDimensions(prhs[5]) != 4)))
    {
      usage();
      return;
    }

  L = (float) *mxGetPr(prhs[1]);
  if (nrhs > 2) verbose = (int) *mxGetPr(prhs[2]);
  if (nrhs > 3) hW      = (int) *mxGetPr(prhs[3]);
  if (nrhs > 4) hP      = (int) *mxGetPr(prhs[4]);

  // Load noise image
  if (!(nrhs > 5))
    sar_noise = NULL;
  else
    {
      sizes = mxGetDimensions(prhs[5]);
      if (sizes[0] != sizes[1])
	{
	  usage();
	  return;
	}
      M = sizes[3];
      N = sizes[2];
      D = sizes[1];
      if (!(sar_noise = sardata_alloc_size(M, N, D)))
	{
	  sarerror_msg_msg("Cannot allocate memory");
	  mexErrMsgTxt(sarerror);
	  return;
	}
      Pr = mxGetData(prhs[5]);
      Pi = mxGetImagData(prhs[5]);
      for (i = 0; i < M; ++i)
	for (j = 0; j < N; ++j)
	  for (k = 0; k < D; ++k)
	    for (l = 0; l < D; ++l)
	      if (mxIsComplex(prhs[5]))
		SARDATA_ACCESS(sar_noise, i, j, k, l) =
		  Pr[(((i * N) + j) * D + k) * D + l]
		  + I * Pi[(((i * N) + j) * D + k) * D + l];
	      else
		SARDATA_ACCESS(sar_noise, i, j, k, l) =
		  Pr[(((i * N) + j) * D + k) * D + l];
    }
  // Load input image
  sizes = mxGetDimensions(prhs[0]);
  if (sizes[0] != sizes[1])
    {
      usage();
      return;
    }
  M = sizes[3];
  N = sizes[2];
  D = sizes[1];
  input = sardata_alloc_size(M, N, D);
  if (!input)
    {
      sarerror_msg_msg("Cannot allocate memory");
      mexErrMsgTxt(sarerror);
      return;
    }
  Pr = mxGetData(prhs[0]);
  Pi = mxGetImagData(prhs[0]);
  for (i = 0; i < M; ++i)
    for (j = 0; j < N; ++j)
      for (k = 0; k < D; ++k)
	for (l = 0; l < D; ++l)
	  if (mxIsComplex(prhs[0]))
	    SARDATA_ACCESS(input, i, j, k, l) =
	      Pr[(((i * N) + j) * D + k) * D + l]
	      + I * Pi[(((i * N) + j) * D + k) * D + l];
	  else
	    SARDATA_ACCESS(input, i, j, k, l) =
	      Pr[(((i * N) + j) * D + k) * D + l];
  // Filter image
  output = sardata_alloc();
  if (nlhs < 2)
    output = sarnlsar(input, output, L, 4, verbose, hW, hP, sar_noise);
  else
    output = sarnlsar(input, output, L, 5, verbose, hW, hP, sar_noise, &look);
  if (!output)
    {
      mexErrMsgTxt(sarerror);
      return;
    }
  // Output
  if (nlhs > 0)
    {
      mwSize sizes[4];
      M = output->M;
      N = output->N;
      D = output->D;
      sizes[0] = D;
      sizes[1] = D;
      sizes[2] = N;
      sizes[3] = M;
      plhs[0] = mxCreateNumericArray(4, sizes, mxSINGLE_CLASS, mxCOMPLEX);
      float* Pr = mxGetData(plhs[0]);
      float* Pi = mxGetImagData(plhs[0]);
      for (i = 0; i < M; ++i)
	for (j = 0; j < N; ++j)
	  for (k = 0; k < D; ++k)
	    for (l = 0; l < D; ++l)
	      {
		Pr[(((i * N) + j) * D + k) * D + l] = crealf(SARDATA_ACCESS(output, i, j, k, l));
		Pi[(((i * N) + j) * D + k) * D + l] = cimagf(SARDATA_ACCESS(output, i, j, k, l));
	      }
    }
  if (nlhs > 1)
    {
      mwSize sizes[2];
      M = look->M;
      N = look->N;
      sizes[0] = N;
      sizes[1] = M;
      plhs[1] = mxCreateNumericArray(2, sizes, mxSINGLE_CLASS, mxREAL);
      float* Pr = mxGetData(plhs[1]);
      for (i = 0; i < M; ++i)
	for (j = 0; j < N; ++j)
	  Pr[i * N + j] = FLTDATA_ACCESS(look, i, j, 0);
    }
  if (sar_noise)
    sardata_free(sar_noise);
  if (look)
    fltdata_free(look);
  sardata_free(input);
  sardata_free(output);
}
