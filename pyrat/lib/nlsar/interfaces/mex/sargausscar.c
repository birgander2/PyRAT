/*
** sargausscar.c: Matlab interface of the gausscar filter
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
** Started on  Wed Jul 24 16:04:59 2013 Charles-Alban Deledalle
** Last update Wed Jul 24 16:05:02 2013 Charles-Alban Deledalle
*/

#include <mex.h>
#include <string.h>
#include "data/sardata.h"
#include "tools/sarerror.h"
#include "algos/carfilter/carfilter.h"

static void usage()
{
  char str[1024];
  sprintf(str, "usage: sarout = sargausscar(sarin [, hW])\n");
  mexErrMsgTxt(str);
}

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  if (nrhs < 1 || nlhs > 1)
    {
      usage();
      return;
    }
  if (mxIsSingle(prhs[0]) != 1 || mxIsComplex(prhs[0]) != 1 ||
      mxGetNumberOfDimensions(prhs[0]) != 4 ||
      (nrhs > 1  && (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) ||
		     mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 1)))
    {
      usage();
      return;
    }
  int hW = nrhs > 1 ? (int) *mxGetPr(prhs[1]) : 1;

  const mwSize* sizes = mxGetDimensions(prhs[0]);
  if (sizes[0] != sizes[1])
    {
      usage();
      return;
    }
  int M = sizes[3];
  int N = sizes[2];
  int D = sizes[1];
  sardata* input = sardata_alloc_size(M, N, D);
  if (!input)
    {
      sarerror_msg_msg("Cannot allocate memory");
      mexErrMsgTxt(sarerror);
      return;
    }
  float* Pr = mxGetData(prhs[0]);
  float* Pi = mxGetImagData(prhs[0]);
  int i, j, k, l;
  for (i = 0; i < M; ++i)
    for (j = 0; j < N; ++j)
      for (k = 0; k < D; ++k)
	for (l = 0; l < D; ++l)
	  SARDATA_ACCESS(input, i, j, k, l) =
	    Pr[(((i * N) + j) * D + k) * D + l]
	    + I * Pi[(((i * N) + j) * D + k) * D + l];

  sardata* output  = sardata_alloc_size(M, N, D);
  if (!output)
    {
      sarerror_msg_msg("Cannot allocate memory");
      mexErrMsgTxt(sarerror);
      return;
    }

  output = sargausscar(input, output, hW);

  if (nlhs >= 1)
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

  sardata_free(input);
  sardata_free(output);
}
