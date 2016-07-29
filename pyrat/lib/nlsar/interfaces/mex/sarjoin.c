/*
** sarjoin.c: Matlab interface to join mono-dimensional complex images
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
** Started on  Wed Jul 24 16:03:58 2013 Charles-Alban Deledalle
** Last update Wed Jul 24 16:04:01 2013 Charles-Alban Deledalle
*/

#include <mex.h>
#include <string.h>
#include "data/sardata.h"
#include "tools/sarerror.h"
#include "algos/join/join.h"

static void usage()
{
  char str[1024];
  sprintf(str,
	  "usage: sarout = sarjoin(sarin1 [sarin2, ..., sarinN], sarout)\n");
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
  int c;
  int D = nrhs;
  sardata** sar_in = malloc(D * sizeof(sardata *));
  if (!sar_in)
    {
      sarerror_msg_perror("Cannot allocate memory");
      mexErrMsgTxt(sarerror);
      return;
    }
  for (c = 0; c < D; ++c)
    {
      if (mxIsSingle(prhs[c]) != 1 || mxIsComplex(prhs[c]) != 1 ||
	  mxGetNumberOfDimensions(prhs[c]) != 4)
	{
	  usage();
	  return;
	}
      const mwSize* sizes = mxGetDimensions(prhs[c]);
      if (sizes[0] != sizes[1])
	{
	  usage();
	  return;
	}
      int M = sizes[3];
      int N = sizes[2];
      if (sizes[0] != 1)
	{
	  sarerror_msg("At least one of the input is not a mono-dimensional complex SAR images");
	  mexErrMsgTxt(sarerror);
	  return;
	}
      sar_in[c] = sardata_alloc_size(M, N, 1);
      if (!sar_in)
	{
	  sarerror_msg_msg("Cannot allocate memory");
	  mexErrMsgTxt(sarerror);
	  return;
	}
      float* Pr = mxGetData(prhs[c]);
      float* Pi = mxGetImagData(prhs[c]);
      int i, j;
      for (i = 0; i < M; ++i)
	for (j = 0; j < N; ++j)
	  SARDATA_ACCESS(sar_in[c], i, j, 0, 0) =
	    Pr[(i * N) + j]
	    + I * Pi[(i * N) + j];
    }
  sardata* sar_out = sardata_alloc();
  if (!sar_out)
    {
      sarerror_msg_msg("Cannot allocate memory");
      mexErrMsgTxt(sarerror);
      return;
    }
  if (!(sar_out = sarjoin(sar_in, D, sar_out)))
    {
      sarerror_msg_msg("Cannot join images");
      mexErrMsgTxt(sarerror);
      return;
    }
  if (nlhs >= 1)
    {
      mwSize sizes[4];
      int M = sar_out->M;
      int N = sar_out->N;
      int D = sar_out->D;
      sizes[0] = D;
      sizes[1] = D;
      sizes[2] = N;
      sizes[3] = M;
      plhs[0] = mxCreateNumericArray(4, sizes, mxSINGLE_CLASS, mxCOMPLEX);
      float* Pr = mxGetData(plhs[0]);
      float* Pi = mxGetImagData(plhs[0]);
      int i, j, k, l;
      for (i = 0; i < M; ++i)
	for (j = 0; j < N; ++j)
	  for (k = 0; k < D; ++k)
	    for (l = 0; l < D; ++l)
	      {
		Pr[(((i * N) + j) * D + k) * D + l] = crealf(SARDATA_ACCESS(sar_out, i, j, k, l));
		Pi[(((i * N) + j) * D + k) * D + l] = cimagf(SARDATA_ACCESS(sar_out, i, j, k, l));
	      }
    }
  for (c = 0; c < D; ++c)
    sardata_free(sar_in[c]);
  free(sar_in);
  sardata_free(sar_out);
}
