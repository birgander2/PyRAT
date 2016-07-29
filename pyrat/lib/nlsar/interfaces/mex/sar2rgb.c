/*
** sar2rgb.c: Matlab interface to RGB export of SAR data
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
** Started on  Wed Jul 24 16:05:38 2013 Charles-Alban Deledalle
** Last update Fri Aug 23 17:10:10 2013 Charles-Alban Deledalle
*/

#include <mex.h>
#include <stdio.h>
#include <string.h>
#include "tools/sarerror.h"
#include "data/sardata.h"
#include "data/iosar.h"

static void usage()
{
  char str[1024];
  sprintf(str, "usage: sar2rgb(sarin, alpha = 3, gamma = 0.7)\n");
  mexErrMsgTxt(str);
}

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  if (nrhs < 1 || nrhs > 3 || nlhs > 1)
    {
      usage();
      return;
    }
  if (!mxIsSingle(prhs[0]) ||
      mxGetNumberOfDimensions(prhs[0]) != 4)
    {
      usage();
      return;
    }
  float alpha = 3;
  if (nrhs >= 2)
    {
      if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) ||
	  mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 1)
	{
	  usage();
	  return;
	}
      alpha = *mxGetPr(prhs[1]);
    }
  float gamma = 0.7;
  if (nrhs >= 3)
    {
      if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) ||
	  mxGetM(prhs[2]) != 1 || mxGetN(prhs[2]) != 1)
	{
	  usage();
	  return;
	}
      gamma = *mxGetPr(prhs[2]);
    }

  const mwSize* sizes = mxGetDimensions(prhs[0]);
  if (sizes[0] != sizes[1])
    {
      usage();
      return;
    }
  int M = sizes[3];
  int N = sizes[2];
  int D = sizes[1];

  sardata* sar_in = sardata_alloc_size(M, N, D);
  float* Pr = mxGetData(prhs[0]);
  float* Pi = mxGetImagData(prhs[0]);

  int i, j, k, l;
  for (i = 0; i < M; ++i)
    for (j = 0; j < N; ++j)
      for (k = 0; k < D; ++k)
	for (l = 0; l < D; ++l)
	  if (mxIsComplex(prhs[0]))
	    SARDATA_ACCESS(sar_in, i, j, k, l) =
	      Pr[(((i * N) + j) * D + k) * D + l]
	      + I * Pi[(((i * N) + j) * D + k) * D + l];
	  else
	    SARDATA_ACCESS(sar_in, i, j, k, l) =
	      Pr[(((i * N) + j) * D + k) * D + l];
  rgbdata* rgb = rgbdata_alloc();
  if (!(rgb = sar2rgb(sar_in, rgb, alpha, gamma)))
    {
      sarerror_msg_msg("Cannot create an RGB representation");
      mexErrMsgTxt(sarerror);
      return;
    }

  M = rgb->M;
  N = rgb->N;
  if (nlhs >= 1)
    {
      mwSize sizes[3];
      sizes[0] = M;
      sizes[1] = N;
      sizes[2] = 3;
      plhs[0] = mxCreateNumericArray(3, sizes, mxUINT8_CLASS, mxREAL);
      unsigned char* Pr = mxGetData(plhs[0]);
      for (i = 0; i < M; ++i)
	for (j = 0; j < N; ++j)
	  for (k = 0; k < 3; ++k)
	    Pr[(k * N + j) * M + i] =
	      RGBDATA_ACCESS(rgb, i, j, k);;
    }
  rgbdata_free(rgb);
}
