/*
** sarwrite.c: Matlab interface to write SAR data to disk
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
** Started on  Wed Jul 24 16:02:34 2013 Charles-Alban Deledalle
** Last update Fri Aug 23 14:59:12 2013 Charles-Alban Deledalle
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
  sprintf(str, "usage: sarwrite(sarin, fileout)\n");
  mexErrMsgTxt(str);
}

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  plhs = plhs;
  if (nrhs != 2 || nlhs > 0)
    {
      usage();
      return;
    }
  if (mxIsSingle(prhs[0]) != 1 || mxIsNumeric(prhs[0]) != 1 ||
      mxGetNumberOfDimensions(prhs[0]) != 4 ||
      !mxIsChar(prhs[1]) || mxGetM(prhs[1]) != 1)
    {
      usage();
      return;
    }
  int fn_out_len = mxGetM(prhs[1]) * mxGetN(prhs[1]) + 1;
  char* fn_out  = mxCalloc(fn_out_len, sizeof(char));
  mxGetString(prhs[1], fn_out, fn_out_len);

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
  int i, j, k, l;
  if (mxIsComplex(prhs[0]))
    {
      float* Pi = mxGetImagData(prhs[0]);
      for (i = 0; i < M; ++i)
	for (j = 0; j < N; ++j)
	  for (k = 0; k < D; ++k)
	    for (l = 0; l < D; ++l)
	      SARDATA_ACCESS(sar_in, i, j, k, l) =
		Pr[(((i * N) + j) * D + k) * D + l]
		+ I * Pi[(((i * N) + j) * D + k) * D + l];
    }
  else
    for (i = 0; i < M; ++i)
      for (j = 0; j < N; ++j)
	for (k = 0; k < D; ++k)
	  for (l = 0; l < D; ++l)
	    SARDATA_ACCESS(sar_in, i, j, k, l) =
	      Pr[(((i * N) + j) * D + k) * D + l];
  if (!(sarwrite(sar_in, fn_out)))
    {
      sarerror_msg_msg("Cannot create file %s", fn_out);
      mexErrMsgTxt(sarerror);
      return;
    }
  sardata_free(sar_in);
}
