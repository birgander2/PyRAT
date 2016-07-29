/*
** sarinfo.c: Matlab interface to extract informations of SAR data files
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
** Started on  Wed Jul 24 16:04:30 2013 Charles-Alban Deledalle
** Last update Wed Jul 24 16:04:33 2013 Charles-Alban Deledalle
*/

#include <mex.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "tools/sarerror.h"
#include "data/sardata.h"
#include "data/iosar.h"

static void usage()
{
  char str[1024];
  sprintf(str, "usage: [M, N, D] = sarinfo(filename)\n");
  mexErrMsgTxt(str);
}

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  if (nrhs != 1 || nlhs > 3)
    {
      usage();
      return;
    }
  if (mxIsChar(prhs[0]) != 1 || mxGetM(prhs[0]) != 1)
    {
      usage();
      return;
    }

  int fn_in_len = mxGetM(prhs[0]) * mxGetN(prhs[0]) + 1;
  char* fn_in  = mxCalloc(fn_in_len, sizeof(char));
  mxGetString(prhs[0], fn_in, fn_in_len);

  sardata* sar_out = sardata_alloc();
  if (!(sar_out = sarread_header(fn_in, sar_out)))
    {
      sarerror_msg_msg("Cannot open file %s", fn_in);
      mexErrMsgTxt(sarerror);
      return;
    }

  if (nlhs >= 1)
    {
      plhs[0] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
      int* Pr = mxGetData(plhs[0]);
      *Pr = sar_out->M;
    }

  if (nlhs >= 2)
    {
      plhs[1] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
      int* Pr = mxGetData(plhs[1]);
      *Pr = sar_out->N;
    }

  if (nlhs >= 3)
    {
      plhs[2] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
      int* Pr = mxGetData(plhs[2]);
      *Pr = sar_out->D;
    }
  sardata_free(sar_out);
}
