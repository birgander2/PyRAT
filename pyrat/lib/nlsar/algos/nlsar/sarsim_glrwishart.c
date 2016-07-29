/*
** sarsim_glrwishart.c: GLR between data following a Wishart distribution
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
** Started on  Wed Jul 24 15:46:27 2013 Charles-Alban Deledalle
** Last update Mon Feb 24 09:20:04 2014 Charles-Alban Deledalle
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "tools/sarerror.h"
#include "tools/mathtools.h"
#include "tools/matrixtools.h"
#include "sarsim.h"

#define DATA_WISHART_C(data, k, l)		(((float complex*) (data))[(k) * D + (l)])
#define DATA_WISHART_DET(data)			((float*) ((float complex*) (data) + D * D))

#define GLRDATA_WISHART_C(stats, i, j, k, l)	(DATA_WISHART_C(SARSIMDATA_ACCESS((stats), (i), (j)), (k), (l)))
#define GLRDATA_WISHART_DET(stats, i, j)	(DATA_WISHART_DET(SARSIMDATA_ACCESS((stats), (i), (j))))

sarsimdata* glrdata_wishart_create(int L, const sardata* sar)
{
  int M = sar->M;
  int N = sar->N;
  int D = sar->D;
  int i, j, k, l;
  sardata* sar_norm = sardata_dup(sar);
  float norm = 0;
  for (i = 0; i < M; ++i)
    for (j = 0; j < N; ++j)
      for (k = 0; k < D; ++k)
	norm += SARDATA_ACCESS(sar, i, j, k, k);
  norm /= M * N;
  for (i = 0; i < M; ++i)
    for (j = 0; j < N; ++j)
      for (k = 0; k < D; ++k)
	for (l = 0; l < D; ++l)
	  SARDATA_ACCESS(sar_norm, i, j, k, l) /= norm;
  sarsimdata* res = malloc(sizeof(sarsimdata));
  if (!res)
    {
      sarerror_perror();
      return NULL;
    }
  res->M = M;
  res->N = N;
  res->size = D * D * sizeof(float complex) + sizeof(float);
  if (!(res->data = malloc(M * N * res->size)))
    {
      sarerror_perror();
      return NULL;
    }
  if (L >= D)
    {
      for (i = 0; i < M; ++i)
	for (j = 0; j < N; ++j)
	  {
	    memcpy(&GLRDATA_WISHART_C(res, i, j, 0, 0),
		   &SARDATA_ACCESS(sar_norm, i, j, 0, 0),
		   D * D * sizeof(float complex));
	    *GLRDATA_WISHART_DET(res, i, j) =
	      det(D, &GLRDATA_WISHART_C(res, i, j, 0, 0));
	  }
    }
  else
    {
      for (i = 0; i < M; ++i)
	for (j = 0; j < N; ++j)
	  {
	    memcpy(&GLRDATA_WISHART_C(res, i, j, 0, 0),
		   &SARDATA_ACCESS(sar_norm, i, j, 0, 0),
		   D * D * sizeof(float complex));
	    for (k = 0; k < D; ++k)
	      for (l = 0; l < D; ++l)
		if (k != l)
		  GLRDATA_WISHART_C(res, i, j, k, l) *= powf(L / (float) D, 1.0/3.0);
	    *GLRDATA_WISHART_DET(res, i, j) =
		det(D, &GLRDATA_WISHART_C(res, i, j, 0, 0));
	  }
    }
  sardata_free(sar_norm);
  return res;
}

sarsimdata* glrdata_wishart_free(sarsimdata* sarsimdata)
{
  if (sarsimdata)
    {
      if (sarsimdata->data)
	free(sarsimdata->data);
      free(sarsimdata);
    }
  return NULL;
}

#define LOG2 0.6931471805599452862267639829951804131269
float lglr_wishart(int D, int L, const void* C1, const void* C2)
{
  int D2 = D * D;
  int k;
  float complex* C12 = malloc(D2 * sizeof(float complex));
  for (k = 0; k < D2; ++k)
    C12[k] = ((float complex*) C1)[k] + ((float complex*) C2)[k];
  float detC1 = *DATA_WISHART_DET(C1);
  float detC2 = *DATA_WISHART_DET(C2);
  float detC12 = det(D, C12);
  free(C12);

  return L * logf(detC12 * detC12 / detC1 / detC2) - 2 * L * D * LOG2;
}

const sarsimfuncs sarsim_glrwishart =
  { &glrdata_wishart_create,
    &glrdata_wishart_free,
    &lglr_wishart };
