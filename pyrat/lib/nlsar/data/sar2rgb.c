/*
** sar2rgb.c: RGB representation of a SAR data
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
** Started on  Wed Jul 24 15:51:53 2013 Charles-Alban Deledalle
** Last update Fri Aug 23 17:07:42 2013 Charles-Alban Deledalle
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools/sarerror.h"
#include "sardata.h"
#include "ionetpbm.h"

#define PI 3.14

static float compute_thr(const float* values, int M, int N, float alpha)
{
  int i, j;
  int MN = M * N;
  float mean = 0;
  float var = 0;
  float value;
  for (i = 0; i < M; ++i)
    for (j = 0; j < N; ++j)
    {
      value = values[i * N + j];
      mean += value;
      var += value * value;
      }
  mean = mean / MN;
  var = var / MN - mean * mean;
  return mean + alpha * sqrtf(var);
}

rgbdata*	 sar2rgb(const sardata* sar, rgbdata* rgb, float alpha, float gamma)
{
  int M = sar->M;
  int N = sar->N;
  int i, j, k;
  float* valuesR = malloc(M * N * sizeof(float));
  float* valuesG = NULL;
  float* valuesB = NULL;
  float zeta;
  float v;

  rgb = rgbdata_realloc_size(rgb, M, N);
  if (!rgb)
    {
      sarerror_msg_msg("Cannot allocate RGB image");
      return NULL;
    }
  if (!valuesR)
    {
      sarerror_msg_perror("Cannot allocate temporary memory");
      return NULL;
    }
  switch (sar->D)
    {
      case 1:
      case 2:
	for (i = 0; i < M; ++i)
	  for (j = 0; j < N; ++j)
	    valuesR[i * N + j] = sqrtf(cabsf(SARDATA_ACCESS(sar, i, j, 0, 0)));
	float thr = compute_thr(valuesR, M, N, alpha);
	for (i = 0; i < M; ++i)
	  for (j = 0; j < N; ++j)
	    {
	      valuesR[i * N + j] /= thr;
	      valuesR[i * N + j] = powf(valuesR[i * N + j], gamma);
	      valuesR[i * N + j] *= 255;
	      RGBDATA_ACCESS(rgb, i, j, 0) =
	      RGBDATA_ACCESS(rgb, i, j, 1) =
	      RGBDATA_ACCESS(rgb, i, j, 2) =
		valuesR[i * N + j] > 255 ? 255 : (int) valuesR[i * N + j];
	    }

	if (sar->D == 1)
	  break;

	rgb = rgbdata_realloc_size(rgb, M, 3 * N);
	for (i = 0; i < M; ++i)
	  for (j = 0; j < N; ++j)
	    {
	      zeta = 1;
	      if (cargf(cpowf(SARDATA_ACCESS(sar, i, j, 0, 1),zeta)) > 0)
		v = ((cargf(cpowf(SARDATA_ACCESS(sar, i, j, 0, 1),zeta))) / (PI));
	      else
		v = ((-cargf(cpowf(SARDATA_ACCESS(sar, i, j, 0, 1),zeta))) / (PI));
	      //v = (int) ((cargf(cpowf(SARDATA_ACCESS(sar, i, j, 0, 1),zeta))+PI) / (PI) * 255);
	      RGBDATA_ACCESS(rgb, i, j + N, 0) =
		RGBDATA_ACCESS(rgb, i, j + N, 1) =
		RGBDATA_ACCESS(rgb, i, j + N, 2) =
		(int) (v * 255);
	    }

	for (i = 0; i < M; ++i)
	  for (j = 0; j < N; ++j)
	    RGBDATA_ACCESS(rgb, i, j + 2 * N, 0) =
	      RGBDATA_ACCESS(rgb, i, j + 2 * N, 1) =
	      RGBDATA_ACCESS(rgb, i, j + 2 * N, 2) =
	      (int) ((2 * cabsf(SARDATA_ACCESS(sar, i, j, 0, 1)) /
		      (cabsf(SARDATA_ACCESS(sar, i, j, 0, 0)) +
		       cabsf(SARDATA_ACCESS(sar, i, j, 1, 1)))
		      ) * 255);
  	break;
      case 3:
      case 6:
	valuesG = malloc(M * N * sizeof(float));
	valuesB = malloc(M * N * sizeof(float));
	if (!valuesG || !valuesB)
	  {
	    sarerror_msg_perror("Cannot allocate temporary memory");
	    return NULL;
	  }
	for (i = 0; i < M; ++i)
	  for (j = 0; j < N; ++j)
	    {
	      valuesB[i * N + j] =
		sqrtf(cabsf(SARDATA_ACCESS(sar, i, j, 0, 0) + SARDATA_ACCESS(sar, i, j, 2, 2)
			     + 2 * crealf(SARDATA_ACCESS(sar, i, j, 0, 2))));
	      valuesR[i * N + j] =
		sqrtf(cabsf(SARDATA_ACCESS(sar, i, j, 0, 0) + SARDATA_ACCESS(sar, i, j, 2, 2)
			     - 2 * crealf(SARDATA_ACCESS(sar, i, j, 0, 2))));
	      valuesG[i * N + j] =
		sqrtf(cabsf(SARDATA_ACCESS(sar, i, j, 1, 1)));
  	    }
	float thrR = compute_thr(valuesR, M, N, alpha);
	float thrG = compute_thr(valuesG, M, N, alpha);
	float thrB = compute_thr(valuesB, M, N, alpha);
	thrR = thrB = (thrR + thrB) / 2;
	for (i = 0; i < M; ++i)
	  for (j = 0; j < N; ++j)
	    {
	      valuesR[i * N + j] /= thrR;
	      valuesG[i * N + j] /= thrG;
	      valuesB[i * N + j] /= thrB;
	      valuesR[i * N + j] = powf(valuesR[i * N + j], gamma);
	      valuesG[i * N + j] = powf(valuesG[i * N + j], gamma);
	      valuesB[i * N + j] = powf(valuesB[i * N + j], gamma);
	      valuesR[i * N + j] *= 255;
	      valuesG[i * N + j] *= 255;
	      valuesB[i * N + j] *= 255;
	      RGBDATA_ACCESS(rgb, i, j, 0) =
		valuesR[i * N + j] > 255 ? 255 : (int) valuesR[i * N + j];
	      RGBDATA_ACCESS(rgb, i, j, 1) =
		valuesG[i * N + j] > 255 ? 255 : (int) valuesG[i * N + j];
	      RGBDATA_ACCESS(rgb, i, j, 2) =
		valuesB[i * N + j] > 255 ? 255 : (int) valuesB[i * N + j];
	    }
	free(valuesG);
	free(valuesB);

	if (sar->D == 3)
	  break;

	rgb = rgbdata_realloc_size(rgb, M, 3 * N);
	for (i = 0; i < M; ++i)
	  for (j = 0; j < N; ++j)
		{
		  zeta = 1;
		  v = 0;
		  for (k = 0; k < 3; ++k)
		    {
		      if (cargf(cpowf(SARDATA_ACCESS(sar, i, j, k, 3+k),zeta)) > 0)
			v += ((cargf(cpowf(SARDATA_ACCESS(sar, i, j, k, 3+k),zeta))) / (PI)) / 3;
		      else
			v += ((-cargf(cpowf(SARDATA_ACCESS(sar, i, j, k, 3+k),zeta))) / (PI)) / 3;
		    }
		  //v = (int) ((cargf(cpowf(SARDATA_ACCESS(sar, i, j, 0, 1),zeta))+PI) / (PI) * 255);
		  RGBDATA_ACCESS(rgb, i , j + N, 0) =
		    RGBDATA_ACCESS(rgb, i , j + N, 1) =
		    RGBDATA_ACCESS(rgb, i , j + N, 2) =
		    (int) (v * 255);
		}
	    for (i = 0; i < M; ++i)
	      for (j = 0; j < N; ++j)
		{
		  v = 0;
		  for (k = 0; k < 3; ++k)
		    {
		      v += ((2 * cabsf(SARDATA_ACCESS(sar, i, j, k, 3+k)) /
			     (cabsf(SARDATA_ACCESS(sar, i, j, k, k)) +
			      cabsf(SARDATA_ACCESS(sar, i, j, 3+k, 3+k)))
			     )) / 3;
		    }
		  RGBDATA_ACCESS(rgb, i , j + 2 * N, 0) =
		    RGBDATA_ACCESS(rgb, i , j + 2 * N, 1) =
		    RGBDATA_ACCESS(rgb, i , j + 2 * N, 2) =
		    (int) (v * 255);
	  }
  	break;
      default:
	sarerror_msg("Dimension %d x %d not implemented", sar->D, sar->D);
	return NULL;
    }
  free(valuesR);
  return rgb;
}
