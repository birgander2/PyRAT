/*
** nlsar.c: implementation of the NL-SAR filter
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
** Started on  Wed Jul 24 15:46:54 2013 Charles-Alban Deledalle
** Last update Tue Nov 12 10:20:01 2013 Charles-Alban Deledalle
*/

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <fftw3.h>
#ifdef OMP
# include <omp.h>
#endif //!OMP
#include "tools/sarprintf.h"
#include "tools/sarwaitbar.h"
#include "data/sardata.h"
#include "tools/sarerror.h"
#include "tools/mathtools.h"
#include "tools/matrixtools.h"
#include "algos/carfilter/carfilter.h"
#include "algos/noisegen/noisegen.h"
#include "sarsim.h"
#include "sarsimstats.h"
#include "nlsar.h"
#include "phi.h"

typedef int int2[2];

#ifdef SPIRALESQUARE
static void spirale(int hW, int2* shifts, int* lengths)
{
  int idx, k, x, y;

  idx = 0;
  for (k = 0; k <= hW; ++k)
    {
      lengths[k] = idx;
      for (x = -k; x <= k; ++x)
	{
	  shifts[idx][0] = x;
	  shifts[idx][1] = -k;
	  idx = idx + 1;
	}
      for (y = (-k+1); y <= k; ++y)
	{
	  shifts[idx][0] = k;
	  shifts[idx][1] = y;
	  idx = idx + 1;
	}
      for (x = (k-1); x >= -k; --x)
	{
	  shifts[idx][0] = x;
	  shifts[idx][1] = k;
	  idx = idx + 1;
	}
      for (y = (k-1); y >= (-k+1); --y)
	{
	  shifts[idx][0] = -k;
	  shifts[idx][1] = y;
	  idx = idx + 1;
	}
      lengths[k] = idx - lengths[k];
    }
}
#else
# define PI (3.14159265358979311599796346854419)

static void spirale(int hW, int2* shifts, int* lengths)
{
  float cx, cy;
  int idx, k, s, x, y;
  char* mask;

  mask = calloc((2*hW+3)*(2*hW+3), sizeof(char));

  idx = 0;
  for (k = 0; k <= hW; ++k)
    {
      lengths[k] = idx;
      if (k == 0)
	{
	  mask[(hW+2)*(2*hW+3)+(hW+2)] = 1;
	  shifts[idx][0] = 0;
	  shifts[idx][1] = 0;
	  idx = idx + 1;
	}
      else
	for (s = 0; s < 8*k; ++s)
	  {
	    cx = k * cosf(((float) s) / (8*k) * 2 * PI);
	    cy = k * sinf(((float) s) / (8*k) * 2 * PI);
	    x = roundf(cx);
	    y = roundf(cy);
	    if (!mask[(x+hW+2)*(2*hW+3)+(y+hW+2)] &&
		(k-0.5)*(k-0.5) < x*x+y*y && x*x+y*y <= (k+0.5)*(k+0.5))
	      {
		mask[(x+hW+2)*(2*hW+3)+(y+hW+2)] = 1;
		shifts[idx][0] = x;
		shifts[idx][1] = y;
		idx = idx + 1;
	      }
	    for (x = floorf(cx); x <= ceilf(cx); ++x)
	      for (y = floorf(cy); y <= ceilf(cy); ++y)
		{
		  if (!mask[(x+hW+2)*(2*hW+3)+(y+hW+2)] &&
		      (k-0.5)*(k-0.5) < x*x+y*y && x*x+y*y <= (k+0.5)*(k+0.5))
		    {
		      mask[(x+hW+2)*(2*hW+3)+(y+hW+2)] = 1;
		      shifts[idx][0] = x;
		      shifts[idx][1] = y;
		      idx = idx + 1;
		    }
		}
	  }
      lengths[k] = idx - lengths[k];
    }
  free(mask);
}
#endif //!SPIRALESQUARE

#define OFFSET(s, k, x, y)			((s)*(M*N*K)+(k)*(M*N)+(x)*N+(y))
#define DIFF_OFFSET(s)				((s) * (Mext * Next))

static sardata* nlsar_core(const sardata*		input,
			   float       			L,
			   int				hW,
			   int				hP,
			   int				hR,
			   int				K,
			   int				S,
			   sarsimstats** const		stats,
			   float			eta2,
			   int				verbose,
			   fltdata**       		outlook,
			   int*				gl_cpt,
			   int				nb_windows)
{
  sardata*	binput;
  sardata*	output;
  sardata**	stack;
  sarsimdata**	inputsarsim;
  fltdata*	look;
  int M = input->M;
  int N = input->N;
  int D = input->D;
  int Mext, Next;
  int offset;
  int W, k, l;
  int i, j;
  int cs;
  int x, y, dx, dy, x_dx, y_dy;
  int2* spirale_shifts;
  int*  spirale_lengths;
  float d, w;
  float span;
  float* sum_w;
  float* sum_w2;
  float* sum_wspan;
  float* sum_wspan2;
  float* diff;
  float cmean, cvar, alpha, alpha_new, clook;
  int percent, cpt;
  int hPmin = 1;
  int s;

  { // Allocations
    inputsarsim = malloc(S * sizeof(sarsimdata*));
    binput = sardata_dup(input);
    for (s = 0; s < S; ++s)
      {
#pragma omp critical
	{
	  if (s)
	    binput = sargausscar(input, binput, s);
	}
	if (!(inputsarsim[s] = sarsim_glrwishart.create(L, binput)))
	  return NULL;
      }
    sardata_free(binput);
    if (!(sum_w      = calloc(M * N * S * K, sizeof(float))) ||
	!(sum_w2     = calloc(M * N * S * K, sizeof(float))) ||
	!(sum_wspan  = calloc(M * N * S * K * D, sizeof(float))) ||
	!(sum_wspan2 = calloc(M * N * S * K * D, sizeof(float))))
      {
	sarerror_perror();
	return NULL;
      }
    if (!(stack = malloc(S*K * sizeof(sardata*))))
      {
	sarerror_perror();
	return NULL;
      }
    for (k = 0; k < S*K; ++k)
      if (!(stack[k] = sardata_calloc_size(M, N, D)))
	{
	  sarerror_perror();
	  return NULL;
	}
    output = sardata_calloc_size(M, N, D);
    look = fltdata_calloc_size(M, N, 1);

    Mext = M+2*hP+1;
    Next = N+2*hP+1;
    diff = malloc(S * Mext * Next * sizeof(float));
  }
  { // Core
    W = (2*hW+1)*(2*hW+1);
    spirale_shifts = malloc(W * sizeof(int2));
    spirale_lengths = malloc((hW+1) * sizeof(int));
    spirale(hW, spirale_shifts, spirale_lengths);
    W = 0;
    for (cs = 0; cs < hW + 1; ++cs)
      W += spirale_lengths[cs];

    // Central pixel
    w = 1;
    for (k = 0; k < K; ++k)
      for (s = 0; s < S; ++s)
	for (x = 0; x < M; ++x)
	  for (y = 0; y < N; ++y)
	    {
	      for (i = 0; i < D; ++i)
		for (j = 0; j < D; ++j)
		  SARDATA_ACCESS(stack[s*K+k], x, y, i, j) =
		    w * SARDATA_ACCESS(input, x, y, i, j);
	      offset = OFFSET(s, k, x, y);
	      sum_w[offset] = w;
	      sum_w2[offset] = w * w;
	      for (i = 0; i < D; ++i)
		{
		  span = SARDATA_ACCESS(input, x, y, i, i);
		  sum_wspan[offset * D + i] = w * span;
		  sum_wspan2[offset * D + i] = w * span * span;
		}
	      FLTDATA_ACCESS(look, x, y, 0) = 1;
	      for (i = 0; i < D; ++i)
		for (j = 0; j < D; ++j)
		  SARDATA_ACCESS(output, x, y, i, j) =
		    w * SARDATA_ACCESS(input, x, y, i, j);
	    }

    // Other pixels
    cpt = 0;
    for (cs = 1; cs < hW+1; ++cs)
      {
	for (l = 0; l < spirale_lengths[cs]; l=l+hR*hR)
	  {
#pragma omp atomic
	    (*gl_cpt)++;
	    ++cpt;
	    dx = spirale_shifts[cpt][0];
	    dy = spirale_shifts[cpt][1];
	    // Compute difference
	    for (s = 0; s < S; ++s)
	      {
		for (x = -hP-1; x < M+hP; ++x)
		  for (y = -hP-1; y < N+hP; ++y)
		    {
		      x_dx = MOD(x + dx, M);
		      y_dy = MOD(y + dy, N);
		      diff[DIFF_OFFSET(s) + (x+1+hP) * Next + (y+1+hP)] =
			sarsim_glrwishart.lsarsim(D, L,
						  SARSIMDATA_ACCESS(inputsarsim[s], MOD(x, M), MOD(y, N)),
						  SARSIMDATA_ACCESS(inputsarsim[s], x_dx, y_dy));
		    }
		// Compute commulative sums
		for (y = 1; y < Next; ++y)
		  diff[DIFF_OFFSET(s) + y] += diff[DIFF_OFFSET(s) + y-1];
		for (x = 1; x < Mext; ++x)
		  diff[DIFF_OFFSET(s) + x * Next] += diff[DIFF_OFFSET(s) + (x-1) * Next];
		for (x = 1; x < Mext; ++x)
		  for (y = 1; y < Next; ++y)
		    diff[DIFF_OFFSET(s) + x * Next + y] +=
		      + diff[DIFF_OFFSET(s) + (x-1) * Next + y]
		      + diff[DIFF_OFFSET(s) + x     * Next + (y-1)]
		      - diff[DIFF_OFFSET(s) + (x-1) * Next + (y-1)];

		// Perform K NL-means
		for (k = 0; k < K; ++k)
		  for (x = 0; x < M; ++x)
		    for (y = 0; y < N; ++y)
		      {
			// Candidate pixel
			x_dx = MOD(x + dx, M);
			y_dy = MOD(y + dy, N);

			// Extract weights for current patch size
			d =
			  + diff[DIFF_OFFSET(s) + (x+1+hP+k*hR+hPmin)   * Next + (y+1+hP+k*hR+hPmin)]
			  - diff[DIFF_OFFSET(s) + (x+1+hP+k*hR+hPmin)   * Next + (y+1+hP-k*hR-hPmin-1)]
			  - diff[DIFF_OFFSET(s) + (x+1+hP-k*hR-hPmin-1) * Next + (y+1+hP+k*hR+hPmin)]
			  + diff[DIFF_OFFSET(s) + (x+1+hP-k*hR-hPmin-1) * Next + (y+1+hP-k*hR-hPmin-1)];
			d /= (2*(k*hR+hPmin)+1) * (2*(k*hR+hPmin)+1);
			w = phi(d, stats[s*K+k]);

			// Accumulation
			offset = OFFSET(s, k, x, y);
			for (i = 0; i < D; ++i)
			  for (j = 0; j < D; ++j)
			    SARDATA_ACCESS(stack[s*K+k], x, y, i, j) +=
			      w * SARDATA_ACCESS(input, x_dx, y_dy, i, j);
			sum_w[offset] += w;
			sum_w2[offset] += w * w;
			for (i = 0; i < D; ++i)
			  {
			    span = SARDATA_ACCESS(input, x_dx, y_dy, i, i);
			    sum_wspan[offset * D + i] += w * span;
			    sum_wspan2[offset * D + i] += w * span * span;
			  }
		      }
	      }
	    if (verbose)
#pragma omp critical
	      {
		percent = (100 * *gl_cpt * hR * hR / W) / nb_windows;
		if (percent != (100 * (*gl_cpt-1) * hR * hR / W) / nb_windows)
		  sarwaitbar_update(percent);
	      }
	  }
	// Aggregations
	for (s = 0; s < S; ++s)
	  for (k = 0; k < K; ++k)
	    for (x = 0; x < M; ++x)
	      for (y = 0; y < N; ++y)
		{
		  offset = OFFSET(s, k, x, y);
		  // Statistics
		  alpha = 0;
		  for (i = 0; i < D; ++i)
		    {
		      cmean = sum_wspan[offset * D + i] / sum_w[offset];
		      cvar = sum_wspan2[offset * D + i] / sum_w[offset] - cmean * cmean;

		      //alpha_new = fabsf(cvar - cmean * cmean * eta2) / (1+eta2);
		      //alpha_new = alpha_new / cvar;
		      alpha_new = fabsf(cvar - cmean * cmean * eta2);
		      alpha_new = alpha_new / (alpha_new + cmean * cmean * eta2);

		      if (alpha_new > alpha)
			alpha = alpha_new;
		    }

		  // LLMMSE update
		  clook = 1. /
		    ((1 - alpha) * (1 - alpha) * (sum_w2[offset]-1) / (sum_w[offset] * sum_w[offset]) +
		     ((1 - alpha) / sum_w[offset] + alpha) * ((1 - alpha) / sum_w[offset] + alpha));
		  // Keep if equivalent number of look has increased
		  if (FLTDATA_ACCESS(look, x, y, 0) < clook)
		    {
		      FLTDATA_ACCESS(look, x, y, 0) = clook;
		      for (i = 0; i < D; ++i)
			for (j = 0; j < D; ++j)
			  SARDATA_ACCESS(output, x, y, i, j) =
			    (1-alpha) * SARDATA_ACCESS(stack[s*K+k], x, y, i, j) / sum_w[offset] +
			    alpha * SARDATA_ACCESS(input, x, y, i, j);
		    }
		}
      }
    free(diff);
    free(spirale_shifts);
    free(spirale_lengths);
  }
  for (k = 0; k < S*K; ++k)
    sardata_free(stack[k]);
  free(stack);
  if (outlook)
    *outlook = look;
  else
    fltdata_free(look);
  free(sum_w);
  free(sum_w2);
  free(sum_wspan);
  free(sum_wspan2);
  for (s = 0; s < S; ++s)
    sarsim_glrwishart.free(inputsarsim[s]);
  free(inputsarsim);
  return output;
}

static sardata* suppress_zero(sardata* data)
{
  int i, j, k;
  double min = INFINITY;
  double minp = INFINITY;
  double value;
  for (i = 0; i < data->M; ++i)
    for (j = 0; j < data->N; ++j)
      for (k = 0; k < data->D; ++k)
	{
	  value = cabsf(SARDATA_ACCESS(data, i, j, k, k));
	  if (value < min)
	    min = value;
	  if (0 < value && value < minp)
	    minp = value;
	}
  if (min == 0.0)
    for (i = 0; i < data->M; ++i)
      for (j = 0; j < data->N; ++j)
	for (k = 0; k < data->D; ++k)
	  if (cabsf(SARDATA_ACCESS(data, i, j, k, k)) == min)
	    SARDATA_ACCESS(data, i, j, k, k) = minp;
  return data;
}

static sardata* nlsar(const sardata*		input,
		      sardata*			output,
		      float			L,
		      int			hW,
		      int			hP,
		      const sardata*		nse,
		      int			verbose,
		      fltdata**       		outlook)
{
  int M = input->M;
  int N = input->N;
  int D = input->D;
  int TN = M;
  int TM = N;
  int ext = hP+hW+3;
  int x, y;
  int w, h;
  int x_ext, y_ext;
  int w_ext, h_ext;
  int dx, dy;
  int k, l, m;
  int i, j;
  sardata* bnse;
  fltdata* look;
  sardata* input_tmp;
  sardata* output_tmp = NULL;
  fltdata* look_tmp = NULL;
  sarsimstats** stats;
  int is_corr;
  float span, eta2, corr, sum_span, sum_span2;
  int hR, K, hPmin = 1;
  int S, s;
  int gl_cpt = 0;
  while (TM * TN * D > 1024*1024)
    {
      if (TM > TN)
	TM /= 2;
      else
	TN /= 2;
    }
#ifdef OMP
  int nb_procs;
  float s_nb_procs;
  nb_procs = omp_get_max_threads();
  s_nb_procs = floorf(sqrtf(nb_procs));
  TM = ceilf(MIN(TM, M) / s_nb_procs);
  TN = ceilf(MIN(TN, N) / (nb_procs / s_nb_procs));
#endif //!OMP
  TM = MAX(TM, 32);
  TN = MAX(TN, 32);

  { // Noise analysis
    if (verbose)
      sarprintf("Noise analysis\n");
    sum_span = 0;
    sum_span2 = 0;
    for (x = 0; x < nse->M; ++x)
      for (y = 0; y < nse->N; ++y)
	{
	  span = SARDATA_ACCESS(nse, x, y, 0, 0);
	  sum_span += span;
	  sum_span2 += span * span;
	}
    sum_span /= nse->M * nse->N;
    sum_span2 /= nse->M * nse->N;
    eta2 = sum_span2 - sum_span * sum_span;
    sum_span2 = 0;
    for (x = 0; x < nse->M-1; ++x)
      for (y = 0; y < nse->N; ++y)
	sum_span2 +=
	  SARDATA_ACCESS(nse, x, y, 0, 0) *
	  SARDATA_ACCESS(nse, x+1, y, 0, 0);
    sum_span2 /= nse->M * nse->N;
    corr = fabsf(sum_span2 - sum_span * sum_span) / eta2;
    eta2 = eta2 / sum_span / sum_span;
    if (verbose)
      sarprintf("\teta2=%.2f [theoretical=%.2f] corr=%.2f\n", eta2, 1./L, corr);
    is_corr = corr > 0.25;
  }
  if (is_corr)
    {
      hR = 2;
      S = 3;
      hW *= 2;
      hP *= 2;
    }
  else
    {
      hR = 1;
      S = 3;
    }
  if (verbose)
    sarprintf("\t=> #scales=%d step=%d\n", S, hR);
  if (fabsf(eta2 * L - 1) > 0.1)
    sarprintf("Empirical eta2 deviates more than 10%% w.r.t theoretical value.\n"
		      "\tIs the theoretical number of look correct?\n"
		      "\tDoes the image of noise contain only noise? No structures?\n"
		      "\tIs the image of noise big enough?");
  K = (hP-hPmin+1)/hR;
  stats = malloc(S * K * sizeof(sarsimstats**));

  { // Compute statistics of SARSIM
    if (verbose)
      sarprintf("Compute similarity statistics\n");
    bnse = sardata_dup(nse);
    for (s = 0; s < S; ++s)
      {
	if (s)
	  bnse = sargausscar(nse, bnse, s);
#pragma omp parallel default(shared) private(k)
	{
# pragma omp for schedule(dynamic) nowait
	  for (k = 0; k < K; ++k)
	    stats[s*K+k] = sarsimstats_create(bnse, &sarsim_glrwishart, NULL, L,
					      k*hR+hPmin, 0, 1, 1024);
	}
	if (verbose)
	  for (k = 0; k < K; k += MAX(1, K-1))
	    sarprintf("\ts=%d hP=%2d mean=%f std=%f\n",
		      s+1, k*hR+hPmin,
		      stats[s*K+k]->mean,
		      stats[s*K+k]->std);
      }
    sardata_free(bnse);
  }
  if (verbose)
#ifdef OMP
    sarprintf("Computation (#proc=%d, #windows=%d)\n", nb_procs, ((M-1)/TM+1)*((N-1)/TN+1));
#else
  sarprintf("Computation (#proc=%d, #windows=%d)\n", 1, ((M-1)/TM+1)*((N-1)/TN+1));
#endif
  look = fltdata_calloc_size(M, N, 1);
  output = sardata_realloc_size(output, M, N, D);

  if (verbose)
    sarwaitbar_open();
#pragma omp parallel default(shared) private(k,l,m,x,y,w,h,x_ext,y_ext,w_ext,h_ext,dx,dy,i,j,input_tmp,output_tmp,look_tmp)
  {
# pragma omp for schedule(dynamic) nowait
    for (m = 0; m < ((M-1)/TM+1) * ((N-1)/TN+1); ++m)
      {
	k = m / ((N-1)/TN+1);
	l = m % ((N-1)/TN+1);

	x = MAX(0, k * TM);
	y = MAX(0, l * TN);
	w = MIN(x + TM, M) - x;
	h = MIN(y + TN, N) - y;
	x_ext = MAX(0, x - ext);
	y_ext = MAX(0, y - ext);
	w_ext = MIN(x + w + ext, M) - x_ext;
	h_ext = MIN(y + h + ext, N) - y_ext;
	input_tmp = sardata_calloc_size(TM + 2 * ext, TN + 2 * ext, D);
	input_tmp = sardata_extract(input, input_tmp,
				    x_ext, y_ext,
				    w_ext, h_ext, 1);
	input_tmp = suppress_zero(input_tmp);
	look_tmp = NULL;
	output_tmp = nlsar_core(input_tmp, L, hW, hP, hR, K, S, stats, eta2, verbose, &look_tmp,
				&gl_cpt, ((M-1)/TM+1)*((N-1)/TN+1));
	for (dx = 0; dx < w; ++dx)
	  for (dy = 0; dy < h; ++dy)
	    {
	      FLTDATA_ACCESS(look, x + dx, y + dy, 0) =
		FLTDATA_ACCESS(look_tmp, x - x_ext + dx, y - y_ext + dy, 0);
	      for (i = 0; i < D; ++i)
		for (j = 0; j < D; ++j)
		  SARDATA_ACCESS(output, x + dx, y + dy, i, j) =
		    SARDATA_ACCESS(output_tmp, x - x_ext + dx, y - y_ext + dy, i, j);
	    }
	sardata_free(input_tmp);
	sardata_free(output_tmp);
	fltdata_free(look_tmp);
      }
  }
  if (verbose)
    sarwaitbar_close();
  for (k = 0; k < S*K; ++k)
    sarsimstats_free(stats[k]);
  free(stats);
  if (outlook)
    *outlook = look;
  else
    fltdata_free(look);
  return output;
}

sardata* sarnlsar(const sardata* input,
		  sardata* output,
		  float L,
		  int n_args,
		  ...)
{
  int hW = 12;
  int hP = 5;
  int verbose = 1;
  sardata* noise = NULL;
  fltdata** look = NULL;
  int i;
  va_list ap;

  va_start(ap, n_args);
  for (i = 0; i < n_args; i++)
    switch (i)
      {
	case 0:
	  verbose = va_arg(ap, int);
	  break;
	case 1:
	  hW = va_arg(ap, int);
	  break;
	case 2:
	  hP = va_arg(ap, int);
	  break;
	case 3:
	  noise = va_arg(ap, sardata*);
	  break;
	case 4:
	  look = va_arg(ap, fltdata**);
	  break;
	default:
	  sarerror_msg("Too many arguments");
	  return NULL;
    }
  va_end(ap);
  if (!noise)
    {
      if (L != (int) L)
	{
	  sarerror_msg("Cannot generate iid %2f non integer looks Wishart image", L);
	  return NULL;
	}
      srand(0); // Enter in deterministic mode
      if (!(noise = wishartrnd(256, 256, input->D, (int) L, 0.95)))
	{
	  sarerror_msg_msg("Cannot generate iid %d looks Wishart image", (int) L);
	  return NULL;
	}
      srand(time(NULL)); // Back to stochasitic mode
      output = nlsar(input, output,
		     L, hW, hP, noise,
		     verbose, look);
      sardata_free(noise);
    }
  else
    {
      noise = suppress_zero(noise);
      output = nlsar(input, output,
		     L, hW, hP, noise,
		     verbose, look);
    }
  if (!output)
    {
      sarerror_msg_msg("Cannot filter the image");
      return NULL;
    }

  return output;
}
