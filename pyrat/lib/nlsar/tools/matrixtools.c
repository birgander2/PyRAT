/*
** matrixtools.c: implementation of matrix tools
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
** Started on  Wed Jul 24 16:07:53 2013 Charles-Alban Deledalle
** Last update Wed Jul 24 16:07:57 2013 Charles-Alban Deledalle
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include "mathtools.h"
#include "matrixtools.h"
#include "sarprintf.h"

#ifdef LAPACK
 int cpotrf_(const char *uplo, long int *n, float complex *a,
	     long int *lda, long int *info);
#endif
#ifdef BLAS
 int ctrsm_(char *side, char *uplo, char *transa, char *diag, 
	    long int *m, long int *n,
	    float complex *alpha, float complex *a, long int *lda,
	    float complex *b, long int *ldb);
#endif

float trace(int D, const float complex* C)
{
  float res = 0;
  int k;
  for (k = 0; k < D; ++k)
    res += cabsf(MATRIX_ACCESS(C, D, k, k));
  return res;
}

float det(int D, const float complex* C)
{
#ifdef LAPACK
  int k, l;
  float complex* CHO_C;
  long int info, longD = D;
  float res;
#endif //LAPACK
  switch (D)
    {
      case 1:
	return (MATRIX_ACCESS(C, D, 0, 0));
      case 2:
	return
	  + (MATRIX_ACCESS(C, D, 0, 0)) * (MATRIX_ACCESS(C, D, 1, 1))
	  - CABS2(MATRIX_ACCESS(C, D, 1, 0));
      case 3:
      	return
      	  + (MATRIX_ACCESS(C, D, 0, 0)) * cabsf(MATRIX_ACCESS(C, D, 1, 1)) * cabsf(MATRIX_ACCESS(C, D, 2, 2))
      	  - (MATRIX_ACCESS(C, D, 0, 0)) * CABS2(MATRIX_ACCESS(C, D, 1, 2))
      	  - (MATRIX_ACCESS(C, D, 1, 1)) * CABS2(MATRIX_ACCESS(C, D, 0, 2))
      	  - (MATRIX_ACCESS(C, D, 2, 2)) * CABS2(MATRIX_ACCESS(C, D, 0, 1))
      	  + 2 * crealf(MATRIX_ACCESS(C, D, 0, 1) * MATRIX_ACCESS(C, D, 1, 2) * MATRIX_ACCESS(C, D, 2, 0));
      default:
#ifdef LAPACK
	// det(C) = Prod(k) | L(k, k) |^2
	// with L the cholesky lower matrix of C
	CHO_C = malloc(D * D * sizeof(float complex));
	for (k = 0; k < D; ++k)
	  for (l = k; l < D; ++l)
	    CHO_C[k * D + l] = MATRIX_ACCESS(C, D, k, l);
	cpotrf_("L", &longD, CHO_C, &longD, &info);
	res = 1;
	for (k = 0; k < D; ++k)
	  res *= CABS2(CHO_C[k * D + k]);
	free(CHO_C);
	return res;
#else //LAPACK
	sarprintf_std_error("Error: Determinant for matrix of dimension %d required LAPACK to be installed.\n", D);
	sarprintf_std_error("       Please install CLAPACK, configure and recompile\n");
	exit(3);
	return 1;
#endif //LAPACK
    }
}

float complex trace_prod_invC1_C2(int D,
				  const float complex* C1,
				  const float complex* C2)
{
#if defined(BLAS) && defined(LAPACK)
  int k, l;
  float complex* CHO_C1;
  float complex* INV_L1_L2;
  long int info, longD = D;
  float res;
  float complex alpha = 1;
#endif //BLAS && LAPACK
  switch (D)
    {
      case 1:
	return crealf(MATRIX_ACCESS(C2, D, 0, 0)) / crealf(MATRIX_ACCESS(C1, D, 0, 0));
      case 2:
	return
	  (+ crealf(MATRIX_ACCESS(C1, D, 1, 1)) * crealf(MATRIX_ACCESS(C2, D, 0, 0))
	   - 2 * crealf(MATRIX_ACCESS(C1, D, 0, 1) * MATRIX_ACCESS(C2, D, 1, 0))
	   + crealf(MATRIX_ACCESS(C1, D, 0, 0)) * crealf(MATRIX_ACCESS(C2, D, 1, 1))) /
	  det(D, C1);
      case 3:
      	return
      	  ((crealf(MATRIX_ACCESS(C1, D, 1, 1)) * crealf(MATRIX_ACCESS(C1, D, 2, 2))
      	    - CABS2(MATRIX_ACCESS(C1, D, 1, 2)))
	   * crealf(MATRIX_ACCESS(C2, D, 0, 0))
      	   +
	   2 * crealf((MATRIX_ACCESS(C1, D, 0, 2) * MATRIX_ACCESS(C1, D, 2, 1)
      		       - MATRIX_ACCESS(C1, D, 0, 1) * MATRIX_ACCESS(C1, D, 2, 2))
      		      * MATRIX_ACCESS(C2, D, 1, 0))
      	   +
      	   2 * crealf((MATRIX_ACCESS(C1, D, 0, 1) * MATRIX_ACCESS(C1, D, 1, 2)
      		       - MATRIX_ACCESS(C1, D, 0, 2) * MATRIX_ACCESS(C1, D, 1, 1))
      		      * MATRIX_ACCESS(C2, D, 2, 0))
      	   +
      	   (crealf(MATRIX_ACCESS(C1, D, 0, 0)) * crealf(MATRIX_ACCESS(C1, D, 2, 2))
      	    - CABS2(MATRIX_ACCESS(C1, D, 0, 2)))
      	   * crealf(MATRIX_ACCESS(C2, D, 1, 1))
      	   +
      	   2 * crealf((MATRIX_ACCESS(C1, D, 0, 2) * MATRIX_ACCESS(C1, D, 1, 0)
      		       - MATRIX_ACCESS(C1, D, 0, 0) * MATRIX_ACCESS(C1, D, 1, 2))
      		      * MATRIX_ACCESS(C2, D, 2, 1))
      	   +
      	   (crealf(MATRIX_ACCESS(C1, D, 0, 0)) * crealf(MATRIX_ACCESS(C1, D, 1, 1))
      	    - CABS2(MATRIX_ACCESS(C1, D, 1, 0)))
      	   * crealf(MATRIX_ACCESS(C2, D, 2, 2))) /
      	  det(D, C1);
      default:
#if defined(BLAS) && defined(LAPACK)
	// tr(C1^-1 * C2) = Sum(k,l) | (L1^-1 * L2)(k,l) |^2
	// with L1 and L2 the cholesky lower matrix of C1 and C2
	CHO_C1 = malloc(D * D * sizeof(float complex));
	INV_L1_L2 = malloc(D * D * sizeof(float complex));
	for (k = 0; k < D; ++k)
	  for (l = k; l < D; ++l)
	    {
	      CHO_C1[k * D + l] = MATRIX_ACCESS(C1, D, k, l);
	      INV_L1_L2[k * D + l] = MATRIX_ACCESS(C2, D, k, l);
	    }
	cpotrf_("L", &longD, CHO_C1, &longD, &info);
	cpotrf_("L", &longD, INV_L1_L2, &longD, &info);
	ctrsm_("L", "L", "N", "N",
	       &longD, &longD, &alpha,
	       CHO_C1, &longD, INV_L1_L2, &longD);
	free(CHO_C1);
	res = 0;
	for (k = 0; k < D; ++k)
	  for (l = k; l < D; ++l)
	    res += CABS2(INV_L1_L2[k * D + l]);
	free(INV_L1_L2);
	return res;
#else //BLAS && LAPCK
	sarprintf_std_error("Calculus for matrix of dimension %d required BLAS and LAPACK to be installed.\n", D);
	sarprintf_std_error("Please install LAPACK, CBLAS, configure and recompile\n");
	exit(3);
	return 1;
#endif //BLAS && LAPCK
    }
}
