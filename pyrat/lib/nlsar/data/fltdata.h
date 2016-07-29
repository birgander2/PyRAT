/*
** fltdata.h: declarations for manipulation of float data
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
** Started on  Wed Jul 24 15:50:58 2013 Charles-Alban Deledalle
** Last update Wed Jul 24 15:51:02 2013 Charles-Alban Deledalle
*/

#ifndef FLTDATA_H_
# define FLTDATA_H_

# define FLTDATA_ACCESS(flt, i, j, k) \
  ((flt)->array[((i) * (flt)->N + (j)) * (flt)->D + (k)])

typedef struct
{
  int    M;
  int    N;
  int    D;
  float* array;
} fltdata;

fltdata*	fltdata_alloc();
fltdata*	fltdata_alloc_size(int M, int N, int D);
fltdata*	fltdata_calloc_size(int M, int N, int D);
fltdata*	fltdata_realloc_size(fltdata* flt, int M, int N, int D);
fltdata*	fltdata_free(fltdata* src);
fltdata*	fltdata_copy(const fltdata* src, fltdata* dst);
fltdata*	fltdata_dup(const fltdata* src);
fltdata*	fltdata_extract(const fltdata* src, fltdata* dst,
				long int xoffset, long int yoffset,
				long int width, long int height,
				long int step);

#endif
