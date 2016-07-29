/*
** interface.h: Interface to NL-SAR Toolbox
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
** Started on  Mon Aug 19 10:19:23 2013 Charles-Alban Deledalle
** Last update Fri Aug 23 17:09:02 2013 Charles-Alban Deledalle
*/

#ifndef NLSARTOOLBOX_H_
# define NLSARTOOLBOX_H_

# include <complex.h>

typedef struct
{
  int    M;
  int    N;
  int    D;
  float complex* array;
} sardata;

typedef struct
{
  int    M;
  int    N;
  unsigned char* array;
} rgbdata;

// Functions for output printing on standard input/ouput
int sarprintf_std(const char* format, ...);
int sarprintf_std_ret(const char* format, ...);
int sarprintf_std_warning(const char* format, ...);
int sarprintf_std_error(const char* format, ...);

// Wrappers for output printing on interface dependent input/ouput
extern int (*sarprintf)(const char* format, ...);
extern int (*sarprintf_ret)(const char* format, ...);
extern int (*sarprintf_warning)(const char* format, ...);
extern int (*sarprintf_error)(const char* format, ...);

// Wrappers for waitbar on interface dependent input/ouput
extern int (*sarwaitbar_open)(const char* format, ...);
extern int (*sarwaitbar_update)(const char* format, ...);
extern int (*sarwaitbar_close)(const char* format, ...);

// Functions for managing error messages
extern char sarerror[2048];

int sarerror_msg(const char* format, ...);
int sarerror_perror();
int sarerror_msg_msg(const char* format, ...);
int sarerror_msg_perror(const char* msg, ...);

// SAR Data Manipulation
sardata*	sardata_alloc();
sardata*	sardata_alloc_size(int M, int N, int D);
sardata*	sardata_calloc_size(int M, int N, int D);
sardata*	sardata_realloc_size(sardata* sar, int M, int N, int D);
sardata*	sardata_free(sardata* src);
sardata*	sardata_copy(const sardata* src, sardata* dst);
sardata*	sardata_dup(const sardata* src);
sardata*	sardata_extract(const sardata* src, sardata* dst,
				long int xoffset, long int yoffset,
				long int width, long int height,
				long int step);
sardata*	sardata_cat(sardata* const* src_list, sardata* dst,
			    int nx, int ny);
sardata*	sardata_join(sardata* const* src_list, sardata* dst,
			     int D);

// SAR I/O
sardata*	sarread(const char* filename, sardata* sar);
sardata*	sarread_extract(const char* filename, sardata* sar,
				long int xoffset, long int yoffset,
				long int width, long int height);
sardata*	sarread_header(const char* filename, sardata* sar);
int		sarwrite(const sardata* sar, const char* filename);

// SAR Filters
sardata*	sarboxcar(const sardata*	input, sardata* output, int hW);
sardata*	sardiskcar(const sardata*	input, sardata* output, int hW);
sardata*	sargausscar(const sardata*	input, sardata* output, int hW);

sardata*	sarnlsar(const sardata*		input, sardata* output,
			 float L, int n_args, ...
			 // int verbose = 1,
			 // int hW = 12,
			 // int hP = 5,
			 // sardata* sar_noise = iid L Wishart,
			 // fltdata** look = NULL,
			 );

// SAR -> RGB rendering
rgbdata*	sar2rgb(const sardata* sar, rgbdata* rgb, float alpha, float gamma);

// RGB Data manipulation
rgbdata*	rgbdata_alloc();
rgbdata*	rgbdata_alloc_size(int M, int N);
rgbdata*	rgbdata_calloc_size(int M, int N);
rgbdata*	rgbdata_realloc_size(rgbdata* rgb, int M, int N);
rgbdata*	rgbdata_free(rgbdata* src);
rgbdata*	rgbdata_copy(const rgbdata* src, rgbdata* dst);
rgbdata*	rgbdata_dup(const rgbdata* src);
rgbdata*	rgbdata_extract(const rgbdata* src, rgbdata* dst,
				int xoffset, int yoffset,
				int width, int height,
				int step);

// RGB I/O
int		ppmwrite(const rgbdata* rgb, const char* filename);

#endif //NLSARTOOLBOX_H_

