/*
** sarnlsar.c: CLI for the NL-SAR filter
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
** Started on  Wed Jul 24 14:44:53 2013 Charles-Alban Deledalle
** Last update Tue Oct  8 11:14:34 2013 Charles-Alban Deledalle
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tools/sarprintf.h"
#include "tools/sarerror.h"
#include "data/sardata.h"
#include "data/fltdata.h"
#include "data/iosar.h"
#include "data/ionetpbm.h"
#include "algos/nlsar/nlsar.h"

static int usage(const char* argv0)
{
  fprintf(stderr, "usage: %s filein.rat fileout.rat L [verbose hW hP noise.rat]\n\n", argv0);
  fprintf(stderr, "\t L \t\t input number of look\n");
  fprintf(stderr, "\t verbose \t display progress bars (default 1)\n");
  fprintf(stderr, "\t hW \t\t radius of the largest search window (default 12)\n");
  fprintf(stderr, "\t hP \t\t half width of the largest patches (default 5)\n");
  fprintf(stderr, "\t noise.rat \t image of noise to learn statistics (default iid L looks Wishart)\n");
  return 1;
}

int main(int argc, char* argv[])
{
  if (argc < 4)
    return usage(argv[0]);
  char* fn_in = argv[1];
  char* fn_out = argv[2];
  char* fn_noise = NULL;
  float L;
  int hW = 12, hP = 5, verbose = 1;
  int argi = 3;
  sardata* input;
  sardata* output = NULL;
  sardata* sar_noise = NULL;

  if (sscanf(argv[argi++], "%f", &L) != 1)
    return usage(argv[0]);
  if (argc >= 5)
    if (sscanf(argv[argi++], "%d", &verbose) != 1)
      return usage(argv[0]);
  if (argc >= 6)
    if (sscanf(argv[argi++], "%d", &hW) != 1)
      return usage(argv[0]);
  if (argc >= 7)
    if (sscanf(argv[argi++], "%d", &hP) != 1)
      return usage(argv[0]);
  if (argc >= 8)
    fn_noise = argv[argi++];
  argi++;

  if (verbose)
    sarprintf("Load image\n");
  input = sardata_alloc();
  if (!(input = sarread(fn_in, input)))
    {
      sarerror_msg_msg("Cannot open file %s", fn_in);
      fprintf(stderr, "%s\n", sarerror);
      return 2;
    }
  if (fn_noise)
    {
      sar_noise = sardata_alloc();
      if (!(sar_noise = sarread(fn_noise, sar_noise)))
	{
	  sarerror_msg_msg("Cannot open file %s", fn_noise);
	  fprintf(stderr, "%s\n", sarerror);
	  return 2;
	}
    }
  output = sardata_alloc();
  if (!(output = sarnlsar(input, output, L, 4, verbose, hW, hP, sar_noise)))
    {
      fprintf(stderr, "%s\n", sarerror);
      return 2;
    }
  if (verbose)
    sarprintf("Save image\n");
  if (!(sarwrite(output, fn_out)))
    {
      sarerror_msg_msg("Cannot create file %s", fn_out);
      fprintf(stderr, "%s\n", sarerror);
      return 2;
    }
  sardata_free(input);
  sardata_free(sar_noise);
  sardata_free(output);
  return 0;
}
