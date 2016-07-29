/*
** sardiskcar.c: CLI for the diskcar filter
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
** Started on  Wed Jul 24 14:48:39 2013 Charles-Alban Deledalle
** Last update Wed Jul 24 14:48:42 2013 Charles-Alban Deledalle
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tools/sarerror.h"
#include "data/sardata.h"
#include "data/iosar.h"
#include "algos/carfilter/carfilter.h"

static int usage(const char* argv0)
{
  fprintf(stderr, "usage: %s filein.rat fileout.rat hW\n\n", argv0);
  fprintf(stderr, "\t hW \t half width of the search window \t [=1]\n");
  return 1;
}

int main(int argc, char* argv[])
{
  if (argc < 3)
    return usage(argv[0]);
  char* fn_in = argv[1];
  char* fn_out = argv[2];
  int hW;

 if (argc < 4)
    hW = 1;
  else
    if (sscanf(argv[3], "%d", &hW) != 1)
     return usage(argv[0]);

  sardata* input = sardata_alloc();
  sardata* output  = sardata_alloc();
  if (!(input = sarread(fn_in, input)))
    {
      sarerror_msg_msg("Cannot open file %s", fn_in);
      fprintf(stderr, "%s\n", sarerror);
      return 2;
    }
  output = sardiskcar(input, output, hW);
  if (!(sarwrite(output, fn_out)))
    {
      sarerror_msg_msg("Cannot create file %s", fn_out);
      fprintf(stderr, "%s\n", sarerror);
      return 2;
    }
  sardata_free(input);
  sardata_free(output);
  return 0;
}
