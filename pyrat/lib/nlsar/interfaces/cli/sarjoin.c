/*
** sarjoin.c: CLI to join mono-dimensional complex images
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
** Started on  Wed Jul 24 14:45:12 2013 Charles-Alban Deledalle
** Last update Wed Jul 24 14:47:22 2013 Charles-Alban Deledalle
*/

#include <stdio.h>
#include <stdlib.h>
#include "tools/sarerror.h"
#include "data/sardata.h"
#include "data/iosar.h"

int main(int argc, char* argv[])
{
  if (argc < 3)
    {
      fprintf(stderr, "usage: %s filein1 [filein2 ... fileinN] fileout\n", argv[0]);
      return 1;
    }
  int k;
  char* fn_in;
  int D = argc - 2;
  char* fn_out = argv[argc - 1];

  sardata** sar_in = malloc(D * sizeof(sardata *));
  if (!sar_in)
    {
      sarerror_msg_perror("Cannot allocate memory");
      fprintf(stderr, "%s\n", sarerror);
      return 2;
    }
  for (k = 0; k < D; ++k)
    {
      fn_in = argv[1 + k];
      sar_in[k] = sardata_alloc();
      if (!(sar_in[k] = sarread(fn_in, sar_in[k])))
	{
	  sarerror_msg_msg("Cannot open file %s", fn_in);
	  fprintf(stderr, "%s\n", sarerror);
	  return 2;
	}
    }
  sardata* sar_out  = sardata_alloc();
  if (!sar_out)
    {
      sarerror_msg_msg("Cannot allocate memory");
      fprintf(stderr, "%s\n", sarerror);
      return 2;
    }
  if (!(sar_out = sardata_join(sar_in, sar_out, D)))
    {
      sarerror_msg_msg("Cannot join images");
      fprintf(stderr, "%s\n", sarerror);
      return 2;
    }
  if (!(sarwrite(sar_out, fn_out)))
    {
      sarerror_msg_msg("Cannot create file %s", fn_out);
      fprintf(stderr, "%s\n", sarerror);
      return 2;
    }
  for (k = 0; k < D; ++k)
    sardata_free(sar_in[k]);
  free(sar_in);
  sardata_free(sar_out);
  return 0;
}
