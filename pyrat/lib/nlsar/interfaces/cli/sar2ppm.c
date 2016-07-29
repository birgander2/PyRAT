/*
** sar2ppm.c: CLI for exporting SAR data files to PPM images
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
** Started on  Wed Jul 24 14:50:10 2013 Charles-Alban Deledalle
** Last update Fri Aug 23 17:09:35 2013 Charles-Alban Deledalle
*/

#include <stdio.h>
#include <stdlib.h>
#include "tools/sarerror.h"
#include "data/sardata.h"
#include "data/rgbdata.h"
#include "data/iosar.h"
#include "data/ionetpbm.h"

int usage(const char* argv0)
{
  fprintf(stderr, "usage: %s filein.rat fileout.pgm [alpha = 3, gamma = 0.7]\n", argv0);
  return 1;
}

int main(int argc, char* argv[])
{
  if (argc < 3)
    return usage(argv[0]);
  char* fn_in  = argv[1];
  char* fn_out = argv[2];
  float alpha, gamma;
  if (argc < 4)
    alpha = 3;
  else
    if (sscanf(argv[3], "%f", &alpha) != 1)
      return usage(argv[0]);
  if (argc < 5)
    gamma = 0.7;
  else
    if (sscanf(argv[4], "%f", &gamma) != 1)
      return usage(argv[0]);
  sardata* sar = sardata_alloc();
  if (!(sar = sarread(fn_in, sar)))
    {
      sarerror_msg_msg("Cannot open file %s", fn_in);
      fprintf(stderr, "%s\n", sarerror);
      return 2;
    }
  rgbdata* rgb = rgbdata_alloc();
  if (!(rgb = sar2rgb(sar, rgb, alpha, gamma)))
    {
      sarerror_msg_msg("Cannot create an RGB representation of %s", fn_in);
      fprintf(stderr, "%s\n", sarerror);
      return 2;
    }
  if (!(ppmwrite(rgb, fn_out)))
    {
      free(rgb);
      sarerror_msg_msg("Cannot create file %s", fn_out);
      fprintf(stderr, "%s\n", sarerror);
      return 2;
    }
  rgbdata_free(rgb);
  sardata_free(sar);
  return 0;
}

