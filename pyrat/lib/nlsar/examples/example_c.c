/*
** example_c.c: Example of a C program using NL-SAR Toolbox
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
** Started on  Fri Aug 23 18:11:47 2013 Charles-Alban Deledalle
** Last update Fri Aug 23 18:11:49 2013 Charles-Alban Deledalle
*/

#include <nlsartoolbox.h>
#include <stdlib.h>
#include <stdio.h>
#include <nlsartoolbox.h>

int (*sarprintf)(const char* format, ...) = &sarprintf_std;
int (*sarprintf_ret)(const char* format, ...) = &sarprintf_std_ret;
int (*sarprintf_warning)(const char* format, ...) = &sarprintf_std_warning;
int (*sarprintf_error)(const char* format, ...) = &sarprintf_std_error;

int (*sarwaitbar_open)(void) = &sarwaitbar_std_open;
int (*sarwaitbar_update)(int) = &sarwaitbar_std_update;
int (*sarwaitbar_close)(void) = &sarwaitbar_std_close;

int main()
{
  sardata* sarimage = sardata_alloc();
  sardata* sarimage_box = sardata_alloc();
  sardata* sarimage_disk = sardata_alloc();
  sardata* sarimage_gauss = sardata_alloc();
  sardata* sarimage_nl = sardata_alloc();

  if (!(sarimage = sarread_header("example.rat", sarimage)))
    {
      sarerror_msg_msg("Cannot open file %s", "example.rat");
      sarprintf_error("%s\n", sarerror);
      sarprintf_error("Please run from command line:\n");
      sarprintf_error("       sarmire example.rat 256 256 3 3\n");
      return EXIT_FAILURE;
    }
  printf("M=%d N=%d D=%d &data=%p\n", sarimage->M, sarimage->N, sarimage->D, sarimage->array);

  if (!(sarimage = sarread("example.rat", sarimage)))
    {
      sarerror_msg_msg("Cannot open file %s", "example.rat");
      sarprintf_error("%s\n", sarerror);
      return EXIT_FAILURE;
    }
  printf("M=%d N=%d D=%d &data=%p\n", sarimage->M, sarimage->N, sarimage->D, sarimage->array);

  sarimage_box = sarboxcar(sarimage, sarimage_box, 1);
  sarimage_disk = sardiskcar(sarimage, sarimage_disk, 1);
  sarimage_gauss = sargausscar(sarimage, sarimage_gauss, 1);
  sarimage_nl = sarnlsar(sarimage, sarimage_nl, 3, 0);

  sarwrite(sarimage_box,   "example_box.rat");
  sarwrite(sarimage_disk,  "example_disk.rat");
  sarwrite(sarimage_gauss, "example_gauss.rat");
  sarwrite(sarimage_nl,    "example_nl.rat");

  sardata_free(sarimage);
  sardata_free(sarimage_box);
  sardata_free(sarimage_disk);
  sardata_free(sarimage_gauss);
  sardata_free(sarimage_nl);

  return EXIT_SUCCESS;
}
