/*
** sarprintf.c: wrapper for print outputs in Python interface
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
** Started on  Wed Jul 24 14:44:11 2013 Charles-Alban Deledalle
** Last update Fri Aug 23 18:10:07 2013 Charles-Alban Deledalle
*/

#include <Python.h>
#include <stdio.h>
#include <stdarg.h>
#include "tools/sarprintf.h"


int sarprintf(const char* format, ...)
{
  va_list args;
  va_start(args, format);
  int res = sarvprintf_std(format, args);
  va_end(args);
  return res;
}

int sarprintf_ret(const char* format, ...)
{
  va_list args;
  va_start(args, format);
  int res = sarvprintf_std_ret(format, args);
  va_end(args);
  return res;
}

int sarprintf_warning(const char* format, ...)
{
  char tmp[2048] = { 0, };
  va_list args;
  va_start(args, format);
  int res = vsprintf(tmp, format, args);
  PyErr_SetString(PyExc_Warning, tmp);
  va_end(args);
  return res;
}

int sarprintf_error(const char* format, ...)
{
  char tmp[2048] = { 0, };
  va_list args;
  va_start(args, format);
  int res = vsprintf(tmp, format, args);
  PyErr_SetString(PyExc_Exception, tmp);
  va_end(args);
  return res;
}
