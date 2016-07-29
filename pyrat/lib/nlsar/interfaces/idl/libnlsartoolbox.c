/*
*   libsarnlsar.c: IDL interface to NL-SAR Toolbox
*
*   This file is part of NL-SAR Toolbox version 0.6.
*
*   Copyright Charles-Alban Deledalle (2013)
*   Email charles-alban.deledalle@math.u-bordeaux1.fr
*
*   This software is a computer program whose purpose is to provide a
*   suite of tools to manipulate SAR images.
*
*   This software is governed by the CeCILL license under French law and
*   abiding by the rules of distribution of free software. You can use,
*   modify and/ or redistribute the software under the terms of the CeCILL
*   license as circulated by CEA, CNRS and INRIA at the following URL
*   "http://www.cecill.info".
*
*   As a counterpart to the access to the source code and rights to copy,
*   modify and redistribute granted by the license, users are provided only
*   with a limited warranty and the software's author, the holder of the
*   economic rights, and the successive licensors have only limited
*   liability.
*
*   In this respect, the user's attention is drawn to the risks associated
*   with loading, using, modifying and/or developing or reproducing the
*   software by the user in light of its specific status of free software,
*   that may mean that it is complicated to manipulate, and that also
*   therefore means that it is reserved for developers and experienced
*   professionals having in-depth computer knowledge. Users are therefore
*   encouraged to load and test the software's suitability as regards their
*   requirements in conditions enabling the security of their systems and/or
*   data to be ensured and, more generally, to use and operate it in the
*   same conditions as regards security.
*
*   The fact that you are presently reading this means that you have had
*   knowledge of the CeCILL license and that you accept its terms.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <idl_export.h>
#include "tools/sarerror.h"
#include "data/sardata.h"
#include "data/fltdata.h"
#include "data/iosar.h"
#include "algos/nlsar/nlsar.h"
#include "algos/carfilter/carfilter.h"


IDL_MSG_BLOCK msg_block;

IDL_VPTR IDL_sarinfo(int argc, IDL_VPTR* argv)
{
  argc = argc;

  IDL_ENSURE_STRING(argv[0]);
  const char* fn = IDL_STRING_STR(&(argv[0]->value.str));

  IDL_VPTR res;
  IDL_MEMINT dims[1] = { 3 };
  IDL_INT* res_data = (IDL_INT*) IDL_MakeTempArray(IDL_TYP_INT,
						   1, dims,
						   IDL_BARR_INI_ZERO, &res);
  sardata* output = sardata_alloc();
  if (!(output = sarread_header(fn, output)))
    {
      sarerror_msg_msg("SARINFO: Cannot open file %s", fn);
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }
  if (output)
    {
      res_data[0] = output->M;
      res_data[1] = output->N;
      res_data[2] = output->D;
    }
  else
    {
      res_data[0] = -1;
      res_data[1] = -1;
      res_data[2] = -1;
    }
  free(output);
  return res;
}

IDL_VPTR IDL_sarread(int argc, IDL_VPTR* argv)
{
  argc = argc;

  IDL_ENSURE_STRING(argv[0]);
  const char* fn = IDL_STRING_STR(&(argv[0]->value.str));

  sardata* output = sardata_alloc();
  if (!(output = sarread(fn, output)))
    {
      sarerror_msg_msg("SARREAD: Cannot open file %s", fn);
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }
  IDL_VPTR res;
  IDL_MEMINT dims[4];
  dims[0] = output->D;
  dims[1] = output->D;
  dims[2] = output->N;
  dims[3] = output->M;
  float complex* res_data = (float complex*) IDL_MakeTempArray(IDL_TYP_COMPLEX,
							       4, dims,
							       IDL_BARR_INI_ZERO, &res);
  res_data = memcpy(res_data, output->array, output->D*output->D*output->N*output->M * sizeof(float complex));
  free(output);
  return res;
}

void IDL_sarwrite(int argc, IDL_VPTR* argv)
{
  sardata* input;

  argc = argc;

  IDL_ENSURE_ARRAY(argv[0]);
  if (argv[0]->type != IDL_TYP_COMPLEX)
    {
      sarerror_msg("SARWRITE: First argument is expected to be an array of float complex");
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return;
    }
  if (argv[0]->value.arr->n_dim != 4)
    {
      sarerror_msg("SARWRITE: First argument should have 4 dimensions");
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return;
    }
  input = sardata_alloc();
  input->D = argv[0]->value.arr->dim[0];
  if (argv[0]->value.arr->dim[1] != input->D)
    {
      sarerror_msg("SARWRITE: First argument is expected to be an array with equal dimenions 1 and 2");
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return;
    }
  input->N = argv[0]->value.arr->dim[2];
  input->M = argv[0]->value.arr->dim[3];
  input->array = (float complex*) argv[0]->value.arr->data;

  IDL_ENSURE_STRING(argv[1]);
  const char* fn = IDL_STRING_STR(&(argv[1]->value.str));

  if (!(sarwrite(input, fn)))
    {
      sarerror_msg_msg("SARWRITE: Cannot create file %s", fn);
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return;
    }
  free(input);
}


IDL_VPTR IDL_sar2rgb(int argc, IDL_VPTR* argv)
{
  IDL_VPTR res;
  float alpha = 3;
  float gamma = 0.7;
  rgbdata* rgb;

  // Input

  IDL_ENSURE_ARRAY(argv[0]);
  if (argv[0]->type != IDL_TYP_COMPLEX)
    {
      sarerror_msg("SAR2RGB: First argument is expected to be an array of float complex");
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }
  if (argv[0]->value.arr->n_dim != 4)
    {
      sarerror_msg("SAR2RGB: First argument should have 4 dimensions");
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }
  sardata* input = sardata_alloc();
  input->D = argv[0]->value.arr->dim[0];
  if (argv[0]->value.arr->dim[1] != input->D)
    {
      sarerror_msg("SAR2RGB: First argument is expected to be an array with equal dimenions 1 and 2");
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }
  input->N = argv[0]->value.arr->dim[2];
  input->M = argv[0]->value.arr->dim[3];
  input->array = (float complex*) argv[0]->value.arr->data;
  if (argc >= 2)
    {
      IDL_ENSURE_SCALAR(argv[1]);
      alpha = IDL_DoubleScalar(argv[1]);
    }
  if (argc >= 3)
    {
      IDL_ENSURE_SCALAR(argv[2]);
      gamma = IDL_DoubleScalar(argv[2]);
    }
  // Processing
  rgb = rgbdata_alloc();
  if (!(rgb = sar2rgb(input, rgb, alpha, gamma)))
    {
      sarerror_msg_msg("SAR2RGB: Cannot create an RGB representation");
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }

  // Output
  IDL_MEMINT dims[3];
  dims[0] = 3;
  dims[1] = input->N;
  dims[2] = input->M;
  char* res_data = (char*) IDL_MakeTempArray(IDL_TYP_BYTE,
					     3, dims,
					     IDL_BARR_INI_ZERO, &res);
  res_data = memcpy(res_data, rgb->array, 3*input->N*input->M * sizeof(char));
  rgbdata_free(rgb);
  free(input);
  return res;
}

IDL_VPTR IDL_sarnlsar(int argc, IDL_VPTR* argv)
{
  IDL_VPTR res;
  sardata* sar_noise = NULL;
  float L;
  int hW = 12, hP = 5, verbose = 1;
  fltdata* look = NULL;

  // Input
  IDL_ENSURE_ARRAY(argv[0]);
  if (argv[0]->type != IDL_TYP_COMPLEX)
    {
      sarerror_msg("SARNLSAR: First argument is expected to be an array of float complex");
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }
  if (argv[0]->value.arr->n_dim != 4)
    {
      sarerror_msg("SARNLSAR: First argument should have 4 dimensions");
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }
  sardata* input = sardata_alloc();
  input->D = argv[0]->value.arr->dim[0];
  if (argv[0]->value.arr->dim[1] != input->D)
    {
      sarerror_msg("SARNLSAR: First argument is expected to be an array with equal dimenions 1 and 2");
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }
  input->N = argv[0]->value.arr->dim[2];
  input->M = argv[0]->value.arr->dim[3];
  input->array = (float complex*) argv[0]->value.arr->data;

  IDL_ENSURE_SCALAR(argv[1]);
  L = IDL_DoubleScalar(argv[1]);
  if (argc >= 3)
    {
      IDL_ENSURE_SCALAR(argv[2]);
      verbose = IDL_LongScalar(argv[2]);
    }
  if (argc >= 4)
    {
      IDL_ENSURE_SCALAR(argv[3]);
      hW = IDL_LongScalar(argv[3]);
    }
  if (argc >= 5)
    {
      IDL_ENSURE_SCALAR(argv[4]);
      hP = IDL_LongScalar(argv[4]);
    }
  if (argc >= 6)
    {
      IDL_ENSURE_ARRAY(argv[5]);
      if (argv[5]->type != IDL_TYP_COMPLEX)
	{
	  sarerror_msg("SARNLSAR: Sixth argument is expected to be an array of float complex");
	  IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
	  return IDL_GettmpInt(EXIT_FAILURE);
	}
      if (argv[5]->value.arr->n_dim != 4)
	{
	  sarerror_msg("SARNLSAR: Sixth argument should have 4 dimensions");
	  IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
	  return IDL_GettmpInt(EXIT_FAILURE);
	}
      sar_noise = sardata_alloc();
      sar_noise->D = argv[5]->value.arr->dim[5];
      if (argv[5]->value.arr->dim[1] != sar_noise->D)
	{
	  sarerror_msg("SARNLSAR: Sixth argument is expected to be an array with equal dimenions 1 and 2");
	  IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
	  return IDL_GettmpInt(EXIT_FAILURE);
	}
      sar_noise->N = argv[5]->value.arr->dim[2];
      sar_noise->M = argv[5]->value.arr->dim[3];
      sar_noise->array = (float complex*) argv[5]->value.arr->data;
    }

  // Processing
  sardata* output = sardata_alloc();
  if (!(output = sarnlsar(input, output, L, 5, verbose, hW, hP, sar_noise, &look)))
    {
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }

  // Output
  IDL_MEMINT dims[4];
  dims[0] = input->D;
  dims[1] = input->D;
  dims[2] = input->N;
  dims[3] = input->M;
  float complex* res_data = (float complex*) IDL_MakeTempArray(IDL_TYP_COMPLEX,
							       4, dims,
							       IDL_BARR_INI_ZERO, &res);
  res_data = memcpy(res_data, output->array, input->D*input->D*input->N*input->M * sizeof(float complex));

  if (sar_noise)
    free(sar_noise);
  if (look)
    fltdata_free(look);
  sardata_free(output);
  free(input);
  return res;
}

IDL_VPTR IDL_sarboxcar(int argc, IDL_VPTR* argv)
{
  IDL_VPTR res;
  int hW = 1;

  // Input
  IDL_ENSURE_ARRAY(argv[0]);
  if (argv[0]->type != IDL_TYP_COMPLEX)
    {
      sarerror_msg("SARBOXCAR: First argument is expected to be an array of float complex");
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }
  if (argv[0]->value.arr->n_dim != 4)
    {
      sarerror_msg("SARBOXCAR: First argument should have 4 dimensions");
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }
  sardata* input = sardata_alloc();
  input->D = argv[0]->value.arr->dim[0];
  if (argv[0]->value.arr->dim[1] != input->D)
    {
      sarerror_msg("SARBOXCAR: First argument is expected to be an array with equal dimenions 1 and 2");
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }
  input->N = argv[0]->value.arr->dim[2];
  input->M = argv[0]->value.arr->dim[3];
  input->array = (float complex*) argv[0]->value.arr->data;

  if (argc >= 2)
    {
      IDL_ENSURE_SCALAR(argv[1]);
      hW = IDL_LongScalar(argv[1]);
    }

  // Processing
  sardata* output = sardata_alloc();
  if (!(output = sarboxcar(input, output, hW)))
    {
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }

  // Output
  IDL_MEMINT dims[4];
  dims[0] = input->D;
  dims[1] = input->D;
  dims[2] = input->N;
  dims[3] = input->M;
  float complex* res_data = (float complex*) IDL_MakeTempArray(IDL_TYP_COMPLEX,
							       4, dims,
							       IDL_BARR_INI_ZERO, &res);
  res_data = memcpy(res_data, output->array, input->D*input->D*input->N*input->M * sizeof(float complex));

  sardata_free(output);
  free(input);
  return res;
}

IDL_VPTR IDL_sardiskcar(int argc, IDL_VPTR* argv)
{
  IDL_VPTR res;
  int hW = 1;

  // Input
  IDL_ENSURE_ARRAY(argv[0]);
  if (argv[0]->type != IDL_TYP_COMPLEX)
    {
      sarerror_msg("SARDISKCAR: First argument is expected to be an array of float complex");
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }
  if (argv[0]->value.arr->n_dim != 4)
    {
      sarerror_msg("SARDISKCAR: First argument should have 4 dimensions");
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }
  sardata* input = sardata_alloc();
  input->D = argv[0]->value.arr->dim[0];
  if (argv[0]->value.arr->dim[1] != input->D)
    {
      sarerror_msg("SARDISKCAR: First argument is expected to be an array with equal dimenions 1 and 2");
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }
  input->N = argv[0]->value.arr->dim[2];
  input->M = argv[0]->value.arr->dim[3];
  input->array = (float complex*) argv[0]->value.arr->data;

  if (argc >= 2)
    {
      IDL_ENSURE_SCALAR(argv[1]);
      hW = IDL_LongScalar(argv[1]);
    }

  // Processing
  sardata* output = sardata_alloc();
  if (!(output = sardiskcar(input, output, hW)))
    {
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }

  // Output
  IDL_MEMINT dims[4];
  dims[0] = input->D;
  dims[1] = input->D;
  dims[2] = input->N;
  dims[3] = input->M;
  float complex* res_data = (float complex*) IDL_MakeTempArray(IDL_TYP_COMPLEX,
							       4, dims,
							       IDL_BARR_INI_ZERO, &res);
  res_data = memcpy(res_data, output->array, input->D*input->D*input->N*input->M * sizeof(float complex));

  sardata_free(output);
  free(input);
  return res;
}

IDL_VPTR IDL_sargausscar(int argc, IDL_VPTR* argv)
{
  IDL_VPTR res;
  int hW = 1;

  // Input
  IDL_ENSURE_ARRAY(argv[0]);
  if (argv[0]->type != IDL_TYP_COMPLEX)
    {
      sarerror_msg("SARGAUSSCAR: First argument is expected to be an array of float complex");
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }
  if (argv[0]->value.arr->n_dim != 4)
    {
      sarerror_msg("SARGAUSSCAR: First argument should have 4 dimensions");
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }
  sardata* input = sardata_alloc();
  input->D = argv[0]->value.arr->dim[0];
  if (argv[0]->value.arr->dim[1] != input->D)
    {
      sarerror_msg("SARGAUSSCAR: First argument is expected to be an array with equal dimenions 1 and 2");
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }
  input->N = argv[0]->value.arr->dim[2];
  input->M = argv[0]->value.arr->dim[3];
  input->array = (float complex*) argv[0]->value.arr->data;

  if (argc >= 2)
    {
      IDL_ENSURE_SCALAR(argv[1]);
      hW = IDL_LongScalar(argv[1]);
    }

  // Processing
  sardata* output = sardata_alloc();
  if (!(output = sargausscar(input, output, hW)))
    {
      IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, sarerror);
      return IDL_GettmpInt(EXIT_FAILURE);
    }

  // Output
  IDL_MEMINT dims[4];
  dims[0] = input->D;
  dims[1] = input->D;
  dims[2] = input->N;
  dims[3] = input->M;
  float complex* res_data = (float complex*) IDL_MakeTempArray(IDL_TYP_COMPLEX,
							       4, dims,
							       IDL_BARR_INI_ZERO, &res);
  res_data = memcpy(res_data, output->array, input->D*input->D*input->N*input->M * sizeof(float complex));

  sardata_free(output);
  free(input);
  return res;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-braces"
int IDL_Load(void)
{
  static IDL_SYSFUN_DEF2 sarinfo_addr[] = {
    {(IDL_FUN_RET) IDL_sarinfo, (char *)"SARINFO", 1, 1, 0, 0}, };
  static IDL_SYSFUN_DEF2 sarread_addr[] = {
    {(IDL_FUN_RET) IDL_sarread, (char *)"SARREAD", 1, 1, 0, 0}, };
  static IDL_SYSFUN_DEF2 sarwrite_addr[] = {
    {(IDL_SYSRTN_GENERIC) IDL_sarwrite, (char *)"SARWRITE", 2, 2, 0, 0}, };
  static IDL_SYSFUN_DEF2 sar2rgb_addr[] = {
    {(IDL_FUN_RET) IDL_sar2rgb, (char *)"SAR2RGB", 1, 3, 0, 0}, };
  static IDL_SYSFUN_DEF2 sarnlsar_addr[] = {
    {(IDL_FUN_RET) IDL_sarnlsar, (char *)"SARNLSAR", 2, 6, 0, 0}, };
  static IDL_SYSFUN_DEF2 sarboxcar_addr[] = {
    {(IDL_FUN_RET) IDL_sarboxcar, (char *)"SARBOXCAR", 1, 2, 0, 0}, };
  static IDL_SYSFUN_DEF2 sardiskcar_addr[] = {
    {(IDL_FUN_RET) IDL_sardiskcar, (char *)"SARDISKCAR", 1, 2, 0, 0}, };
  static IDL_SYSFUN_DEF2 sargausscar_addr[] = {
    {(IDL_FUN_RET) IDL_sargausscar, (char *)"SARGAUSSCAR", 1, 2, 0, 0}, };

  return
    IDL_SysRtnAdd(sarinfo_addr, TRUE, IDL_CARRAY_ELTS(sarinfo_addr)) &
    IDL_SysRtnAdd(sarread_addr, TRUE, IDL_CARRAY_ELTS(sarread_addr)) &
    IDL_SysRtnAdd(sarwrite_addr, FALSE, IDL_CARRAY_ELTS(sarwrite_addr)) &
    IDL_SysRtnAdd(sar2rgb_addr, TRUE, IDL_CARRAY_ELTS(sar2rgb_addr)) &
    IDL_SysRtnAdd(sarnlsar_addr, TRUE, IDL_CARRAY_ELTS(sarnlsar_addr)) &
    IDL_SysRtnAdd(sarboxcar_addr, TRUE, IDL_CARRAY_ELTS(sarboxcar_addr)) &
    IDL_SysRtnAdd(sardiskcar_addr, TRUE, IDL_CARRAY_ELTS(sardiskcar_addr)) &
    IDL_SysRtnAdd(sargausscar_addr, TRUE, IDL_CARRAY_ELTS(sargausscar_addr));
}
#pragma GCC diagnostic pop
