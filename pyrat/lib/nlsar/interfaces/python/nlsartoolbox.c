/*
** nlsartoolbox.c: Python interface of NL-SAR Toolbox
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
** Started on  Fri Aug 23 18:10:47 2013 Charles-Alban Deledalle
** Last update Fri Aug 23 18:10:49 2013 Charles-Alban Deledalle
*/

/*
#include "/opt/anaconda/include/python3.4m/Python.h"
#include "/opt/anaconda/lib/python3.4/site-packages/numpy/core/include/numpy/arrayobject.h"
*/
#include <Python.h>
#include <numpy/arrayobject.h>
#include "data/sardata.h"
#include "data/iosar.h"
#include "tools/sarerror.h"
#include "algos/nlsar/nlsar.h"
#include "algos/carfilter/carfilter.h"

// To avoid a segfault and a warning
#define PYINIT {self = self; import_array1(NULL);}
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))

struct module_state {
    PyObject *error;
};



static PyObject* py_sarread(PyObject* self, PyObject* args)
{
  char* fn;
  sardata* output;
  PyObject* py_out;
  npy_intp dims[4];
  int ndims;

  PYINIT;

  // Input
  if (!PyArg_ParseTuple(args, "s", &fn))
    return NULL;

  // Processing
  output = sardata_alloc();
  if (!(output = sarread(fn, output)))
    {
      sarerror_msg_msg("Cannot open file %s", fn);
      PyErr_SetString(PyExc_IOError, sarerror);
      return NULL;
    }

  // Output
  ndims = 4;
  dims[0] = output->M;
  dims[1] = output->N;
  dims[2] = output->D;
  dims[3] = output->D;
  py_out = PyArray_SimpleNewFromData(ndims, dims,
				     NPY_CFLOAT,
				     output->array);
  PyArray_ENABLEFLAGS((PyArrayObject *) py_out, NPY_OWNDATA);

  free(output);

  return Py_BuildValue("O", py_out);
}

static sardata* PyArrayObject2sardata(PyArrayObject* py_in)
{
  sardata* input = NULL;
  int D;
  int i, j, k, l;

  if (!((py_in->nd == 2) || (py_in->nd == 4 && py_in->dimensions[2] == py_in->dimensions[3])))
    {
      PyErr_SetString(PyExc_ValueError, "Unexpected dimensions");
      return NULL;
    }
  if (py_in->nd == 4)
      D = py_in->dimensions[2];
  else
      D = 1;
  input = sardata_alloc_size(py_in->dimensions[0], py_in->dimensions[1], D);
  for (i = 0; i < input->M; ++i)
    for (j = 0; j < input->N; ++j)
      for (k = 0; k < input->D; ++k)
	for (l = 0; l < input->D; ++l)
	  SARDATA_ACCESS(input, i, j, k, l) =
	    *((float complex *)((char *) py_in->data +
				py_in->strides[0] * i +
				py_in->strides[1] * j +
				py_in->strides[2] * k +
				py_in->strides[3] * l));
  return input;
}

static PyObject* py_sarwrite(PyObject* self, PyObject* args)
{
  PyArrayObject* py_in;
  char* fn;
  sardata* input = NULL;

  PYINIT;

  // Input
  if (!PyArg_ParseTuple(args, "Os", (PyObject *) &py_in, &fn))
    return NULL;
  if (py_in == NULL)
    return NULL;
  if (!(input = PyArrayObject2sardata(py_in)))
    return NULL;

  // Processing
  if (!(sarwrite(input, fn)))
    {
      sarerror_msg_msg("Cannot create file %s", fn);
      PyErr_SetString(PyExc_IOError, sarerror);
      return NULL;
    }
  sardata_free(input);

  // Output
  return Py_True;
}

static PyObject* py_sarinfo(PyObject* self, PyObject* args)
{
  char* fn;
  sardata* output;
  int M, N, D;
  PYINIT;

  // Input
  if (!PyArg_ParseTuple(args, "s", &fn))
    return NULL;

  // Processing
  output = sardata_alloc();
  if (!(output = sarread_header(fn, output)))
    {
      sarerror_msg_msg("Cannot open file %s", fn);
      PyErr_SetString(PyExc_IOError, sarerror);
      return NULL;
    }
  M = output->M;
  N = output->N;
  D = output->D;
  sardata_free(output);

  // Output
  return Py_BuildValue("iii", M, N, D);
}

static PyObject* py_sar2rgb(PyObject* self, PyObject* args)
{
  float alpha = 3, gamma = 0.7;
  PyArrayObject* py_in;
  PyObject* py_out;
  sardata* input = NULL;
  rgbdata* rgb;
  npy_intp dims[3];
  int ndims;

  PYINIT;

  // Input
  if (!PyArg_ParseTuple(args, "O|ff", (PyObject *) &py_in, &alpha, &gamma))
    return NULL;
  if (py_in == NULL)
    return NULL;
  if (!(input = PyArrayObject2sardata(py_in)))
    return NULL;

  // Processing
  rgb = rgbdata_alloc();
  if (!(rgb = sar2rgb(input, rgb, alpha, gamma)))
    {
      sarerror_msg_msg("Cannot create an RGB representation");
      PyErr_SetString(PyExc_Exception, sarerror);
      return NULL;
    }

  // Output
  ndims = 3;
  dims[0] = rgb->M;
  dims[1] = rgb->N;
  dims[2] = 3;
  py_out = PyArray_SimpleNewFromData(ndims, dims,
				     NPY_UBYTE,
				     rgb->array);
  PyArray_ENABLEFLAGS((PyArrayObject *) py_out, NPY_OWNDATA);
  free(rgb);
  sardata_free(input);

  return Py_BuildValue("O", py_out);
}

static PyObject* py_sarshow(PyObject* self, PyObject* args)
{
  PyObject* py_rgb;

  PYINIT;

  // Input
  py_rgb = py_sar2rgb(self, args);

  PyObject* module_plt = PyImport_Import(PyBytes_FromString("matplotlib.pyplot"));
  if (!module_plt)
    {
      sarerror_msg("Module matplotlib.pyplot unkown");
      PyErr_SetString(PyExc_ImportError, sarerror);
      return NULL;
    }

  PyObject* imshow = PyObject_GetAttrString(module_plt, "imshow");
  if (!imshow)
    {
      sarerror_msg("Function imshow unkown");
      PyErr_SetString(PyExc_NameError, sarerror);
      return NULL;
    }

  PyObject* keywords_interpolation = PyDict_New();
  PyDict_SetItemString(keywords_interpolation, "interpolation", PyBytes_FromString("nearest"));
  PyDict_SetItemString(keywords_interpolation, "aspect", PyBytes_FromString("equal"));
  PyObject_Call(imshow, PyTuple_Pack(1, py_rgb), keywords_interpolation);

  PyObject* axis = PyObject_GetAttrString(module_plt, "axis");
  if (!axis)
    {
      sarerror_msg("Function axis unkown");
      PyErr_SetString(PyExc_NameError, sarerror);
      return NULL;
    }
  PyObject_CallObject(axis, PyTuple_Pack(1, PyBytes_FromString("off")));

  return Py_True;
}

static PyObject* py_sarnlsar(PyObject* self, PyObject* args)
{
  float L;
  int hW = 12, hP = 5, verbose = 1;
  PyArrayObject* py_in;
  PyArrayObject* py_noise = NULL;
  PyObject* py_out;
  PyObject* py_look;
  npy_intp dims[4];
  int ndims;
  sardata* input;
  sardata* sar_noise = NULL;
  sardata* output;
  fltdata* look = NULL;

  PYINIT;

  // Input
  if (!PyArg_ParseTuple(args, "Of|iiiO", (PyObject *) &py_in, &L, &verbose, &hW, &hP, &py_noise))
    return NULL;
  if (py_in == NULL)
    return NULL;
  if (!(input = PyArrayObject2sardata(py_in)))
    return NULL;
  if (py_noise)
    {
      if (!(sar_noise = PyArrayObject2sardata(py_noise)))
	return NULL;
    }

  // Processing
  output = sardata_alloc();
  if (!(output = sarnlsar(input, output, L, 5, verbose, hW, hP, sar_noise, &look)))
    {
      PyErr_SetString(PyExc_Exception, sarerror);
      return NULL;
    }

  // Output
  ndims = 4;
  dims[0] = output->M;
  dims[1] = output->N;
  dims[2] = output->D;
  dims[3] = output->D;
  py_out = PyArray_SimpleNewFromData(ndims, dims,
				     NPY_CFLOAT,
				     output->array);
  PyArray_ENABLEFLAGS((PyArrayObject *) py_out, NPY_OWNDATA);

  ndims = 2;
  dims[0] = look->M;
  dims[1] = look->N;
  py_look = PyArray_SimpleNewFromData(ndims, dims,
				      NPY_FLOAT,
				      look->array);
  PyArray_ENABLEFLAGS((PyArrayObject *) py_out, NPY_OWNDATA);

  if (sar_noise)
    sardata_free(sar_noise);
  if (look)
    free(look);
  sardata_free(input);
  free(output);

  return Py_BuildValue("OO", py_out, py_look);
}

static PyObject* py_sarboxcar(PyObject* self, PyObject* args)
{
  int hW = 1;
  PyArrayObject* py_in;
  PyObject* py_out;
  npy_intp dims[4];
  int ndims;
  sardata* input;
  sardata* output = NULL;

  PYINIT;

  // Input
  if (!PyArg_ParseTuple(args, "O|i", (PyObject *) &py_in, &hW))
    return NULL;
  if (py_in == NULL)
    return NULL;
  if (!(input = PyArrayObject2sardata(py_in)))
    return NULL;

  // Processing
  output = sardata_alloc();
  if (!(output = sarboxcar(input, output, hW)))
    {
      PyErr_SetString(PyExc_Exception, sarerror);
      return NULL;
    }

  // Output
  ndims = 4;
  dims[0] = output->M;
  dims[1] = output->N;
  dims[2] = output->D;
  dims[3] = output->D;
  py_out = PyArray_SimpleNewFromData(ndims, dims,
				     NPY_CFLOAT,
				     output->array);
  PyArray_ENABLEFLAGS((PyArrayObject *) py_out, NPY_OWNDATA);

  sardata_free(input);
  free(output);

  return Py_BuildValue("O", py_out);
}

static PyObject* py_sardiskcar(PyObject* self, PyObject* args)
{
  int hW = 1;
  PyArrayObject* py_in;
  PyObject* py_out;
  npy_intp dims[4];
  int ndims;
  sardata* input;
  sardata* output = NULL;

  PYINIT;

  // Input
  if (!PyArg_ParseTuple(args, "O|i", (PyObject *) &py_in, &hW))
    return NULL;
  if (py_in == NULL)
    return NULL;
  if (!(input = PyArrayObject2sardata(py_in)))
    return NULL;

  // Processing
  output = sardata_alloc();
  if (!(output = sardiskcar(input, output, hW)))
    {
      PyErr_SetString(PyExc_Exception, sarerror);
      return NULL;
    }

  // Output
  ndims = 4;
  dims[0] = output->M;
  dims[1] = output->N;
  dims[2] = output->D;
  dims[3] = output->D;
  py_out = PyArray_SimpleNewFromData(ndims, dims,
				     NPY_CFLOAT,
				     output->array);

  PyArray_ENABLEFLAGS((PyArrayObject *) py_out, NPY_OWNDATA);

  sardata_free(input);
  free(output);

  return Py_BuildValue("O", py_out);
}

static PyObject* py_sargausscar(PyObject* self, PyObject* args)
{
  int hW = 1;
  PyArrayObject* py_in;
  PyObject* py_out;
  npy_intp dims[4];
  int ndims;
  sardata* input;
  sardata* output = NULL;

  PYINIT;

  // Input
  if (!PyArg_ParseTuple(args, "O|i", (PyObject *) &py_in, &hW))
    return NULL;
  if (py_in == NULL)
    return NULL;
  if (!(input = PyArrayObject2sardata(py_in)))
    return NULL;

  // Processing
  output = sardata_alloc();
  if (!(output = sargausscar(input, output, hW)))
    {
      PyErr_SetString(PyExc_Exception, sarerror);
      return NULL;
    }

  // Output
  ndims = 4;
  dims[0] = output->M;
  dims[1] = output->N;
  dims[2] = output->D;
  dims[3] = output->D;
  py_out = PyArray_SimpleNewFromData(ndims, dims,
				     NPY_CFLOAT,
				     output->array);
  PyArray_ENABLEFLAGS((PyArrayObject *) py_out, NPY_OWNDATA);

  sardata_free(input);
  free(output);

  return Py_BuildValue("O", py_out);
}

static PyMethodDef nlsartoolbox_methods[] = {
  {"sarread",     py_sarread,     METH_VARARGS},
  {"sarwrite",    py_sarwrite,    METH_VARARGS},
  {"sarinfo",     py_sarinfo,     METH_VARARGS},
  {"sar2rgb",     py_sar2rgb,     METH_VARARGS},
  {"sarshow",     py_sarshow,     METH_VARARGS},
  {"sarnlsar",    py_sarnlsar,    METH_VARARGS},
  {"sarboxcar",   py_sarboxcar,   METH_VARARGS},
  {"sardiskcar",  py_sardiskcar,  METH_VARARGS},
  {"sargausscar", py_sargausscar, METH_VARARGS},
  {NULL, NULL}
};

/*
void initnlsartoolbox()
{
  (void) Py_InitModule("nlsartoolbox", nlsartoolbox_methods);
}
*/

static int nlsartoolbox_traverse(PyObject *m, visitproc visit, void *arg) 
{
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int nlsartoolbox_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "nlsartoolbox",
        NULL,
        sizeof(struct module_state),
        nlsartoolbox_methods,
        NULL,
        nlsartoolbox_traverse,
        nlsartoolbox_clear,
        NULL
};

PyObject *PyInit_nlsartoolbox()
{
    PyObject *module = PyModule_Create(&moduledef);   
    return module;
}
