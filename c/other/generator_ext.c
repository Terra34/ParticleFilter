#include <python3.8/Python.h>
#include <numpy/arrayobject.h>
#include <stdlib.h>
#include <math.h>
#include "../headers/generators.h"

/*
    Python module based on generators.h
*/

static PyObject *_MarsagliaPolar(PyObject *self, PyObject *args)
{
    PyArrayObject *py_mat;
    float *mat;
    npy_intp dims[1] = { 0 };
    int n, i;
    int seed;
    struct gaussGenState state;
    if (!PyArg_ParseTuple(args, "ii", &n, &seed))
        return NULL;

    dims[0] = n;
    py_mat = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT32);
    mat = (float *)py_mat->data;

    setSeed(seed);
    initializeGauss(&state);
    for (i = 0; i < n; ++i)
      mat[i] = gauss(&state);

    return PyArray_Return(py_mat);
}

static PyObject *_BoxMuller(PyObject *self, PyObject *args)
{
    PyArrayObject *py_mat;
    float *mat;
    npy_intp dims[1] = { 0 };
    int n, i;
    int seed;
    struct gaussGenState state;
    if (!PyArg_ParseTuple(args, "ii", &n, &seed))
        return NULL;

    dims[0] = n;
    py_mat = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT32);
    mat = (float *)py_mat->data;

    setSeed(seed);
    initializeGauss(&state);
    for (i = 0; i < n; ++i)
      mat[i] = gaussbm(&state);

    return PyArray_Return(py_mat);
}

static PyObject *_Inverse(PyObject *self, PyObject *args)
{
    PyArrayObject *py_mat;
    float *mat;
    npy_intp dims[1] = { 0 };
    int n, i;
    int seed;
    if (!PyArg_ParseTuple(args, "ii", &n, &seed))
        return NULL;

    dims[0] = n;
    py_mat = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT32);
    mat = (float *)py_mat->data;

    setSeed(seed);
    for (i = 0; i < n; ++i)
      mat[i] = gaussInv();

    return PyArray_Return(py_mat);
}

static PyObject *_Uniform(PyObject *self, PyObject *args)
{
    PyArrayObject *py_mat;
    float *mat;
    npy_intp dims[1] = { 0 };
    int n, i;
    int seed;
    if (!PyArg_ParseTuple(args, "ii", &n, &seed))
        return NULL;

    dims[0] = n;
    py_mat = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT32);
    mat = (float *)py_mat->data;

    setSeed(seed);
    for (i = 0; i < n; ++i)
      mat[i] = uniform();

    return PyArray_Return(py_mat);
}

static PyObject *_Ziggurat(PyObject *self, PyObject *args)
{
    PyArrayObject *py_mat;
    float *mat;
    npy_intp dims[1] = { 0 };
    int n, i;
    int seed;
    if (!PyArg_ParseTuple(args, "ii", &n, &seed))
        return NULL;

    dims[0] = n;
    py_mat = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT32);
    mat = (float *)py_mat->data;

    setSeed(seed);
    zigset();
    for (i = 0; i < n; ++i)
      mat[i] = ziggurat();

    return PyArray_Return(py_mat);
}

static PyObject *_LFSR(PyObject *self, PyObject *args)
{
    PyArrayObject *py_mat;
    float *mat;
    npy_intp dims[1] = { 0 };
    int n, i;
    int taps[3] = {0, 3, 31}; // change later to be a parameter
    struct genState state;
    if (!PyArg_ParseTuple(args, "i", &n))
        return NULL;

    dims[0] = n;
    py_mat = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT32);
    mat = (float *)py_mat->data;

    initializeGenerator(&state, 32, 8, taps, 0.f, 1.f);
    seedGenerator(&state, ((long long)1<<32) - 1);

    for (i = 0; i < n; ++i)
      mat[i] = generate(&state);

    return PyArray_Return(py_mat);
}

static PyMethodDef genMethods[] = {
    {"marsaglia", _MarsagliaPolar, METH_VARARGS, "Marsaglia's polar normal distribution random generator implementation."},
    {"boxmuller", _BoxMuller, METH_VARARGS, "Box-Muller normal distribution random generator implementation."},
    {"inverse", _Inverse, METH_VARARGS, "Inverse cumulative distribution approximation with relative error less than 1.15 x 10e-9 in the entire region."},
    {"uniform", _Uniform, METH_VARARGS, "Marsaglia SHR3 method for generating random uniform numbers."},
    {"ziggurat", _Ziggurat, METH_VARARGS, "Ziggurat algorithm for generating random normally distributed numbers."},
    {"lfsr", _LFSR, METH_VARARGS, "Linear feedback shift register implementation for generating gaussian-like distribution."},
    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef genModule = {
    PyModuleDef_HEAD_INIT,
    "generators",
    "Generator module with gaussian random generators - Contains: Inverse, Box-Muller, Marsaglia Polar, Ziggurat and Linear feedback shift register implementations.",
    -1,
    genMethods
};


PyMODINIT_FUNC PyInit_generators(void){
    import_array();
    return PyModule_Create(&genModule);
}
