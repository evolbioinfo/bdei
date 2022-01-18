#pragma once
#include <array>
#include <string>
#include "_python.hpp"

// Docstring for the module
static char module_docstring[] =
"pybdei: this module provides an interface for BDEI, for fast and accurate epidemiological parameter estimation from phylogenetic trees with the Birth-Death Exposed-Infectious (BDEI) model.";

// Docstring for the Solution::solve() function
static char infer_docstring[] =
"Runs the parameter inference";
static char likelihood_docstring[] =
"Runs the likelihood calculation";

// Available functions in the pybdei module
static PyObject *_pybdei_infer(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *_pybdei_likelihood(PyObject *self, PyObject *args, PyObject *kwargs);

// Module interface
static PyMethodDef module_methods[] = {
    {"infer", (PyCFunction) _pybdei_infer, METH_VARARGS | METH_KEYWORDS, infer_docstring},
    {"likelihood", (PyCFunction) _pybdei_likelihood, METH_VARARGS | METH_KEYWORDS, likelihood_docstring},
    {NULL, NULL, 0, NULL}
};

#ifdef PYTHON3
static struct PyModuleDef _pybdeimodule = {
    PyModuleDef_HEAD_INIT,
    "_pybdei",
    module_docstring,
    -1,
    module_methods
};
#endif