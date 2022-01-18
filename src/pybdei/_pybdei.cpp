#define PY_SSIZE_T_CLEAN
#include "_pybdei.hpp"
#include "pool.hpp"
#include "queue.hpp"
#include "semaphore.hpp"
#include "BDEI.hpp"
#include <string>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#ifdef PYTHON3
PyMODINIT_FUNC PyInit__pybdei(void){
    import_array();
    return PyModule_Create(&_pybdeimodule);
}
#else
PyMODINIT_FUNC init_pybdei(void){
    PyObject *m = Py_InitModule3("_pybdei", module_methods, module_docstring);
    if(m==NULL)
        return;
    import_array();
}
#endif

static PyObject *_pybdei_infer(PyObject *self, PyObject *args, PyObject *kwargs) {
    char* treename;
    int nbiter = 0; //
    double p = -1;
    double la = -1;
    double mu = -1;
    double psi = -1;
    double T = 0;
    int u = 0;
    int nt = 1;

    PyObject *startobj, *ubobj;

    // Define keywords
    static const char *kwlist[] = {"f", "start", "ub", "mu", "la", "psi", "p", "T", "u", "nt", "nbiter", NULL};

    // Interpret input arguments.
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "sOO|dddddiii",
        const_cast<char**>(kwlist), &treename, &startobj, &ubobj, &mu, &la, &psi, &p, &T, &u, &nt, &nbiter)) {
        PyErr_Format(PyExc_ValueError, "Could not cast the input arguments.");
        return NULL;
    }
    // Interpret input objects as numpy arrays
    PyObject *startarray = PyArray_FROM_OTF(startobj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyObject *ubarray = PyArray_FROM_OTF(ubobj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    // If that didn't work, throw an exception
    if(startarray==NULL or ubarray==NULL){
        Py_XDECREF(startarray);
        Py_XDECREF(ubarray);
        PyErr_Format(PyExc_ValueError, "Could not convert the start/upper bound values.");
        return NULL;
    }

    // Get pointers to the data as c++-types
    PyArrayObject *startarray_arr = reinterpret_cast<PyArrayObject*>(startarray);
    PyArrayObject *ubarray_arr = reinterpret_cast<PyArrayObject*>(ubarray);
    double *ts = (double*)PyArray_DATA(startarray_arr);
    double *ubs = (double*)PyArray_DATA(ubarray_arr);

    // Call the C++ functions to infer parameters
    Solution sol = *inferParameters(treename, ts, ubs, mu, la, psi, p, T, u, nbiter);

    PyObject *pysol=PyList_New(13);
    PyList_SetItem(pysol, 0, Py_BuildValue("d", sol.mu));
    PyList_SetItem(pysol, 1, Py_BuildValue("d", sol.la));
    PyList_SetItem(pysol, 2, Py_BuildValue("d", sol.psi));
    PyList_SetItem(pysol, 3, Py_BuildValue("d", sol.p));

    PyList_SetItem(pysol, 4, Py_BuildValue("d", sol.mu_min));
    PyList_SetItem(pysol, 5, Py_BuildValue("d", sol.mu_max));

    PyList_SetItem(pysol, 6, Py_BuildValue("d", sol.la_min));
    PyList_SetItem(pysol, 7, Py_BuildValue("d", sol.la_max));

    PyList_SetItem(pysol, 8, Py_BuildValue("d", sol.psi_min));
    PyList_SetItem(pysol, 9, Py_BuildValue("d", sol.psi_max));

    PyList_SetItem(pysol, 10, Py_BuildValue("d", sol.p_min));
    PyList_SetItem(pysol, 11, Py_BuildValue("d", sol.p_max));

    PyList_SetItem(pysol, 12, Py_BuildValue("d", sol.likelihood));

    // Clean up
    Py_DECREF(startarray);
    Py_DECREF(ubarray);

    return pysol;
}


static PyObject *_pybdei_likelihood(PyObject *self, PyObject *args, PyObject *kwargs) {
    char* treename;
    double p = -1;
    double la = -1;
    double mu = -1;
    double psi = -1;
    double T = 0;
    int u = 0;

    // Define keywords
    static const char *kwlist[] = {"f", "mu", "la", "psi", "p", "T", "u", NULL};

    // Interpret input arguments.
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "sdddddi",
        const_cast<char**>(kwlist), &treename, &mu, &la, &psi, &p, &T, &u)) {
        PyErr_Format(PyExc_ValueError, "Could not cast the input arguments.");
        return NULL;
    }

    // Call the C++ functions to infer parameters
    double res = calculateLikelihood(treename, mu, la, psi, p, T, u);
    return Py_BuildValue("d", res);
}
