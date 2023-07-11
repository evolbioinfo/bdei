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
    double pie = -1;
    double ut = 0;
    int u = -1;
    int nt = 0;
    int nstarts = 1;
    int debug = 2;

    PyObject *startobj, *ubobj;

    // Define keywords
    static const char *kwlist[] = {"f", "start", "ub", "pie", "mu", "la", "psi", "p", "u", "ut", "nt", "nbiter", "debug", "nstarts", NULL};

    // Interpret input arguments.
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "sOO|dddddidiiii",
        const_cast<char**>(kwlist), &treename, &startobj, &ubobj, &pie, &mu, &la, &psi, &p, &u, &ut, &nt, &nbiter, &debug, &nstarts)) {
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

    PyObject *pysol = NULL;
    // Call the C++ functions to infer parameters
    try {
        Solution sol = *inferParameters(treename, ts, ubs, pie, mu, la, psi, p, u, ut, nbiter, nt, debug, nstarts);
        pysol = PyList_New(15);
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

        PyList_SetItem(pysol, 13, Py_BuildValue("d", sol.cpu_time));

        PyList_SetItem(pysol, 14, Py_BuildValue("i", sol.nb_iter));
    } catch(const std::invalid_argument& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
    }

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
    double pie = -1;
    double ut = 0;
    int u = -1;
    int nt = 0;
    int debug = 2;

    // Define keywords
    static const char *kwlist[] = {"f", "mu", "la", "psi", "p", "pie", "u", "ut", "nt", "debug", NULL};

    // Interpret input arguments.
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "sdddddidii",
        const_cast<char**>(kwlist), &treename, &mu, &la, &psi, &p, &pie, &u, &ut, &nt, &debug)) {
        PyErr_Format(PyExc_ValueError, "Could not cast the input arguments.");
        return NULL;
    }

    // Call the C++ functions to infer parameters
    try {
        double res = calculateLikelihood(treename, mu, la, psi, p, pie, u, ut, nt, debug);
        return Py_BuildValue("d", res);
    } catch(const std::invalid_argument& e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }
}
