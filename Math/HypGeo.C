#include "Math/HypGeo.H"
#include "ATOOLS/Org/Exception.H"

// using namespace RESUM;

double _HypGeo_3F2_Py(double a1, double a2, double a3, double b1, double b2, double x);
double _HypGeo_3F2_Arb(double a1, double a2, double a3, double b1, double b2, double x,
                         int regularized=0, int prec=64);

namespace RESUM {

  double HypGeo_3F2(double a1, double a2, double a3, double b1, double b2, double x) {
#ifdef USING_ARBLIB
    return _HypGeo_3F2_Arb(a1,  a2,  a3,  b1,  b2,  x);
#elif defined USING_PYTHON
    return _HypGeo_3F2_Py(a1,  a2,  a3,  b1,  b2,  x);
#else
    THROW(fatal_error,"No hypergeometric function available. Use --enable-arblib or --enable-python.");
    return 0;
#endif
  }
}


#ifdef USING_PYTHON
#include <python2.7/Python.h>
double _HypGeo_3F2_Py(double a1, double a2, double a3, double b1, double b2, double x){
    PyObject *pName, *pModule, *pDict, *pFunc;
    PyObject *pArgs, *pValue;
    double result;
//     std::cout << "testing" << std::endl;
    Py_Initialize();
    pName = PyString_FromString("mpmath");
    /* Error checking of pName left out */

    pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule != NULL) {
        pFunc = PyObject_GetAttrString(pModule, "hyp3f2");
        /* pFunc is a new reference */

        if (pFunc && PyCallable_Check(pFunc)) {
            pArgs = PyTuple_New(6);
            pValue = PyFloat_FromDouble(a1);
            if (!pValue) {
                Py_DECREF(pArgs);
                Py_DECREF(pModule);
                fprintf(stderr, "Cannot convert argument\n");
                return 1;
            }
            PyTuple_SetItem(pArgs, 0, pValue);
            
            pValue = PyFloat_FromDouble(a2);
            if (!pValue) {
                Py_DECREF(pArgs);
                Py_DECREF(pModule);
                fprintf(stderr, "Cannot convert argument\n");
                return 1;
            }
            PyTuple_SetItem(pArgs, 1, pValue);
            
            pValue = PyFloat_FromDouble(a3);
            if (!pValue) {
                Py_DECREF(pArgs);
                Py_DECREF(pModule);
                fprintf(stderr, "Cannot convert argument\n");
                return 1;
            }
            PyTuple_SetItem(pArgs, 2, pValue);
            
            pValue = PyFloat_FromDouble(b1);
            if (!pValue) {
                Py_DECREF(pArgs);
                Py_DECREF(pModule);
                fprintf(stderr, "Cannot convert argument\n");
                return 1;
            }
            PyTuple_SetItem(pArgs, 3, pValue);
            
            pValue = PyFloat_FromDouble(b2);
            if (!pValue) {
                Py_DECREF(pArgs);
                Py_DECREF(pModule);
                fprintf(stderr, "Cannot convert argument\n");
                return 1;
            }
            PyTuple_SetItem(pArgs, 4, pValue);
            
            pValue = PyFloat_FromDouble(x);
            if (!pValue) {
                Py_DECREF(pArgs);
                Py_DECREF(pModule);
                fprintf(stderr, "Cannot convert argument\n");
                return 1;
            }
            PyTuple_SetItem(pArgs, 5, pValue);
            
            pValue = PyObject_CallObject(pFunc, pArgs);
            Py_DECREF(pArgs);
            if (pValue != NULL) {
                result = PyFloat_AsDouble(pValue);
                Py_DECREF(pValue);
            }
            else {
                Py_DECREF(pFunc);
                Py_DECREF(pModule);
                PyErr_Print();
                fprintf(stderr,"Call failed\n");
                return 1;
            }
        }
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \"%s\"\n", "hyp3f2");
        }
        Py_XDECREF(pFunc);
        Py_DECREF(pModule);
    }
    else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", "mpmath");
        return 1;
    }
    Py_Finalize();
    return(result);
}
#endif

#ifdef USING_ARBLIB
#include "arb.h"
#include "arb_hypgeom.h"

double _HypGeo_3F2_Arb(double a1, double a2, double a3, double b1, double b2, double z, int regularized, int prec){
  // init arb variables
  arb_t result; arb_init(result);
  arb_t arg; arb_init(arg); arb_set_d(arg,z);
  arb_ptr a = _arb_vec_init(3);
  arb_set_d(a,a1); arb_set_d(a+1,a2); arb_set_d(a+2,a3);
  arb_ptr b = _arb_vec_init(2);
  arb_set_d(b,b1); arb_set_d(b+1,b2);

  // actual calculation
  arb_hypgeom_pfq(result, a, 3, b, 2, arg, regularized, prec);

  // clean up
  arb_clear(arg);
  arb_clear(a+2); arb_clear(a+1); arb_clear(a);
  delete a;
  arb_clear(b+1); arb_clear(b);
  delete b;
  const double ret = arf_get_d(arb_midref(result), ARF_RND_NEAR);
  arb_clear(result);

  return ret;
}
#endif
