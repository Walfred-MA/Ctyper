//
//  Regression.cpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#include "Regression.hpp"


double Regression::Regess( double *kernal_vec, double *weightnorm, uint size, double * &coefs)
{
    Py_Initialize();
    import_array();
    
    NumpyMatrix matrix(weightnorm, size);
    NumpyVector vector(kernal_vec, size);
    
    cout<<"check8.01"<<endl;
    
    PyObject *pArgs = PyTuple_New(2);
    PyTuple_SetItem(pArgs, 0, reinterpret_cast<PyObject*>(matrix.data));
    PyTuple_SetItem(pArgs, 1, reinterpret_cast<PyObject*>(vector.data));
    
    cout<<"check8.1"<<endl;
    
    PyRun_SimpleString("from scipy.optimize import lsq_linear,nnls\n");
    PyObject *pName = PyUnicode_FromString("scipy.optimize");
    PyObject *pModule = PyImport_Import(pName);
    PyObject *pFunc = PyObject_GetAttrString(pModule, "nnls");
    
    cout<<"check8.2"<<endl;
    
    PyObject * Result = PyObject_CallObject(pFunc, pArgs);
    
    coefs = reinterpret_cast<double*>(PyArray_DATA(reinterpret_cast<PyArrayObject*>(PyTuple_GetItem(Result, 0))));
    double residuel = PyFloat_AsDouble(PyTuple_GetItem(Result, 1));
    
    cout<<"check8.3"<<endl;
    
    Py_DECREF (pName);
    Py_DECREF (pModule);
    Py_DECREF (pFunc);
    Py_DECREF (Result);
    Py_DECREF (pArgs);
    Py_DECREF (matrix.data);
    Py_DECREF (vector.data);
    
    Py_Finalize();
    
    return residuel;
    
}

pair<double, double*> Regression::Call(double *kernal_vec, double *weightnorm, ull kmer_counts, double *unweightnorm, uint size)
{
    
    double * regression = NULL;
    
    double residuel = Regess(kernal_vec, weightnorm, size, regression);
    
    cout<<"check8.1"<<endl;
    
    float regressed_kmer = 0.0;
    
    for (int i = 0; i < size; ++i)
    {
        regressed_kmer += unweightnorm[i*size + i] * regression[i];
    }
    
    cout<<"check8.2"<<endl;
    
    float correction = kmer_counts/regressed_kmer;
    
    for (int i = 0; i < size; ++i)
    {
        regression[i] *= correction;
    }
     
    return make_pair( residuel, regression);
}
