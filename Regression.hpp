//
//  Regression.hpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#ifndef Regression_hpp
#define Regression_hpp

#include <stdio.h>
#include <Python.h>
#include <string>
#include <iostream>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#include "config.hpp"

using namespace std;

class NumpyMatrix
{
    
public:
    NumpyMatrix(double* C_Array, uint size)
    {
        npy_intp dims[2]{size, size};
        data = reinterpret_cast<PyArrayObject*>(PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, reinterpret_cast<void*>(C_Array)));
    }

    PyArrayObject* data;
    
};

class NumpyVector
{
    
public:
    NumpyVector(double* C_Array, uint size)
    {
        npy_intp dims[1]{size};
        data = reinterpret_cast<PyArrayObject*>(PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, reinterpret_cast<void*>(C_Array)));
    }

    PyArrayObject* data;
    
};

class Regression
{
public:
    Regression(){};
    ~Regression(){};
    
    double Regess(double *kernal_vec, double *weightnorm, uint size, double * &coefs);
    pair<double, double*> Call(double *kernal_vec, double *weightnorm, ull kmer_counts, double *unweightnorm, uint size);
    
    //double * coefs = NULL;
    //double residuel = 0.0;
};



#endif /* Regression_hpp */
