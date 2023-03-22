//
//  Regression.hpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#ifndef Regression_hpp
#define Regression_hpp

#include <stdio.h>
#include <chrono>
#include <Python.h>
#include <string>
#include <iostream>
#include <cstdlib>
#include <vector>
#include "config.hpp"

using namespace std;

/*
class NumpyMatrix
{
    
public:
    NumpyMatrix(float* C_Array, uint size)
    {
        npy_intp dims[2]{size, size};
        data = reinterpret_cast<PyArrayObject*>(PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, reinterpret_cast<void*>(C_Array)));
    }
    
    PyArrayObject* data;
    
};

class NumpyVector
{
    
public:
    NumpyVector(float* C_Array, uint size)
    {
        npy_intp dims[1]{size};
        data = reinterpret_cast<PyArrayObject*>(PyArray_SimpleNewFromData(1, dims, NPY_FLOAT, reinterpret_cast<void*>(C_Array)));
    }
    
    PyArrayObject* data;
    
};
 
 */

class Regression
{
public:
    Regression()
    {
        //Py_Initialize();
    };
    ~Regression()
    {
        //Py_Finalize();
    };
    
    int lawson_hanson_nnls(const float *kernal_vec, const float *weightnorm, uint16 size, float *coefs, float *residuel);
    //int Regess3(float *kernal_vec, float *weightnorm, uint size, float * &coefs, float &residuel);
    //int Regess2(float *kernal_vec, float *weightnorm, uint size, float * &coefs, float &residuel);
    //int Regess(float *kernal_vec, float *weightnorm, uint size, float * &coefs, float &residuel);
    void Call(uint size , float *kernal_vec, float *weightnorm, float total_lambda, float *unweightnorm, float * coefs, float * residuel);
    //int test();
    
    vector<bool> active_or_passive = vector<bool> (MAX_UINT16);
        
    vector<uint16> passive_set = vector<uint16> (MAX_UINT16);
    
    vector<uint16> passive_set_old = vector<uint16> (MAX_UINT16);
        
    vector<float> x_trial_vec = vector<float> (MAX_UINT16);

};



#endif /* Regression_hpp */
