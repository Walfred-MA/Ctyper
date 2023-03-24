//
//  Regression.hpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#ifndef Regression_hpp
#define Regression_hpp

#include <stdio.h>
#include <string>
#include <iostream>
#include <cstdlib>

#include "config.hpp"

using namespace std;

class Regression
{
public:
    Regression()
    {};
    ~Regression()
    {};
    
    int lawson_hanson_nnls(const FLOAT_T *kernal_vec, const FLOAT_T *weightnorm, uint16 size, FLOAT_T *coefs, FLOAT_T *residuel);

    void Call(uint size , FLOAT_T *kernal_vec, FLOAT_T *weightnorm, float total_lambda, const uint *kmercounts, FLOAT_T * coefs, FLOAT_T * residuel);

private:
    
    vector<bool> active_or_passive = vector<bool> (MAX_UINT16);
        
    vector<uint16> passive_set = vector<uint16> (MAX_UINT16);
    
    vector<uint16> passive_set_old = vector<uint16> (MAX_UINT16);
        
};



#endif /* Regression_hpp */
