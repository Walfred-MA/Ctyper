//
//  KmerMatrix.hpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#ifndef KmerMatrix_hpp
#define KmerMatrix_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <tuple>
#include <mutex>

using namespace std;


#include "config.hpp"

class KmerMatrix
{
    
public:
    
    KmerMatrix()
    {};
    ~KmerMatrix()
    {};
    
    void getNorm(const uint16* kmervec, const uint16* kmermatrix, const float depth, const uint16 gnum, const uint knum, FLOAT_T* norm_vec, FLOAT_T* norm_matrix, double &total_lambda);
    
private:
        
    unique_ptr<double > row_offsites = unique_ptr<double >( new double [MAX_UINT16] ) ;
    unique_ptr<double > diag_offsites = unique_ptr<double >( new double [MAX_UINT16] ) ;
    
    
    
};

#endif /* KmerMatrix_hpp */
