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

using namespace std;


#include "config.hpp"

class KmerMatrix
{
    
public:
    
    KmerMatrix()
    {};
    ~KmerMatrix()
    {};
    
    void getNorm(const uint16* kmervec, const uint16* kmermatrix, const float depth, const uint16 gnum, const uint knum, float* norm_vec, float* norm_matrix, float &total_lambda);
    
private:
        
    unique_ptr<float> row_offsites = unique_ptr<float>( new float[MAX_UINT16] ) ;
    unique_ptr<float> diag_offsites = unique_ptr<float>( new float[MAX_UINT16] ) ;
    
    
    
};

#endif /* KmerMatrix_hpp */
