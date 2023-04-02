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
#include <memory>
using namespace std;


#include "config.hpp"

class KmerMatrix
{
    
public:
    
    KmerMatrix()
    {};
    ~KmerMatrix()
    {};
    
    void getNorm(const uint16* kmervec, const uint16* kmermatrix, const FLOAT_T depth, const uint16 gnum, const uint knum, FLOAT_T* norm_vec, FLOAT_T* norm_matrix, FLOAT_T &total_lambda);
    
private:
        
    unique_ptr<FLOAT_T> row_offsites = unique_ptr<FLOAT_T >( new FLOAT_T [MAX_UINT16] ) ;
    unique_ptr<FLOAT_T> diag_offsites = unique_ptr<FLOAT_T>( new FLOAT_T [MAX_UINT16] ) ;
    
    
    
};

#endif /* KmerMatrix_hpp */
