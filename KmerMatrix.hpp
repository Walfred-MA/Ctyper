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

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "config.hpp"
#include "PriorData.hpp"

class KmerMatrix
{
    
public:
    
    KmerMatrix()
    {
        weightnorm = (double*) realloc(weightnorm,sizeof(double*) * (genenum * genenum + 1));
        
        kernal_vec = (double*) realloc(kernal_vec,sizeof(double*) * (genenum+1) );
        
        norm_offsites = (double*) realloc(norm_offsites, sizeof (double) * (genenum+1) );
    };
    ~KmerMatrix()
    {
        free(weightnorm);
        
        free(kernal_vec);
        
        free(norm_offsites);
    };
    
    void InitiateMatrix();
    void FinishMatrix(const double * prior);
        
    void getNormEachRow(const int count, const uint16 *row, const uint16 rowsize, const uint16 sign);
    void getNorm(const PriorChunk* Data, const uint16* allkmervec, const float depth);
    
    float depth = 14.0;
    
    double* weightnorm = NULL;
    double* kernal_vec = NULL;
    
    double* norm_offsites = NULL ;
    double norm_offsite = 0;
    
    double vec_offsite = 0;
    double kmer_counts = 0;
        
    uint16 genenum = 0;
    size_t kmernum = 0;
    
    const uint16* kmervec;
    
    size_t alloc_size;
    
};

#endif /* KmerMatrix_hpp */
