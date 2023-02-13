//
//  Kmatrix.hpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 1/31/23.
//  Copyright © 2023 USC_Mark. All rights reserved.
//

#ifndef Kmatrix_hpp
#define Kmatrix_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <tuple>

using namespace std;

typedef unsigned short uint16;
typedef unsigned int uint;
typedef unsigned __int128 u128;
typedef unsigned long long ull;
typedef unsigned char  uint8;

#include "ktable.hpp"


class Kmatrix
{
    
public:
    
    Kmatrix(const char* mfile, float dep, ull msize)
    {
        StrLine.reserve(10000);
        row.reserve(10000);
        depth = dep;
        kmervec = (uint16* )malloc(sizeof(uint16) * msize);
        memset(kmervec, 0, sizeof(uint16) * msize);
        tablefile = new ktable(mfile);
    };
    ~Kmatrix()
    {
        for (int i = 0; i<matrixinfo.size()-1;++i)
        {
            free(weightnorms[i]);
            free(kernal);
        }
        
        
        free(kernal);
        free(weightnorm);
        free(kmervec);
        
        delete tablefile;
    }
        
    
    void LoadRow(const char* line, int linesize);
    int LoadHeader(const char* line, int linesize);
    void OpenFile(const char* inputfile);
    void ProcessMatrix();
    void eachRowNorm(int count);
    void FinishMatrix();
    void InitiateMatrix();
    void Process();
    void write(const char* outputfile, const char* prefix );
    ull kmerindex = 0;
    
    vector<double*> weightnorms;
    vector<double*> kernals;
    vector<double>  totalkmers;
    vector<tuple<string,ull,ull,uint16>> matrixinfo;
    uint16* kmervec=NULL;
    
    tuple<string,ull,ull,uint16>* curr_info;
    double* weightnorm = NULL;
    double* kernal = NULL;
    double* totalkmer= NULL;
    
    int sign = 1;
    
    float depth = 14.0;
    int klen = 31;
    
    double vec_offsite = 0;
    double* norm_offsites ;
    double norm_offsite = 0;
    uint genenum = 0;
    
private:
    
    ktable *tablefile;
    std::string StrLine;
    vector<uint16> row;
    
};

#endif /* Kmatrix_hpp */
