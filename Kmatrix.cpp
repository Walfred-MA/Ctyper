//
//  Kmatrix.cpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 1/31/23.
//  Copyright © 2023 USC_Mark. All rights reserved.
//

#include "Kmatrix.hpp"


void Kmatrix::eachRowNorm(int count)
{
    
    if (count <= 3) return;
    
    double ori_weight = 1.0;
    if (sign==1 && row.size() == 1)
    {
        ori_weight = 0.05 ;
    }
    
    double count_f = count/depth;
    float count_i = (int(count_f+ 0.5) > 1.0) ? int(count_f+ 0.5) : 1.0;
    
    double weight_change = (1.00/(count_i*count_i)-ori_weight);
    uint16 genenum = get<3>(*curr_info);
    (*totalkmer) += count_f;
    
    
    double value = weight_change;
    double vec_value = 1.00/count_i;
    
    if (sign < 0)
    {
        vec_offsite += vec_value;
        for (uint16 gindex: row)
        {
            kernal[gindex] -= vec_value;
        }
    }
    else
    {
        for (uint16 gindex: row)
        {
            kernal[gindex] += vec_value;
        }
    }
    
    if (weight_change == 0)
    {
        return;
    }

        
    if (sign < 0)               //this kmer is a common kmer and only sample missing this kmer included
    {
        norm_offsite += value;            //common kmer, default is add this weight change to all samples
        
        for (int i = 0; i < row.size(); ++i)            //sample missing, ignored from the chance weight
        {
            weightnorm[row[i]*genenum+row[i]] -= value/2;   //will be multiple 2 later
            norm_offsites[row[i]] -= value;
            
            for (int j = i+1; j < row.size(); ++j)
            {
                weightnorm[row[i]*genenum+row[j]] += value;           //correct for dounle counted index
            }
        }
    }
    
    else
    {
        for (int i = 0; i < row.size(); ++i)
        {
            weightnorm[row[i]*genenum+row[i]] += value/2;
            
            for (int j = i+1; j < row.size(); ++j)
            {
                weightnorm[row[i]*genenum+row[j]] += value;
            }
        }

    }
    
    
}

void Kmatrix::InitiateMatrix()
{
    weightnorm = (double*) malloc(sizeof(double*) * (genenum * genenum + 1));
    memset(weightnorm, 0, sizeof (double) *  (genenum * genenum + 1) );
    weightnorms.push_back(weightnorm);
    
    kernal = (double*) malloc(sizeof(double*) * (genenum+1) );
    memset(kernal, 0, sizeof (double) * (genenum+1) );
    kernals.push_back(kernal);
    
    totalkmers.push_back(0.0);
    totalkmer = &totalkmers.back();
    
    
    norm_offsites = (double*) malloc(sizeof (double) * (genenum+1) );
    memset(norm_offsites, 0, sizeof (double) * (genenum+1) );
    norm_offsite = 0;
    vec_offsite = 0;
}

void Kmatrix::FinishMatrix()
{
    
    for (int i = 0; i < genenum; ++i)
    {
        kernal[i] += vec_offsite;
        weightnorm[i*genenum+i] *= 2;
        weightnorm[i*genenum+i] += norm_offsite;
        
        for (int j = i + 1; j < genenum; ++j)
        {
            weightnorm[i*genenum+j] += norm_offsites[i];
            weightnorm[i*genenum+j] += norm_offsites[j];
            weightnorm[i*genenum+j] += norm_offsite;
            
            weightnorm[j*genenum + i] = weightnorm[i*genenum+j];
        }
    }
    
    free(norm_offsites);
    
}

void Kmatrix::ProcessMatrix()
{

    int linenum = 0;
    
    while (tablefile->nextLine(StrLine))
    {
        if (StrLine[0] == '#')
        {
            linenum = LoadHeader(StrLine.c_str(), (int)StrLine.length());
            InitiateMatrix();
            
            auto end = kmerindex + linenum;
                        
            for (; kmerindex < end; ++kmerindex)
            {
                if (!tablefile->nextLine(StrLine))
                {
                    std::cerr << "ERROR: error in kmer matrix file "<<std::endl;
                    std::_Exit(EXIT_FAILURE);
                    break;
                }
                
                LoadRow(StrLine.c_str(), (int) StrLine.size());
                
                eachRowNorm (kmervec[kmerindex]);
            }
                        
            FinishMatrix();
        }
    }
    
};

void Kmatrix::Process()
{
    
    ProcessMatrix();
}

int Kmatrix::LoadHeader(const char* line, int linesize)
{
    
    matrixinfo.emplace_back();
    curr_info = &matrixinfo.back();
    
    string& name = get<0>(*curr_info);
    auto& kmerstart = get<1>(*curr_info);
    auto& numkmers = get<2>(*curr_info);
    auto& numgenes = get<3>(*curr_info);
    
    kmerstart = kmerindex;
    
    numgenes = 0;
    numkmers = 0;
    int startpos = 0;
    char c;
    for (; startpos < linesize; ++startpos)
    {
        if (line[startpos] == '\t') break;
    }
    
    name = std::string(line, line + startpos);
    
    for (startpos = startpos + 1; startpos < linesize; ++startpos)
    {
        c = line[startpos];
        if (c == '\t') break;
        
        numkmers *= 10;
        numkmers += c - '0';
    }
    
    for (startpos = startpos + 1; startpos < linesize; ++startpos)
    {
        c = line[startpos];
        if (c == '\t' || c=='\0' || c =='\n') break;
        
        numgenes *= 10;
        numgenes += c - '0';
    }
    
    
    genenum = numgenes;
    
    return (int)numkmers;
}


void Kmatrix::LoadRow(const char* line, int linesize)
{
    row.clear();
    
    int startpos = 0;
    int currnum = 0;
    for (; startpos < linesize; ++startpos)
    {
        if (line[startpos] == '\t') break;
        if (line[startpos] == ',')
        {
            row.push_back(currnum);
            currnum = 0;
        }
        
        if (line[startpos] >= '0' && line[startpos] <= '9')
        {
            currnum *= 10;
            currnum += line[startpos] - '0';
        }
    }
    
    if (line[startpos+1] == '0')
    {
        sign = -1;
    }
    else
    {
        sign = 1;
    }

}

void Kmatrix::write(const char* outputfile, const char* prefix )
{
    FILE *fwrite=fopen(outputfile, "w");
    
    if (fwrite==NULL)
    {
        std::cerr << "ERROR: Cannot write file: " << outputfile << endl;
        
        std::_Exit(EXIT_FAILURE);
    }
    
    fprintf(fwrite,"@%s\n", prefix);
    for (int i = 0; i < matrixinfo.size(); ++i)
    {
        string name = get<0>(matrixinfo[i]);
        auto totalkmer = totalkmers[i];
        uint genenum = (uint)get<3>(matrixinfo[i]);
        double* kernal = kernals[i];
        double* weightnorm = weightnorms[i];
        
        
        fprintf(fwrite,">%s\t%.4lf\n", name.c_str(), totalkmer);
        
        for (int j = 0; j < genenum ; j++)
        {
            fprintf(fwrite,"%.4lf,", kernal[j]);
        }
        
        for (int j = 0; j < genenum * genenum ; j++)
        {
            if (j % genenum == 0) fprintf(fwrite,"\n");
            fprintf(fwrite,"%.4lf,", weightnorm[j]);
        }
        fprintf(fwrite,"\n");
        
    }
    
    fclose(fwrite);
    
    return ;
}

/*
void getSquare(vector<int>& row, vector<vector<float>> &kernel, vector<vector<float>> &kernel_unw, vector<float> &kmervector, vector<float> &kmervector_unw,float weight,float count_f)
{
    for (int i = 0; i < row.size(); ++i)
    {
        uint16 index1 = row[i];
        kernel[index1][index1] += weight/2;
        kernel_unw[index1][index1] += 0.5;
        kmervector[i] += weight*count_f;
        kmervector_unw[i] += count_f;
        
        for (int j = i+1; j < row.size(); ++i)
        {
            uint16 index2 = row[j];
            kernel[index1][index2] += weight;
            kernel_unw[index1][index2] += 1;
        }
    }
}
*/




