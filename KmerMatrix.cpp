//
//  KmerMatrix.cpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#include "KmerMatrix.hpp"

void KmerMatrix::getNormEachRow(const int count, const uint16 *row, const uint16 rowsize, const uint16 sign)
{
    if (count <= 3) return;
    
    double ori_weight = 1.0;
    if (sign==1 && rowsize == 1)
    {
        ori_weight = 0.05 ;
    }
    
    double count_f = count/depth;
    float count_i = (int(count_f+ 0.5) > 1.0) ? int(count_f+ 0.5) : 1.0;
    
    double weight_change = (1.00/(count_i*count_i)-ori_weight);
    kmer_counts += count_f;
    
    
    double value = weight_change;
    double vec_value = 1.00/count_i;
    
    if (sign == 0)
    {
        vec_offsite += vec_value;
        for (int i = 0; i < rowsize; ++i)
        {
            kernal_vec[row[i]] -= vec_value;
        }

    }
    else
    {
        for (int i = 0; i < rowsize; ++i)
        {
            kernal_vec[row[i]] += vec_value;
        }
    }
    
    if (weight_change == 0)
    {
        return;
    }

        
    if (sign == 0)               //this kmer is a common kmer and only sample missing this kmer included
    {
        norm_offsite += value;            //common kmer, default is add this weight change to all samples
        
        for (int i = 0; i < rowsize; ++i)            //sample missing, ignored from the chance weight
        {
            weightnorm[row[i]*genenum+row[i]] -= value/2;   //will be multiple 2 later
            norm_offsites[row[i]] -= value;
            
            for (int j = i+1; j < rowsize; ++j)
            {
                weightnorm[row[i]*genenum+row[j]] += value;           //correct for dounle counted index
            }
        }
    }
    
    else
    {
        for (int i = 0; i < rowsize; ++i)
        {
            weightnorm[row[i]*genenum+row[i]] += value/2;
            
            for (int j = i+1; j < rowsize; ++j)
            {
                weightnorm[row[i]*genenum+row[j]] += value;
            }
        }
    }
}

void KmerMatrix::InitiateMatrix()
{
    
    if (genenum > alloc_size || 2 * genenum < alloc_size)
    {
        weightnorm = (double*) realloc(weightnorm,sizeof(double*) * (genenum * genenum + 1));
        
        kernal_vec = (double*) realloc(kernal_vec,sizeof(double*) * (genenum+1) );
        
        norm_offsites = (double*) realloc(norm_offsites, sizeof (double) * (genenum+1) );
        
        alloc_size = genenum;
    }
    
    memset(weightnorm, 0, sizeof (double) *  (genenum * genenum + 1) );
    
    memset(kernal_vec, 0, sizeof (double) * (genenum+1) );
    
    memset(norm_offsites, 0, sizeof (double) * (genenum+1) );
    
    norm_offsite = 0;
    
    vec_offsite = 0;
    
    kmer_counts = 0 ;
    
}

void KmerMatrix::FinishMatrix(const double * prior)
{
    
    for (int i = 0; i < genenum; ++i)
    {
        kernal_vec[i] += vec_offsite;
        
        weightnorm[i*genenum+i] = prior[i*genenum+i];
        weightnorm[i*genenum+i] *= 2;
        weightnorm[i*genenum+i] += norm_offsite;
        
        for (int j = i + 1; j < genenum; ++j)
        {
            weightnorm[i*genenum+j] = prior[i*genenum+j];
            
            weightnorm[i*genenum+j] += norm_offsites[i];
            weightnorm[i*genenum+j] += norm_offsites[j];
            weightnorm[i*genenum+j] += norm_offsite;
            
            weightnorm[j*genenum + i] = weightnorm[i*genenum+j];
        }
    }
    
    
}

void KmerMatrix::getNorm(const PriorChunk* Data, const uint16* allkmervec, const float dep)
{
    
    genenum = (uint16) Data-> genenum;
    kmernum =  Data -> kmernum;
    kmervec = &allkmervec[Data ->kmervec_start];
    depth = dep;
    
    InitiateMatrix();
    
    const uint16* rowdata = Data -> kmer_matrix;
    for (size_t i = 0; i < kmernum; ++i)
    {
        const uint16 rowsize = rowdata[0] ;
        
        getNormEachRow(kmervec[i], &rowdata[2], rowsize - 1, rowdata[1]);
        
        rowdata = &rowdata[rowsize];
    }
    
    FinishMatrix(Data -> prior_norm);
    
}





/*
void KmerMatrix::RegressMatrix()
{
    Regression regress;
    
    regress.Call(weightnorm, kernal, genenum);
    RawResults = regress.coefs;
    residuel = regress.residuel;
    
}


void KmerMatrix::Round()
{
    treefile->nextLine(Phylotree);
    TreeRound rounder;
    
    rounder.Run(Phylotree,RawResults);
}

void KmerMatrix::ProcessEachMatrix()
{
    RegressMatrix();
    Round();
    write();
}


void KmerMatrix::ProcessMatrix()
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
            ProcessEachMatrix();
        }
    }
    
};



void KmerMatrix::Process()
{
    
    ProcessMatrix();
}




int KmerMatrix::LoadHeader(const char* line, int linesize)
{
    
    matrixinfo.emplace_back();
    curr_info = &matrixinfo.back();
    
    string& prefix = get<0>(*curr_info);
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
    
    prefix = std::string(line, line + startpos);
    
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


void KmerMatrix::LoadRow(const char* line, int linesize)
{
    rowsize = 0;
    
    int startpos = 0;
    int currnum = 0;
    for (; startpos < linesize; ++startpos)
    {
        if (line[startpos] == '\t') break;
        if (line[startpos] == ',')
        {
            row[rowsize++]=currnum;
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


void KmerMatrix::write()
{
    FILE *fwrite=fopen(outputfile, "w");
    
    if (fwrite==NULL)
    {
        std::cerr << "ERROR: Cannot write file: " << outputfile << endl;
        
        std::_Exit(EXIT_FAILURE);
    }
    
    for (int i = 0; i < matrixinfo.size(); ++i)
    {
        string name = get<0>(matrixinfo[i]);
        auto totalkmer = totalkmers[i];
        uint genenum = (uint)get<3>(matrixinfo[i]);
        double* kernal = kernals[i];
        double* weightnorm = weightnorms[i];
        
        
        fprintf(fwrite,">%s\t%.4lf\n", prefix.c_str(), totalkmer);
        
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

*/
