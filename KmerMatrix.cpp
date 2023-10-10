//
//  KmerMatrix.cpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#include "KmerMatrix.hpp"

#include <fstream>
#include <string>
#include <unordered_set>

extern bool optioncorr;

inline void getEachRowValue(const FLOAT_T depth, const int count, const char sign, const uint16 uniqcounter,const uint16 flag,const FLOAT_T mean_repeat,FLOAT_T &total_lambda, FLOAT_T &norm_value, FLOAT_T &weight_value)
{
    float ori_weight = 1.0;
    if (sign==1 && uniqcounter == 1)
    {
        ori_weight = 0.05 ;
    }
    
    if (count < 3 )
    {
        if (mean_repeat > 1.0) weight_value += (ori_weight/mean_repeat - ori_weight);
        return ;
    }
    
    
    float count_f;           //float copy number value
    float count_i;   //estimate number of copy and at least one copy
    float new_weight;
    
    count_f = 1.0 * count/depth  ;
    
    if (optioncorr && (flag & 0x3F) > errorcutoff1 )
    {
        float corr = 0.01 * (flag & 0xFFC0)/64;
        count_f *= corr;
    }
    
    count_i = (int(count_f+ 0.5) > 1.0) ? int(count_f+ 0.5) : 1.0;   //estimate number of copy and at least one copy
    
    new_weight =  1.00/(count_i*count_i * mean_repeat);
    
    norm_value += count_i * new_weight ;
    weight_value += new_weight - ori_weight;
    //weight_value += (new_weight - ori_weight);  //weight is reversely proportion to square of estimated copy number, we calculate offsite to the original weight
               //norm vector value of this kmer = count_i * weight = 1.00/count_i
        
    total_lambda += count_f;
}

//A function to change values of norm vec and norm matrix for a list of kmers found in exactly the same list of samples
//This version is for binaray major kmers (frequencies > 50% and no sample has more than 1)
//We reverse calculate values for samples that missing this kmer
inline void getEachRowNorm_major(const uint16 rowsize, const uint16 *row, const uint16 gnum, FLOAT_T* norm_vec, FLOAT_T* norm_matrix, const FLOAT_T norm_value, const FLOAT_T  weight_value, FLOAT_T  &vec_offsite, FLOAT_T&matrix_offsite, FLOAT_T* row_offsites)
{
    
    vec_offsite += norm_value;            //major kmer, default is add this weight change to whole norm vector
    for (int i = 0; i < rowsize; ++i)
    {
        norm_vec[row[i]] -= norm_value;
    }
    
    if (weight_value == 0)
    {
        return;
    }
    
    matrix_offsite += weight_value;            //major kmer, default is add this weight change to whole matrix
    
    for (int i = 0; i < rowsize; ++i)            //sample missing, ignored its row/col from the weight change
    {
        norm_matrix[row[i]*gnum+row[i]] -= weight_value * 0.5;   //we will double diagonal elements after
        row_offsites[row[i]] -= weight_value;               //row[i] sample is missing this kmer, should not be affected, use offsite to correct back
        
        for (int j = i+1; j < rowsize; ++j)
        {
            norm_matrix[row[i]*gnum+row[j]] += weight_value;           //we use add to replace matrix multiplication because most kmers only found once in a locus, when row[i] == row[j], it should be doubled, so we will double diagonal elements after
        }
    }
    
}

//This version is for minor or non-binaray major kmers
inline void getEachRowNorm_minor(const uint16 rowsize, const uint16 *row, const uint16 gnum, FLOAT_T* norm_vec, FLOAT_T* norm_matrix, const FLOAT_T norm_value, const FLOAT_T  weight_value)
{
    for (int i = 0; i < rowsize; ++i)
    {
        norm_vec[row[i]] += norm_value;
    }
    
    if (weight_value == 0)
    {
        return;
    }

    for (int i = 0; i < rowsize; ++i)
    {
        norm_matrix[row[i]*gnum+row[i]] += weight_value/2;
        
        for (int j = i+1; j < rowsize; ++j)
        {
            norm_matrix[row[i]*gnum+row[j]] += weight_value;
        }
    }

}


//Add offsites back to norm vector and norm matrix
inline void AddOffsites(FLOAT_T *norm_vec, FLOAT_T *norm_matrix, const FLOAT_T  vec_offsite, const FLOAT_T  matrix_offsite, const FLOAT_T  *row_offsites, const FLOAT_T *diag_offsites, uint16 gnum)
{
    
    for (int i = 0; i < gnum; ++i)
    {
        norm_vec[i] += vec_offsite;                     //offsite for norm vector
        
        norm_matrix[i*gnum+i] *= 2;                     //this is doubling diagonal mentioned above
        norm_matrix[i*gnum+i] -= diag_offsites[i];       //we only want to double offsites, but prior values also doubled, should be changed back
        norm_matrix[i*gnum+i] += matrix_offsite;         //offsite for every cell
        
        for (int j = i + 1; j < gnum; ++j)
        {
            
            norm_matrix[i*gnum+j] += row_offsites[i];     //offsite for ith sample
            norm_matrix[i*gnum+j] += row_offsites[j];      //offsite for jth sample
            norm_matrix[i*gnum+j] += matrix_offsite;        //offsite for every cell
            
            norm_matrix[j*gnum + i] = norm_matrix[i*gnum+j];       //square matrix is symmetric
        }
    }
        
}


void KmerMatrix::getNorm(const uint16* kmervec, const uint16* kmermatrix, const FLOAT_T depth, const uint16 gnum, const uint knum, FLOAT_T* norm_vec, FLOAT_T* norm_matrix, FLOAT_T  &total_lambda)
{
    //memset(norm_matrix , 0,sizeof(FLOAT_T)*gnum*gnum);
    
    for (size_t i = 0; i < gnum; ++i)
    {
        row_offsites.get()[i] = 0;
        diag_offsites.get()[i] = norm_matrix[i * gnum + i];
    }
    
    FLOAT_T matrix_offsite = 0, vec_offsite = 0;
    
    FLOAT_T norm_value = 0.0, weight_value = 0.0, mean_repeat = 1.0, mean_repeat_this = 1.0;
    
    uint16 uniqcounter=0, j =0;
    
    float weight = 1.0;
    uint16 lastsize = 0;
    uint16 lastsign = 0;
    const uint16* lastnorm = NULL;
    const uint16* rowdata = kmermatrix;
    
    for (size_t i = 0; i < knum ; ++i)
    {
       
        switch (rowdata[0])
        {
            case '_':
                
                getEachRowValue(depth, kmervec[i], 0, 2,  rowdata[5],1.0 ,total_lambda, norm_value, weight_value);
                break;
            case '=':
                
                if (rowdata[2] && (rowdata[3] || rowdata[4] >= 30 ))
                {
                    mean_repeat_this = 1.0;
                }
                else
                {
                    mean_repeat_this = mean_repeat;
                }
                
                
                getEachRowValue(depth, kmervec[i], 1, uniqcounter, rowdata[5], mean_repeat_this, total_lambda, norm_value, weight_value);
                break;
                
            case '-':
                
                if (lastsign == '-')
                {
                    getEachRowNorm_major(lastsize, lastnorm, gnum, norm_vec,  norm_matrix, norm_value, weight_value,
                                         vec_offsite, matrix_offsite, row_offsites.get());
                }
                else if (lastsign == '+')
                {
                    getEachRowNorm_minor(lastsize, lastnorm, gnum, norm_vec,  norm_matrix, norm_value, weight_value);
                }
                
                lastsize = rowdata[1];
                lastnorm = &rowdata[FIXCOL];
                lastsign = rowdata[0];
                norm_value = 0.0;
                weight_value = 0.0;
                
                getEachRowValue(depth, kmervec[i], 0,  2, rowdata[5], 1.0, total_lambda, norm_value, weight_value);
                break;
            case'+':
                
                if (lastsign == '-')
                {
                    getEachRowNorm_major(lastsize, lastnorm, gnum, norm_vec,  norm_matrix, norm_value, weight_value,
                                         vec_offsite, matrix_offsite, row_offsites.get());
                }
                else if (lastsign == '+')
                {
                    getEachRowNorm_minor(lastsize, lastnorm, gnum, norm_vec,  norm_matrix, norm_value, weight_value);
                }
                
                lastsize = rowdata[1];
                lastnorm = &rowdata[FIXCOL];
                lastsign = rowdata[0];
                norm_value = 0.0;
                weight_value = 0.0;
                
                uniqcounter = 1;
                for (j = FIXCOL + 1 ; j < FIXCOL + rowdata[1]  ; ++j)
                {
                    if ( rowdata[j] != rowdata[j-1] ) ++ uniqcounter;
                }
                mean_repeat = MAX(1.0, 1.0 * rowdata[1])/uniqcounter;
                mean_repeat *= mean_repeat  ;
                if (rowdata[2] && (rowdata[3] || rowdata[4] >= 30 ))
                {
                    mean_repeat_this = 1.0;
                }
                else
                {
                    mean_repeat_this = mean_repeat;
                }
                
                
                getEachRowValue(depth, kmervec[i], 1, uniqcounter, rowdata[5], mean_repeat_this, total_lambda, norm_value, weight_value);
                break;
            default:
                break;
        }
                    
        rowdata = &rowdata[rowdata[1] + FIXCOL];
        
    }
    
    if (lastsign == '-')
    {
        getEachRowNorm_major(lastsize, lastnorm, gnum, norm_vec,  norm_matrix, norm_value, weight_value,
                             vec_offsite, matrix_offsite, row_offsites.get());
    }
    else if (lastsign == '+')
    {
        getEachRowNorm_minor(lastsize, lastnorm, gnum, norm_vec,  norm_matrix, norm_value, weight_value);
    }
    
    
    AddOffsites(norm_vec, norm_matrix, vec_offsite, matrix_offsite, row_offsites.get(), diag_offsites.get(), gnum);
    
}

