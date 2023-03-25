//
//  KmerMatrix.cpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#include "KmerMatrix.hpp"

inline void getEachRowValue(const float depth, const int count, const char sign, const uint16 rowsize, double &total_lambda, double &norm_value, double &weight_value)
{
    
    float ori_weight = 1.0;
    if (sign==1 && rowsize == 1)
    {
        ori_weight = 0.05 ;
    }
    
    float count_f;           //float copy number value
    float count_i;   //estimate number of copy and at least one copy
    float new_weight;
    if (count <= 3)
    {
        count_f = 0.0;
        count_i = 1.0;
        new_weight = ori_weight;
        
    }
    else
    {
        count_f = 1.0 * count/depth;
        count_i = (int(count_f+ 0.5) > 1.0) ? int(count_f+ 0.5) : 1.0;   //estimate number of copy and at least one copy
        new_weight =  1.00/(count_i*count_i);
        norm_value += 1.00/count_i;
    }
    
    weight_value += new_weight - ori_weight;
    //weight_value += (new_weight - ori_weight);  //weight is reversely proportion to square of estimated copy number, we calculate offsite to the original weight
               //norm vector value of this kmer = count_i * weight = 1.00/count_i
    
    total_lambda += count_f;
}

//A function to change values of norm vec and norm matrix for a list of kmers found in exactly the same list of samples
//This version is for binaray major kmers (frequencies > 50% and no sample has more than 1)
//We reverse calculate values for samples that missing this kmer
inline void getEachRowNorm_major(const uint16 rowsize, const uint16 *row, const uint16 gnum, FLOAT_T* norm_vec, FLOAT_T* norm_matrix, const double norm_value, const double  weight_value, double  &vec_offsite, double &matrix_offsite, double* row_offsites)
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
inline void getEachRowNorm_minor(const uint16 rowsize, const uint16 *row, const uint16 gnum, FLOAT_T* norm_vec, FLOAT_T* norm_matrix, const double norm_value, const double  weight_value)
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
inline void AddOffsites(FLOAT_T *norm_vec, FLOAT_T *norm_matrix, const double  vec_offsite, const double  matrix_offsite, const double  *row_offsites, const double *diag_offsites, uint16 gnum)
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


void KmerMatrix::getNorm(const uint16* kmervec, const uint16* kmermatrix, const float depth, const uint16 gnum, const uint knum, FLOAT_T* norm_vec, FLOAT_T* norm_matrix, double  &total_lambda)
{
        
    memset(row_offsites.get(), 0, sizeof(double) * gnum);
    
    for (size_t i = 0; i < gnum; ++i)
    {
        diag_offsites.get()[i] = norm_matrix[i * gnum + i];
    }
    
    double matrix_offsite = 0, vec_offsite = 0;
    
    double norm_value = 0.0, weight_value = 0.0;
    
    uint16* rowdata = (uint16*) kmermatrix;
    for (size_t i = 0; i < knum; ++i)
    {
        
        switch (rowdata[0])
        {
            case '_':
                getEachRowValue(depth, kmervec[i], 0, rowdata[1], total_lambda, norm_value, weight_value);
                break;
            case '=':
                getEachRowValue(depth, kmervec[i], 1, rowdata[1], total_lambda, norm_value, weight_value);
                break;
            case '-':
                getEachRowValue(depth, kmervec[i], 0, rowdata[1], total_lambda, norm_value, weight_value);
                
                getEachRowNorm_major(rowdata[1], &rowdata[2], gnum, norm_vec,  norm_matrix, norm_value, weight_value,
                                     vec_offsite, matrix_offsite, row_offsites.get());
                norm_value = 0.0;
                weight_value = 0.0;
                break;
            case '+':
                getEachRowValue(depth, kmervec[i], 1, rowdata[1], total_lambda, norm_value, weight_value);
                
                getEachRowNorm_minor(rowdata[1], &rowdata[2], gnum, norm_vec,  norm_matrix, norm_value, weight_value);
                
                norm_value = 0.0;
                weight_value = 0.0;
                break;
            default:
                
                break;
        }
        
        rowdata = &rowdata[rowdata[1] + 2];
        
    }
    
    AddOffsites(norm_vec, norm_matrix, vec_offsite, matrix_offsite, row_offsites.get(), diag_offsites.get(), gnum);
    
    
}
