//
//  KmerMatrix.cpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#include "KmerMatrix.hpp"

inline void getEachRowValue(const float depth, const int count, const char sign, const uint16 rowsize, float &total_lambda, float &norm_value, float &weight_correct)
{
    if (count <= 3) return;
    
    float ori_weight = 1.0;
    if (sign==1 && rowsize == 1)
    {
        ori_weight = 0.05 ;
    }
    
    float count_f = count/depth;           //float copy number value
    float count_i = (int(count_f+ 0.5) > 1.0) ? int(count_f+ 0.5) : 1.0;   //estimate number of copy and at least one copy
    float new_weight = 1.00/(count_i*count_i);
    
    weight_correct += (0.00 - ori_weight);
    //weight_correct += (new_weight - ori_weight);  //weight is reversely proportion to estimated copy number, we calculate offsite to the original weight
    norm_value += 1.00/count_i;             //norm vector value of this kmer = count_i * weight = 1.00/count_i
    
    total_lambda += count_f;
}

//A function to change values of norm vec and norm matrix for a list of kmers found in exactly the same list of samples
//This version is for binaray major kmers (frequencies > 50% and no sample has more than 1)
//We reverse calculate values for samples that missing this kmer
inline void getEachRowNorm_major(const uint16 rowsize, const uint16 *row, const uint16 gnum, float* &norm_vec, float* &norm_matrix,float &vec_offsite, float &matrix_offsite, float* row_offsites, const float norm_value, const float weight_correct)
{
    
    vec_offsite += norm_value;            //major kmer, default is add this weight change to whole norm vector
    for (int i = 0; i < rowsize; ++i)
    {
        norm_vec[row[i]] -= norm_value;
    }
    
    if (weight_correct == 0)
    {
        return;
    }
    
    matrix_offsite += weight_correct;            //major kmer, default is add this weight change to whole matrix
    
    for (int i = 0; i < rowsize; ++i)            //sample missing, ignored its row/col from the weight change
    {
        norm_matrix[row[i]*gnum+row[i]] -= weight_correct/2;   //we will double diagonal elements after
        row_offsites[row[i]] -= weight_correct;               //row[i] sample is missing this kmer, should not be affected, use offsite to correct back
        
        for (int j = i+1; j < rowsize; ++j)
        {
            norm_matrix[row[i]*gnum+row[j]] += weight_correct;           //we use add to replace matrix multiplication because most kmers only found once in a locus, when row[i] == row[j], it should be doubled, so we will double diagonal elements after
        }
    }
    
}

//This version is for minor or non-binaray major kmers
inline void getEachRowNorm_minor(const uint16 rowsize, const uint16 *row, const uint16 gnum, float* &norm_vec, float* &norm_matrix, const float norm_value, const float weight_correct)
{
    for (int i = 0; i < rowsize; ++i)
    {
        norm_vec[row[i]] += norm_value;
    }
    
    if (weight_correct == 0)
    {
        return;
    }

    for (int i = 0; i < rowsize; ++i)
    {
        norm_matrix[row[i]*gnum+row[i]] += weight_correct/2;
        
        for (int j = i+1; j < rowsize; ++j)
        {
            norm_matrix[row[i]*gnum+row[j]] += weight_correct;
        }
    }

}


//Add offsites back to norm vector and norm matrix
inline void AddOffsites(float *norm_vec, float *norm_matrix, const float vec_offsite, const float matrix_offsite, const float *diag_offsites, const float *row_offsites, uint16 gnum)
{
    
    for (int i = 0; i < gnum; ++i)
    {
        norm_vec[i] += vec_offsite;                     //offsite for norm vector
        
        norm_matrix[i*gnum+i] -= diag_offsites[i];       //we only want to double offsites, but prior values also doubled, should be changed back
        norm_matrix[i*gnum+i] *= 2;                     //this is doubling diagonal mentioned above
        norm_matrix[i*gnum+i] += diag_offsites[i];       //we only want to double offsites, but prior values also doubled, should be changed back
        norm_matrix[i*gnum+i] += matrix_offsite;         //offsite for every cell
        
        for (int j = i + 1; j < gnum; ++j)
        {
            
            norm_matrix[i*gnum+j] += row_offsites[i];     //offsite for ith sample
            //norm_matrix[i*gnum+j] += row_offsites[j];      //offsite for jth sample
            norm_matrix[i*gnum+j] += matrix_offsite;        //offsite for every cell
            
            norm_matrix[j*gnum + i] = norm_matrix[i*gnum+j];       //square matrix is symmetric
        }
    }
    
    
}


void KmerMatrix::getNorm(const uint16* kmervec, const uint16* kmermatrix, const float depth, const uint16 gnum, const uint knum, float* norm_vec, float* norm_matrix, float &total_lambda)
{
    
    cout<<"weight"<<endl;
    
    memset(row_offsites.get(), 0, sizeof(float) * MAX_UINT16);
    
    for (size_t i = 0; i < gnum; ++i)
    {
        diag_offsites.get()[i] = norm_matrix[i * gnum + i];
    }
    
    float matrix_offsite = 0, vec_offsite = 0;
    
    float norm_value = 0.0, weight_correct = 0.0;
    
    uint16* rowdata = (uint16*) kmermatrix;
    for (size_t i = 0; i < knum; ++i)
    {
        switch (rowdata[0])
        {
            case '_':
                getEachRowValue(depth, kmervec[i], 0, rowdata[1], total_lambda, norm_value, weight_correct);
                break;
            case '=':
                getEachRowValue(depth, kmervec[i], 1, rowdata[1], total_lambda, norm_value, weight_correct);
                break;
            case '-':
                getEachRowValue(depth, kmervec[i], 0, rowdata[1], total_lambda, norm_value, weight_correct);
                
                getEachRowNorm_major(rowdata[1], &rowdata[2], gnum, norm_vec,  norm_matrix, vec_offsite, matrix_offsite, row_offsites.get(), norm_value, weight_correct);
                norm_value = 0.0;
                weight_correct = 0.0;
                break;
            case '+':
                getEachRowValue(depth, kmervec[i], 1, rowdata[1], total_lambda, norm_value, weight_correct);
                
                getEachRowNorm_minor(rowdata[1], &rowdata[2], gnum, norm_vec,  norm_matrix, norm_value, weight_correct);
                
                norm_value = 0.0;
                weight_correct = 0.0;
                break;
            default:
                
                break;
        }
        
        rowdata = &rowdata[rowdata[1] + 2];
        
    }
    
    AddOffsites(norm_vec, norm_matrix, vec_offsite, matrix_offsite, row_offsites.get(), diag_offsites.get(), gnum);
    
}
