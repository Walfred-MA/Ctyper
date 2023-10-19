//
//  Regression.cpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//
#include <Eigen/Dense>
#include <vector>
#include <limits>
#include <algorithm>
#include "Regression.hpp"
#include <iomanip>

#define max_repetitions 5
extern bool optioncorr;

void trial_solution(vector<uint16>& passive_set, const uint16 passive_num, const Vector_T &y, const Matrix_T &A, const uint16 size, Vector_T &s)
{
    assert( passive_num <= size );
    
    s.head(size).setZero();
    
    
    if (passive_num == 1)
    {
        s(passive_set[0]) = (FLOAT_T)y(passive_set[0]) / (FLOAT_T)A(passive_set[0], passive_set[0]);
        return ;
    }
     
    
    Matrix_T sub_matrix = Matrix_T(passive_num, passive_num);
    Vector_T sub_vector = Vector_T(passive_num);
    
    for (uint16 i = 0; i < passive_num; ++i)
    {
        sub_vector(i) = y(passive_set[i]);
        
        for (uint16 j = 0; j < passive_num; ++j)
        {
            sub_matrix(i, j) = A(passive_set[i], passive_set[j]);
        }
    }
    
    Vector_T s_ = sub_matrix.colPivHouseholderQr().solve(sub_vector);
        
    for (size_t i = 0; i < passive_num; ++i)
    {
        s(passive_set[i]) = s_(i);
    }
    
}

int Regression::lawson_hanson_nnls(const FLOAT_T *kernal_vec, const FLOAT_T *weightnorm, uint16 size, FLOAT_T *coefs, FLOAT_T *residuel)
{
    int max_iterations = MIN(2*size, 1000);
    
    const FLOAT_T tol = size*numeric_limits<FLOAT_T>::epsilon();
    
    const Vector_T y = Eigen::Map<Vector_T> ((FLOAT_T*) kernal_vec, size);
    
    const Matrix_T A = Eigen::Map<Matrix_T> ((FLOAT_T*) weightnorm, size, size);
        
    Vector_T x = Vector_T::Zero(size), x_trial(size) , r (size);
    
    uint16 no_update = 0;
    
    active_or_passive.assign(size,0); //0 for active, 1 for passive
    uint16 passive_num = 0;
    uint16 passive_num_old = 0;
    
    r = y - A * x;
    
    int j = 0;
    while (passive_num < size && j++ < max_iterations)
    {
        // Find the max residual max_r in active set and its index;
        
        int max_r_index = 0;
        FLOAT_T max_r = std::numeric_limits<FLOAT_T>::min();

        for (uint16 i = 0; i < size; ++i)
        {
            if (active_or_passive[i]) continue; //in passive, ignore
            
            if (r[i] > max_r)
            {
                max_r = r[i];
                max_r_index = i;
            }
        }
                
        if (max_r <= tol)  break; //non-negative achieved

        //make a copy for passive_set for late comparison
        passive_num_old = passive_num;
        passive_set_old.assign(passive_set.begin(), passive_set.begin() + passive_num_old);
        
        // Move the variable j from the passive set to active set .
        passive_set[passive_num++] = max_r_index;
        active_or_passive[max_r_index] = 1;
        
        trial_solution(passive_set, passive_num, y, A, size, x_trial); //solve || y - A * x || on passive vectors
        
        int k = 0;
        while (passive_num && k++ < max_iterations)
        {

            // Find the index i in the passive set P that has the smallest alpha value.
            FLOAT_T alpha = std::numeric_limits<FLOAT_T>::max();
            int alpha_index = -1;
            for (size_t i = 0; i < passive_num; ++i)
            {
                int idx = passive_set[i];
                FLOAT_T alpha_i = x_trial[idx] ;
                if (alpha_i < alpha)
                {
                    alpha = alpha_i;
                    alpha_index = idx;
                }
            }
            
                        
            if (alpha > tol)
            {
                break;
            }
                        
            int passive_shrink_index = 0;
            for (size_t i = 0; i < passive_num; ++i)
            {
                int idx = passive_set[i];
                if (x_trial[idx] > tol)
                {
                    passive_set[passive_shrink_index++] = passive_set[i];
                }
                else
                {
                    active_or_passive[max_r_index] = 0;
                }
            }
            
            passive_num = passive_shrink_index;
            
            // Update x.
            x = x + alpha * (x_trial - x);
            
            if (passive_num) trial_solution(passive_set, passive_num, y, A, size, x_trial);
            
        }

        x = x_trial.eval();
        // Update the residual.
        r = y - A * x;
        
        if(passive_num == passive_num_old  && equal(passive_set_old.begin(), passive_set_old.begin() + passive_num, passive_set.begin() ) )
        {
            no_update ++;
        }
        else
        {
            no_update = 0;
        }
        
        
        if (no_update > max_repetitions) break;
        
    }
    
    r = y - A * x;
    
    for (int i = 0; i < size; ++i)
    {
        coefs[i] = x(i);
        residuel[i] = r(i);
    }
    
    
    return 0;
}


inline void getlocalnum_mul(const FLOAT_T* coefs, const uint16* rowdata,const FLOAT_T totalnum, FLOAT_T &localnum, const vector<FLOAT_T>& grouptotalnums, vector<FLOAT_T>& grouplocalnums, const vector<uint16> &groups)
{
    
    uint16 geneindex = 0;
    
    switch (rowdata[0])
    {
        case '_': case '=':
            break;

        case '-':
                        
            localnum = totalnum;
            grouplocalnums = grouptotalnums;
            
            for (int i = 0 ; i < rowdata[1] ; i++)
            {
                geneindex = rowdata[ FIXCOL + i];
                localnum -= coefs[geneindex];
                grouplocalnums[groups[geneindex]] -= coefs[geneindex];
            }
            
            break;
            
        case '+':
            
            localnum = 0;
            std::fill(grouplocalnums.begin(), grouplocalnums.end(), 0);

            for (int i = 0 ; i < rowdata[1] ; i++)
            {
                geneindex = rowdata[ FIXCOL + i];
                localnum += coefs[geneindex];
                grouplocalnums[groups[geneindex]] += coefs[geneindex];
            }
            
            break;
        default:
            
            break;
    }
}

inline void GetMedianNofil_mul(const FLOAT_T* coefs, const uint16* kmervec, const uint16* rowdata, const uint knum, const FLOAT_T &totalnum, const vector<FLOAT_T> &grouptotalnums,const vector<uint16> &groups, const uint16 groupnum, vector<size_t> &grouptotalobs, vector<vector<FLOAT_T>> &allratios)
{
    
    FLOAT_T localnum = 0.0;
    vector<FLOAT_T> grouplocalnums (groupnum + 1, 0.0);
    
    for (size_t i = 0; i < knum; ++i)
    {
        getlocalnum_mul(coefs, rowdata, totalnum, localnum, grouptotalnums, grouplocalnums, groups);
        
        if (localnum >= 0.5 && kmervec[i] > 2)
        {
            grouplocalnums[groupnum] = localnum;
            
            FLOAT_T ratio = kmervec[i] / localnum;

            for (int i = 0 ; i < groupnum + 1; ++i)
            {
                if (grouplocalnums[i] >= 0.5 && grouptotalnums[i] >= 0.2)
                {
                    allratios[i][grouptotalobs[i]++] = ratio;
                }
            }
        }
        
        rowdata = &rowdata[rowdata[1] + FIXCOL];
    }
}


inline void GetMedianAttemp1_mul(const FLOAT_T* coefs, const uint16* kmervec, const uint16* rowdata, const uint knum, const FLOAT_T &totalnum, const vector<FLOAT_T> &grouptotalnums,const vector<uint16> &groups, const uint16 groupnum, vector<size_t> &grouptotalobs, vector<vector<FLOAT_T>> &allratios)
{
    
    FLOAT_T localnum = 0.0;
    vector<FLOAT_T> grouplocalnums (groupnum + 1, 0.0);
    
    for (size_t i = 0; i < knum; ++i)
    {
        getlocalnum_mul(coefs, rowdata, totalnum, localnum, grouptotalnums, grouplocalnums, groups);
        
        if (localnum >= 0.5 && kmervec[i] > 2)
        {
            grouplocalnums[groupnum] = localnum;
            
            if ( (rowdata[5] & 0x3F) < errorcutoff1 )
            {
                FLOAT_T ratio = kmervec[i] / localnum;
                
                if (optioncorr && (rowdata[5] & 0x3F) >= errorcutoff1 )
                {
                    float corr = 0.01 * (rowdata[5] & 0xFFC0)/64;
                    ratio *= corr;
                }

                for (int i = 0 ; i < groupnum + 1; ++i)
                {
                    if (grouplocalnums[i] >= 0.5 && grouptotalnums[i] >= 0.2)
                    {
                        allratios[i][grouptotalobs[i]++] = ratio;
                    }
                }
            }
        }
        
        rowdata = &rowdata[rowdata[1] + FIXCOL];
    }
}

inline void GetMedianAttemp2_mul(const FLOAT_T* coefs, const uint16* kmervec, const uint16* rowdata, const uint knum, const FLOAT_T &totalnum, const vector<FLOAT_T> &grouptotalnums,const vector<uint16> &groups, const uint16 groupnum, vector<size_t> &grouptotalobs, vector<bool>& oldtotalobs, vector<vector<FLOAT_T>> &allratios)
{
    
    FLOAT_T localnum = 0.0;
    vector<FLOAT_T> grouplocalnums (groupnum + 1, 0.0);
    
    for (size_t i = 0; i < knum; ++i)
    {
        
        getlocalnum_mul(coefs, rowdata, totalnum, localnum, grouptotalnums, grouplocalnums, groups);

        if ( localnum >= 0.5 && kmervec[i] > 2)
        {
            FLOAT_T ratio = kmervec[i] / localnum;
            grouplocalnums[groupnum] = localnum;
            
            if (optioncorr && (rowdata[5] & 0x3F) >= errorcutoff1 )
            {
                float corr = 0.01 * (rowdata[5] & 0xFFC0)/64;
                ratio *= corr;
            }
            
            if ( (rowdata[5] & 0x3F) < errorcutoff1)
            {
                
                for (int i = 0 ; i < groupnum + 1; ++i)
                {
                    if (oldtotalobs[i] == 0 && grouplocalnums[i] >= 0.5 && grouptotalnums[i] >= 0.2)
                    {
                        fill(allratios[i].begin() + grouptotalobs[i], allratios[i].begin() + grouptotalobs[i]+3, ratio);
                        grouptotalobs[i] += 3;
                    }
                }
            }
            
            else if ((rowdata[5] & 0x3F) < errorcutoff2)
            {
                for (int i = 0 ; i < groupnum + 1; ++i)
                {
                    if (oldtotalobs[i] == 0 && grouplocalnums[i] >= 0.5 && grouptotalnums[i] >= 0.2)
                    {
                        allratios[i][grouptotalobs[i]++] = ratio;
                    }
                }
            }
        }
        
        rowdata = &rowdata[rowdata[1] + FIXCOL];
        
    }
}

inline void GetMedianAttemp3_mul(const FLOAT_T* coefs, const uint16* kmervec, const uint16* rowdata, const uint knum, const FLOAT_T &totalnum, const vector<FLOAT_T> &grouptotalnums,const vector<uint16> &groups, const uint16 groupnum, vector<size_t> &grouptotalobs, vector<bool>& oldtotalobs, vector<vector<FLOAT_T>> &allratios)
{
    FLOAT_T localnum = 0.0;
    vector<FLOAT_T> grouplocalnums (groupnum + 1, 0.0);
    
        
    for (size_t i = 0; i < knum; ++i)
    {
        
        getlocalnum_mul(coefs, rowdata, totalnum, localnum, grouptotalnums, grouplocalnums, groups);

        if ( localnum >= 0.5 && kmervec[i] > 2)
        {
            FLOAT_T ratio = kmervec[i] / localnum;
            grouplocalnums[groupnum] = localnum;
            
            if (optioncorr && (rowdata[5] & 0x3F) >= errorcutoff1 )
            {
                float corr = 0.01 * (rowdata[5] & 0xFFC0)/64;
                ratio *= corr;
            }
            
            if ( (rowdata[5] & 0x3F) < errorcutoff1)
            {
                for (int i = 0 ; i < groupnum + 1; ++i)
                {
                    if (oldtotalobs[i] == 0 && grouplocalnums[i] >= 0.5 && grouptotalnums[i] >= 0.2)
                    {
                        fill(allratios[i].begin() + grouptotalobs[i], allratios[i].begin() + grouptotalobs[i]+6, ratio);
                        grouptotalobs[i] += 6;
                    }
                }
            }
            
            else if ((rowdata[5] & 0x3F) < errorcutoff2)
            {
                for (int i = 0 ; i < groupnum + 1; ++i)
                {
                    if (oldtotalobs[i] == 0 && grouplocalnums[i] >= 0.5 && grouptotalnums[i] >= 0.2)
                    {
                        fill(allratios[i].begin() + grouptotalobs[i], allratios[i].begin() + grouptotalobs[i]+2, ratio);
                        grouptotalobs[i] += 2;
                    }
                }
            
            }
            
            else
            {
                for (int i = 0 ; i < groupnum + 1; ++i)
                {
                    if (oldtotalobs[i] == 0 && grouplocalnums[i] >= 0.5 && grouptotalnums[i] >= 0.2)
                    {
                        allratios[i][grouptotalobs[i]++] = ratio;
                    }
                }
            }
        }
        
        rowdata = &rowdata[rowdata[1] + FIXCOL];
    }
}


void aggregateCorr_mul(FLOAT_T * coefs, const uint16* kmervec, const uint16* kmermatrix, const uint16 gnum, const uint knum, const vector<uint16> &groups, const uint16 groupnum, vector<FLOAT_T> &corrections)
{
    
    vector<size_t> grouptotalobs (groupnum + 1, 0);
    vector<FLOAT_T> grouptotalnums (groupnum + 1, 0.0);
    vector<vector<FLOAT_T>> allratios(groupnum + 1);
    
    FLOAT_T totalnum = 0.0;
    for (int i = 0; i < gnum; ++i)
    {
        totalnum += coefs[i];
        grouptotalnums[groups[i]] += coefs[i];
    }
    grouptotalnums[groupnum] = MAX(totalnum, 1.0);
    
    
    for (uint16 index = 0; index < groupnum + 1;++index)
    {
        if ( grouptotalnums[index] >= 0.2) allratios[index].resize(knum,0);
    }
    
    if (!optioncorr)
    {
        GetMedianNofil_mul(coefs, kmervec,kmermatrix,knum, totalnum,grouptotalnums,groups, groupnum, grouptotalobs, allratios);
    }
    else
    {
        GetMedianAttemp1_mul(coefs, kmervec,kmermatrix,knum, totalnum,grouptotalnums,groups, groupnum, grouptotalobs, allratios);
        
        vector<bool> oldtotalobs ( grouptotalobs.size() , 0 );
        size_t totalobs = sufficient;
        
        for (uint16 index = 0; index < groupnum + 1;++index)
        {
            if ( grouptotalnums[index] >= 0.2 && grouptotalobs[index] < sufficient)
            {
                totalobs = grouptotalobs[index];
                allratios[index].resize(10 * sufficient + knum,0);
            }
            else if (grouptotalobs[index] >= sufficient)
            {
                oldtotalobs[index] = 1;
            }
        }
        
        
        if ( totalobs < sufficient)
        {
            
            
            GetMedianAttemp2_mul(coefs, kmervec,kmermatrix,knum, totalnum,grouptotalnums,groups, groupnum, grouptotalobs, oldtotalobs, allratios);
            
            
            totalobs = 3 * sufficient;
            for (uint16 index = 0; index < groupnum + 1;++index)
            {
                if ( grouptotalnums[index] >= 0.2 && oldtotalobs[index] ==0 && grouptotalobs[index] < 3 * sufficient)
                {
                    totalobs = grouptotalobs[index];
                }
                else if (grouptotalobs[index] >= 3 * sufficient)
                {
                    oldtotalobs[index] = 1;
                }
            }
            
            if ( totalobs < 3 * sufficient)
            {
                GetMedianAttemp3_mul(coefs, kmervec,kmermatrix,knum, totalnum,grouptotalnums,groups, groupnum, grouptotalobs, oldtotalobs, allratios);
            }
        }
        
    }
    
    for (uint16 i = 0; i < groupnum + 1; ++i)
    {
        if (grouptotalnums[i] >= 0.2)
        {
            std::sort(allratios[i].begin(), allratios[i].begin() + grouptotalobs[i]);
            
            corrections[i] = allratios[i][grouptotalobs[i]/2];
        }
    }
}

inline void getlocalnum(const FLOAT_T* coefs, const uint16* rowdata,const FLOAT_T totalnum, FLOAT_T &localnum)
{
    switch (rowdata[0])
    {
        case '_': case '=':
            break;

        case '-':
            
            localnum = totalnum;
            for (int i = 0 ; i < rowdata[1] ; i++)
            {
                localnum -= coefs[rowdata[ FIXCOL + i]];
            }
            
            break;
            
        case '+':
            
            localnum = 0;
            for (int i = 0 ; i < rowdata[1] ; i++)
            {
                localnum += coefs[rowdata[ FIXCOL + i]];
            }
            
            break;
        default:
            
            break;
    }
}


inline void GetMedianNofil(const FLOAT_T* coefs, const uint16* kmervec, const uint16* rowdata, const uint knum, FLOAT_T &totalnum,size_t &totalobs, vector<FLOAT_T>& ratios)
{
    FLOAT_T localnum = 0.0;
    for (size_t i = 0; i < knum; ++i)
    {
        
        getlocalnum(coefs, rowdata, totalnum, localnum);
        
        if ( localnum >= 0.5 && kmervec[i] > 2)
        {
            FLOAT_T ratio = kmervec[i] / localnum;
                
            ratios[totalobs] = ratio;
            totalobs ++;
        }
        
        rowdata = &rowdata[rowdata[1] + FIXCOL];
    }
}


inline void GetMedianAttemp1(const FLOAT_T* coefs, const uint16* kmervec, const uint16* rowdata, const uint knum, FLOAT_T &totalnum,size_t &totalobs, vector<FLOAT_T>& ratios)
{
    FLOAT_T localnum = 0.0;
    for (size_t i = 0; i < knum; ++i)
    {
        
        getlocalnum(coefs, rowdata, totalnum, localnum);
        
        if ( localnum >= 0.5 && kmervec[i] > 2)
        {
            FLOAT_T ratio = kmervec[i] / localnum;
                
            if ( (rowdata[5] & 0x3F) < errorcutoff1)
            {
                ratios[totalobs] = ratio;
                totalobs ++;
            }
        }
        
        rowdata = &rowdata[rowdata[1] + FIXCOL];
    }
    
    
}

inline void GetMedianAttemp2(const FLOAT_T* coefs, const uint16* kmervec, const uint16* rowdata, const uint knum, FLOAT_T &totalnum, size_t &totalobs, vector<FLOAT_T>& ratios)
{
    FLOAT_T localnum = 0.0;
    
    for (size_t i = 0; i < knum; ++i)
    {
        
        getlocalnum(coefs, rowdata, totalnum, localnum);
        
        if ( localnum >= 0.5 && kmervec[i] > 2)
        {
            FLOAT_T ratio = kmervec[i] / localnum;
                
            if (optioncorr && (rowdata[5] & 0x3F) >= errorcutoff1 )
            {
                float corr = 0.01 * (rowdata[5] & 0xFFC0) / 64;
                ratio *= corr;
            }
                
            if ((rowdata[5] & 0x3F) < errorcutoff1)
            {
                fill(ratios.begin() + totalobs, ratios.begin() + totalobs+3, ratio);
                totalobs += 3;
            }
            else if ((rowdata[5] & 0x3F) < errorcutoff2)
            {
                ratios[totalobs] = ratio;
                totalobs ++;
            }
        }
        
        rowdata = &rowdata[rowdata[1] + FIXCOL];
    }
}

inline void GetMedianAttemp3(const FLOAT_T* coefs, const uint16* kmervec, const uint16* rowdata, const uint knum, FLOAT_T &totalnum,size_t &totalobs, vector<FLOAT_T>& ratios)
{
    FLOAT_T localnum = 0.0;
    for (size_t i = 0; i < knum; ++i)
    {
        getlocalnum(coefs, rowdata, totalnum, localnum);

        if ( localnum >= 0.5 && kmervec[i] > 2)
        {
            FLOAT_T ratio = kmervec[i] / localnum;
                
            if (optioncorr && (rowdata[5] & 0x3F) >= errorcutoff1 )
            {
                float corr = 0.01 * (rowdata[5] & 0xFFC0) / 64;
                ratio *= corr;
            }
                
            if ((rowdata[5] & 0x3F) < errorcutoff1)
            {
                fill(ratios.begin() + totalobs, ratios.begin() + totalobs+6, ratio);
                totalobs += 6;
            }
            else if ((rowdata[5] & 0x3F) < errorcutoff2)
            {
                fill(ratios.begin() + totalobs, ratios.begin() + totalobs+2, ratio);
                totalobs += 2;
            }
            else
            {
                ratios[totalobs] = ratio;
                totalobs ++;
            }
        }
        
        rowdata = &rowdata[rowdata[1] + FIXCOL];
    }
}

FLOAT_T aggregateCorr(const FLOAT_T * coefs, const uint16* kmervec, const uint16* kmermatrix, const uint16 gnum, const uint knum)
{
    
    vector<FLOAT_T> ratios(knum, 0.0);
        
    size_t totalobs  = 0;
    FLOAT_T totalnum = 0.0;
    for (int i = 0; i < gnum; ++i)
    {
        totalnum += coefs[i];
    }
    
    if (!optioncorr)
    {
        GetMedianNofil(coefs, kmervec, kmermatrix, knum, totalnum, totalobs, ratios);
    }
    else
    {
        GetMedianAttemp1(coefs, kmervec, kmermatrix, knum, totalnum, totalobs, ratios);
        
        if ( totalobs < sufficient)
        {
            ratios.resize(10 * sufficient + knum);
            GetMedianAttemp2(coefs, kmervec, kmermatrix, knum,totalnum, totalobs, ratios);
            
            if ( totalobs < 3 * sufficient)
            {
                GetMedianAttemp3(coefs, kmervec, kmermatrix, knum,totalnum, totalobs, ratios);
            }
        }
        
    }
    std::sort(ratios.begin(), ratios.begin() + totalobs);
    
    
    FLOAT_T median = ratios[totalobs/2];
    

    return median;
}


void Regression::Call(const uint16* kmervec, const uint16* kmermatrix, const FLOAT_T depth, const uint16 gnum, const uint knum, FLOAT_T* norm_vec, FLOAT_T* norm_matrix, FLOAT_T &total_lambda, const uint *kmercounts, FLOAT_T * coefs, FLOAT_T * residuels, const size_t numgroup, const vector<uint16> &groups)
{
    
    lawson_hanson_nnls(norm_vec, norm_matrix, gnum, coefs, residuels);
    
    /*
    float regressed_kmer = 0.0;
    
    for (int i = 0; i < gnum; ++i)
    {
        regressed_kmer += kmercounts[i] * coefs[i];
    }
    regressed_kmer = (regressed_kmer > 0 ) ? regressed_kmer : 1;
     
    float correction = total_lambda/regressed_kmer;
    */
    if (numgroup <= 1)
    {
        FLOAT_T correction = aggregateCorr(coefs, kmervec, kmermatrix, gnum, knum) / depth;
        
        for (int i = 0; i < gnum; ++i)
        {
            
            coefs[i] *= correction;
            
            residuels[i] *= correction;
            
        }
    }
    else
    {
        vector<FLOAT_T> corrections(numgroup + 1, 0.0);
        
        aggregateCorr_mul(coefs, kmervec, kmermatrix, gnum, knum, groups, numgroup, corrections);
        
        for (int index = 0; index < numgroup ; ++index)
        {
            corrections[index] = corrections[index]/depth;
        }

        FLOAT_T correction = corrections[numgroup]/depth;
        
        for (int i = 0; i < gnum; ++i)
        {
            
            if (corrections[groups[i]] > 0)
            {
                coefs[i] *= corrections[groups[i]];
                
                residuels[i] *= corrections[groups[i]];
            }
            else
            {
                coefs[i] *= correction;
                
                residuels[i] *= correction;
            }
        }
    }
}

