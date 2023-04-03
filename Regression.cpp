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
    int max_iterations = MIN(size, 1000);
    
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


void Regression::Call(uint size,  FLOAT_T *kernal_vec, FLOAT_T *weightnorm, float total_lambda, const uint *kmercounts, FLOAT_T * coefs, FLOAT_T * residuels)
{
    
    lawson_hanson_nnls(kernal_vec, weightnorm, size, coefs, residuels);
        
    float regressed_kmer = 0.0;
    
    for (int i = 0; i < size; ++i)
    {
        regressed_kmer += kmercounts[i] * coefs[i];
    }
    
    float correction = total_lambda/regressed_kmer;
    
    for (int i = 0; i < size; ++i)
    {
        coefs[i] *= correction;
    }
}

