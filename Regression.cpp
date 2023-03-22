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
        s(passive_set[0]) = (double)y(passive_set[0]) / (double)A(passive_set[0], passive_set[0]);
        cout<<"divide:"<<std::setprecision(15)<< (double)y(passive_set[0])<<"," << (double)A(passive_set[0], passive_set[0])<<endl;
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
    
    //Vector_T s_ = A.inverse() * sub_vector;
    
    for (size_t i = 0; i < passive_num; ++i)
    {
        s(passive_set[i]) = s_(i);
    }
    
    /*
    cout <<"A: "<<",";
    for (uint16 i = 0; i < passive_num; ++i)
    {
        for (size_t j = 0; j < passive_num; ++j)
        {
            cout<<A(passive_set[i],passive_set[j] )<<",";
        }
        cout << endl;
    }
    
    cout <<"y: "<<",";
    for (uint16 i = 0; i < passive_num; ++i)
    {
        cout <<  y(passive_set[i]) <<",";
    }
    cout<<endl;
    
    cout <<"s: "<<",";
    for (uint16 i = 0; i < passive_num; ++i)
    {
        cout <<  s(passive_set[i]) <<",";
    }
    cout<<endl;
    
    */
}

int Regression::lawson_hanson_nnls(const FLOAT_T *kernal_vec, const FLOAT_T *weightnorm, uint16 size, FLOAT_T *coefs, FLOAT_T *residuel)
{
    int max_iterations = MIN(size, 1000);
    
    const FLOAT_T tol = size*numeric_limits<FLOAT_T>::epsilon();
    
    const Vector_T y = Eigen::Map<Vector_T> ((FLOAT_T*) kernal_vec, size);
    
    const Matrix_T A = Eigen::Map<Matrix_T> ((FLOAT_T*) weightnorm, size, size);
        
    memset(coefs, 0, sizeof(FLOAT_T) * size);
    Vector_T x = Eigen::Map<Vector_T>(coefs, size);
    
    memset(x_trial_vec.data(), 0, sizeof(FLOAT_T) * size);
    Vector_T x_trial = Eigen::Map<Vector_T>(x_trial_vec.data(), size);

    memset(residuel, 0, sizeof(FLOAT_T) * size);
    Vector_T r = Eigen::Map<Vector_T>(residuel, size);
    
    cout<<"y matrix"<<endl;
    for (uint16 i = 0; i < size; ++i)
    {
        cout << y(i) <<",";
    }
    cout<<endl;
    
    cout<<"A matrix"<<endl;
    for (uint16 i = 0; i < size; ++i)
    {
        for (uint16 j = 0; j < size; ++j)
        {
            cout << A(i, j) <<",";
        }
        cout<<endl;
    }
    
    
    
    
    
    uint16 no_update = 0;
    
    active_or_passive.assign(size,0); //0 for active, 1 for passive
    uint16 passive_num = 0;
    uint16 passive_num_old = 0;
    
    r = y - A * x;
    
    
    while (passive_num < size)
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
        
        cout<<"max: "<<max_r<<endl;
        
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
            
            cout<<"min: ";
            for (int i = 0; i < passive_num; ++i)
            {
                cout<<passive_set[i]<<",";
            }
            
            cout<<alpha<<endl;
            
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
            
            trial_solution(passive_set, passive_num, y, A, size, x_trial);
            
            cout<<"passive_set: ";
            for (int i = 0; i < passive_num; ++i)
            {
                cout<<passive_set[i]<<",";
            }
            cout<<endl;
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
    
    for (int i = 0; i < size; ++i)
    {
        cout<<x(i)<<",";
    }
    cout<<endl;

    
    
    cout<<"checklawson"<<endl;
    return 0;
}

 
/*
// Solve the least squares subproblem for the variables in the passive set P.
Matrix_T A_P(size, passive_num);
for (size_t i = 0; i < passive_num; ++i)
{
    A_P.col(i) = A.col(passive_set[i]);
}
Vector_T x_P_new = A_P.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(y);
Vector_T x_P_old = x;

for (size_t i = 0; i < passive_num; ++i)
{
    x[passive_set[i]] = x_P_new[i];
}

if ((x_P_new.array() >= 0).all())
{
    break;
}
*/

/*
int Regression::Regess2( float *kernal_vec, float *weightnorm, uint size, float * &coefs, float &residuel)
{
    //Threads_lock.lock();
    
    Py_Initialize();
    
    _import_array();
    
    //Threads_lock.unlock();
    
    NumpyMatrix *matrix = new NumpyMatrix(weightnorm, size);
    NumpyVector *vector = new NumpyVector(kernal_vec, size);
    
    
    
    
    PyObject *pArgs = PyTuple_New(2);
    PyTuple_SetItem(pArgs, 0, reinterpret_cast<PyObject*>(matrix->data));
    PyTuple_SetItem(pArgs, 1, reinterpret_cast<PyObject*>(vector->data));
    
    
    //PyRun_SimpleString("from scipy.optimize import lsq_linear,nnls\n");
    PyObject *pName = PyUnicode_FromString("scipy.optimize");
    PyObject *pModule = PyImport_Import(pName);
    PyObject *pFunc = PyObject_GetAttrString(pModule, "nnls");
    
    
    PyObject * Result = PyObject_CallObject(pFunc, pArgs);
    
    coefs = reinterpret_cast<float*>(PyArray_DATA(reinterpret_cast<PyArrayObject*>(PyTuple_GetItem(Result, 0))));
    residuel = PyFloat_AsDouble (PyTuple_GetItem(Result, 1));
    
    
    Py_DECREF (pName);
    Py_DECREF (pModule);
    Py_DECREF (pFunc);
    Py_DECREF (Result);
    Py_DECREF (pArgs);
    Py_DECREF (matrix);
    Py_DECREF (vector);
    
    Py_Finalize();
    
    return 0;
    
}


int Regression::Regess3( float *kernal_vec, float *weightnorm, uint size, float * &coefs, float &residuel)
{
    
    auto *mat = new nsNNLS::denseMatrix(size, size, weightnorm);
    
    auto *vec = new nsNNLS::vector (size, kernal_vec);
    
    auto *solver = new nsNNLS::nnls(mat, vec, MIN(size/2, 1000));
    
    solver->optimize();
    nsNNLS::vector* sol = solver->getSolution();
    //nsNNLS::vector* res = solver->getResidual();
    
    coefs = sol->getData();
    
    residuel = 0.0;
    unique_ptr<float> residuels = unique_ptr<float>(new float[size]);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j <= i; ++j)
        {
            residuels.get()[j] += weightnorm[i*size + j] * coefs[i];
        }
    }
    
    for (int i = 0; i < size; ++i)
    {
        residuel += ( kernal_vec[i] - residuels.get()[i] )  * ( kernal_vec[i] - residuels.get()[i] ) ;
    }
    
    return 0;
}


int Regression::Regess( float *kernal_vec, float *weightnorm, uint size, float * &coefs, float &residuel)
{
    
    auto *mat = new nsNNLS::denseMatrix(size, size, weightnorm);
    
    auto *vec = new nsNNLS::vector (size, kernal_vec);
    
    auto *solver = new nsNNLS::nnls(mat, vec, MIN(size/2, 1000));
    
    solver->optimize();
    nsNNLS::vector* sol = solver->getSolution();
    //nsNNLS::vector* res = solver->getResidual();
    
    coefs = sol->getData();
    
    residuel = 0.0;
    unique_ptr<float> residuels = unique_ptr<float>(new float[size]);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j <= i; ++j)
        {
            residuels.get()[j] += weightnorm[i*size + j] * coefs[i];
        }
    }
    
    for (int i = 0; i < size; ++i)
    {
        residuel += ( kernal_vec[i] - residuels.get()[i] )  * ( kernal_vec[i] - residuels.get()[i] ) ;
    }
    
    return 0;
}
 
 */

void Regression::Call(uint size,  FLOAT_T *kernal_vec, FLOAT_T *weightnorm, float total_lambda, FLOAT_T *unweightnorm, FLOAT_T * coefs, FLOAT_T * residuels)
{
    
    float *v1 = new float[size];
    float *m1 = new float[size*size];
    
    
    for (int i = 0; i < size; ++i)
    {
        v1[i] = kernal_vec[i];
    }
    
    for (int i = 0; i < size * size; ++i)
    {
        m1[i] = weightnorm[i];
    }
    
    cout<<size<<endl;
    
    lawson_hanson_nnls(kernal_vec, weightnorm, size, coefs, residuels);
    
    cout<<"check8.1"<<endl;
    
    float regressed_kmer = 0.0;
    
    for (int i = 0; i < size; ++i)
    {
        regressed_kmer += unweightnorm[i*size + i] * coefs[i];
    }
    
    cout<<"check8.2"<<endl;
    
    float correction = total_lambda/regressed_kmer;
    
    for (int i = 0; i < size; ++i)
    {
        coefs[i] *= correction;
        cout<<coefs[i]<<",";
    }
    cout<<endl;
    cout<<correction<<endl;
    
    cout<<"check8.3"<<endl;
     
}




/*
int Regression::test()
{
    uint size = 7000;
    float *v1 = new float[size];
    float *m1 = new float[size*size];
    float *m2 = new float[size*size];
    //memcpy(v1, v, 3*sizeof(float));
    //memcpy(&m1[9], m, 9*sizeof(float));
    
    
    for (int i = 0; i < size; ++i)
    {
        float x = MAX(rand()%10000 * 1+ rand()%10000 * 0.000002 * 0.5 * (rand()%5 - 2) * 0.5, 0 );
        float y = MAX(0, x - rand()%10000 * 0.02 + rand()%10000 * 0.0000002 * 0.5 * (rand()%5 - 2) * 0.5 );
        for (int j = 0; j <= i; ++j)
        {
            if (i==0 || j ==0)
            {
                m1[i*size + j ] = x ;
                m1[j*size + i ] = x ;
            }
            else
            {
                m1[i*size + j ] = y;
                m1[j*size + i ] = y;
            }
            //m1[i*200 + j ] += rand()%100 * 0.005 * (rand()%5 - 2);
            //m1[j*200 + i ] += rand()%100 * 0.005 * (rand()%5 - 2);
        }
        v1[i] = x ;
    }
    
    
    
    float * regression = NULL;
    float residuel = 0.0;
    
    cout<<"checkbegin"<<endl;
    auto begin = std::chrono::high_resolution_clock::now();
    Regess( v1, m1, size, regression, residuel);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    
    cout<<"time: "<<elapsed.count()* 1e-9 <<endl;
    cout<<"res: "<<residuel<<endl;
    
    float *sum = new float[size];
    memset(sum, 0, sizeof(double)*size);
    
    for (int i = 0; i < size; i++)
    {
        if (regression[i] > 0.001)
        {
            for (int j = 0; j <= i; ++j)
            {
                sum[j] += m1[i*size + j]* regression[i];
            }
        }
    }
    double allsum = 0.0;
    for (int i = 0; i < size; i++)
    {
        double sum0 = (v1[i]-sum[i]);
        sum0 *= (v1[i] - sum[i])/10000;
        if (!isnan(sum[i])) allsum += sum0;
    }
    
    cout<<"sum of resi:"<<allsum<<endl;
    
    for (int i = 0; i < size; i++)
    {
        if (regression[i] > 0.001)cout <<i <<":"<< regression[i] << ' ';
    }
       
    cout << endl;
    
    
    regression = NULL;
    residuel = 0.0;
    
    cout<<"checkbegin"<<endl;
    begin = std::chrono::high_resolution_clock::now();
    Regess( v1, m1, size, regression, residuel);
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    
    cout<<"time: "<<elapsed.count()* 1e-9 <<endl;
    cout<<"res: "<<residuel<<endl;
    
    memset(sum, 0, sizeof(double)*size);
    
    for (int i = 0; i < size; i++)
    {
        if (regression[i] > 0.001)
        {
            for (int j = 0; j <= i; ++j)
            {
                sum[j] += m1[i*size + j]* regression[i];
            }
        }
    }
    allsum = 0.0;
    for (int i = 0; i < size; i++)
    {
        double sum0 = (v1[i]-sum[i]);
        sum0 *= (v1[i] - sum[i])/10000;
        if (!isnan(sum[i])) allsum += sum0;
    }
    cout<<endl;
    
    
    for (int i = 0; i < size; i++)
    {
        if (regression[i] > 0.001)cout <<i <<":"<< regression[i] << ' ';
    }
       
    cout << endl;
    
    return 0;
}

*/
