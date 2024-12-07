#ifndef config_hpp
#define config_hpp

#include <stdio.h>
#include <Eigen/Dense>
#include <atomic>

typedef unsigned int uint;
typedef unsigned __int128 u128;
typedef unsigned long long ull;
typedef unsigned short  uint16;
typedef unsigned char  uint8;
typedef unsigned long long ull_atom;

//large_prime will determine the hash space and memory usage: recommend use 536870909, 1073741827, or 2147483647
#define large_prime 2147483647
#define prime10M 9999991
#define MAX_UINT16 65535
#define MAX_LINE 10000000
#define Comb2( size ) (size + 1) * size / 2
#define MIN( A , B ) ( A <= B ) ? A : B
#define MAX( A , B ) ( A >= B ) ? A : B
#define MAXABS( A , B ) ( abs(A) >= abs(B) ) ? A : B
#define spair std::pair<std::string,std::string>
typedef double FLOAT_T;
typedef Eigen::MatrixXd Matrix_T ;
typedef Eigen::VectorXd Vector_T ;

#define FLOAT_T double
#define FIXCOL 6

#define genomesize_mean 6320012150.0
#define genomesize_male 6270605410.0
#define genomesize_female 6320012150.0

extern bool optioncorr;

#define errorcutoff1 7
#define errorcutoff2 11


#define sufficient 1000
#define corrstartpoint 0.3

#define windowmerge 15
#define minwindowcutoff 3
#define windowunit 30

#endif /* config_hpp */
