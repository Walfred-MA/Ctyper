//
//  TreeRound.hpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#ifndef TreeRound_hpp
#define TreeRound_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <tuple>
#include <mutex>
#include <memory>

using namespace std;

#include "config.hpp"

class node
{
public:
    node(node* p=NULL, int n=0, float d =0.0):index(n), dist(d), parent(p)
    {};
    
    FLOAT_T *coef(FLOAT_T* unround_coefs, FLOAT_T* unround_coefs_2);
    int *round(int* round_coefs, int* round_coefs_2);
    
    float leaveto_root(FLOAT_T* unround_coefs, int * round_coefs, FLOAT_T *non_leaves_unrounds, int *non_leaves_rounds);
    void rootto_leave(FLOAT_T* unround_coefs, int * round_coefs, FLOAT_T *non_leaves_unrounds, int *non_leaves_rounds);
    
    void clear();
    node* add(node * child);
  
    int index = 0;
    float dist = 0.0;
    
    node* parent;
    node* children[2];
    int8_t numchildren = 0;
};

class TreeRound
{
public:
    TreeRound()
    {};
    ~TreeRound()
    {};
    void Run(const node* nodes, FLOAT_T* ori_coefs, size_t size, int * round_coefs);
    
private:
    
    unique_ptr<FLOAT_T > leaves_unrounds = unique_ptr<FLOAT_T >( new FLOAT_T [MAX_UINT16] ) ;
    unique_ptr<FLOAT_T > non_leaves_unrounds = unique_ptr<FLOAT_T >( new FLOAT_T [MAX_UINT16] ) ;
    unique_ptr<int > non_leaves_rounds = unique_ptr<int >( new int [MAX_UINT16] ) ;
    
};



#endif /* TreeRound_hpp */
