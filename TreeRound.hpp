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
    
    pair<const node*, float> leaveto_root(FLOAT_T* reminder, int * rounds) const;
    //void rootto_leave(FLOAT_T* unround_coefs, int * round_coefs, FLOAT_T *non_leaves_unrounds, int *non_leaves_rounds);
    
    void clear();
    node* add(node * child);
  
    int index = 0;
    float dist = 0.0;
    float size = 0.0;
    
    node* parent=NULL;
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
    
    unique_ptr<FLOAT_T > reminder = unique_ptr<FLOAT_T >( new FLOAT_T [MAX_UINT16] ) ;
    
};



#endif /* TreeRound_hpp */
