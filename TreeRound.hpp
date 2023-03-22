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

using namespace std;

#include "config.hpp"

class node
{
public:
    node(node* p=NULL, int n=0, float d =0.0, float c = 0.0):index(n), dist(d), coef(c), parent(p)
    {};
    
    node* push(int current_index);
    void build(string &newick,vector<node*>& allnodes);
    void leaveto_root();
    void rootto_leave();
    void clear();
    node* add(node * child);
  
    int index = 0;
    float dist = 0.0;
    float coef = 0.0;
    float total_coef = 0.0;
    float round = 0.0;
    float total_round = 0.0;
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
    void Run(node* nodes, float* ori_coefs, size_t size, int * round_coefs);
    
   
    
};



#endif /* TreeRound_hpp */
