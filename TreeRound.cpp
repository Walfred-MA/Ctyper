//
//  TreeRound.cpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#include "TreeRound.hpp"

void node::clear()
{
    index = 0;
    dist = 0.0;
    parent = NULL;
    numchildren = 0;
    
}

FLOAT_T *node::coef(FLOAT_T* unround_coefs, FLOAT_T* unround_coefs_2)
{
    if (index > 0 )
    {
        return &unround_coefs[abs(index) - 1];
    }
    
    else
    {
        return &unround_coefs_2[abs(index) - 1];
    }
}

int *node::round(int* unround_coefs, int* unround_coefs_2)
{
    if (index > 0 )
    {
        return &unround_coefs[abs(index) - 1];
    }
    
    else
    {
        return &unround_coefs_2[abs(index) - 1];
    }
}




node* node::add(node * child)
{
    child->clear();
    children[numchildren ++] = child;
    child->parent = this;
    
    return child;
}


float node::leaveto_root(FLOAT_T* unround_coefs, int * round_coefs, FLOAT_T *non_leaves_unrounds, int *non_leaves_rounds)
{
    FLOAT_T &coef = *(this->coef(unround_coefs, non_leaves_unrounds));
    int &round  = *(this->round(round_coefs, non_leaves_rounds));
    
    if (numchildren)
    {
        float total_dist = (1-(children[0]->dist+children[1]->dist));
        
        coef += total_dist * (children[0]->leaveto_root(unround_coefs, round_coefs, non_leaves_unrounds, non_leaves_rounds) + children[1]->leaveto_root(unround_coefs, round_coefs, non_leaves_unrounds, non_leaves_rounds));
        
        //total_coef +=  children[0]->total_coef + children[1]->total_coef;
        
        //total_round += children[0]->total_round + children[1]->total_round;
    }

    
    round = int(coef + 0.5);
    
    //total_round += round;
        
    coef -= round;
    
    return coef;
}

void node::rootto_leave(FLOAT_T* unround_coefs, int * round_coefs, FLOAT_T *non_leaves_unrounds, int *non_leaves_rounds)
{
    
    int &round  = *(this->round(round_coefs, non_leaves_rounds));
    
    if (round && numchildren)
    {
        node *prefer_leave;
        if (round > 0)
        {
            prefer_leave = (*children[1]->coef(unround_coefs, non_leaves_unrounds) - *children[1]->round(round_coefs, non_leaves_rounds) > *children[0]->coef(unround_coefs, non_leaves_unrounds) - *children[0]->round(round_coefs, non_leaves_rounds)) ?  children[1]: children[0];
        }
        else
        {
            prefer_leave = (*children[1]->coef(unround_coefs, non_leaves_unrounds) - *children[1]->round(round_coefs,non_leaves_rounds) < *children[0]->coef(unround_coefs, non_leaves_unrounds) - *children[0]->round(round_coefs,non_leaves_rounds)) ?  children[1]: children[0];
        }
        
        *prefer_leave->round(round_coefs, non_leaves_rounds) += round;
        //prefer_leave->total_round += round;
    }
    
    if (numchildren)
    {
        children[0]->rootto_leave(unround_coefs, round_coefs, non_leaves_unrounds, non_leaves_rounds);
        children[1]->rootto_leave(unround_coefs, round_coefs, non_leaves_unrounds, non_leaves_rounds);
    }
}

void TreeRound::Run(const node* nodes, FLOAT_T* unround_coefs, size_t size, int * round_coefs)
{
    memcpy(leaves_unrounds.get(), unround_coefs, sizeof (FLOAT_T) *  size);
    memset(non_leaves_unrounds.get(), 0, sizeof(FLOAT_T) * size );
    memset(non_leaves_rounds.get(), 0, sizeof(int) * size );

    
    node root = nodes[0];
    
    root.leaveto_root(leaves_unrounds.get(), round_coefs, non_leaves_unrounds.get(),  non_leaves_rounds.get());
    
    root.rootto_leave(leaves_unrounds.get(), round_coefs, non_leaves_unrounds.get(),  non_leaves_rounds.get());
    
    
}
