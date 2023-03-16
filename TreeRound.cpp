//
//  TreeRound.cpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#include "TreeRound.hpp"

node * node::push(int current_index)
{
    node * newnode = new node(this, current_index);
    children[numchildren++] = newnode;
    
    return newnode;
}

void node::clear()
{
    index = 0;
    dist = 0.0;
    coef = 0.0;
    total_coef = 0.0;
    round = 0.0;
    total_round = 0.0;
    parent = NULL;
    numchildren = 0;
}

node* node::add(node * child)
{
    child->clear();
    children[numchildren ++] = child;
    child->parent = this;
    
    return child;
}

void node::build(string & newick,vector<node*>& allnodes)
{
    node* current_node = this;
    double current_num = 0.0;
    int current_index = 0;
    float ifdeci = 0;
    
    if (newick.size()==0) return ;
    
    int notuselast = (newick[newick.size()-1] == ';') ;
    
    char c;
    for (int index=0; index < newick.size() - notuselast; ++index)
    {
        c = newick[index];
        switch(c)
        {
            case ' ': case '\n':
                break;
            case '(':
                current_node = current_node->push(current_index ++);
                allnodes.push_back(current_node);
                break;
            case ')': case ';':
                current_node->dist = current_num;
                ifdeci = 0;
                current_node = current_node->parent;
                break;
            case ',':
                current_node->dist = current_num;
                ifdeci = 0;
                current_node = current_node->parent->push(current_index ++);
                break;
            case ':':
                ifdeci = 1;
                break;
            case '.':
                ifdeci *= 0.1;
            default:
                if (ifdeci>0)
                {
                    if (ifdeci==1) current_num *= 10;
                    current_num += ifdeci * (c - '0');
                    if (ifdeci<1) ifdeci *= 0.1;
                }
                break;
        }
    }
}

void node::leaveto_root()
{
    
    if (numchildren)
    {
        float total_dist = (1-(children[0]->dist+children[1]->dist));
        
        coef += total_dist * (children[0]->coef + children[1]->coef);
        
        //total_coef +=  children[0]->total_coef + children[1]->total_coef;
        
        //total_round += children[0]->total_round + children[1]->total_round;
    }
    
    round = int(coef + 0.5);
    
    //total_round += round;
    
    coef -= round;
}

void node::rootto_leave()
{
    
    if (round && numchildren)
    {
        node *prefer_leave;
        if (round > 0)
        {
            prefer_leave = (children[1]->coef - children[1]->round > children[0]->coef - children[0]->round) ?  children[1]: children[0];
        }
        else
        {
            prefer_leave = (children[1]->coef - children[1]->round < children[0]->coef - children[0]->round) ?  children[1]: children[0];
        }
        
        prefer_leave->round += round;
        //prefer_leave->total_round += round;
    }
    
    if (numchildren)
    {
        children[0]->rootto_leave();
        children[1]->rootto_leave();
    }
}

void TreeRound::Run(node** nodes, double* unround_coefs, size_t size, int * round_coefs)
{
    
    
    int leaveindex =0;
    for (int i=0; i< size; ++i)
    {
        node * node = nodes[i];
        if (!node->numchildren) node->coef = unround_coefs[leaveindex++];
    }
    
    node * root = nodes[0];
    
    root->leaveto_root();
    
    root->rootto_leave();
    
    leaveindex = 0;
    for (int i=0; i< size; ++i)
    {
        node * node = nodes[i];
        if (!node->numchildren) round_coefs[leaveindex++] = node->coef;
    }
    
}
