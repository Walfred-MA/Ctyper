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


node* node::add(node * child)
{
    child->clear();
    children[numchildren ++] = child;
    child->parent = this;
    
    return child;
}


void reproject(pair<const node*, float>& target, pair<const node*, float>& source, FLOAT_T* reminders, int * rounds)
{
    FLOAT_T &reminder_source = reminders[source.first->index];
    FLOAT_T &reminder_target = reminders[target.first->index];
        
    float total_dist = target.second + source.second;
    
    float proj_val = reminder_source * (1-total_dist) * source.first->size / target.first->size;
    
    reminder_target += proj_val;
    
    float round = ((reminder_target) > 0) ? int(reminder_target + 0.5) : int(reminder_target - 0.5);
    
    reminder_target -= round;
    rounds[target.first->index] += round;
    
}


void whichiscloser(pair<const node*, float>& left, pair<const node*, float>&right, pair<const node*, float>& closer, pair<const node*, float>& further,  FLOAT_T* reminders)
{
    
    float left_reminder = left.first->size * reminders[left.first->index];
    float right_reminder = right.first->size * reminders[right.first->index];
    
    float total_reminder = left_reminder + right_reminder ;
    
    if ( ( right_reminder > left_reminder ) == ( total_reminder >= 0 ) )  //when total reminder is positive, pick the bigger postive reminder
    {
        closer =  right ;
        
        further = left;
    }
    
    else  //when total reminder is negative, pick the bigger negative reminder
    {
        closer =  left ;
        
        further = right ;
    }
    
    
}


pair<const node*, float> node::leaveto_root(FLOAT_T* reminders, int * rounds) const
{
    pair<const node*, float> closer, further;
    
    if (numchildren)
    {
        
        pair<const node*, float> left = children[0]->leaveto_root(reminders, rounds);
        pair<const node*, float> right = children[1]->leaveto_root(reminders, rounds);
        
        if (reminders[right.first->index] == 0)
        {
            closer = left;
        }
        
        else if (reminders[left.first->index] == 0)
        {
            closer = right;
        }
        
        else
        {
            whichiscloser(left, right, closer, further, reminders);
            
            reproject(closer, further, reminders, rounds);
        }
        
        closer.second += this->dist;
        
    }
    else
    {
        FLOAT_T &reminder = reminders[this->index];
        
        if (reminder == 0)
        {
            closer = make_pair(this, this->dist);
        }
        
        else
        {
            float round = ((reminder) > 0) ? int(reminder + 0.5) : int(reminder - 0.5);
            
            reminder -= round;
            rounds[this->index] += round;
            
            closer = make_pair(this, this->dist);
        }
    }
    return closer;
}


void TreeRound::Run(const node* nodes, FLOAT_T* unround_coefs, size_t size, int * round_coefs)
{
    memcpy(reminder.get(), unround_coefs, sizeof (FLOAT_T) *  size+10);
    
    const node& root = nodes[0];

    root.leaveto_root(reminder.get(), round_coefs);

}
