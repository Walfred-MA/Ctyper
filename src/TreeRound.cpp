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


inline pair<const node*, float> reproject(pair<const node*, float>& left, pair<const node*, float>& right, vector<FLOAT_T>& rounding_residuels, FLOAT_T* reminders, int * rounds, const FLOAT_T *residuels, const FLOAT_T *weightnorm, int size)
{
    FLOAT_T &reminder_left = reminders[left.first->index];
    FLOAT_T &reminder_right = reminders[right.first->index];
    
    int index_left = left.first->index;
    int index_right = right.first->index;
    
    FLOAT_T norm_value_left = 0.0, norm_value_right = 0.0, norm_dot_residuel_left = 0.0, norm_dot_residuel_right = 0.0, residuel = 0.0 ;
    
    for (int i =0; i < size; ++i)
    {
        residuel = residuels[i] + rounding_residuels[i] -  reminder_left * weightnorm[size*i+index_left] - reminder_right * weightnorm[size*i+index_right] ;
        
        norm_dot_residuel_left += residuel * weightnorm[size*i+index_left];
        
        norm_dot_residuel_right += residuel * weightnorm[size*i+index_right];
        
        norm_value_left += weightnorm[size*i+index_left] * weightnorm[size*i+index_left];
        
        norm_value_right += weightnorm[size*i+index_right] * weightnorm[size*i+index_right];
    }
    
    FLOAT_T new_reminder_left = - norm_dot_residuel_left / norm_value_left;
    
    FLOAT_T new_reminder_right = - norm_dot_residuel_right / norm_value_right;
    
    norm_dot_residuel_left = 0.0; norm_dot_residuel_right = 0.0; residuel = 0.0 ;
    
    for (int i =0; i < size; ++i)
    {
        residuel = residuels[i] + rounding_residuels[i] +  (new_reminder_right - reminder_right) * weightnorm[size*i+index_right] ;
        
        norm_dot_residuel_left += residuel * weightnorm[size*i+index_left];
        
        
        residuel = residuels[i] + rounding_residuels[i] +  (new_reminder_left - reminder_left) * weightnorm[size*i+index_left] ;
        
        norm_dot_residuel_right += residuel * weightnorm[size*i+index_right];
    }
    
    FLOAT_T legacy_reminder_left = - norm_dot_residuel_left / norm_value_left;
    
    FLOAT_T legacy_reminder_right = - norm_dot_residuel_right / norm_value_right;
    
    
    FLOAT_T residuel_left = 0.0, residuel_right = 0.0;
    for (int i =0; i < size; ++i)
    {
        residuel = residuels[i] + rounding_residuels[i] + (new_reminder_left - reminder_left) * weightnorm[size*i+index_left] - reminder_right * weightnorm[size*i+index_right] ;
        
        residuel_left += residuel * residuel;
        
        residuel = residuels[i] + rounding_residuels[i] + (new_reminder_right - reminder_right) * weightnorm[size*i+index_right] - reminder_left * weightnorm[size*i+index_left] ;
        
        residuel_right += residuel * residuel;
    }
    
    
    if (residuel_right < residuel_left)
    {
        /*
        cout<< "distance: " << left.second + right.second << endl;
        cout<< "source: "<< reminder_left << endl;
        cout<< "old: "<<reminder_right << endl;
        cout<< "new: "<<new_reminder_right << endl;
        cout<< "reser: "<<index_left<< ","<<reminder_left + legacy_reminder_left << endl;
        */
        
        for (int i =0; i < size; ++i)
        {
            rounding_residuels[i] += (new_reminder_right - reminder_right) * weightnorm[size*i+index_right] - reminder_left * weightnorm[size*i+index_left]  ;
        }
        
        reminders[left.first->index] += legacy_reminder_left;
        
        float round = ((new_reminder_right) > 0) ? int(new_reminder_right + 0.5) : int(new_reminder_right - 0.5);
        
        reminders[index_right] = new_reminder_right - round;
        rounds[index_right] += round;
        
        return right;
    }
    else
    {
        /*
        cout<< "distance: " << left.second + right.second << endl;
        cout<< "source: "<< reminder_right << endl;
        cout<< "old: "<<reminder_left << endl;
        cout<< "new: "<<new_reminder_left << endl;
        cout<< "reser: "<<index_right<< ","<<reminder_right + legacy_reminder_right << endl;
    	*/    
        for (int i =0; i < size; ++i)
        {
            rounding_residuels[i] += (new_reminder_left - reminder_left) * weightnorm[size*i+index_left] - reminder_right * weightnorm[size*i+index_right] ;
        }
        
        
        reminders[right.first->index] += legacy_reminder_right;
        
        float round = ((new_reminder_left) > 0) ? int(new_reminder_left + 0.5) : int(new_reminder_left - 0.5);
        
        reminders[index_left] = new_reminder_left - round;
        rounds[index_left] += round;
        
        return left;
    }
    
    
}

/*
inline void whichiscloser(pair<const node*, float>& left, pair<const node*, float>&right, pair<const node*, float>& closer, pair<const node*, float>& further,  FLOAT_T* reminders, const FLOAT_T *residuels, const FLOAT_T *weightnorm, int size)
{
    int left_index = left.first->index;
    int right_index = right.first->index;
    
    float left_reminder = left.first->size * reminders[left_index];
    float right_reminder = right.first->size * reminders[right_index];
    
    float total_reminder = left_reminder + right_reminder ;
    
    FLOAT_T residuel_change = 0.0;
    if ( abs(reminders[left_index] + reminders[right_index]) > 0.01 or left_reminder/right_reminder > 0.3 or right_reminder/left_reminder > 0.3 )
    {
        FLOAT_T left_roundresiduel = 0.0, right_roundresiduel = 0.0;
        
        FLOAT_T left_roundint = (left_reminder > 0 ) ? 1: -1;
        FLOAT_T right_roundint = (right_reminder > 0 ) ? 1: -1;
        
        
        for (int i =0; i < size; ++i)
        {
            residuel_change = ( residuels[i] + (left_roundint - left_reminder) * weightnorm[size*i+left_index] )  ;
            
            left_roundresiduel += residuel_change * residuel_change;
            
            residuel_change = ( residuels[i] + (right_roundint - right_reminder) * weightnorm[size*i+right_index] )  ;
            
            right_roundresiduel += residuel_change * residuel_change;
            
        }
        
        if ( right_roundresiduel < left_roundresiduel )
        {
            closer =  right ;
            
            further = left;
        }
        
        else
        {
            closer =  left ;
            
            further = right ;
        }

        
    }
    
    else if ( ( right_reminder > left_reminder ) == ( total_reminder >= 0 ) )  //when total reminder is positive, pick the bigger postive reminder
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
*/

pair<const node*, float> node::leaveto_root(FLOAT_T* reminders ,int * rounds, vector<FLOAT_T> &rounding_residuels, const FLOAT_T *residuels, const FLOAT_T *weightnorm, const int size) const
{
    pair<const node*, float> closer, further;
    
    if (numchildren)
    {
        pair<const node*, float> left, right;
        
        if (children[0]->numchildren > 0 && children[1]->numchildren > 0)
        {
            vector<FLOAT_T> rounding_residuels_second = rounding_residuels;
            
            left = children[0]->leaveto_root(reminders,  rounds, rounding_residuels, residuels, weightnorm, size);
            
            right = children[1]->leaveto_root(reminders, rounds, rounding_residuels_second, residuels, weightnorm, size);
            
            for (size_t i = 0; i < size; ++i)
            {
                 rounding_residuels[i] += rounding_residuels_second[i];
            }
        }
        else
        {
            left = children[0]->leaveto_root(reminders,  rounds, rounding_residuels, residuels, weightnorm, size);
            
            right = children[1]->leaveto_root(reminders, rounds, rounding_residuels, residuels, weightnorm, size);
        }
        
        
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
            closer = reproject(left , right, rounding_residuels, reminders, rounds, residuels, weightnorm, size);
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
            float round = ((reminder) >= 0) ? int(reminder + 0.5) : int(reminder - 0.5);
            
            reminder -= round;
            
            rounds[this->index] += round;
            
            closer = make_pair(this, this->dist);
        }
    }
    return closer;
}


void TreeRound::Run(const node* nodes, FLOAT_T* unround_coefs, size_t size, int * round_coefs, FLOAT_T * reminder_coefs, const FLOAT_T *residuels, const FLOAT_T *weightnorm)
{
    memcpy(reminder_coefs, unround_coefs, sizeof (FLOAT_T) *  size+10);
    
    const node& root = nodes[0];
    
    vector<FLOAT_T> rounding_residuels (size, 0.0);
    
    root.leaveto_root(reminder_coefs, round_coefs, rounding_residuels, residuels, weightnorm, size);
}
