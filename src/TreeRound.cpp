//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License. 
//  If you use this code, please cite our work.
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
    
    FLOAT_T norm_value_left = 0.0, norm_value_right = 0.0, norm_dot_residuel_left = 0.0, norm_dot_residuel_right = 0.0, residuel = 0.0;
    
    
    //here we calculate the coeficients of left and right side as distinations after reprojection
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

    
    
    //calculate residuel after reprojection
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
        
        norm_dot_residuel_left = 0.0; norm_dot_residuel_right = 0.0; residuel = 0.0 ;
        for (int i =0; i < size; ++i)
        {
            residuel =  reminder_left * weightnorm[size*i+index_left] ;
                        
            norm_dot_residuel_right += residuel * weightnorm[size*i+index_right];
        }
        
        new_reminder_right = norm_dot_residuel_right / norm_value_right + reminder_right;
        
        
        //calculate the left coeficients as sources after reprojection
        norm_dot_residuel_left = 0.0; norm_dot_residuel_right = 0.0; residuel = 0.0 ;
        
        for (int i =0; i < size; ++i)
        {
            residuel = residuels[i] + rounding_residuels[i] +  (new_reminder_right - reminder_right) * weightnorm[size*i+index_right] ;
            
            norm_dot_residuel_left += residuel * weightnorm[size*i+index_left];
            
            residuel = residuels[i] + rounding_residuels[i] +  (new_reminder_left - reminder_left) * weightnorm[size*i+index_left] ;
            
            norm_dot_residuel_right += residuel * weightnorm[size*i+index_right];
        }
        
        FLOAT_T legacy_reminder_left = - norm_dot_residuel_left / norm_value_left;
        
        
        
        reminders[left.first->index] = legacy_reminder_left;
        
        float round = ((new_reminder_right) > 0) ? int(new_reminder_right + 0.5) : int(new_reminder_right - 0.5);
        
        reminders[index_right] = new_reminder_right - round;
        rounds[index_right] += round;
        
        return right;
    }
    else
    {
        norm_dot_residuel_left = 0.0; norm_dot_residuel_right = 0.0; residuel = 0.0 ;
        for (int i =0; i < size; ++i)
        {
            residuel = reminder_right * weightnorm[size*i+index_right] ;
            
            norm_dot_residuel_left += residuel * weightnorm[size*i+index_left];
        }
        
        
        new_reminder_left = norm_dot_residuel_left / norm_value_left + reminder_left;
        
        
        //calculate the left coeficients as sources after reprojection
        norm_dot_residuel_left = 0.0; norm_dot_residuel_right = 0.0; residuel = 0.0 ;
        
        for (int i =0; i < size; ++i)
        {
            residuel = residuels[i] + rounding_residuels[i] +  (new_reminder_right - reminder_right) * weightnorm[size*i+index_right] ;
            
            norm_dot_residuel_left += residuel * weightnorm[size*i+index_left];
            
            residuel = residuels[i] + rounding_residuels[i] +  (new_reminder_left - reminder_left) * weightnorm[size*i+index_left] ;
            
            norm_dot_residuel_right += residuel * weightnorm[size*i+index_right];
        }
        
        //FLOAT_T legacy_reminder_left = - norm_dot_residuel_left / norm_value_left;
        
        FLOAT_T legacy_reminder_right = - norm_dot_residuel_right / norm_value_right;
        
        reminders[right.first->index] = legacy_reminder_right;
        
        float round = ((new_reminder_left) > 0) ? int(new_reminder_left + 0.5) : int(new_reminder_left - 0.5);
        
        reminders[index_left] = new_reminder_left - round;
        rounds[index_left] += round;
        
        
        return left;
    }
    
    
}


pair<const node*, float> node::leaveto_root(FLOAT_T* reminders ,int * rounds, vector<FLOAT_T> &rounding_residuels, const FLOAT_T *residuels, const FLOAT_T *weightnorm, const int size) const
{
    pair<const node*, float> closer, further;
    
    if (numchildren)
    {
        pair<const node*, float> left, right;
        
        if (children[0]->numchildren > 0 && children[1]->numchildren > 0)
        {
            
            //left = children[0]->leaveto_root(reminders,  rounds, rounding_residuels, residuels, weightnorm, size);
            //right = children[1]->leaveto_root(reminders, rounds, rounding_residuels, residuels, weightnorm, size);
            
            vector<FLOAT_T> &rounding_residuels_second = rounding_residuels;
            left = children[0]->leaveto_root(reminders,  rounds, rounding_residuels, residuels, weightnorm, size);
            right = children[1]->leaveto_root(reminders, rounds, rounding_residuels_second, residuels, weightnorm, size);
            for (size_t i = 0; i < size; ++i)
            {
                 rounding_residuels[i] += rounding_residuels_second[i] - rounding_residuels[i];
            }
            
        }
        else
        {
            //left = children[0]->leaveto_root(reminders,  rounds, rounding_residuels, residuels, weightnorm, size);
            //right = children[1]->leaveto_root(reminders, rounds, rounding_residuels, residuels, weightnorm, size);
            
            vector<FLOAT_T> &rounding_residuels_second = rounding_residuels;
            right = children[1]->leaveto_root(reminders, rounds, rounding_residuels, residuels, weightnorm, size);
            left = children[0]->leaveto_root(reminders,  rounds, rounding_residuels_second, residuels, weightnorm, size);
            
            for (size_t i = 0; i < size; ++i)
            {
                 rounding_residuels[i] += rounding_residuels_second[i] - rounding_residuels[i];
            }
            
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
    memcpy(reminder_coefs, unround_coefs, sizeof (FLOAT_T) *  size);
    
    const node& root = nodes[0];
    
    vector<FLOAT_T> rounding_residuels (size, 0.0);
    
    root.leaveto_root(reminder_coefs, round_coefs, rounding_residuels, residuels, weightnorm, size);
}
