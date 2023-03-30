//
//  KmerHash.cpp
//  CTyper
//
//  Created by Wangfei MA on 3/29/23.
//

#include "KmerHash.hpp"


static int Search(const uint40 *arr, const uint size, const uint40 x)
{
    for (size_t i = 0 ; i < size; ++i)
    {
        if (arr[i] == x) return i;
    }
    return -1;
}


uint* Kmer32_hash::add(const ull key, const uint val)
{
    uint hash_ = key % modsize;
    uint40 reminder_ = uint40(key / modsize);
    
    uint &size = key_sizes[hash_];
    uint40 *&bucket_key = keys[hash_];
    uint *&bucket_val = values[hash_];
    
    int loc = Search(bucket_key, size, reminder_);
    
    if ( loc == -1 )
    {
        size ++;
        
        bucket_key = (uint40*) realloc(bucket_key, sizeof(uint40)*(size));
        bucket_val = (uint*) realloc(bucket_val, sizeof(uint)*(size));
        
        bucket_key[size - 1] = reminder_;
        bucket_val[size - 1] = val;
                    
        return &bucket_val[size - 1];
    }
    
    return &bucket_val[loc];

}

uint* Kmer32_hash::find(const ull kmer_int)
{
    uint hash_ = kmer_int % modsize;
    uint40 reminder_ = uint40(kmer_int / modsize);
    
    uint size = key_sizes[hash_];
    uint40 *bucket_key = keys[hash_];
    uint *bucket_val = values[hash_];
    
    int loc = Search(bucket_key, size, reminder_);
    
    if (loc < 0)
    {
        return NULL;
    }
    else
    {
        return &bucket_val[loc];
    }
}



uint *Kmer64_hash::add(const u128 kmer_int, const uint index)
{
    auto map_find = data.find(kmer_int);
    
    if (map_find == data.end())
    {
        data[kmer_int] = kmer_int;
        
        return &data[kmer_int];
    }
    
    return &map_find->second;
}
uint *Kmer64_hash::find(const u128 kmer_int)
{
    auto map_find = data.find(kmer_int);
    
    if (map_find == data.end())
    {
        return NULL;
    }
    
    return &map_find->second;
}
