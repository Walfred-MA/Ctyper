//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License. 
//  If you use this code, please cite our work.
//   

#include "KmerHash.hpp"
#define uint unsigned int

/*
static int CompareItem40(const void *a, const void *b) 
{
    const uint40 *key = (const uint40 *)a;
    const item40_t *item = (const item40_t *)b;
    if (key->first != item->first.first)
        return (key->first > item->first.first) ? 1 : -1;
    if (key->second != item->first.second)
        return (key->second > item->first.second) ? 1 : -1;
    return 0;
}

static int BSearch(const item40_t *arr, const uint size, const uint40 x) 
{
    item40_t *item = bsearch(&x, arr, size, sizeof(item40_t), CompareItem40);
    if (item != NULL) {
        return (int)(item - arr); // The index is the difference between pointers
    }
    return -1;
}
*/
static int Search(const item40_t *arr, const uint size, const uint40 x)
{
    for (uint i = 0 ; i < size; ++i)
    {
        if (arr[i].first.first == x.first && arr[i].first.second == x.second) return i;
    }
    return -1;
}



uint* Kmer32_hash::add(const ull kmer_int, const uint val)
{
    std::hash<std::uint64_t> hasher;
    ull key = hasher(kmer_int);
    uint hash_ = key % modsize;
    uint40 reminder_ = uint40(key / modsize);

    uint &size = key_sizes[hash_];
    item40_t *&bucket = items[hash_];

    int loc = Search(bucket, size, reminder_);

    if ( loc == -1 )
    {
        size ++;
        
        bucket = (item40_t *) realloc(bucket, sizeof(item40_t)*(size));

        bucket[size - 1] = make_pair(reminder_, val);
                    
        return &bucket[size - 1].second;
    }
    
    return &bucket[loc].second;

}


uint* Kmer32_hash::find(const ull kmer_int)
{
    std::hash<std::uint64_t> hasher;
    ull key = hasher(kmer_int);
    uint hash_ = key % modsize;
    uint40 reminder_ = uint40(key / modsize);
    
    uint size = key_sizes[hash_];
    item40_t *bucket = items[hash_];
        
    int loc = Search(bucket, size, reminder_);
    
    if (loc < 0)
    {
        return NULL;
    }
    else
    {
        return &bucket[loc].second;
    }
}


uint* Kmer32_hash::find(const ull kmer_int, uint &hash_)
{
    std::hash<std::uint64_t> hasher;
    ull key = hasher(kmer_int);
    hash_ = key % modsize;
    uint40 reminder_ = uint40(key / modsize);
    
    uint size = key_sizes[hash_];
    item40_t *bucket = items[hash_];
        
    int loc = Search(bucket, size, reminder_);
    
    if (loc < 0)
    {
        return NULL;
    }
    else
    {
        return &bucket[loc].second;
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
