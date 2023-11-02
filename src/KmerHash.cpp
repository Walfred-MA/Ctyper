//
//  KmerHash.cpp
//  CTyper
//
//  Created by Wangfei MA on 3/29/23.
//

#include "KmerHash.hpp"
#define uint unsigned int


static int CompareItem40(const void *a, const void *b) {
    const item40_t *itemA = (const item40_t *)a;
    const item40_t *itemB = (const item40_t *)b;
    if (itemA->first.first != itemB->first.first)
        return (itemA->first.first > itemB->first.first) ? 1 : -1;
    if (itemA->first.second != itemB->first.second)
        return (itemA->first.second > itemB->first.second) ? 1 : -1;
    return 0;
}

static int Search(const item40_t *arr, const uint size, const uint40 x) {
    item40_t key = { .first = { .first = x.first, .second = x.second } };
    item40_t *item = bsearch(&key, arr, size, sizeof(item40_t), CompareItem40);
    if (item != NULL) {
        return (int)(item - arr); // The index is the difference between pointers
    }
    return -1;
}

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
