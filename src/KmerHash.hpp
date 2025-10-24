//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License. 
//  If you use this code, please cite our work.
//   

#ifndef KmerHash_hpp
#define KmerHash_hpp

#include <stdio.h>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <functional>
#include <cstdint>
#include "config.hpp"

#define MAX_UINT24 16777215
#define item40_t std::pair<uint40,uint>
#define uint unsigned int

using namespace std;


struct hash_128
{
    size_t operator()(const u128& num128) const
    {
        auto hash1 = hash<ull>{}((ull)(num128%large_prime));
        return hash1;
    }
};

using kmer64_dict = std::unordered_map<u128, uint, hash_128> ;
using kmer32_dict = std::unordered_map<ull, uint > ;

using kmer64_dict_nt = std::unordered_map<u128, uint8, hash_128> ;
using kmer32_dict_nt = std::unordered_map<ull, uint8>;

using kmer64_dict_mul = std::unordered_map<u128, uint*, hash_128> ;
using kmer32_dict_mul = std::unordered_map<ull, uint* > ;

struct uint40
{
    uint40 (const ull x): first( x % UINT_MAX ) , second( x / UINT_MAX )
    {};
    
    auto operator==(const uint40& other) const
    {
        return (this->first == other.first && this->second == other.second);
    }
 
    uint first;
    uint8 second;
};


class Kmer32_hash
{
public:
    Kmer32_hash(int size): modsize(MAX(MAX_UINT24,size))
    {
        key_sizes = new uint[modsize]();              // Initialize all to 0
        items = new item40_t*[modsize]();             // Initialize all to nullptr
    };
    
    ~Kmer32_hash()
    {
        for (int i = 0 ; i < modsize; ++i)
        {
            if (key_sizes[i])
            {
                free(items[i]);
            }
        }
        free(items);
    }
    
    uint* add(const ull kmer_int, const uint index);
    uint* find(const ull kmer_int);
    uint* find(const ull kmer_int, uint &hash_);
    const size_t modsize;

private:
    item40_t** items = new item40_t*[modsize];
    uint* key_sizes = new uint[modsize];
};


class Kmer64_hash
{
public:
    Kmer64_hash(int size): modsize(size)
    {};
    
    uint* add(const u128 kmer_int, const uint index);
    uint* find(const u128 kmer_int);
    
private:
    const size_t modsize;
    kmer64_dict data;
    
};

/*
class Group_hash
{
public:
    Group_hash(size_t t):size(t)
    {
        groupindex.resize(size,0);
        redirect.resize(size,0);
        multiindex.resize(100000000);
        
    };
    
    ~Group_hash()
    {}
    
    uint* add(const uint index, const uint groupindex);
    uint find(const ull kmer_int);
    
    vector<uint16> groupindex;
    vector<bool> redirect;
    vector<vector<uint>> multiindex;
    size_t size;
    size_t multisize = 100000000;
};
 */
#endif /* KmerHash_hpp */
