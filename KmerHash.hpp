//
//  KmerHash.hpp
//  CTyper
//
//  Created by Wangfei MA on 3/29/23.
//

#ifndef KmerHash_hpp
#define KmerHash_hpp

#include <stdio.h>
#include <vector>
#include <unordered_map>
#include <tuple>
#include "config.hpp"

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
    
    auto operator>(const uint40& other) const
    {
        return (this->first > other.first || ( this->first == other.first && this->second > other.second) );
    }

    uint first;
    uint8 second;
};

class Kmer32_hash
{
public:
    Kmer32_hash(int size): modsize(MAX(4194304,size))
    {};
    
    ~Kmer32_hash()
    {
        for (int i = 0 ; i < modsize; ++i)
        {
            if (key_sizes[i])
            {
                free(keys[i]);
                free(values[i]);
            }
        }
        free(keys);
        free(values);
    }
    
    uint* add(const ull kmer_int, const uint index);
    uint* find(const ull kmer_int);

private:
    const size_t modsize;
    uint40** keys = new uint40*[modsize];
    uint** values = new uint*[modsize];
    int* key_sizes = new int[modsize] ();
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


#endif /* KmerHash_hpp */
