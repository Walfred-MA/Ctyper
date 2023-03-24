//
//  KmerCounter.hpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#ifndef KmerCounter_hpp
#define KmerCounter_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstring>
#include <unordered_set>
#include <algorithm>
#include <atomic>
#include <mutex>
#include <thread>
#include <atomic>
#include <unordered_map>
#include <tuple>

#include "FastaReader.hpp"
#include "FastqReader.hpp"
//#include "CramReader.hpp"
#include "KtableReader.hpp"

using namespace std;

struct hash_128
{
    size_t operator()(const u128& num128) const
    {
        auto hash1 = hash<ull>{}((ull)(num128%large_prime));
        return hash1;
    }
};

template <typename T1>
string int_to_kmer(T1 value, int size = 31)
{
    string kmer (31,' ');
    for (int i = 0; i<31 ; ++i)
    {
        kmer[30-i] = "ACGT"[value % 4];
        value /= 4;
    }
    return kmer;
}

static inline bool base_to_int(const char base, int &converted)
{

    switch (base)
    {
        case 'A' : case 'a':
            converted=0b00;
            break;
        case 'C' : case 'c':
            converted=0b01;
            break;
        case 'G' : case 'g':
            converted=0b10;
            break;
        case 'T' : case 't':
            converted=0b11;
            break;
        default:
            return 0;
    }
    
    return 1;
}



using kmer64_dict = std::unordered_map<u128, uint, hash_128> ;
using kmer32_dict = std::unordered_map<ull, uint > ;

using kmer64_dict_nt = std::unordered_map<u128, uint8, hash_128> ;
using kmer32_dict_nt = std::unordered_map<ull, uint8>;

using kmer64_dict_mul = std::unordered_map<u128, uint*, hash_128> ;
using kmer32_dict_mul = std::unordered_map<ull, uint* > ;

template <int dictsize>
class KmerCounter
{
    using kmer_int = typename std::conditional<(dictsize>32), u128, ull>::type;
    using kmer_hash_type = typename std::conditional<(dictsize>32), kmer64_dict, kmer32_dict>::type;
    using kmer_hash_type_mul = typename std::conditional<(dictsize>32), kmer64_dict_mul, kmer32_dict_mul>::type;
 
public:
    
    KmerCounter ()
    {};
    ~KmerCounter()
    {};
    
    ull read_target(KtableReader &fastafile);
    
    ull read_target(FastaReader &fastafile);
    
    ull read_target(const char* infile);
        
    template <class typefile>
    void count_kmer(typefile &file, uint16* samplevecs);
    
    void count_kmer_(char* infile, uint16* samplevecs);
    
    void read_files(std::vector<std::string>& inputfiles, std::vector<std::string>& outputfiles, std::vector<std::string>& prefixes,std::vector<float>& deps,int numthread);
    
    void Call(const char* infile, uint16* samplevecs, const int nthreads);
    
    kmer_hash_type kmer_hash;
    kmer_hash_type_mul kmer_multi_hash;
    
    uint totalkmers = 0;
    
    const int klen = 31;

};

template <typename T1, typename T2>
static bool initiate_counter(T2 &kmer_hash, T1 &larger_kmer, uint &kindex)
{
    
    typename T2::iterator map_find = kmer_hash.find(larger_kmer);

    if (map_find == kmer_hash.end())
    {
        kmer_hash[larger_kmer] = (int) kindex++;
        
        return false;
    }
    
    return true;
}

template <typename T1, typename T2, typename T3>
static bool initiate_counter_mul(T2 &kmer_hash, T3 &kmer_multi_hash, T1 &larger_kmer, uint &kindex)
{
    typename T2::iterator map_find = kmer_hash.find(larger_kmer);

    if (map_find == kmer_hash.end())
    {
        kmer_hash[larger_kmer] = (int) kindex++;
        
        return false;
    }
    
    else if (map_find->second == MAX_UINT16 )
    {
        uint* &data = kmer_multi_hash.find(larger_kmer)->second;
        if (data[0]%5 == 2)
        {
            data = (uint*) realloc(data, sizeof(uint)*(data[0]+1 + 5));
        }
            
        data[0]++;
        data[data[0]] = kindex++;
        
        return false;
    }
    
    else
    {
        kmer_hash[larger_kmer] = MAX_UINT16;
        
        uint* newarray = (uint*) malloc(sizeof(uint)*(3));
        newarray[0] = 2;
        newarray[1] = map_find->second;
        newarray[2] = kindex++;
        
        kmer_multi_hash[larger_kmer] = newarray;
        
        return true;
    }
    
    return true;
}


template <typename T1, typename T2, typename T3>
static void update_counter(T2 &kmer_hash, T3 &kmer_multi_hash, T1 &larger_kmer, uint16* vec)
{
    typename T2::iterator map_find = kmer_hash.find(larger_kmer);

    if (map_find != kmer_hash.end())
    {
        
        uint index = map_find->second;
        
        if (index < MAX_UINT16 )
        {
            if ( vec[index] < MAX_UINT16 - 1) vec[index] ++;
        }
        else
        {
            uint* &data = kmer_multi_hash.find(larger_kmer)->second;
            
            uint num_num = data[0];
            
            for (uint i = 1 ; i < num_num + 1; ++i)
            {
                if ( vec[data[i]] < MAX_UINT16 - 1 )  vec[data[i] ] ++;
            }
        }
    }
}





template <typename T>
static void kmer_deconpress(const char * kmer_compress, T &larger_kmer)
{
    larger_kmer = 0;
    for (int pos = 0; pos < 30; ++pos)
    {
        if (kmer_compress[pos] == '\t') break;
        larger_kmer <<= 6;
        larger_kmer += kmer_compress[pos] - '0';
    }
}


template <typename T>
static void kmer_read_c(char base, const int klen, std::size_t &current_size, T &current_kmer, T &reverse_kmer)
{
    int converted = 0;
    T reverse_converted;
    
    if (base == '\n' || base == ' ') return;
    
    if (base_to_int(base, converted))
    {
        
        current_kmer <<= ( 8*sizeof(current_kmer) - 2*klen + 2 );
        current_kmer >>=  ( 8*sizeof(current_kmer) - 2*klen );
        current_kmer += converted;
        
        reverse_kmer >>= 2;
        reverse_converted = 0b11-converted;
        reverse_converted <<= (2*klen-2);
        reverse_kmer += reverse_converted;
        
    }
    
    else
    {
        current_size = -1;
        current_kmer = 0;
        reverse_kmer = 0;
    }
}

template <int dictsize>
ull KmerCounter<dictsize>::read_target(FastaReader &fastafile)
{
    
    std::size_t current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
        
    std::string StrLine;
    StrLine.resize(MAX_LINE);
    
    while (fastafile.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '@':  case '+': case '>':
                current_size = 0;
                continue;
            case ' ': case '\n': case '\t':
                continue;
            default:
                break;
        }
        
        for (auto base: StrLine)
        {
            if (base == '\0') break;
            
            if (base == '\n' || base == ' ') continue;

            kmer_read_c(base, klen, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen || base >= 'a') continue;
                                
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
            if (initiate_counter(kmer_hash, larger_kmer, totalkmers)) continue;
                
        }
    }
    
    
    fastafile.Close();
    
    return totalkmers;

};

template <int dictsize>
ull KmerCounter<dictsize>::read_target(KtableReader &ktablefile)
{
    
    kmer_int current_kmer = 0;
        
    std::string StrLine;
    StrLine.resize(MAX_LINE);
    
    char base;
    int pos;
    while (ktablefile.nextLine_kmer(StrLine))
    {
        for (pos = 0; pos < klen; ++pos)
        {
            base = StrLine[pos];
            if (base == '\t') break;
        }
        ++pos;
        
        kmer_deconpress(&StrLine.c_str()[pos], current_kmer);
        
        initiate_counter_mul(kmer_hash, kmer_multi_hash, current_kmer, totalkmers);
        
    }
    
    ktablefile.Close();
    
    return totalkmers;
        
};
 
template <int dictsize>
ull KmerCounter<dictsize>::read_target(const char* inputfile)
{
    int pathlen = (int)strlen(inputfile);
    
    if ( (pathlen > 2 && strcmp(inputfile+(pathlen-3),".fa")==0) || (pathlen > 6 && strcmp(inputfile+(pathlen-6),".fasta") == 0 ))
    {
        FastaReader readsfile(inputfile);
        read_target(readsfile);
    }
    else
    {
        
        KtableReader readsfile(inputfile);
        read_target(readsfile);
    }
    
    return totalkmers;
}



template <int dictsize>
template <class typefile>
void KmerCounter<dictsize>::count_kmer(typefile &file, uint16* samplevecs)
{
    ull totalsamplekmers = 0;
    std::size_t current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    //uint64_t ifmasked = 0;
    //int num_masked = 0 ;
    std::string StrLine;
    std::string Header;
    
    while (file.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '@':  case '+': case '>':
                current_size = 0;
                current_kmer = 0;
                reverse_kmer = 0;
                continue;
            case ' ': case '\n': case '\t':
                continue;
            default:
                break;
        }
        for (auto base: StrLine)
        {
            
            if (base == '\0') break;
                        
            if (base == '\n' || base == ' ') continue;
 
            kmer_read_c(base, klen, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen) continue;
            
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
            
            update_counter(kmer_hash, kmer_multi_hash, larger_kmer, samplevecs);
            
        }
    }
        
    file.Close();
    return ;
};


template <int dictsize>
void KmerCounter<dictsize>::count_kmer_(char* inputfile, uint16* samplevecs)
{
    
    int pathlen = (int)strlen(inputfile);
    
    if ( pathlen > 2 && strcmp(inputfile+(pathlen-3),".gz") == 0 )
    {
        FastqReader readsfile(inputfile);
        count_kmer(readsfile, samplevecs);
    }
        
    else if ( (pathlen > 2 && strcmp(inputfile+(pathlen-3),".fa")==0) || (pathlen > 6 && strcmp(inputfile+(pathlen-6),".fasta") == 0 ))
    {
        FastaReader readsfile(inputfile);
        
        count_kmer(readsfile, samplevecs);
    }
    else if (pathlen > 5 && ( strcmp(inputfile+(pathlen-3),".cram") == 0 ))
    {
        //CramReader readsfile(inputfile);
        //totalsamplekmers = count_kmer(readsfile, samplevecs);
    }
    
};

template <int dictsize>
void KmerCounter<dictsize>::Call(const char* infile, uint16* samplevecs, const int nthreads)
{
    std::vector<std::thread> threads;
        
    for(int i=0; i< nthreads; ++i)
    {
        threads.push_back(std::thread( std::bind( &KmerCounter<dictsize>::count_kmer_, this, (char*)infile, samplevecs ) ));
    }
    
    
    for(int i=0; i< nthreads; ++i)
    {
        threads[i].join();
    }
    
    
}


#endif /* KmerCounter_hpp */
