//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License. 
//  If you use this code, please cite our work.
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
#include <filesystem>

#include "FastaReader.hpp"
#include "FastqReader.hpp"
#include "CramReader.hpp"
#include "KtableReader.hpp"
#include "KmerHash.hpp"
#include "config.hpp"

extern int ifprofile;
extern int forceunmap;
extern string refpath;
extern int inputdep;

struct Hash10M
{
    std::size_t operator()(ull key) const
    {
        return static_cast<std::size_t>(key % prime10M);
    }
};

/*
static inline ull reverse_compliment(uint64_t x)
{
    ull x0  = x;
    x = __builtin_bitreverse64(x);
    uint64_t y = ((x & 0x5555555555555555ULL) << 1) | ((x & 0xAAAAAAAAAAAAAAAAULL) >> 1);
        
    return (~y) >>2;
}
*/

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
        case 'A' :
            converted=0b00;
            break;
        case 'C' :
            converted=0b01;
            break;
        case 'G' :
            converted=0b10;
            break;
        case 'T' :
            converted=0b11;
            break;
        default:
            return 0;
    }
    
    return 1;
}


template <int dictsize>
class KmerCounter
{
public:
    
    using kmer_int = typename std::conditional<(dictsize>32), u128, ull>::type;
    using kmer_hash_type = typename std::conditional<(dictsize>32), Kmer64_hash, Kmer32_hash>::type;
    using kmer_hash_type_mul = typename std::conditional<(dictsize>32), kmer64_dict_mul, kmer32_dict_mul>::type;
    
    KmerCounter (size_t size):kmer_hash(size)
    {};
    ~KmerCounter()
    {};
        
    ull read_target(KtableReader &fastafile);
    
    ull read_target(KtableReader &fastafile, vector<pair<size_t,size_t>>& poses);
    
    ull read_target(FastaReader &fastafile);
    
    ull read_target(const char* infile);
    ull read_target(const char* infile, vector<pair<size_t,size_t>>& poses);
	
    void load_backgrounds(const char * backfile);
    
    void LoadRegion(std::vector<char *> &r, int hla, int unmap);
        
    template <class typefile>
    void count_kmer_(typefile &file, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, vector<uint16> &nBg, const int nthreads);
    
    template <bool EnableLock>
    void count_kmer(CramReader &file, uint16* samplevecs , ull_atom &nBases, ull_atom &nReads ,vector<uint16> &nBg);
    
    template <bool EnableLock>
    void count_kmer(FastaReader &file, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, vector<uint16> &nBg);
   
    template <bool EnableLock>
    void count_kmer(FastqReader &file, uint16* samplevecs , ull_atom &nBases, ull_atom &nReads ,vector<uint16> &nBg);
    
    template <bool EnableLock>
    void count_kmer(FileReader &file, uint16* samplevecs , ull_atom &nBases, ull_atom &nReads ,vector<uint16> &nBg);
 
    void QuickCoverage(const char* file, vector<uint16>& nBg, const int nthreads);
    
    void read_files(std::vector<std::string>& inputfiles, std::vector<std::string>& outputfiles, std::vector<std::string>& prefixes,std::vector<float>& deps,int numthread);
        
    void Call(const char* infile, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, vector<uint16>& nBg, const int nthreads);
    
    kmer_hash_type kmer_hash;
    kmer_hash_type_mul kmer_multi_hash;

    vector<ull> backgrounds = vector<ull>(background_prime, 0);
    
    uint totalkmers = 1;
    std::vector<char *>* regions = NULL;
    
    int ifhla = 0;
    int ifunmap = 0;
    std::mutex counting_lock;
    std::mutex lock;

};

template <typename T1, typename T2>
static void initiate_counter(T2 &kmer_hash, T1 &larger_kmer, uint &kindex)
{
    
    auto map_find = kmer_hash.add(larger_kmer, kindex);

    if (*map_find == kindex)
    {
        kindex++;
    }
}

template <typename T1, typename T2, typename T3>
static void initiate_counter_mul(T2 &kmer_hash, T3 &kmer_multi_hash, T1 &larger_kmer, uint &kindex)
{
    auto map_find = kmer_hash.add(larger_kmer, kindex);

    if (*map_find == kindex)
    {
        kindex++;
    }
    
    else if (*map_find == UINT_MAX )
    {
        uint* &data = kmer_multi_hash.find(larger_kmer)->second;
        if (data[0]%5 == 2)
        {
            data = (uint*) realloc(data, sizeof(uint)*(data[0]+1 + FIXCOL));
        }
            
        data[0]++;
        data[data[0]] = kindex++;
    }
    
    else
    {
        uint* newarray = (uint*) malloc(sizeof(uint)*(3));
        newarray[0] = 2;
        newarray[1] = *map_find;
        newarray[2] = kindex++;
        
        kmer_multi_hash[larger_kmer] = newarray;
        
        *map_find = UINT_MAX ;
        
    }

}



template <bool EnableLock,typename T1, typename T2, typename T3>
static void update_counter(T2 &kmer_hash, T3 &kmer_multi_hash, T1 &larger_kmer, uint16* vec, uint &hash_ , mutex &lock)
{
    auto map_find = kmer_hash.find(larger_kmer, hash_);
    
    if (map_find != NULL)
    {
        uint index = *map_find;
        
        if (index < UINT_MAX )
        {
            //temp = vec[index] + 1;
            //vec[index] = temp & -(temp < MAX_UINT16);
            if constexpr (EnableLock)
            {
                std::lock_guard<std::mutex> guard(lock);
                if (__builtin_expect(vec[index] < MAX_COUNT, 0)) vec[index]++;
                                
            }
            else
            {
                if (__builtin_expect(vec[index] < MAX_COUNT, 0)) vec[index]++;
            }

        }
        else
        {
            uint* &data = kmer_multi_hash.find(larger_kmer)->second;
            
            uint num_num = data[0];
            
            if constexpr (EnableLock)
            {
                std::lock_guard<std::mutex> guard(lock);
                for (index = 1 ; index < num_num + 1; ++index)
                {
                    if ( __builtin_expect(vec[data[index]] < MAX_COUNT, 0) ) vec[data[index] ] ++;
                }
            }
            else
            {
                for (index = 1 ; index < num_num + 1; ++index)
                {
                    if ( __builtin_expect(vec[data[index]] < MAX_COUNT, 0) ) vec[data[index] ] ++;
                }
            }
            
        }
    }
}


template <bool EnableLock, typename T1, typename T2, typename T3>
static void update_counter(T2 &kmer_hash, T3 &kmer_multi_hash, T1 &larger_kmer, uint16* vec, uint &hash_ , mutex &lock, const bool strd)
{
    auto map_find = kmer_hash.find(larger_kmer, hash_);
    
    if (map_find != NULL)
    {
        uint index = *map_find;
        
        if (index < UINT_MAX )
        {
            //temp = vec[index] + 1;
            //vec[index] = temp & -(temp < MAX_UINT16);
            if constexpr (EnableLock)
            {
                std::lock_guard<std::mutex> guard(lock);
                if (__builtin_expect((vec[index] & 0x3fff) < MAX_COUNT, 0)) vec[index]++;
                if (strd)
                {
                    vec[index] |= 0x8000;
                }
                else
                {
                    vec[index] |= 0x4000;
                }
                                
            }
            else
            {
                if (__builtin_expect((vec[index] & 0x3fff) < MAX_COUNT, 0)) vec[index]++;
                if (strd)
                {
                    vec[index] |= 0x8000;
                }
                else
                {
                    vec[index] |= 0x4000;
                }
            }

        }
        else
        {
            uint* &data = kmer_multi_hash.find(larger_kmer)->second;
            
            uint num_num = data[0];
            
            if constexpr (EnableLock)
            {
                std::lock_guard<std::mutex> guard(lock);
                for (index = 1 ; index < num_num + 1; ++index)
                {
                    if ( __builtin_expect(( vec[data[index]] & 0x3fff) < MAX_COUNT, 0) ) vec[data[index] ] ++;
                    if (strd)
                    {
                        vec[data[index]] |= 0x8000;
                    }
                    else
                    {
                        vec[data[index]] |= 0x4000;
                    }
                }
            }
            else
            {
                for (index = 1 ; index < num_num + 1; ++index)
                {
                    if ( __builtin_expect(( vec[data[index]] & 0x3fff) < MAX_COUNT, 0) ) vec[data[index] ] ++;
                    if (strd)
                    {
                        vec[data[index]] |= 0x8000;
                    }
                    else
                    {
                        vec[data[index]] |= 0x4000;
                    }
                }
            }
            
        }
    }
}


template <typename T1, typename T2, typename T3>
static void update_counter_val(T2 &kmer_hash, T3 &kmer_multi_hash, T1 &larger_kmer, uint16* vec, uint &hash_ , uint val, mutex &lock)
{
    auto map_find = kmer_hash.find(larger_kmer, hash_);
    
    if (map_find != NULL)
    {
        uint index = *map_find;
        
        if (index < UINT_MAX )
        {
            //temp = vec[index] + 1;
            //vec[index] = temp & -(temp < MAX_UINT16);
            //lock.lock();
            if ( __builtin_expect(vec[index] < MAX_COUNT, 0) )
            {
                vec[index] = MIN( vec[index] + val,  MAX_COUNT );
            }
            //lock.unlock();
        }
        else
        {
            uint* &data = kmer_multi_hash.find(larger_kmer)->second;
            
            uint num_num = data[0];
            
            //lock.lock();
            for (index = 1 ; index < num_num + 1; ++index)
            {
                
                if ( __builtin_expect(vec[data[index]] < MAX_COUNT, 0) )
                {
                    vec[data[index] ] = MIN( vec[index] + val,  MAX_COUNT );
                }
                
                //temp = vec[data[index]] + 1;
                //vec[data[index]] = temp & -(temp < MAX_UINT16);

            }
            //lock.unlock();
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
static void kmer_read_31(char base,  std::size_t &current_size, T &current_kmer, T &reverse_kmer)
{
    int converted = 0;
    T reverse_converted;
    
    const T mask = ((T)1 << (klen*2)) - 1;
    const int shift = 2*klen - 2;
    
    //if (base == '\n' || base == ' ') return;
    
    if (base_to_int(base, converted))
    {
        current_kmer = ((current_kmer << 2) & mask) | converted;
        reverse_kmer = (reverse_kmer >> 2) | ((T)(0b11 ^ converted) << shift);
        
    }
    
    else
    {
        current_size = 0;
        current_kmer = 0;
        reverse_kmer = 0;
    }
}


template <int dictsize>
void KmerCounter<dictsize>::LoadRegion(std::vector<char *> &r, int hla , int unmap )
{
    regions = &r;
    ifhla = hla;
    ifunmap = unmap;
}

template <int dictsize>
void KmerCounter<dictsize>::load_backgrounds(const char * backgroudfile)
{
    FastaReader fastafile(backgroudfile);
    
    std::size_t current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
        
    std::hash<std::uint64_t> hasher;
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

            kmer_read_31(base, ++current_size, current_kmer, reverse_kmer);
            
            if (current_size < klen || base >= 'a') continue;
            
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
            
            backgrounds[larger_kmer%backgrounds.size()] = larger_kmer;
                            
        }
    }
    
    
    fastafile.Close();
    
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

            kmer_read_31(base, ++current_size, current_kmer, reverse_kmer);
            
            if (current_size < klen || base >= 'a') continue;
            
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
            
            initiate_counter(kmer_hash, larger_kmer, totalkmers);
                
        }
    }
    
    
    fastafile.Close();
    
    return totalkmers;

};

template <int dictsize>
ull KmerCounter<dictsize>::read_target(KtableReader &ktablefile)
{
    
    if (ktablefile.file == NULL and ktablefile.file_bgzf == NULL) return 0;
    
    kmer_int current_kmer = 0;
    
    std::string StrLine;
    StrLine.resize(MAX_LINE);
    
    char base;
    int pos;
    int elecounter;
    
    while (ktablefile.nextLine_kmer(StrLine))
    {
        kmer_deconpress(&StrLine.c_str()[THREECOLLENGTH], current_kmer);
        initiate_counter_mul(kmer_hash, kmer_multi_hash, current_kmer, totalkmers);
        
    }
    
    ktablefile.Close();
    
    return totalkmers;
        
};

template <int dictsize>
ull KmerCounter<dictsize>::read_target(KtableReader &ktablefile, vector<pair<size_t,size_t>>& poses)
{
    
    if (ktablefile.file == NULL and ktablefile.file_bgzf == NULL) return 0;
    
    kmer_int current_kmer = 0;
    
    std::string StrLine;
    StrLine.resize(MAX_LINE);
    
    char base;
    int pos;
    int elecounter;
    
    for (const auto& [fpos, knum]: poses)
    {
        ktablefile.Seek(fpos);
        
        int kindex = 0;
        while (ktablefile.nextLine_kmer(StrLine))
        {
            if (kindex++ >= knum) break;

            kmer_deconpress(&StrLine.c_str()[THREECOLLENGTH], current_kmer);
            
            initiate_counter_mul(kmer_hash, kmer_multi_hash, current_kmer, totalkmers);
        }
    }
    
    
    ktablefile.Close();
    
    return totalkmers;
        
};
 
template <int dictsize>
ull KmerCounter<dictsize>::read_target(const char* inputfile, vector<pair<size_t,size_t>>& poses)
{
    int pathlen = (int)strlen(inputfile);
    
    KtableReader readsfile(inputfile);
    read_target(readsfile, poses);
    
    return totalkmers;
}

template <int dictsize>
ull KmerCounter<dictsize>::read_target(const char* inputfile)
{
    KtableReader readsfile(inputfile);
    read_target(readsfile);

    return totalkmers;
}

template <int dictsize>
template <bool EnableLock>
void KmerCounter<dictsize>::count_kmer(FileReader &file, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, vector<uint16> &nBg)
{
    std::size_t current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    //uint64_t ifmasked = 0;
    //int num_masked = 0 ;
    uint hash_;
    std::string StrLine;
    std::string Header;
    StrLine.resize(MAX_LINE);
    
    uint cindex = 0;
    uint count = 0;
    bool lastline  = 0;
    bool getnumber = 0;
    
    file.Load();
    while (file.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '>':
                lastline = 1;
                getnumber = 0;
                break;
            case ' ': case '\n': case '\t': case '#' : case '@':
                continue;
            default:
                getnumber = lastline;
                lastline = 0;
                break;
        }
        
        if (getnumber)
        {
            count = 0;
            for (auto base: StrLine)
            {
                if (base == '\0') break;
                
                if (base == '\n' || base == ' ' || base == '\t') break;
                
                count = 10* count + base - '0';
            }
                        
            update_counter_val(kmer_hash, kmer_multi_hash, current_kmer, samplevecs, hash_,count,counting_lock);
            
            if ( backgrounds[current_kmer%background_prime] == current_kmer  )
            {
                nBg[current_kmer%background_prime] += count;
            }
        }
        else
        {
            current_size = 0;
            for (cindex =0; cindex < StrLine.size(); ++cindex)
            {
                char base = StrLine[cindex];
                
                if (base == '\0') break;
                
                if (base == '\n' || base == ' ' || base == '\t') break;
                
                kmer_read_31(base, ++current_size, current_kmer, reverse_kmer);
                
                if (current_size < klen) continue;
            }
            
            current_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
            count = 0;
            cindex ++;
            for (; cindex < StrLine.size(); ++cindex)
            {
                char base = StrLine[cindex];
                
                if (base == '\0') break;
                
                if (base == '\n' || base == ' ' || base == '\t') break;
                
                count = 10* count + base - '0';
            }
            
            update_counter_val(kmer_hash, kmer_multi_hash, current_kmer, samplevecs, hash_,count,counting_lock);
            
            if ( backgrounds[current_kmer%background_prime] == current_kmer  )
            {
                std::lock_guard<std::mutex> guard(counting_lock);
                nBg[current_kmer%background_prime] += count;
            }
            
        }
    
        nReads+=1;
        if (nReads % 10000000 == 0)
        {
            cerr << "processed " << nReads / 1000000 << "M reads." << endl;
        }
    }
        
    
    return ;
};


template <int dictsize>
template <bool EnableLock>
void KmerCounter<dictsize>::count_kmer(FastaReader &file, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, vector<uint16> &nBg)
{
    std::size_t current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    //uint64_t ifmasked = 0;
    //int num_masked = 0 ;
    uint hash_;
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
        for (int i = 0; i < StrLine.size() ; ++i)
        {
            auto base = StrLine[i];
            
            if (base == '\0') break;
                        
            if (base == '\n' || base == ' ') continue;
 
            kmer_read_31(base, ++current_size, current_kmer, reverse_kmer);
            
            if (current_size < klen) continue;
            
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
            
            update_counter<EnableLock>(kmer_hash, kmer_multi_hash, larger_kmer, samplevecs, hash_,counting_lock);
            
            if ( backgrounds[larger_kmer%background_prime] == larger_kmer  )
            {
                std::lock_guard<std::mutex> guard(lock);
                nBg[larger_kmer%background_prime] ++;
            }
        }
        
        //nBases+=MAX(0, StrLine.length() - klen + 1);
        nReads+=1;
        if (nReads % 10000000 == 0)
        {
            cerr << "processed " << nReads / 1000000 << "M reads." << endl;
        }
    }
        
    return ;
};



template <int dictsize>
template <bool EnableLock>
void KmerCounter<dictsize>::count_kmer(FastqReader &file, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, vector<uint16> &nBg)
{
    std::size_t current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    //uint64_t ifmasked = 0;
    //int num_masked = 0 ;

    char lastheader = '0';
    uint hash_ = 0;
    size_t rstart =0 , rlen = 0;
    char *StrLine = NULL;
    char *Header = NULL;
    while (file.nextLine(StrLine,rlen,Header))
    {
        current_size = 0;
        for (int pos = 0; pos < klen - 1; ++pos)
        {
            char base = StrLine[pos];
                
            kmer_read_31(base, ++current_size, current_kmer, reverse_kmer);
        }
        for (int pos = 0; pos < rlen; ++pos)
        {
            char base = StrLine[pos];
                        
            kmer_read_31(base, ++current_size, current_kmer, reverse_kmer);
            
            //if (current_size < klen) continue;
            
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
                        
            bool kmer_strd = ((current_kmer >= reverse_kmer) ? 1 : 0 );  // read is on forward strand
            
            update_counter<EnableLock>(kmer_hash, kmer_multi_hash, larger_kmer, samplevecs, hash_, counting_lock, kmer_strd);
            
            if ( __builtin_expect(backgrounds[larger_kmer%background_prime] == larger_kmer , 0))
            {
                std::lock_guard<std::mutex> guard(lock);
                nBg[larger_kmer%background_prime] ++;
            }

        }
        
        //nBases+=MAX(0, StrLine.length() - klen + 1);
        nReads+=1;
        if (nReads % 10000000 == 0)
        {
            cerr << "processed " << nReads / 1000000 << "M reads." << endl;
        }
    }
        
    return ;
};

template <int dictsize>
template <bool EnableLock>
void KmerCounter<dictsize>::count_kmer(CramReader &file, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, vector<uint16> &nBg)
{
    std::size_t current_size = 0;
    uint hash_ = 0;
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    uint8_t* StrLine = NULL;
    std::string Header;
    size_t rlen = 0;
    char base;
    if (samplevecs == NULL)
    {
        auto SRread = bam_init1();
        while (file.nextLine(StrLine, rlen, SRread))
        {
            current_size = 0;
            for (int pos = 0; pos < rlen; ++pos)
            {
                base = seq_nt16_str[bam_seqi(StrLine,pos)];

                kmer_read_31(base, ++current_size, current_kmer, reverse_kmer);
                
                if (current_size < klen) continue;
                
                auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
                
                if ( backgrounds[larger_kmer%background_prime] == larger_kmer  )
                {
                    std::lock_guard<std::mutex> guard(lock);
                    nBg[larger_kmer%background_prime] ++;
                }
            }
            nReads+=1;
        }
        
        bam_destroy1(SRread);
        
        return;
    }
    
    bool ifbackground = 0;
    auto SRread = bam_init1();
    while (file.nextLine(StrLine, rlen, SRread))
    {
        current_size = 0;
        bool current_strd = !(SRread->core.flag & BAM_FREVERSE);
        for (int pos = 0; pos < rlen; ++pos)
        {
            base = seq_nt16_str[bam_seqi(StrLine,pos)];
            //if (base=='\n') break;
            kmer_read_31(base, ++current_size, current_kmer, reverse_kmer);
            
            if (current_size < klen) continue;
                    
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:  reverse_kmer;
           
            bool kmer_strd = (current_strd)
                ? ((current_kmer >= reverse_kmer) ? 1 : 0 )  // read is on forward strand
                : ((current_kmer >= reverse_kmer) ? 0 : 1 ); // read is on reverse strand
            
            update_counter<EnableLock>(kmer_hash, kmer_multi_hash, larger_kmer, samplevecs, hash_, counting_lock, kmer_strd);
            
            if ( backgrounds[larger_kmer%background_prime] == larger_kmer  )
            {
                std::lock_guard<std::mutex> guard(lock);
                nBg[larger_kmer%background_prime] ++;
            }
        }
    	//nBases+=MAX(0, StrLine.length() - klen + 1);
        nReads+=1;
        if (nReads % 10000000 == 0)
        {
            cerr << "processed " << nReads / 1000000 << "M reads." << endl;
        }
    }
    
    bam_destroy1(SRread);
    
    for (size_t i = 0; i < totalkmers; ++i)
    {
        if (samplevecs[i] < 0xC000 && (samplevecs[i] & 0x3fff) >= 7)
        {
            samplevecs[i] = 0;
        }
        else
        {
            samplevecs[i] &= 0x3fff;  // Clears the two highest bits (15 and 14)
        }
    }
    
    return ;
};


template <int dictsize>
template <class typefile>
void KmerCounter<dictsize>::count_kmer_(typefile &file, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, vector<uint16> &nBg, const int nthreads)
{
    if (nthreads > 1 )
    {
        std::vector<std::thread> threads;
        
        for (int i = 0; i < nthreads; ++i)
        {
            threads.emplace_back([this, &file, &samplevecs, &nBases, &nReads, &nBg]() {
                this->template count_kmer<true>(file, samplevecs, nBases, nReads, nBg);
            });
        }
        
        for (auto& t : threads) t.join();
    }
    else
    {
        this->template count_kmer<false>(file, samplevecs, nBases, nReads, nBg);
    }
}

template <int dictsize>
void KmerCounter<dictsize>::QuickCoverage(const char* file, vector<uint16>& nBg, const int nthreads)
{
    CramReader readsfile(file, refpath, nthreads);
    std::vector<char*> depth_region_dummy = DEPTH_REGION;
    
    readsfile.useunmap = 2;
    readsfile.LoadRegion(depth_region_dummy, 0, 0);
    ull_atom nBases =0, nReads=0;
    count_kmer_(readsfile, NULL, nBases, nReads, nBg, nthreads);
    readsfile.Close();
}
template <int dictsize>
void KmerCounter<dictsize>::Call(const char* inputfile, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, vector<uint16> &nBg, const int nthreads)
{
    
    vector<string> files;
    if (std::filesystem::is_directory(inputfile))
    {
        for (const auto& entry : std::filesystem::directory_iterator(inputfile))
        {
            if (std::filesystem::is_regular_file(entry.path()))
            {
                files.push_back(entry.path().string());
            }
        }
    }
    else
    {
        files.push_back(inputfile);
    }
    
    for (string& filestring: files)
    {
        const char* file = filestring.c_str();
        
        int pathlen = (int)strlen(file);
        
        if ( pathlen >= 3 && strcmp(file+(pathlen-3),".gz") == 0 || (pathlen >= 6 && strcmp(file+(pathlen-6),".fastq") == 0 ))
        {
            FastqReader readsfile(file);
            count_kmer_(readsfile, samplevecs, nBases, nReads, nBg, nthreads);
            readsfile.Close();
        }
            
        else if ( (pathlen >= 3 && strcmp(file+(pathlen-3),".fa")==0) || (pathlen >= 6 && strcmp(file+(pathlen-6),".fasta") == 0 ))
        {
            FastaReader readsfile(file);
            
            count_kmer_(readsfile, samplevecs, nBases, nReads, nBg, 1);
            readsfile.Close();
        }
        else if ( ( pathlen >= 5 && ( strcmp(file+(pathlen-5),".cram") == 0 )) || (pathlen >= 4 && strcmp(file+(pathlen-4),".bam") == 0 ) ||  (pathlen >= 4 && strcmp(file+(pathlen-4),".sam") == 0 )    )
        {
                        
            if ( regions != NULL && regions->size() && inputdep == 0)
            {
                QuickCoverage(file, nBg, nthreads);
            }
            
            CramReader readsfile(file, refpath, nthreads);
            
            readsfile.LoadRegion(*regions, ifhla, ifunmap);
            
            if ( regions != NULL && regions->size() && inputdep == 0 )
            {
                std::vector<uint16> nBg_ (background_prime, 0);
                
                count_kmer_(readsfile, samplevecs, nBases, nReads, nBg_, nthreads);
            }
            else
            {
                count_kmer_(readsfile, samplevecs, nBases, nReads, nBg, nthreads);
            }	
            
            readsfile.Close();

            for (size_t i = 0; i < totalkmers; ++i)
            {
                if (samplevecs[i] < 0xC000 && (samplevecs[i] & 0x3fff) >= 7)
                {
                    samplevecs[i] = 0;
                }
                else
                {
                    samplevecs[i] &= 0x3fff;  // Clears the two highest bits (15 and 14)
                }
            }
    	
        }
        
        else if ( ( pathlen >= 4 && ( strcmp(file+(pathlen-4),".txt") == 0 )) || (pathlen >= 3 && strcmp(file+(pathlen-3),".jy") == 0 )  )
        {
            FileReader readsfile(file);
            count_kmer_(readsfile, samplevecs, nBases, nReads, nBg, nthreads);
            readsfile.Close();
        }
    }
    
}




#endif /* KmerCounter_hpp */
