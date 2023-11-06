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
#include <filesystem>

#include "FastaReader.hpp"
#include "FastqReader.hpp"
#include "CramReader.hpp"
#include "KtableReader.hpp"
#include "KmerHash.hpp"

struct Hash10M {
    std::size_t operator()(ull key) const {
        return static_cast<std::size_t>(key % prime10M);
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
    using kmer_int = typename std::conditional<(dictsize>32), u128, ull>::type;
    using kmer_hash_type = typename std::conditional<(dictsize>32), Kmer64_hash, Kmer32_hash>::type;
    using kmer_hash_type_mul = typename std::conditional<(dictsize>32), kmer64_dict_mul, kmer32_dict_mul>::type;
 
public:
    
    KmerCounter (size_t size):kmer_hash(4*size)
    {ifbackgrounds.resize(kmer_hash.modsize, 0);};
    ~KmerCounter()
    {};
    
    ull read_target(KtableReader &fastafile);
    
    ull read_target(FastaReader &fastafile);
    
    ull read_target(const char* infile);
	
    void load_backgrounds(const char * backfile);
        
    template <class typefile>
    void count_kmer_(typefile &file, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg, const int nthreads);
    
    void count_kmer(CramReader &file, uint16* samplevecs , ull_atom &nBases, ull_atom &nReads ,ull_atom &nBg);
    
    void count_kmer(FastaReader &file, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg);
   
    void count_kmer(FastqReader &file, uint16* samplevecs , ull_atom &nBases, ull_atom &nReads ,ull_atom &nBg);
 
    void read_files(std::vector<std::string>& inputfiles, std::vector<std::string>& outputfiles, std::vector<std::string>& prefixes,std::vector<float>& deps,int numthread);
    
    void Call(const char* infile, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg, const int nthreads);
    
    kmer_hash_type kmer_hash;
    kmer_hash_type_mul kmer_multi_hash;
    
    vector<bool> ifbackgrounds;
    unordered_set<ull, Hash10M> backgrounds;
    
    uint totalkmers = 0;
    std::vector<char *> regions;
    
    const int klen = 31;
    
    std::mutex counting_lock;

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


template <typename T1, typename T2, typename T3>
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
            //lock.lock();
            if ( __builtin_expect(vec[index] < MAX_UINT16 - 1, 0) )
            {
                vec[index] ++;
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
                
                if ( __builtin_expect(vec[data[index]] < MAX_UINT16 - 1, 0) )
                {
                    vec[data[index] ] ++;
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
    const int klen = 31;
    
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
        current_size = -1;
        current_kmer = 0;
        reverse_kmer = 0;
    }
}

/*
template <int dictsize>
void KmerCounter<dictsize>::LoadRegion(std::vector<char *> &r)
{
    this->regions = r;
}
*/

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

            kmer_read_31(base, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen || base >= 'a') continue;
                                
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
            
            backgrounds.insert(larger_kmer);
            
            ull key = hasher(larger_kmer);
            uint hash_ = key % kmer_hash.modsize;
            
            ifbackgrounds[hash_] = 1;
                
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
   
    uint hash_; 
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

            kmer_read_31(base, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen || base >= 'a') continue;
                                
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
    
    kmer_int current_kmer = 0;
        
    std::string StrLine;
    StrLine.resize(MAX_LINE);
    
    char base;
    int pos;
    int elecounter;
    while (ktablefile.nextLine_kmer(StrLine))
    {
        pos = 0;
        
        elecounter = 0;
        for (; pos < klen; ++pos)
        {
            base = StrLine[pos];
            if (base == '\t' && ++elecounter > FIXCOL-3) break;
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
void KmerCounter<dictsize>::count_kmer(FastaReader &file, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg)
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
        for (auto base: StrLine)
        {
            
            if (base == '\0') break;
                        
            if (base == '\n' || base == ' ') continue;
 
            kmer_read_31(base, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen) continue;
            
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
            
            update_counter(kmer_hash, kmer_multi_hash, larger_kmer, samplevecs, hash_,counting_lock);
            
	    if ( __builtin_expect( ifbackgrounds[hash_], 0 ) &&  __builtin_expect(backgrounds.find(larger_kmer)!=backgrounds.end() , 1)  )
	    {
		nBg ++;
            }
            
        }
        
    //nBases+=MAX(0, StrLine.length() - klen + 1);
    nReads+=1;
    if (nReads % 10000000 == 0) {
      cerr << "processed " << nReads / 1000000 << "M reads." << endl;
    }
    }
        
    
    return ;
};



template <int dictsize>
void KmerCounter<dictsize>::count_kmer(FastqReader &file, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg)
{
    std::size_t current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    //uint64_t ifmasked = 0;
    //int num_masked = 0 ;
    std::string Header;

    uint hash_ = 0;
    size_t rlen = 0;
    const char *StrLine = NULL;
    
    while (file.nextLine_prt(StrLine, rlen))
    {
        for (int pos = 0; pos < rlen; ++pos)
        {
            char base = StrLine[pos];
            
            //if (base=='\n') break;
            
            kmer_read_31(base, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen) continue;
            
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
            
            update_counter(kmer_hash, kmer_multi_hash, larger_kmer, samplevecs, hash_, counting_lock);
            
	    if ( __builtin_expect( ifbackgrounds[hash_], 0 ) &&  __builtin_expect(backgrounds.find(larger_kmer)!=backgrounds.end() , 1)  )
            {
            	nBg ++;
            }
        }
        
        //nBases+=MAX(0, StrLine.length() - klen + 1);
    nReads+=1;
    if (nReads % 10000000 == 0) {
      cerr << "processed " << nReads / 1000000 << "M reads." << endl;
    }
    }
        
    return ;
};

template <int dictsize>
void KmerCounter<dictsize>::count_kmer(CramReader &file, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg)
{
    std::size_t current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    //uint64_t ifmasked = 0;
    //int num_masked = 0 ;
    std::string Header;
    
    uint hash_ ;
    size_t rlen = 0;
    uint8_t *StrLine = NULL;
    bool ifbackground = 0;
    char base; 
    while (file.nextLine_prt(StrLine, rlen))
    {
        for (int pos = 0; pos < rlen; ++pos)
        {
	    base = seq_nt16_str[bam_seqi(StrLine,pos)];
            //if (base=='\n') break;
            kmer_read_31(base, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen) continue;
            
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;

            update_counter(kmer_hash, kmer_multi_hash, larger_kmer, samplevecs, hash_, counting_lock);
           
            if ( __builtin_expect( ifbackgrounds[hash_], 0 ) &&  __builtin_expect(backgrounds.find(larger_kmer)!=backgrounds.end() , 1)  )
            {
                nBg ++;
            }

        }
        
    	//nBases+=MAX(0, StrLine.length() - klen + 1);
	nReads+=1;
	if (nReads % 10000000 == 0) {
	  cerr << "processed " << nReads / 1000000 << "M reads." << endl;
	}
    }
        
    return ;
};


template <int dictsize>
template <class typefile>
void KmerCounter<dictsize>::count_kmer_(typefile &file, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg, const int nthreads)
{
    
    std::vector<std::thread> threads;
        
    for(int i=0; i< nthreads; ++i)
    {
        threads.push_back(std::thread([this, &file, &samplevecs, &nBases, &nReads, &nBg]()
        {
            this->count_kmer(file, samplevecs, nBases, nReads, nBg);
        }));
    }
    
    for(int i=0; i< nthreads; ++i)
    {
        threads[i].join();
    }
    
};

template <int dictsize>
void KmerCounter<dictsize>::Call(const char* inputfile, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg, const int nthreads)
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
        
        if ( pathlen > 2 && strcmp(file+(pathlen-3),".gz") == 0 )
        {
            
            FastqReader readsfile(file);
            count_kmer_(readsfile, samplevecs, nBases, nReads, nBg, nthreads);
            readsfile.Close();
        }
            
        else if ( (pathlen > 2 && strcmp(file+(pathlen-3),".fa")==0) || (pathlen > 6 && strcmp(file+(pathlen-6),".fasta") == 0 ))
        {
            FastaReader readsfile(file);
            
            count_kmer_(readsfile, samplevecs, nBases, nReads, nBg, 1);
            readsfile.Close();
        }
        else if (pathlen > 5 && ( strcmp(file+(pathlen-5),".cram") == 0 ))
        {
            CramReader readsfile(file);
            readsfile.LoadRegion(regions);
            count_kmer_(readsfile, samplevecs, nBases, nReads, nBg, nthreads);
            readsfile.Close();
        }
    }
    
}


#endif /* KmerCounter_hpp */
