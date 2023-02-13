//
//  KmerCounter.hpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 2/2/23.
//  Copyright © 2023 USC_Mark. All rights reserved.
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
#include <thread>
#include <atomic>
#include <mutex>
#include <thread>
#include <atomic>
#include <unordered_map>
#include <tuple>

#include "fasta.hpp"
#include "fastq.hpp"
//#include "CramReader.hpp"
#include "ktable.hpp"
#include "Kmatrix.hpp"

using namespace std;

typedef unsigned int uint;
typedef unsigned __int128 u128;
typedef unsigned long long ull;
typedef unsigned short  uint16;
typedef unsigned char  uint8;

#define large_prime 2147483647

struct hash_128
{
    
    size_t operator()(const u128& num128) const
    {
        auto hash1 = hash<ull>{}((ull)(num128%large_prime));
        return hash1;
    }
};

static bool base_to_int(char base, int &converted)
{

    switch (base)
    {
        case 'A' : case 'a':
            converted=0b00;
            break;
        case 'T' : case 't':
            converted=0b11;
            break;
        case 'C' : case 'c':
            converted=0b01;
            break;
        case 'G' : case 'g':
            converted=0b10;
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

template <int dictsize>
class kmer_counter
{
    using kmer_int = typename std::conditional<(dictsize>32), u128, ull>::type;
    using kmer_dict_type = typename std::conditional<(dictsize>32), kmer64_dict, kmer32_dict>::type;
    using kmer_dict_type_nt = typename std::conditional<(dictsize>32), kmer64_dict_nt, kmer32_dict_nt>::type;
    
    kmer_dict_type target_hash;
    
    ull totalkmers = 0;

    
    int klen = 31 , knum = 0;
    bool iftarget = 0;
    uint targetindex = 0;
    std::atomic_uint restfileindex ;
    std::mutex Threads_lock;
    
    std::vector<std::thread*> threads;
    std::vector<std::string> inputfiles;
    std::vector<std::string> outputfiles;
    std::vector<std::string> prefixes;
    std::vector<float> depths;
    
    const char* mfile = "";
public:
    
    kmer_counter (int kmersize): klen(31)
    {
    };
    ~kmer_counter()
    {
    };
    
    void read_counttarget(ktable &fastafile);
    
    void read_counttarget(fasta &fastafile);
    
    void read_target(const char* infile, const char* mfile);
        
    template <class typefile>
    ull count_kmer(typefile &fastafile, uint16* samplevecs);
    
    void read_files(std::vector<std::string>& inputfiles, std::vector<std::string>& outputfiles, std::vector<std::string>& prefixes,std::vector<float>& deps,int numthread);
    
    void read_file();
            
    void write(const char * outputfile, const char* input, Kmatrix* matrix);

};

template <typename T1, typename T2>
static bool initiate_counter(T2 &target_hash, T1 &larger_kmer, ull &kindex)
{
    
    typename T2::iterator map_find = target_hash.find(larger_kmer);

    if (map_find == target_hash.end())
    {
    
        target_hash[larger_kmer] = (int) kindex++;
        
        return false;
    }
    
    return true;
}


template <typename T1, typename T2>
static void update_counter(T2 &target_hash, T1 &larger_kmer, uint16* vec)
{
    
    typename T2::iterator map_find = target_hash.find(larger_kmer);

    if (map_find != target_hash.end())
    {
        uint index = map_find->second;
        
        if ( vec[index] < 65535) vec[index] ++;
    }

}


template <typename T>
static void kmer_read_c(char base, int klen, std::size_t &current_size, T &current_kmer, T &reverse_kmer)
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
void kmer_counter<dictsize>::read_counttarget(fasta &fastafile)
{
    
    std::size_t current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
        
    std::string StrLine;
    
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
            if (initiate_counter(target_hash, larger_kmer, totalkmers)) continue;
                
        }
    }
    
    
    fastafile.Close();

};

template <int dictsize>
void kmer_counter<dictsize>::read_counttarget(ktable &fastafile)
{
        
    std::size_t current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    iftarget = 1;
    
    std::string StrLine;
    char base;
    int pos;
    while (fastafile.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '#':
                //matrix.LoadHeader(StrLine, StrLine.c_str(), StrLine.length());
                continue;
            case ' ': case '\n': case '\t':
                continue;
            default:
                break;
        }
        
        for (pos = 0; pos < klen; ++pos)
        {
            base = StrLine[pos];
            if (base == '\0' || base == '\n' || base == ' ' || base == '\t') continue;
            kmer_read_c(base, klen, current_size, current_kmer, reverse_kmer);
        }
        
        if (pos < klen) continue;
        
        auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
        
        initiate_counter(target_hash, larger_kmer, totalkmers);
        
        pos++;
        //matrix.LoadRow(StrLine.c_str() + pos , StrLine.length() - pos );
        
    }
    
    fastafile.Close();
        
};
 
template <int dictsize>
void kmer_counter<dictsize>::read_target(const char* inputfile, const char* matrixfile)
{
    int pathlen = (int)strlen(inputfile);
    
    if ( (pathlen > 2 && strcmp(inputfile+(pathlen-3),".fa")==0) || (pathlen > 6 && strcmp(inputfile+(pathlen-6),".fasta") == 0 ))
    {
        fasta readsfile(inputfile);
        read_counttarget(readsfile);
    }
    else
    {
        ktable readsfile(inputfile);
        read_counttarget(readsfile);
    }
    mfile = matrixfile;
    
}



template <int dictsize>
template <class typefile>
ull kmer_counter<dictsize>::count_kmer(typefile &fastafile, uint16* samplevecs)
{
    ull totalsamplekmers = 0;
    std::size_t current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    //uint64_t ifmasked = 0;
    //int num_masked = 0 ;
    std::string StrLine;
    std::string Header;
    
    while (fastafile.nextLine(StrLine))
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
            
            update_counter(target_hash, larger_kmer, samplevecs);

        }
    }
        
    fastafile.Close();
    return totalsamplekmers;
};

template <int dictsize>
void kmer_counter<dictsize>::read_files(std::vector<std::string>& inputs, std::vector<std::string>& outputs, std::vector<std::string>& prefs, std::vector<float>& deps, int nthreads)
{
    
    inputfiles = inputs;
    outputfiles = outputs;
    prefixes = prefs;
    restfileindex = 0;
    depths = deps;
    
    std::vector<std::thread*> threads;
    
    for(int i=0; i< nthreads; ++i)
    {
        std::thread *newthread_ = new std::thread(&kmer_counter<dictsize>::read_file, this);
        threads.push_back(newthread_);
    }
    
    
    for(int i=0; i< nthreads; ++i)
    {
        threads[i]->join();
    }
    
}


template <int dictsize>
void kmer_counter<dictsize>::read_file()
{
    
    while (restfileindex < inputfiles.size())
    {
        
        Threads_lock.lock();
        
        int inputindex = restfileindex++ ;
        
        Threads_lock.unlock();
        
        if (inputindex >= inputfiles.size()) break;
        
        float depth;
        if (inputindex <= depths.size())
        {
            depth= depths[inputindex];
        }
        else
        {
            depth = 14.0;
        }
        
        Kmatrix* matrix = new Kmatrix(mfile, depth, totalkmers);
        uint16* samplevecs = matrix->kmervec;
                
        const char* inputfile = inputfiles[inputindex ].c_str();
        const char* outputfile = outputfiles[inputindex ].c_str();
        string prefix;
        if (prefixes.size()>inputindex)
        {
            prefix = string(prefixes[inputindex]);
        }
        else
        {
            prefix = string(inputfile);
        }
        
        int pathlen = (int)strlen(inputfile);
        
        if ( pathlen > 2 && strcmp(inputfile+(pathlen-3),".gz") == 0 )
        {
            fastq readsfile(inputfile);
            count_kmer(readsfile, samplevecs);
        }
            
        else if ( (pathlen > 2 && strcmp(inputfile+(pathlen-3),".fa")==0) || (pathlen > 6 && strcmp(inputfile+(pathlen-6),".fasta") == 0 ))
        {
            fasta readsfile(inputfile);
            count_kmer(readsfile, samplevecs);
        }
        else if (pathlen > 5 && ( strcmp(inputfile+(pathlen-3),".cram") == 0 ))
        {
            //CramReader readsfile(inputfile);
            //count_kmer(readsfile, samplevecs);
        }
        
        matrix->Process();
        
        matrix->write(outputfile, prefix.c_str());
        
        delete matrix;
        
    }
        
}



/*
template <int dictsize>
void kmer_counter<dictsize>::write(const char * outputfile, const char* input, Kmatrix* matrix)
{
    
    FILE *fwrite=fopen(outputfile, "w");
    
    
    if (fwrite==NULL)
    {
        std::cerr << "ERROR: Cannot write file: " << outputfile << endl;
        
        std::_Exit(EXIT_FAILURE);
    }
        
    fprintf(fwrite,"@%s\n", input);
    for (int i = 1; i < (matrix->matrixinfo).size(); ++i)
    {
        string name = get<0>(matrix->matrixinfo[i]);
        auto totalkmer = matrix->totalkmers[i];
        uint genenum = (uint)get<3>(matrix->matrixinfo[i]);
        double* kernal = matrix->kernals[i];
        double* weightnorm = matrix->weightnorms[i];
        
        
        fprintf(fwrite,">%s\t%.4lf\n", name.c_str(), totalkmer);
        
        for (int j = 0; j < genenum ; j++)
        {
            fprintf(fwrite,"%.4lf,", kernal[j]);
        }
        
        for (int j = 0; j < genenum * genenum ; j++)
        {
            if (j % genenum == 0) fprintf(fwrite,"\n");
            fprintf(fwrite,"%.4lf,", weightnorm[j]);
        }
        fprintf(fwrite,"\n");
        
    }
    
    fclose(fwrite);
    
    return ;
}
*/


#endif /* KmerCounter_hpp */
