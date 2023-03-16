//
//  Processor.hpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#ifndef Processor_hpp
#define Processor_hpp

#include <atomic>
#include <mutex>
#include <thread>
#include <atomic>
#include <memory>
#include <stdio.h>
#include <string>
#include <vector>

#include "config.hpp"
#include "FastaReader.hpp"
#include "FastqReader.hpp"
//#include "CramReader.hpp"
#include "KtableReader.hpp"
#include "KmerMatrix.hpp"
#include "KmerCounter.hpp"
#include "Regression.hpp"


using namespace std;


template <int ksize>
class Genotyper
{
    using kmer_int = typename std::conditional<(ksize>32), u128, ull>::type;
    using kmer_hash_type = typename std::conditional<(ksize>32), kmer64_dict, kmer32_dict>::type;
    
public:
    Genotyper(size_t k, size_t p, KmerCounter<ksize> &c, PriorData &priordata):knum(k), pnum(p), counter(c), priordata_manager(priordata), kmer_counts(new uint16[k+1]), results(new int[MAX_UINT16]), finished_group(p)
    {};
    
    void counting(const std::string& inputfile)
    {
        cout<<"check2:"<<inputfile<<endl;
        counter.count_kmer(inputfile.c_str(), kmer_counts.get());
        finishcounting = 1;
    };
    
    void runOneGroup(const PriorChunk* priorData, const std::string& inputfile, const std::string& outputfile, const float depth)
    {
        
        matrix.getNorm(priorData, kmer_counts.get(), depth);
        
        cout<<"check8,"<<gnum<<endl;
        
        for (int i = 0 ; i < gnum; ++i)
        {
            cout<<matrix.kernal_vec[i]<<",";
        }
        cout<<endl;
        
        pair<double, double*> regress = regresser.Call(matrix.kernal_vec, matrix.weightnorm, matrix.kmer_counts, priorData->prior_norm,gnum);
        
        cout<<"check9"<<endl;
        
        tree.Run(priorData->phylo_tree, regress.second, gnum, results.get());
        
        cout<<"check10"<<endl;
        
        write(outputfile, inputfile, priorData->prefix ,regress, results.get());
        
        cout<<"check11"<<endl;
    };
    
    void write(const std::string& outputfile, const string &sample, const string &genename, const pair<double, double*> regress, const int* results)
    {
        FILE *fwrite=fopen(outputfile.c_str(), "w");
        
        if (fwrite==NULL)
        {
            std::cerr << "ERROR: Cannot write file: " << outputfile << endl;
            
            std::_Exit(EXIT_FAILURE);
        }
        
        fprintf(fwrite,">%s\t%s\n", genename.c_str(), sample.c_str());
        fprintf(fwrite,"rsdl: %.4lf\n", regress.first);
        
        const float cutoff = 1.0 / (gnum + 1);
        for (int i = 0; i < gnum; ++i)
        {
            if (regress.second[i] > cutoff) fprintf(fwrite,"%d:%.4lf,", i,regress.second[i]);
        }
        fprintf(fwrite,"\n");
        
        for (int i = 0; i < gnum; ++i)
        {
            if (results[i] > 0) fprintf(fwrite,"%d:%d,", i,results[i]);
        }
        
        fprintf(fwrite,"\n");
        
        fclose(fwrite);
        
        return ;
    };

    
    void clear()
    {
        memset(kmer_counts.get(), 0, sizeof(uint16) * knum);
        finished_group.assign(pnum , 0);
        finishcounting = 0;
    };
    
    void run(const std::string& inputfile, const std::string& outputfile, const float depth)
    {
        clear();
        counting(inputfile);
        
        for (int i = 0 ; i < 100000; ++i)
        {
            if (kmer_counts.get()[i]) cout<<kmer_counts.get()[i]<<",";
        }
        cout<<endl;
                
        for (int i = 0; i < pnum; ++i)
        {
            
            PriorChunk* priorData = priordata_manager.getNextChunk(finished_group);
            
            gnum = priorData->genenum;
            
            memset(results.get(), 0, sizeof(int) * gnum);
            
            cout<<"checkx:"<<gnum<<endl;
            
            
            
            cout<<"check7"<<endl;
            
            runOneGroup (priorData, inputfile, outputfile, depth);
        }
                
    };
private:
    unique_ptr<uint16> kmer_counts;
    unique_ptr<int> results;
    vector<bool> finished_group;
    
    const size_t knum, pnum;
    size_t gnum;
    bool finishcounting = 0;
    
    PriorData &priordata_manager;
    KmerCounter<ksize> &counter;
    KmerMatrix matrix;
    Regression regresser;
    TreeRound tree;
};

template <int ksize>
class Processor
{
    
public:
    Processor(std::vector<std::string>& infiles, std::vector<std::string>& outfiles, std::vector<float> &d, std::string &mfile, std::unordered_set<std::string> &g, std::vector<char *> &r,  const int n):inputfiles(infiles), outputfiles(outfiles), depths(d), genes(g), matrixfile(mfile), priordata_manager(mfile) ,regions(r),  nthreads(n)
    {};
    
    
    void Run();
    void Load();
    void Onethread();
    
    const std::vector<std::string>& inputfiles;
    const std::vector<std::string>& outputfiles;
    const std::string &matrixfile;
    const std::vector<char *> &regions;
    const std::vector<float> &depths;
    const std::unordered_set<std::string> &genes;
    const int nthreads;
  
private:
    uint totalkmers, totalgroups;
    std::atomic_uint restfileindex{0};    
    std::mutex Threads_lock;
    
    KmerCounter<ksize> Counter;
    PriorData priordata_manager;
};


template <int ksize>
void Processor<ksize>::Run()
{
    if (genes.size() > 0)
    {
        totalgroups = priordata_manager.LoadIndex(genes);
    }
    else
    {
        totalgroups = priordata_manager.LoadIndex();
    }
    
    totalkmers = Counter.read_target(matrixfile.c_str());
    
    cout<<"check3:"<<Counter.kmer_hash.size()<<endl;
    std::vector<std::thread*> threads;
    
    for(int i=0; i< nthreads; ++i)
    {
        std::thread *newthread_ = new std::thread(&Processor<ksize>::Onethread, this);
        threads.push_back(newthread_);
    }
    
    
    for(int i=0; i< nthreads; ++i)
    {
        threads[i]->join();
    }
}


template <int ksize>
void Processor<ksize>::Onethread()
{
    
    Genotyper<ksize> genotyper(totalkmers, totalgroups, Counter, priordata_manager);
    
    
    while (restfileindex < inputfiles.size() )
    {
        cout<<"check4"<<endl;
                
        Threads_lock.lock();
        
        int inputindex = restfileindex++ ;
        
        Threads_lock.unlock();
        
        if (inputindex >= inputfiles.size() ) break;
        
        float depth;
        if (inputindex < depths.size())
        {
            depth= depths[inputindex];
        }
        else if (depths.size())
        {
            depth = depths[depths.size() - 1];
        }
        else
        {
            depth = 14.0;
        }
        
        cout<<"check0:"<<inputindex<<endl;
        
        genotyper.run(inputfiles[inputindex], outputfiles[inputindex], depth);
    }
    
}


#endif /* Processor_hpp */
