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
#include "PriorData.hpp"

#define DefaultSize 4000

using namespace std;


template <int ksize>
class Genotyper
{
    using kmer_int = typename std::conditional<(ksize>32), u128, ull>::type;
    using kmer_hash_type = typename std::conditional<(ksize>32), kmer64_dict, kmer32_dict>::type;
    
public:
    
    unique_ptr<int> results;
    unique_ptr<FLOAT_T> coefs;
    unique_ptr<FLOAT_T> residuels;
    const size_t knum, pnum;
    
    Genotyper(size_t k, size_t p, KmerCounter<ksize> &c, PriorData &priordata):
    knum(k),
    pnum(p),
    counter(c),
    priordata_manager(priordata),
    kmer_counts(new uint16[k+1]),
    
    norm_vec(new FLOAT_T[DefaultSize]),
    norm_matrix(new FLOAT_T[DefaultSize*DefaultSize]),
    
    coefs(new FLOAT_T[MAX_UINT16]),
    residuels(new FLOAT_T[MAX_UINT16]),
    
    results(new int[MAX_UINT16]),
    finished_group(p)
    {};
    
    void counting(const std::string& inputfile)
    {
        cout<<"check2:"<<inputfile<<endl;
        counter.Call(inputfile.c_str(), kmer_counts.get(), 1);
        finishcounting = 1;
    };
    
    void runOneGroup(const PriorChunk* priorData, const std::string& inputfile, const std::string& outputfile, const float depth)
    {
        /*
        cout<<priorData->prefix<<endl;
        for (int i = 0; i < 100; ++i)
        {
            cout<<kmer_counts.get()[priorData ->kmervec_start+i]<<",";
        }
        cout<<endl;
        
        for (int i = 0; i < priorData->genenum * priorData->genenum; ++i)
        {
            cout<<norm_matrix.get()[i]<<",";
            if (i%priorData->genenum == 0) cout<<endl;
        }
        cout<<endl;
        */
        
        matrix.getNorm(&kmer_counts.get()[priorData ->kmervec_start], priorData->kmer_matrix, depth, priorData->genenum, priorData->kmervec_size,
                       norm_vec.get(), norm_matrix.get(), total_lambda);
        
        
        for (int i = 0; i < priorData->genenum; ++i)
        {
            cout<<norm_vec.get()[i]<<",";
        }
        cout<<endl;
        
        for (int i = 0; i < priorData->genenum * priorData->genenum; ++i)
        {
            if (i%priorData->genenum == 0) cout<<endl;
            cout<<norm_matrix.get()[i]<<",";
        }
        cout<<endl;
                
        regresser.Call(priorData->genenum, norm_vec.get(), norm_matrix.get(),  total_lambda, priorData->prior_norm, coefs.get(), residuels.get());
        
        cout<<"check9"<<endl;
        
        tree.Run(priorData->phylo_tree, coefs.get(), gnum, results.get());
        
        cout<<"check10"<<endl;
        
        write(outputfile, inputfile, priorData->prefix );
        
        cout<<"check11"<<endl;
    };
    
    void write(const std::string& outputfile, const string &sample, const string &genename)
    {
        FILE *fwrite=fopen(outputfile.c_str(), "w");
        
        if (fwrite==NULL)
        {
            std::cerr << "ERROR: Cannot write file: " << outputfile << endl;
            
            std::_Exit(EXIT_FAILURE);
        }
        
        fprintf(fwrite,">%s\t%s\n", genename.c_str(), sample.c_str());
        //fprintf(fwrite,"rsdl: %.4lf\n", regress.first);
        
        const float cutoff = 1.0 / (gnum + 1);
        for (int i = 0; i < gnum; ++i)
        {
            if (coefs.get()[i] > cutoff) fprintf(fwrite,"%d:%.4lf,", i,coefs.get()[i]);
        }
        fprintf(fwrite,"\n");
        
        for (int i = 0; i < gnum; ++i)
        {
            if (results.get()[i] > 0) fprintf(fwrite,"%d:%d,", i,results.get()[i]);
        }
        
        fprintf(fwrite,"\n");
        
        fclose(fwrite);
        
        return ;
    };

    void newsample()
    {
        memset(kmer_counts.get(), 0, sizeof(uint16) * knum);
        finished_group.assign(pnum , 0);
        finishcounting = 0;
    }
    
    void newgroup(const PriorChunk* priorData)
    {
        size_t newalloc_size = alloc_size;
        
        if(  (  gnum > alloc_size &&
               (newalloc_size = gnum)
             )
           ||
             ( alloc_size > DefaultSize &&
               gnum <= MIN(alloc_size/2, DefaultSize) &&
               ( newalloc_size = MAX(DefaultSize , gnum) )
             )
          )
        {
            norm_vec.reset(new FLOAT_T[newalloc_size]),
            
            norm_matrix.reset(new FLOAT_T[newalloc_size*newalloc_size]),
            
            results.reset(new int[newalloc_size]),
            
            alloc_size = newalloc_size;
        }
        
        memcpy(norm_matrix.get(), priorData->prior_norm, sizeof (FLOAT_T) *  gnum * gnum);
                        
        memset(norm_vec.get(), 0, sizeof (FLOAT_T) * gnum );
        
        memset(results.get(), 0, sizeof(int) * gnum);
        
        total_lambda = 0;

                
    };
    
    void run(const std::string& inputfile, const std::string& outputfile, const float depth)
    {
        newsample();
        
        counting(inputfile);
        
        for (int i = 0; i < pnum; ++i)
        {
            
            PriorChunk* priorData = priordata_manager.getNextChunk(finished_group);
            
            gnum = priorData->genenum;
            
            newgroup(priorData);
            
            cout<<"checkx:"<<gnum<<endl;
            
            cout<<"check7"<<endl;
            
            runOneGroup (priorData, inputfile, outputfile, depth);
        }
                
    };
    
    

private:
    
    unique_ptr<uint16> kmer_counts;
    unique_ptr<FLOAT_T> norm_vec;
    unique_ptr<FLOAT_T> norm_matrix;
    
    double total_lambda;
    
    PriorData &priordata_manager;
    KmerCounter<ksize> &counter;
    KmerMatrix matrix;
    Regression regresser;
    TreeRound tree;
    
    size_t gnum;
    size_t alloc_size = 1000;
    bool finishcounting = 0;
    vector<bool> finished_group;
    
};

template <int ksize>
class Processor
{
    
public:
    Processor(
              std::vector<std::string>& infiles,
              std::vector<std::string>& outfiles,
              std::vector<float> &d,
              std::string &mfile,
              std::unordered_set<std::string> &g,
              std::vector<char *> &r,
              const int n
              ):
    inputfiles(infiles),
    outputfiles(outfiles),
    depths(d),
    genes(g),
    matrixfile(mfile),
    priordata_manager(mfile),
    regions(r),
    nthreads(n)
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
    std::atomic_uint restfileindex = 0;
    std::mutex Threads_lock;
    
    KmerCounter<ksize> Counter;
    Regression regresser;
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
    
    Genotyper genotyper(totalkmers, totalgroups, Counter,  priordata_manager);
    
    
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
