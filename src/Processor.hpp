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
#include <chrono>

#include "config.hpp"
#include "FastaReader.hpp"
#include "FastqReader.hpp"
#include "KtableReader.hpp"
#include "KmerMatrix.hpp"
#include "KmerCounter.hpp"
#include "KmerWindow.hpp"
#include "Regression.hpp"
#include "PriorData.hpp"


#define DefaultSize 2000


using namespace std;



template <int ksize>
class Genotyper
{
    using kmer_int = typename std::conditional<(ksize>32), u128, ull>::type;
    using kmer_hash_type = typename std::conditional<(ksize>32), Kmer64_hash, Kmer32_hash>::type;
    using kmer_hash_type_mul = typename std::conditional<(ksize>32), kmer64_dict_mul, kmer32_dict_mul>::type;
    
public:
    
    unique_ptr<int> results;
    unique_ptr<FLOAT_T> reminders;
    unique_ptr<FLOAT_T> coefs;
    unique_ptr<FLOAT_T> residuels;
    const size_t knum, pnum;
    const uint window;
    const int Nsubthreads;
    Genotyper(size_t k, size_t p, KmerCounter<ksize> &c ,PriorData &priordata, const int w, const int N):
    knum(k),
    pnum(p),
    window(w),
    counter(c),
    priordata_manager(priordata),
    kmer_counts(new uint16[k+1]),
    Nsubthreads(N),
    
    norm_vec(new FLOAT_T[DefaultSize]),
    norm_matrix(new FLOAT_T[DefaultSize*DefaultSize]),
    
    coefs(new FLOAT_T[MAX_UINT16]),
    residuels(new FLOAT_T[MAX_UINT16]),
    reminders(new FLOAT_T[MAX_UINT16]),
    results(new int[MAX_UINT16]),
    
    finished_group(p)
    {};
    
    void counting(const std::string& inputfile)
    {

        cerr << "counting kmers for sample: " << inputfile<<endl;

        auto begin = std::chrono::high_resolution_clock::now();
        
        //ull_atom totalbases_atom = 0, totalreads_atom = 0, totalbgs_atom = 0;
		
        counter.Call(inputfile.c_str(), kmer_counts.get(), totalbases, totalreads, totalbgs, Nsubthreads);
        
        //totalbases= totalbases_atom ;
        //totalreads = totalreads_atom ;
        //totalbgs = totalbgs_atom ;

        finishcounting = 1;

	auto end = std::chrono::high_resolution_clock::now();
        
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            
        cerr<<"finished counting "<< inputfile <<" at time: "<<elapsed.count()* 1e-9 <<endl;
        
    };
    
    void runOneGroup(const PriorChunk* priorData, const std::string& inputfile, const std::string& outputfile, const float depth, std::mutex& Threads_lock)
    {
        
        cout << "generating kmer matrix for sample: " << inputfile<<endl;
        matrix.getNorm(&kmer_counts.get()[priorData ->kmervec_start], priorData->kmer_matrix, depth, priorData->genenum, priorData->kmervec_size,
                       norm_vec.get(), norm_matrix.get(), total_lambda);
        
        cout << "regressing to references for sample: " << inputfile<<endl;
        regresser.Call(&kmer_counts.get()[priorData ->kmervec_start], priorData->kmer_matrix, depth, priorData->genenum, priorData->kmervec_size, norm_vec.get(), norm_matrix.get(),  total_lambda, priorData->gene_kmercounts, coefs.get(), residuels.get(), priorData->numgroups, priorData->genegroups,priorData->groupkmernums);

        cout << "rounding for sample: " << inputfile<<endl;
        
        tree.Run(priorData->phylo_tree, coefs.get(), gnum, &results.get()[0], &reminders.get()[0], residuels.get(), norm_matrix.get());

        cout << "determine window residuels: " << inputfile<<endl;

    
        KmerWindow kmerwindow(window);
        kmerwindow.resize(priorData->pathsizes);
        kmerwindow.WindowCovers(&kmer_counts.get()[priorData ->kmervec_start], priorData->kmer_matrix, depth, priorData->genenum, priorData->kmervec_size, priorData->genenum, &results.get()[0], total_obs, total_exp);
        
        vector<vector<tuple<int, int, float, string>>> PatialCopies(priorData->pathnames.size()+1);
        kmerwindow.PartialCopy(PatialCopies, &reminders.get()[0], priorData->genenames, priorData->pathnames, depth);
        
        write(priorData, outputfile, inputfile, priorData->prefix, priorData->genenames, PatialCopies, kmerwindow.windowcovers, depth, Threads_lock);

        cout<<"finish run"<<endl;
    };
    
    void write(const PriorChunk* priorData, const std::string& outputfile, const string &sample, const string &prefix, const vector<string>&genenames_ori, const vector<vector<tuple<int, int, float, string>>>& PatialCopies, const vector<vector<tuple<int,int,int>>>& windowcovers, const float depth, std::mutex& Threads_lock)
    {
        
        auto genenames(genenames_ori);
        for (auto &genename: genenames)
        {
            genename = genename.substr(0,genename.find('\t', 0));
        }
        
        std::unique_lock<std::mutex> lck(Threads_lock);
        
        FILE *fwrite;
        
        if (outputfile != "stdout")
        {
            fwrite=fopen(outputfile.c_str(), "a");
        }
        else
        {
            fwrite=stdout;
        }
        
        if (fwrite==NULL)
        {
            std::cerr << "ERROR: Cannot write file: " << outputfile << endl;
            
            std::_Exit(EXIT_FAILURE);
        }
        
        fprintf(fwrite,">%s\t%s\n", prefix.c_str(), sample.c_str());
        //fprintf(fwrite,"rsdl: %.4lf\n", regress.first);
        fprintf(fwrite,"lambda: %llu/%llu\n",total_obs, total_exp);
        
        fprintf(fwrite,"regress: ");
        const float cutoff = 0.5 / (gnum + 1);
        for (int i = 0; i < gnum; ++i)
        {
            if (coefs.get()[i] > cutoff) fprintf(fwrite,"%s:%.2lf,", genenames[i].c_str(),coefs.get()[i]);
        }
        fprintf(fwrite,"\n");

	/*
        fprintf(fwrite,"reproject: ");
        for (int i = 0; i < gnum; ++i)
        {
            auto total =reminders.get()[i] + results.get()[i];
            
            string info = "";
            if (priorData->gene_kmercounts[i] < 1000) info = "(aux)";
            
            if ( total > 0.01) fprintf(fwrite,"%s:%.2lf%s,", genenames[i].c_str(),total, info.c_str());
        }
        fprintf(fwrite,"\n");
	*/
        
        fprintf(fwrite,"round: ");
        for (int i = 0; i < gnum; ++i)
        {
            int result = results.get()[i];
            
            if (result > 0)
            {
                
                for (int j = 0 ; j < result  ; ++j)
                {
                    string info = "";
                    if (priorData->gene_kmercounts[i] < 1000) info = "(aux)";
        
                    fprintf(fwrite,"%s%s,", genenames[i].c_str(), info.c_str());
                }
            }
        }
        fprintf(fwrite,"\n");
        
        fprintf(fwrite,"result: ");
        for (int i = 0; i < gnum; ++i)
        {
            int result = results.get()[i];
            
            if (result  > 0)
            {
                for (int j = 0 ; j < result ; ++j)
                {
                    if (priorData->gene_kmercounts[i] >= 1000) fprintf(fwrite,"%s,", genenames[i].c_str());
                }
            }
        }
        fprintf(fwrite,"\n");
        
        
        for (int path = 1; path < PatialCopies.size(); ++path)
        {
            auto& patials = PatialCopies[path];
            
            std::string pathname = "";
            if (priorData->pathnames.size() > path)
            {
                pathname = priorData->pathnames[path];
            }
            else
            {
                pathname = string("Path") + to_string(path);
            }
            
            for (auto &patial : patials)
            {
                
                if (get<1>(patial) - get<0>(patial) > 100) fprintf(fwrite,"Partial\t%s\t%s\t%d\t%d\t%d\n", get<3>(patial).c_str(), pathname.c_str(), 30*(get<0>(patial)-1), 30*(get<1>(patial)-1),  (int) get<2>(patial));
            }
        }
        
        
        
        
        for (int path = 1; path < windowcovers.size(); ++path )
        {
            std::string pathname = "";
            if (priorData->pathnames.size() >= path)
            {
                pathname = string("Path") + to_string(path) + string("\t") + priorData->pathnames[path];
            }
            else
            {
                pathname = string("Path") + to_string(path);
            }
            
            auto &windowcover = windowcovers[path];
            
            ull totalwindow = 0;
            for (auto &thepair: windowcover)
            {
                totalwindow += get<1>(thepair);
            }
            if (totalwindow < 100) continue;
            
            fprintf(fwrite,"Coverage %s: ", pathname.c_str());
            
            int lastcpnum = -1;
            for (auto &thepair: windowcover)
            {
                int query = get<0>(thepair);
                int ref = get<1>(thepair);
                int cpnum = get<2>(thepair);
                
                if (ref > 0)
                {
                    if (cpnum != lastcpnum)
                    {
                        fprintf(fwrite,"%d|", cpnum);
                        lastcpnum = cpnum;
                    }
                }
                
                if (query > 0 or ref >0)
                {
                    
                    fprintf(fwrite,"%d/%d,", query, ref);
                }
                else
                {
                    fprintf(fwrite,",", query, ref);
                }
            }
            fprintf(fwrite,"\n");
        }
        
        fclose(fwrite);
        
        return ;
    };

    void newsample()
    {
        memset(kmer_counts.get(), 0, sizeof(uint16) * knum);
        finished_group.assign(pnum , 0);
        finishcounting = 0;
        totalbases = 0;
        totalreads = 0;
        totalbgs = 0;
    }
    
    void newgroup(const PriorChunk* priorData)
    {
	size_t newalloc_size = MAX( DefaultSize, gnum + 10 );
	
	if (MAX(newalloc_size, alloc_size) > DefaultSize)
        {
            norm_vec.reset(new FLOAT_T[newalloc_size]);
            
            try_allocate_unique(norm_matrix, newalloc_size*newalloc_size, newalloc_size*newalloc_size);
            
            coefs.reset(new FLOAT_T[newalloc_size]);
            
            residuels.reset(new FLOAT_T[newalloc_size]);
            
            results.reset(new int[newalloc_size]);
            
            alloc_size = newalloc_size;
        }

	/*
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
            
            coefs.reset(new FLOAT_T[newalloc_size]),
            
            residuels.reset(new FLOAT_T[newalloc_size]),
            
            results.reset(new int[newalloc_size]),
            
            alloc_size = newalloc_size;
        }
        */
        //memset(norm_matrix.get(), 0, sizeof (FLOAT_T) *  gnum * gnum);
        
        memcpy(norm_matrix.get(), priorData->prior_norm, sizeof (FLOAT_T) *  gnum * gnum);
                        
        memset(norm_vec.get(), 0, sizeof (FLOAT_T) * gnum );
        
        memset(coefs.get(), 0, sizeof(FLOAT_T) * gnum);
        
        memset(residuels.get(), 0, sizeof(FLOAT_T) * gnum);
        
        memset(results.get(), 0, sizeof(int) * gnum);
        
        //kmerwindow.resize(priorData->pathsizes);
                        
        total_lambda = 0;
        total_obs = 0;
        total_exp = 0;
                
    };
    
    void run(const std::string& inputfile, const std::string& outputfile, float depth,std::mutex& Threads_lock)
    {
        cerr<<"running for sample: "<<inputfile << endl;

        newsample();
        
        counting(inputfile);
        
        if (depth <= 0)
        {
            depth = ( 0.5 * totalbgs )/counter.backgrounds.size();
        }
        
        FILE *fwrite;
        
        if (outputfile != "stdout")
        {
            fwrite=fopen(outputfile.c_str(), "a");
        }
        else
        {
            fwrite=stdout;
        }
        
        if (fwrite==NULL)
        {
            std::cerr << "ERROR: Cannot write file: " << outputfile << endl;
            std::_Exit(EXIT_FAILURE);
        }
                
        fprintf(fwrite,"@totalreads: %llu, totalbackgrounds: %llu/%llu \n", totalreads, totalbgs, counter.backgrounds.size());
        
        fclose(fwrite);
        
        for (int i = 0; i < pnum; ++i)
        {
            
            auto begin = std::chrono::high_resolution_clock::now();
            
            PriorChunk* priorData = priordata_manager.getNextChunk(finished_group);
            
            cout << "running gene " << priorData->prefix << " for sample " << inputfile << endl ;

            gnum = priorData->genenum;
            
            newgroup(priorData);
                        
            runOneGroup (priorData, inputfile, outputfile, depth, Threads_lock);
    
            finished_group[priorData->index] = 1;
            
            priordata_manager.FinishChunk(priorData);

            auto end = std::chrono::high_resolution_clock::now();

	    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            
	    cout<<"finished sample:" << inputfile << "for gene:"<< priorData->prefix << " for :" << elapsed.count()* 1e-9 <<endl;
        }

	cerr<<"finished sample: "<<inputfile << endl;
    };
    
    

private:
    
    unique_ptr<uint16> kmer_counts;
    unique_ptr<FLOAT_T> norm_vec;
    unique_ptr<FLOAT_T> norm_matrix;
    
    FLOAT_T total_lambda =0;
    ull total_exp = 0, total_obs = 0;
    
    PriorData &priordata_manager;
    KmerCounter<ksize> &counter;
    KmerMatrix matrix;
    Regression regresser;
    TreeRound tree;
    //KmerWindow kmerwindow;
    
    ull totalbases = 0, totalreads = 0, totalbgs = 0;
    size_t gnum;
    size_t alloc_size = DefaultSize;
    bool finishcounting = 0;
    vector<bool> finished_group;
    
};

template <int ksize>
class Processor
{
    using kmer_int = typename std::conditional<(ksize>32), u128, ull>::type;
    using kmer_hash_type = typename std::conditional<(ksize>32), Kmer64_hash, Kmer32_hash>::type;
    using kmer_hash_type_mul = typename std::conditional<(ksize>32), kmer64_dict_mul, kmer32_dict_mul>::type;
    
public:
    Processor
    (
              std::vector<std::string>& infiles,
              std::vector<std::string>& outfiles,
              std::vector<float> &d,
              std::string &mfile,
              std::string &bfile,
              std::unordered_set<std::string> &g,
              std::vector<char *> &r,
              const int w,
              const int n,
              const int N
    ):
    inputfiles(infiles),
    outputfiles(outfiles),
    depths(d),
    genes(g),
    matrixfile(mfile),
    backgroundfile(bfile),
    priordata_manager(mfile, 2*n),
    regions(r),
    window(w),
    nthreads(n),
    Nsubthreads(N)
    {};
    
    
    void Run();
    void Load();
    void Onethread();
    
    const std::vector<std::string>& inputfiles;
    const std::vector<std::string>& outputfiles;
    const std::string &matrixfile;
    const std::string &backgroundfile;
    std::vector<char *> &regions;
    const std::vector<float> &depths;
    const std::unordered_set<std::string> &genes;
    const int window;
    const int nthreads;
    const int Nsubthreads;
    
private:
    uint totalkmers, totalgroups;
    std::atomic_uint restfileindex = {0};
    std::mutex Threads_lock;
    std::mutex Threads_lock2;
    
    //kmer_hash_type kmer_hash;
    //kmer_hash_type_mul kmer_multi_hash;
    //unordered_set<ull, Hash10M> backgroud;
    KmerCounter<ksize> *Counter =NULL;
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
   
	cout<<"reading all kmer targets"<<endl;
    
    //kmer_hash.initiate(4 * priordata_manager.totalkmers);
 
    Counter = new KmerCounter<ksize>(2 * priordata_manager.totalkmers);
    if (backgroundfile.length() > 0) Counter->load_backgrounds(backgroundfile.c_str());
        
    //Counter->LoadRegion(regions);
    totalkmers = Counter->read_target(matrixfile.c_str());
    
    cout<<"finishing reading targets, start genotyping"<<endl;

    std::vector<std::unique_ptr<std::thread>> threads;
    
    for(int i=0; i< nthreads; ++i)
    {
        threads.push_back(std::unique_ptr<std::thread>(new std::thread(&Processor<ksize>::Onethread, this)));
    }
    
    
    for(int i=0; i< nthreads; ++i)
    {
        threads[i].get()->join();
    }
    
    
    
}


template <int ksize>
void Processor<ksize>::Onethread()
{
    
    unique_ptr<Genotyper<ksize>> genotyper = unique_ptr<Genotyper<ksize>>(new Genotyper<ksize>(totalkmers, totalgroups,  *Counter, priordata_manager, window, Nsubthreads));
    
    
    while (restfileindex < inputfiles.size() )
    {
        while(Threads_lock.try_lock());
        
        int inputindex = restfileindex++ ;
        
        Threads_lock.unlock();
        
        if (inputindex >= inputfiles.size() ) break;
        
    
        string outputfile;
        if (inputindex < outputfiles.size())
        {
            outputfile = outputfiles[inputindex];
        }
        else if (outputfiles.size())
        {
            outputfile = outputfiles[outputfiles.size()-1];
        }
        else
        {
            outputfile = "stdout";
        }
                
        float depth = -1;
        if (inputindex < depths.size())
        {
            depth= depths[inputindex];
        }
        else if (depths.size() > 0)
        {
            depth= depths[depths.size() -1];
        }
        
        genotyper.get()->run(inputfiles[inputindex], outputfile, depth, Threads_lock2);
    }
    
    
}


#endif /* Processor_hpp */

