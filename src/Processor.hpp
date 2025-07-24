//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License. 
//  If you use this code, please cite our work.
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
#include "KmerProfile.hpp"
#include "KmerWindow.hpp"
#include "Regression.hpp"
#include "PriorData.hpp"
#include "BedReader.hpp"

#define DefaultSize 1000


extern bool ifgeneregion;
extern int ifprofile;
extern string profilesummaryfile;

using namespace std;


template <int ksize>
class Genotyper
{
    using kmer_int = typename std::conditional<(ksize>32), u128, ull>::type;
    using kmer_hash_type = typename std::conditional<(ksize>32), Kmer64_hash, Kmer32_hash>::type;
    using kmer_hash_type_mul = typename std::conditional<(ksize>32), kmer64_dict_mul, kmer32_dict_mul>::type;
    
public:
    
    unique_ptr<int[]> results;
    unique_ptr<FLOAT_T[]> reminders;
    unique_ptr<FLOAT_T[]> coefs;
    unique_ptr<FLOAT_T[]> residuels;
    size_t alloc_size = DefaultSize;
    unique_ptr<FLOAT_T[]> norm_vec;
    unique_ptr<FLOAT_T[]> norm_matrix;
    const size_t knum, pnum;
    const uint window;
    Genotyper(string inputfile, string outputfile, size_t k, size_t p, unique_ptr<uint16[]> &c, const int w):
    knum(k),
    pnum(p),
    window(w),
    kmer_counts(c),
    norm_vec(std::make_unique<FLOAT_T[]>(DefaultSize)),
    norm_matrix(std::make_unique<FLOAT_T[]>(DefaultSize*DefaultSize)),
    coefs(std::make_unique<FLOAT_T[]>(MAX_UINT16)),
    residuels(std::make_unique<FLOAT_T[]>(MAX_UINT16)),
    reminders(std::make_unique<FLOAT_T[]>(MAX_UINT16)),
    results(std::make_unique<int[]>(MAX_UINT16))
    {};
    
    
    void runOneGroup(const PriorChunk* priorData, const std::string& inputfile, const std::string& outputfile, const float depth, std::mutex& Threads_lock)
    {
        
        auto gnum = priorData->genenum;
        
        
        //cerr << "generating kmer matrix for sample: " << inputfile<<endl;
        matrix.getNorm(&kmer_counts.get()[priorData ->kmervec_start], priorData->kmer_matrix, depth, priorData->genenum, priorData->kmervec_size,
                       norm_vec.get(), norm_matrix.get(), total_lambda);

        //cerr << "regressing to references for sample: " << inputfile<<endl;
        regresser.Call(&kmer_counts.get()[priorData ->kmervec_start], priorData->kmer_matrix, depth, priorData->genenum, priorData->kmervec_size, norm_vec.get(), norm_matrix.get(),  total_lambda, priorData->gene_kmercounts, coefs.get(), residuels.get(), priorData->numgroups, priorData->genegroups,priorData->groupkmernums);
        
        
        tree.Run(priorData->phylo_tree, coefs.get(), gnum, &results.get()[0], &reminders.get()[0], residuels.get(), norm_matrix.get());

        //cerr << "determine window residuels: " << inputfile<<endl;

        vector<double> likelihoods(gnum,0.0);
        KmerWindow kmerwindow(window);
        kmerwindow.KmerMatch(&kmer_counts.get()[priorData ->kmervec_start], priorData->kmer_matrix, depth, priorData->genenum, priorData->kmervec_size, &results.get()[0], likelihoods);
        
        write(priorData, outputfile, inputfile, priorData->prefix, priorData->genenames, likelihoods, kmerwindow.windowcovers, depth, Threads_lock);

    };
    
    void write(const PriorChunk* priorData, const std::string& outputfile, const string &sample, const string &prefix, const vector<string>&genenames_ori, const vector<double>& likelihoods, const vector<vector<tuple<int,int,int>>>& windowcovers, const float depth, std::mutex& Threads_lock)
    {
        
        std::unique_lock<std::mutex> lck(Threads_lock);
        
        auto gnum = priorData->genenum;
        
        auto genenames(genenames_ori);
        for (auto &genename: genenames)
        {
            genename = genename.substr(0,genename.find('\t', 0));
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
        
        fprintf(fwrite,">%s\t%s\n", prefix.c_str(), sample.c_str());
        
        /*
        fprintf(fwrite,"regress: ");
        const float cutoff = 0.5 / (gnum + 1);
        for (int i = 0; i < gnum; ++i)
        {
            if (coefs.get()[i] > cutoff) fprintf(fwrite,"%s:%.2lf,", genenames[i].c_str(),coefs.get()[i]);
        }
        fprintf(fwrite,"\n");
        */
        
        fprintf(fwrite,"regress: ");
        for (int i = 0; i < gnum; ++i)
        {
            float result = ((float) results.get()[i]) + reminders.get()[i] ;
            
            if (results.get()[i] > 0)
            {
                if (priorData->gene_kmercounts[i] >= 1000) fprintf(fwrite,"%s:%.2f,", genenames[i].c_str(), result);
            }
        }
        fprintf(fwrite,"\n");
        
        fprintf(fwrite,"score: ");
        for (int i = 0; i < gnum; ++i)
        {
            float result = ((float) results.get()[i]) + reminders.get()[i] ;
            
            if (results.get()[i] > 0)
            {
                if (priorData->gene_kmercounts[i] >= 1000) fprintf(fwrite,"%s:%.2f,", genenames[i].c_str(), 100*likelihoods[i]);
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
        
        if (outputfile != "stdout") fclose(fwrite);
        
        return ;
    };

    void newgroup(const PriorChunk* priorData)
    {
        auto gnum = priorData->genenum;
        
        size_t newalloc_size = MAX( DefaultSize, gnum + 10 );
	
        if (MAX(newalloc_size, alloc_size) > DefaultSize)
        {

            norm_vec = std::make_unique<FLOAT_T[]>(newalloc_size);
            try_allocate_unique(norm_matrix, newalloc_size*newalloc_size, newalloc_size*newalloc_size);
            coefs = std::make_unique<FLOAT_T[]>(newalloc_size);
            residuels = std::make_unique<FLOAT_T[]>(newalloc_size);
            results = std::make_unique<int[]>(newalloc_size);
                    
            alloc_size = newalloc_size;
        }
        
        //memcpy(norm_matrix.get(), priorData->prior_norm, sizeof (FLOAT_T) *  gnum * gnum);
        
        memset(norm_matrix.get(), 0, sizeof (FLOAT_T) *  gnum * gnum);
        
        memset(norm_vec.get(), 0, sizeof (FLOAT_T) * gnum );
        
        memset(coefs.get(), 0, sizeof(FLOAT_T) * gnum);
        
        memset(residuels.get(), 0, sizeof(FLOAT_T) * gnum);
        
        memset(results.get(), 0, sizeof(int) * gnum);
        
        //kmerwindow.resize(priorData->pathsizes);
                        
        total_lambda = 0;
        total_obs = 0;
        total_exp = 0;
                
    };
    
    void run(PriorChunk* priorData, const std::string& inputfile, const std::string& outputfile, float depth, std::mutex& Threads_lock)
    {
        auto begin = std::chrono::high_resolution_clock::now();
        
        //std::cerr << "Running sample: " << inputfile << " for gene "<< priorData->prefix << endl ;
        
        auto gnum = priorData->genenum;
        
        newgroup(priorData);
                    
        runOneGroup (priorData, inputfile, outputfile, depth, Threads_lock);

        auto end = std::chrono::high_resolution_clock::now();

	    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
            
        std::cerr<<"finished sample:" << inputfile << " for gene:"<< priorData->prefix << " for :" << elapsed.count()* 1e-9 << "s"<<endl;

    };

private:
    
    FLOAT_T total_lambda =0;
    ull total_exp = 0, total_obs = 0;
    
    unique_ptr<uint16[]>&  kmer_counts;
    KmerMatrix matrix;
    Regression regresser;
    TreeRound tree;
        
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
              std::string &bedfile,
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
    Bedfile(bedfile),
    priordata_manager(mfile, N*n),
    regions(r),
    window(w),
    nthreads(n),
    Nsubthreads(N)
    {};
    
    
    void Run();
    void Load();
    void LoadBed();
    void MergeRegions();
    void Onethread();
    void Count(string &inputfile,string &outputfile, float &depth,unique_ptr<uint16[]>& kmer_counts);
    void Genotype(string &inputfile,string &outputfile, float depth,unique_ptr<uint16[]>& kmer_counts,ull_atom &numfinished,std::vector<bool>& finished_group, std::mutex&  Threads_lock2);
    
    const std::vector<std::string>& inputfiles;
    const std::vector<std::string>& outputfiles;
    const std::string &matrixfile;
    const std::string &Bedfile;
    const std::string &backgroundfile;
    std::vector<char *> &regions;
    int ifhla=0, ifunmap = 0;
    const std::vector<float> &depths;
    std::unordered_set<std::string> &genes;
    const int window;
    const int nthreads;
    const int Nsubthreads;
    
private:
    uint totalkmers, totalgroups;
    std::atomic_uint restfileindex = {0};
    std::mutex Threads_lock;
    std::mutex Threads_lock2;
    
    KmerCounter<ksize> *Counter =NULL;
    PriorData priordata_manager;
};


template <int ksize>
void Processor<ksize>::Run()
{
    if (genes.size() > 0)
    {
        if (ifgeneregion)
        {
            totalgroups = priordata_manager.LoadIndex(genes, regions);
            
            if (!regions.size())
            {
                std::cerr << "ERROR: Missing gene in the record " << endl;
                std::_Exit(EXIT_FAILURE);
            }
            
            for (auto &gene: genes)
            {
                if ( strncmp(gene.c_str(), "HLA-", 4) || gene == "HLA") //flag of loading all HLA decoys
                {
                    regions.push_back("HLA");
                    break;
                }
            }
        }
        else
        {
            totalgroups = priordata_manager.LoadIndex(genes);
        }
        
    }
    else
    {
        totalgroups = priordata_manager.LoadIndex();
    }
   
    if (!Bedfile.empty())
    {
        LoadBed();
    }
    
    Counter = new KmerCounter<ksize>(2*priordata_manager.totalkmers);
    
    if (regions.size())
    {
        std::cerr<<"Note: Target run has been used, making sure you are using correct gene/matrix/prefix names and your BEDfile match with your reference MD5 value.\n"<<endl;
        
        auto size_before = regions.size();
        regions.erase(
            std::remove_if(regions.begin(), regions.end(),
                [](char* region) { return std::string(region) == "HLA"; }),
            regions.end()
        );
        
        if (regions.size() < size_before) ifhla = 1;
        
        size_before = regions.size();
        regions.erase(
            std::remove_if(regions.begin(), regions.end(),
                [](char* region) { return strncmp(region, "Unmap", 5) == 0; }),
            regions.end()
        );
        if (regions.size() < size_before) ifunmap = 1;
        
        MergeRegions();
        
        if (forceunmap) ifunmap = forceunmap;   //forceunamp = 0,1,2, 0: unset, 1: useunmap, 2:notuse
    }
    else
    {
        std::cerr<<"Note: Global run has been used, all NGS reads will be used and all genes in the database will be genotyped. \n"<<endl;
    }
    
    Counter->LoadRegion(regions, ifhla, ifunmap);
    
    if (backgroundfile.length() > 0)
    {
        std::cerr<<"Note: NGS coverage not provided, backgrouds kmers will be used to determine NGS coverage, using file: " <<backgroundfile << endl <<endl;
        Counter->load_backgrounds(backgroundfile.c_str());
    }
    else
    {
        std::cerr<<"Note: NGS coverage provided. MAKING SURE THIS: 31-mer depth = (1 - 30/read_length) × sequencing_depth = 0.8 × sequencing_depth for 150bps NGS." << endl <<endl;
    }
    
    std::cerr<<"reading all kmer targets"<<endl;

    if (genes.size())
    {
        
        totalkmers = Counter->read_target(matrixfile.c_str(), priordata_manager.file_pos);
    }
    else
    {
        totalkmers = Counter->read_target(matrixfile.c_str());
    }
    
    
    std::cerr<<"finishing reading targets, start genotyping"<<endl;

    std::vector<std::unique_ptr<std::thread>> threads;
    
    for(int i=0; i< nthreads; ++i)
    {
        threads.push_back(std::unique_ptr<std::thread>(new std::thread(&Processor<ksize>::Onethread, this)));
    }
    
    
    for(int i=0; i< nthreads; ++i)
    {
        threads[i].get()->join();
    }
    
    if (!profilesummaryfile.empty())
    {
        BedReader breader;
        breader.MergeFiles(outputfiles, profilesummaryfile);
    }
    
}


template <int ksize>
void Processor<ksize>::Onethread()
{
    
    while (restfileindex < inputfiles.size() )
    {
        Threads_lock.lock();
        
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
        
        string inputfile = inputfiles[inputindex];
        auto begin = std::chrono::high_resolution_clock::now();
        
        if (ifprofile == 0)
        {
            unique_ptr<uint16[]> kmer_counts(new uint16[totalkmers+1]);
            std::fill_n(kmer_counts.get(), totalkmers + 1, 0);
            
            Count(inputfile,outputfile, depth, kmer_counts);
            
            std::vector<std::thread> subthreads;
            
            std::vector<bool> finished_group(totalgroups, 0);
            ull_atom numfinished = 0;
            
            ull total = 0;
            for (int i = 1; i < totalkmers+1; ++i)
            {
                total += kmer_counts[i] ;
            }
            
            subthreads.emplace_back([&]() {
                Genotype(inputfile, outputfile, depth, kmer_counts, numfinished, finished_group, Threads_lock2);
            });

            
            
            for (auto& t : subthreads)
            {
                t.join();
            }
        }
        
        else
        {
            unique_ptr<uint16_t[]> kmer_counts(new uint16_t[1]()); 
            
            KmerProfile<ksize> Profiler(Counter->kmer_hash,Counter->kmer_multi_hash,priordata_manager.prefixes,outputfile.c_str());
            
            
            Profiler.group_names = priordata_manager.prefixes;
            Profiler.kmer_ranges = priordata_manager.kmer_ranges;
            ull_atom totalbgs =0, totalbases = 0, totalreads = 0;
            Profiler.LoadRegion(regions, ifhla, ifunmap);
            Profiler.Call(inputfile.c_str(), kmer_counts.get(), totalbases, totalreads, totalbgs, Nsubthreads);
            return;
        }
        
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin);
        std::cerr<<"Fnished sample:" << inputfile << " for :" << elapsed.count()* 1e-9 <<"s"<<endl;
    }
}

template <int ksize>
void Processor<ksize>::Genotype(string &inputfile,string &outputfile, float depth,unique_ptr<uint16[]>& kmer_counts, ull_atom &numfinished, std::vector<bool>& finished_group, std::mutex&  Threads_lock2)
{
    
    Genotyper<ksize> genotyper( inputfile, outputfile, totalkmers, totalgroups, kmer_counts, window);
        
    for(int i=0; i < finished_group.size(); ++i)
    {
        Threads_lock2.lock();
        ull numfinished_ = numfinished++;

        if (numfinished_ >= finished_group.size())
        {
            Threads_lock2.unlock();
            break;
        }

        PriorChunk* priorData = priordata_manager.getNextChunk(finished_group);
        if (priorData)
        {
            finished_group[priorData->index] = 1;
        }
        Threads_lock2.unlock();

        if (!priorData) break;

        genotyper.run(priorData, inputfile, outputfile, depth, Threads_lock2);
        priordata_manager.FinishChunk(priorData);
    }
    
}


template <int ksize>
void Processor<ksize>::Count(string &inputfile,string &outputfile, float &depth,unique_ptr<uint16[]>& kmer_counts)
{
    std::cerr << "counting kmers for sample: " << inputfile<<endl;

    auto begin = std::chrono::high_resolution_clock::now();
    
    ull_atom totalbases = 0, totalreads = 0;
    std::vector<uint16> totalbgs (background_prime, 0);

    
    Counter->Call(inputfile.c_str(), kmer_counts.get(), totalbases, totalreads, totalbgs, Nsubthreads);

    auto end = std::chrono::high_resolution_clock::now();
    
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        
    std::cerr<<"finished counting "<< inputfile <<" at time: "<<elapsed.count()* 1e-9 <<endl;
    
    size_t uniqbgs = std::count_if(totalbgs.begin(), totalbgs.end(), [](uint16 x) { return x > 2; });
    size_t sumbgs = std::accumulate(totalbgs.begin(), totalbgs.end(), 0);
    
    if (depth <= 0)
    {
        depth = ( 0.5 * sumbgs )/uniqbgs;
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
    
            
    fprintf(fwrite,"@totalreads: %llu, totalbackgrounds: %llu/%llu \n", (ull)totalreads, sumbgs, uniqbgs);
    
    if (outputfile != "stdout") fclose(fwrite);
}

template <int ksize>
void Processor<ksize>::LoadBed()
{
    BedReader breader;
    
    breader.ReaderFile(Bedfile, genes, regions);
}

template <int ksize>
void Processor<ksize>::MergeRegions()
{
    BedReader breader;
    
    breader.MergeRegions(regions);
}


#endif /* Processor_hpp */

