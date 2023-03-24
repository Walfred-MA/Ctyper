//
//  PriorData.hpp
//  CTyper
//
//  Created by Wangfei MA on 3/8/23.
//

#ifndef PriorData_hpp
#define PriorData_hpp

#include <stdio.h>
#include <vector>
#include <unordered_set>

#include "config.hpp"
#include "KmerCounter.hpp"
#include "Regression.hpp"
#include "TreeRound.hpp"
#include "KtableReader.hpp"

#define buffer_size 100

using namespace std;

struct PriorChunk
{
    
    ~PriorChunk()
    {
        free(gene_kmercounts);
        free(kmer_matrix);
        free(prior_norm);
        free(phylo_tree);
    }
    
    string prefix = "";
    size_t genenum = 0;
    
    vector<string> genenames;

    FLOAT_T* prior_norm= NULL;
    size_t prior_norm_allocsize = 0;
    
    uint16* kmer_matrix = NULL;
    size_t kmer_matrix_allocsize = 0;
    
    node* phylo_tree = NULL;
    size_t phylo_tree_allocsize = 0;
    
    uint* gene_kmercounts = NULL;
    size_t gene_kmercounts_allocsize = 0;
    
    size_t kmervec_start= 0, kmervec_size = 0;

};



class PriorData
{
    
public:
    
    PriorData(const string &path): datapath(path), file(path.c_str())
    {};
    ~PriorData()
    {}
    
    size_t LoadIndex(const unordered_set<string>& geneset);
    size_t LoadIndex();
    
    void LoadFile();
    void CloseFile();

    PriorChunk* getChunkData(const size_t Chunkindex);
    PriorChunk* getNextChunk(const vector<bool>& finished);
    

private:
    
    void LoadHeader(PriorChunk &Chunk);
    void LoadAlleles(PriorChunk &Chunk);
    void LoadCounts(PriorChunk &Chunk);
    void LoadMatrix(PriorChunk &Chunk, size_t index_size);
    void LoadNorm(PriorChunk &Chunk);
    void LoadTree(PriorChunk &Chunk);
    
    
    void FinishChunk(PriorChunk* Chunk_prt);
    
    vector<string> prefixes;
    vector<pair<size_t,size_t>> kmervec_pos;
    vector<pair<size_t,size_t>> file_pos;
    vector<size_t> indexed_matrix_sizes;
    
    vector<PriorChunk> Buffers = vector<PriorChunk>( buffer_size );
    vector<size_t> Buffer_indexes = vector<size_t>( buffer_size );
    vector<size_t> Buffer_working_counts = vector<size_t>( buffer_size ) ;
    
    size_t buff_index = 0, total_buff;
    const string datapath;
    
    size_t LoadRow(uint16* matrix, size_t index, string &StrLine);
    KtableReader file;
    std::mutex IO_lock;
    
};


#endif /* PriorData_hpp */
