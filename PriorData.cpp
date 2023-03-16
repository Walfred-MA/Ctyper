//
//  PriorData.cpp
//  CTyper
//
//  Created by Wangfei MA on 3/8/23.
//

#include "PriorData.hpp"



void strsplit(string& str, vector<string>& eles, char deli)
{
    string::size_type start = 0 , end = 0 ;
    size_t len = strlen(str.c_str());
    while (start < len)
    {
        end = str.find(deli, start);
        if (end == std::string::npos) end = len ;
        eles.push_back(str.substr(start, end - start));
        start = end + 1;
    }
}


size_t PriorData::LoadIndex(const unordered_set<string>& geneset)
{
    std::ifstream pathfile(datapath + ".index");

    if(!pathfile)
    {
        std::cout<<"Error opening index file"<<std::endl;
        return 0;
    }
    std::string line;
    
    while (std::getline(pathfile, line))
    {
        string::size_type pos = line.find('\t');
        string genename = line.substr(0, pos);
        if (geneset.size() == 0 or geneset.find(genename) != geneset.end())
        {
            vector<string> eles;
            
            strsplit(line, eles, '\t');
    
            file_pos.push_back(make_pair(stoi(eles[1]), stoi(eles[2])));
        }
    }
    
    return file_pos.size();
}

size_t PriorData::LoadIndex()
{
    
    cout<<datapath<<endl;
    std::ifstream pathfile(datapath + ".index");

    if(!pathfile)
    {
        std::cout<<"Error opening index file"<<std::endl;
        return 0;
    }
    std::string line;
    
    size_t total_kmers = 0 ;
    size_t readsize = 0;
    while ( std::getline(pathfile, line) )
    {
        string::size_type pos = line.find('\t');
        
        string genename = line.substr(0, pos);
        
        vector<string> eles;
        
        strsplit(line, eles, '\t');
        
        prefixes.push_back(eles[0]);

        file_pos.push_back(make_pair(stoi(eles[1]), stoi(eles[2])));
        
        kmervec_pos.push_back(make_pair(total_kmers , total_kmers  + stoi(eles[3])));
        
        indexed_matrix_sizes.push_back(stoi(eles[4]));
        
        total_kmers += stoi(eles[3]);
        
    }
    
    return file_pos.size();
}

void PriorData::LoadHeader(PriorChunk &Chunk)
{
    string StrLine;
    StrLine.resize(MAX_LINE);
    
    if (!file.nextLine(StrLine))
    {
        std::cerr << "ERROR: error in kmer matrix file "<<std::endl;
        std::_Exit(EXIT_FAILURE);
        return;
    }
    
    const size_t len = strlen(StrLine.c_str());
    
    tuple<string,size_t,uint16> curr_info;
        
    int startpos = 0;
    char c;
    for (; startpos < len ; ++startpos)
    {
        if (StrLine[startpos] == '\t') break;
    }
    
    cout<<"header:"<<StrLine<<endl;
    
    Chunk.prefix = StrLine.substr(0, startpos);
    
    auto &curr_genenum = Chunk.genenum;
    auto &curr_kmernum = Chunk.kmernum;
    
    curr_kmernum = 0;
    for (startpos = startpos + 1; startpos < len ; ++startpos)
    {
        c = StrLine[startpos];
        if (c == '\t' || c=='\0' || c =='\n') break;
        
        curr_kmernum *= 10;
        curr_kmernum += c - '0';
    }
    
    curr_genenum = 0;
    for (startpos = startpos + 1; startpos < len ; ++startpos)
    {
        c = StrLine[startpos];
        if (c == '\t') break;
        
        curr_genenum *= 10;
        curr_genenum += c - '0';
    }
    
}

void PriorData::LoadNorm(PriorChunk &Chunk)
{
    string StrLine;
    StrLine.resize(MAX_LINE);
    
    if (!file.nextLine_norm(StrLine))
    {
        std::cerr << "ERROR: error in kmer matrix file "<<std::endl;
        std::_Exit(EXIT_FAILURE);
        return;
    }
    
    const size_t len = strlen(StrLine.c_str());
    
    size_t& curr_genenum = Chunk.genenum;
    double*& prior_norm = Chunk.prior_norm;
    size_t& prior_norm_size = Chunk.prior_norm_size;
    
    if (curr_genenum *curr_genenum > prior_norm_size || curr_genenum *curr_genenum  < prior_norm_size)
    {
        prior_norm = (double *) realloc(prior_norm, sizeof(double) * curr_genenum *curr_genenum  );
        prior_norm_size = curr_genenum * curr_genenum  ;
    }
    
    double element = 0.0;
    double decimal = 1.0;
    bool ifdecimal = 0;
    uint16 rowindex = 0, colindex = 0;
    
    char c;
    for (int startpos = 1; startpos < len; ++startpos)
    {
        c = StrLine[startpos];
        switch (c)
        {
            case '\t': case '\n':
                if (++colindex >= curr_genenum)
                {
                    rowindex++;
                    colindex -= curr_genenum - rowindex;
                }
                prior_norm[curr_genenum * rowindex + colindex ] = element;
                ifdecimal = 0;
                element = 0.0;
                decimal = 1.0;
                break;
            case '.':
                ifdecimal = 1;
                break;
            default:
                if (ifdecimal)
                {
                    decimal *= 0.1;
                    element += (c - '0') * decimal;
                }
                else
                {
                    element *= 10;
                    element += c - '0';
                }
        }
    }
    
    cout<<"checknorm"<<endl;
}

size_t PriorData::LoadRow(uint16* matrix, size_t rindex)
{
    string StrLine;
    StrLine.resize(MAX_LINE);
    
    if (!file.nextLine(StrLine))
    {
        std::cerr << "ERROR: error in kmer matrix file "<<std::endl;
        std::_Exit(EXIT_FAILURE);
        return NULL;
    }
    
    const size_t len = strlen(StrLine.c_str());
    //size_t count = std::count_if( StrLine.begin(), StrLine.end(), []( char c ){return c ==',';}) + 3;
    uint16 &sign = matrix[1];
    if (StrLine[1] == '-')
    {
        sign = 0;
    }
    else
    {
        sign = 1;
    }
    
    uint16 rownum = 2;
    uint16 element = 0;
    char c;
    
    size_t startpos = 2;
    for (; startpos < len ; ++startpos)
    {
        if (StrLine[startpos] == '\t') break;
    }
    
    for (; startpos < len; ++startpos)
    {
        c = StrLine[startpos];
        switch (c)
        {
            case ',': case '\n':
                matrix[rownum ++] = element;
                element = 0;
                break;
            default:
                element <<= 6;
                element += c - '0';
        }
    }
    
    if (element) matrix[rownum ++] = element;
    
    matrix[0] = rownum - 2;
    
    return rownum;
    
}

void PriorData::LoadMatrix(PriorChunk &Chunk)
{
    string StrLine;
    StrLine.resize(MAX_LINE);
    
    file.nextLine(StrLine);
    
    uint16*& kmer_matrix = Chunk.kmer_matrix;
    const size_t kmernum = Chunk.kmernum;
    size_t& kmer_matrix_allocsize = Chunk.kmer_matrix_allocsize;
    const size_t& indexsize = Chunk.kmer_matrix_indexsize;
    
    if (indexsize > kmer_matrix_allocsize || indexsize < kmer_matrix_allocsize)
    {
        kmer_matrix = (uint16 *) realloc(kmer_matrix, sizeof(uint16) * indexsize );
                
        kmer_matrix_allocsize = indexsize ;
    }
    
    uint16* matrix = kmer_matrix;
    size_t total = 0;
    for (size_t rindex =0; rindex < kmernum; ++ rindex)
    {
        size_t rsize = LoadRow(matrix , rindex);
        matrix = &matrix[rsize];
        total += rsize;
    }
    
    cout<<"checkmatrix"<<endl;
}

void PriorData::LoadTree(PriorChunk &Chunk)
{
    string StrLine;
    StrLine.resize(MAX_LINE);
    
    if (!file.nextLine(StrLine))
    {
        std::cerr << "ERROR: error in kmer matrix file "<<std::endl;
        std::_Exit(EXIT_FAILURE);
        return;
    }
        
    size_t len = strlen(StrLine.c_str());
    
    if (len==0) return ;
    
    size_t count = std::count_if( StrLine.begin(), StrLine.begin()+len, []( char c ){return c ==':';}) + 1;
    
    size_t& phylo_tree_size = Chunk.phylo_tree_size;
    node**& phylo_tree = Chunk.phylo_tree;
    
    
    if (count > phylo_tree_size || 2*count < phylo_tree_size)
    {
        phylo_tree = (node**) realloc(phylo_tree, sizeof(node*) * count); 
        
        for (size_t index =phylo_tree_size; index < count; ++index)
        {
            phylo_tree[index] = new node();
        }
    }
        
    for (size_t index =0; index < count; ++index)
    {
        phylo_tree[index]->clear();
    }
    
    node* current_node = phylo_tree[0];
    double current_num = 0.0;
    uint16 current_index = 1;
    float ifdeci = 0;
    
    cout<<"check19,"<<phylo_tree_size<<endl;
    cout<<"check20,"<<Chunk.genenum<<endl;
    
    int notuselast = (StrLine[len-2] == ';') ;
    
    char c;
    for (int pos=1; pos < len - notuselast; ++pos)
    {
        c = StrLine[pos];
        switch(c)
        {
            case ' ': case '\n':
                break;
            case '(':
                current_node = current_node->add(phylo_tree[current_index++]);
                break;
            case ')': case ';':
                current_node->dist = current_num;
                ifdeci = 0;
                current_node = current_node->parent;
                break;
            case ',':
                current_node->dist = current_num;
                ifdeci = 0;
                current_node = current_node->parent->add(phylo_tree[current_index++]);
                break;
            case ':':
                ifdeci = 1;
                break;
            case '.':
                ifdeci *= 0.1;
            default:
                if (ifdeci>0)
                {
                    if (ifdeci==1) current_num *= 10;
                    current_num += ifdeci * (c - '0');
                    if (ifdeci<1) ifdeci *= 0.1;
                }
                break;
        }
    }
    
    cout<<"check21"<<endl;
}


PriorChunk* PriorData::getChunkData(size_t Chunkindex)
{
    auto chunk_region = file_pos[Chunkindex];
    size_t chunk_start = chunk_region.first;
    
    file.Seek(chunk_start);
    cout<<"chunkstart:"<<chunk_start<<endl;
    
    size_t i = 0;
    for (; i < buffer_size; ++i)
    {
        if (Buffer_working_counts[i] == 0 ) break;
    }
    
    Buffer_working_counts[i]++;
    PriorChunk &Chunk = Buffers[i];
    
    auto &kmervec_range = kmervec_pos[Buffer_indexes[i]];
    
    Chunk.kmervec_start = kmervec_range.first;
    Chunk.kmer_matrix_indexsize = indexed_matrix_sizes[Buffer_indexes[i]] + 10;
        
    LoadHeader(Chunk);
    
    LoadTree(Chunk);
    cout<<"check13"<<endl;
    LoadNorm(Chunk);
    
    cout<<"check14"<<endl;
    
    LoadMatrix(Chunk);
    
    cout<<"check16"<<endl;
    return &Chunk;
    
}

void PriorData::FinishChunk(PriorChunk* Chunk_prt)
{
    lock_guard<mutex> IO(IO_lock);
    
    for (int i = 0; i < buffer_size; ++i)
    {
        if (&Buffers[i] == Chunk_prt )
        {
            Buffer_working_counts[i]--;
            break;
        }
    }
}

PriorChunk* PriorData::getNextChunk(const vector<bool>& finished)
{
    lock_guard<mutex> IO(IO_lock);
    
    cout<<"check9"<<endl;
    
    for (size_t i = 0 ; i < Buffer_indexes.size(); ++i)
    {
        auto buffer_index = Buffer_indexes[i];
        
        cout<<"check10,"<<buffer_index<<endl;
        
        if (buffer_index && finished[buffer_index -1] == 0)
        {
            return &Buffers[i];
        }
    }
    
    
    size_t i = 0;
    for (; i < finished.size(); ++i)
    {
        if (finished[i] == 0) break;
    }
    
    
    
    return getChunkData(i);
}
