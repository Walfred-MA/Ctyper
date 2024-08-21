//
//  PriorData.cpp
//  CTyper
//
//  Created by Wangfei MA on 3/8/23.
//

#include "PriorData.hpp"

using namespace std;

void strsplit(string& str, vector<uint16>& eles, char deli, string::size_type start = 0 )
{
    size_t end = str.find(",");
    size_t len = strlen(str.c_str());
    
    while (end != std::string::npos && end < len)
    {
        eles.push_back(std::stoi(str.substr(start, end - start)));
        start = end + 1;
        end = str.find(",", start);
    }
}

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
    
    totalkmers = 0 ;
    while (std::getline(pathfile, line))
    {
        string::size_type pos = line.find('\t');
        
        string genename = line.substr(0, pos);
        
        if (geneset.size() == 0 or geneset.find(genename) != geneset.end())
        {
            vector<string> eles;
            
            strsplit(line, eles, '\t');
            
            prefixes.push_back(eles[0]);

            file_pos.push_back(make_pair(stol(eles[1]), stol(eles[2])));
            
            kmervec_pos.push_back(make_pair(totalkmers , totalkmers  + stol(eles[3])));
            
            indexed_matrix_sizes.push_back(stol(eles[4]));
            
            totalkmers += stol(eles[3]);
        }
    }
    
    return file_pos.size();
}

size_t PriorData::LoadIndex()
{
    
    std::ifstream pathfile(datapath + ".index");

    if(!pathfile)
    {
        std::cout<<"Error opening index file"<<std::endl;
        return 0;
    }
    std::string line;
    
    totalkmers = 0 ;
    while ( std::getline(pathfile, line) )
    {
        string::size_type pos = line.find('\t');
        
        string genename = line.substr(0, pos);
        
        vector<string> eles;
        
        strsplit(line, eles, '\t');
        
        prefixes.push_back(eles[0]);

        file_pos.push_back(make_pair(stol(eles[1]), stol(eles[2])));
        
        kmervec_pos.push_back(make_pair(totalkmers , totalkmers  + stol(eles[3])));
        
        indexed_matrix_sizes.push_back(stol(eles[4]));
        
        totalkmers += stol(eles[3]);
        
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
        
    Chunk.prefix = StrLine.substr(0, startpos);
    
    auto &curr_genenum = Chunk.genenum;
    auto &curr_kmernum = Chunk.kmervec_size;
    
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

void PriorData::LoadSizes(PriorChunk &Chunk)
{
    string StrLine;
    StrLine.resize(MAX_LINE);
    
    if (!file.nextLine(StrLine))
    {
        std::cerr << "ERROR: error in kmer matrix file "<<std::endl;
        std::_Exit(EXIT_FAILURE);
        return;
    }
    
    uint* &gene_kmercounts = Chunk.gene_kmercounts;
    auto &curr_genenum = Chunk.genenum;
    auto &allocsize = Chunk.gene_kmercounts_allocsize;
    
    if (curr_genenum > allocsize || 2 * curr_genenum < allocsize)
    {
        gene_kmercounts = (uint *) realloc(gene_kmercounts, sizeof(uint) * curr_genenum  );
        
        allocsize = curr_genenum  ;
    }
    
    const size_t len = strlen(StrLine.c_str());
    
    int index =0;
    int ele = 0;
    char c;
    for (int startpos = 1; startpos < len; ++startpos)
    {
        c = StrLine[startpos];
        switch (c)
        {
            case ' ': case '\n':
                gene_kmercounts[index++] = ele;
                ele = 0;
                break;
            default:
                ele *= 10;
                ele += c - '0';
        }
    }
    
    if (ele) gene_kmercounts[index++] = ele;
        
}


void PriorData::LoadAlleles(PriorChunk &Chunk)
{
    
    const size_t& curr_genenum = Chunk.genenum;
    vector<string>& genenames = Chunk.genenames;
    vector<string>& pathnames = Chunk.pathnames;
    
    if (curr_genenum > genenames.size() || 2 * curr_genenum  < genenames.size())
    {
        genenames.resize(curr_genenum);
    }
    
    string StrLine;
    StrLine.resize(MAX_LINE);
    
    pathnames.clear();
    pathnames.push_back("");
    while (file.nextLine_start(StrLine, '+'))
    {
        string line = StrLine.substr(1, StrLine.find('\n') - 1);
        size_t second_underscore = line.find('_', line.find('_') + 1);
        if (second_underscore != std::string::npos) line = line.substr(second_underscore + 1);
        
        size_t first_tab = line.find('\t');
        size_t second_tab = line.find('\t', first_tab + 1);

        if (first_tab != std::string::npos && second_tab != std::string::npos)
        {
            for (size_t i = first_tab + 1; i < second_tab; ++i) 
            {
                if (line[i] == '-') line[i] = '_';
            }
            line[first_tab] = '_';
        }

        //pathnames.push_back(StrLine.substr(StrLine.find_last_of('\t') + 1, StrLine.find('\n') - StrLine.find_last_of('\t') - 1));
        pathnames.push_back(line);
    }
    
    int index = 0;
    for (int i =0 ; i<  curr_genenum ; ++i)
    {
        if (!file.nextLine_start(StrLine, '>'))
        {
            std::cerr << "ERROR: error in kmer matrix file "<<std::endl;
            std::_Exit(EXIT_FAILURE);
            return;
        }
        genenames[index++] =  StrLine.substr(1,MIN( StrLine.find('\n', 0)-1  , StrLine.find(';', 0) -1) );
    }

}


void PriorData::LoadNorm(PriorChunk &Chunk)
{
    /*
    size_t& curr_genenum = Chunk.genenum;
    string StrLine;
    for ( int i = 0 ; i < curr_genenum ; ++i)
    {
        
        if (!file.nextLine_norm(StrLine))
        {
            std::cerr << "ERROR: error in kmer matrix file "<<std::endl;
            std::_Exit(EXIT_FAILURE);
            return;
        }
    }
    
    return;
    */
    
    size_t& curr_genenum = Chunk.genenum;
    
    
    FLOAT_T*& prior_norm = Chunk.prior_norm;
    size_t& prior_norm_allocsize = Chunk.prior_norm_allocsize;
     
    
    if (curr_genenum *curr_genenum != prior_norm_allocsize)
    {

	assert(curr_genenum < MAX_UINT16);
	    
        try_allocate(prior_norm, curr_genenum *curr_genenum, curr_genenum *curr_genenum);
        
        prior_norm_allocsize = curr_genenum * curr_genenum  ;
    }
    
    
    float element = 0.0;
    float decimal = 1.0;
    bool ifdecimal = 0;
    uint16 rowindex = 0, colindex = 0;
    
    string StrLine;
    
    
    for ( int i = 0 ; i < curr_genenum ; ++i)
    {
        
        if (!file.nextLine_norm(StrLine))
        {
            std::cerr << "ERROR: error in kmer matrix file "<<std::endl;
            std::_Exit(EXIT_FAILURE);
            return;
        }
        
        size_t len = strlen(StrLine.c_str());
                
        int norm_c = 0;
        char c;
        for (int startpos = 1; startpos < len; ++startpos)
        {
            c = StrLine[startpos];
            switch (c)
            {
                case '\t': case ' ': case '\n':
                    prior_norm[curr_genenum * rowindex + colindex ] = element;
                    prior_norm[curr_genenum * colindex + rowindex ] = element;
                    ifdecimal = 0;
                    element = 0.0;
                    decimal = 1.0;
                    norm_c++;
                    
                    if (++colindex >= curr_genenum )
                    {
                        rowindex++;
                        colindex = rowindex;
                    }
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
    }
    
    
    
    node*& phylo_tree = Chunk.phylo_tree;
    
    int leaveindex = 0;
    for (size_t index =0; index < Chunk.nodenum ; ++index)
    {
        if (phylo_tree[index].numchildren == 0)
        {
            int leaveindex_ = leaveindex++;
            phylo_tree[index].size = prior_norm[curr_genenum*leaveindex_+leaveindex_];
        }
    }
    
}


void PriorData::LoadGroups(PriorChunk &Chunk)
{
   
    size_t& curr_genenum = Chunk.genenum;
    string StrLine;
   
    Chunk.genegroups.resize(curr_genenum, 0);
    std::fill(Chunk.genegroups.begin(), Chunk.genegroups.end(), 0);

    Chunk.numgroups = 0;
    StrLine.resize(MAX_LINE);
    for ( int i = 0 ; i < curr_genenum ; ++i)
    {
        if (!file.nextLine_genegroup(StrLine) )
        {
            break;
        }
        
        Chunk.numgroups ++ ;
        
        vector<uint16> eles;
        
        strsplit(StrLine, eles, ',', 1);
        
        for (uint16 ele: eles)
        {
            Chunk.genegroups[ele] = i ;
        }
    }
    
    Chunk.groupkmernums.resize(Chunk.numgroups, 0);;
    std::fill(Chunk.groupkmernums.begin(), Chunk.groupkmernums.end(), 0);    


}
size_t PriorData::LoadRow(uint16* matrix, size_t rindex, string &StrLine, vector<uint> &pathsizes)
{
        
    if (!file.nextLine(StrLine))
    {
        std::cerr << "ERROR: error in kmer matrix file "<<std::endl;
        std::_Exit(EXIT_FAILURE);
        return NULL;
    }
        
    const size_t len = strlen(StrLine.c_str());
    //size_t count = std::count_if( StrLine.begin(), StrLine.end(), []( char c ){return c ==',';}) + 3;
    matrix[0] = StrLine[1];
    
    uint16 rownum = FIXCOL;
    uint16 element = 0;
    char c;
        
    size_t startpos = 3;
    
    uint16 tag1=0, tag2 = 0;
    for (; startpos < len ; ++startpos)
    {
        if (StrLine[startpos] == '|') break;
        tag1 <<= 6;
        tag1 += StrLine[startpos] - '0';
    }
    ++startpos;
    
    for (; startpos < len ; ++startpos)
    {
        if (StrLine[startpos] == '\t') break;
        tag2 <<= 6;
        tag2 += StrLine[startpos] - '0';
    }
    tag2 <<= 6;
    ++startpos;
    
    matrix[5] = tag1 + tag2;

    uint16 path=0;
    for (; startpos < len ; ++startpos)
    {
        if (StrLine[startpos] == '\t') break;
        path <<= 6;
        path += StrLine[startpos] - '0';
    }
    ++startpos;
    
    matrix[2] = path;
    
    if (pathsizes.size() <= path) pathsizes.resize(path + 1, 0);

    
    uint loc=0;
    for (; startpos < len ; ++startpos)
    {
        if (StrLine[startpos] == '\t') break;
        loc <<= 6;
        loc += StrLine[startpos] - '0';
    }
    ++startpos;
    
    
    pathsizes[path] = MAX(pathsizes[path], loc);
    
    uint16 part2 = loc & 0xFFFF;
    uint16 part1 = (loc >> 16) & 0xFFFF;
    
    matrix[3] = part1;
    matrix[4] = part2;

    for (; startpos < len ; ++startpos)
    {
        if (StrLine[startpos] == '\t') break;
    }
    ++startpos;
    
    for (; startpos < len; ++startpos)
    {
        
        c = StrLine[startpos];
        switch (c)
        {
            case ',':
                matrix[rownum ++] = element;
                element = 0;
                break;
            default:
                element <<= 6;
                element += c - '0';
        }
    }
    
    matrix[1] = rownum - FIXCOL;
    
    return rownum;
    
}

void GetGroupKmerNum(PriorChunk &Chunk, uint16* matrix, uint &ifingroup_counter, vector<bool> &ifingroup, size_t gnum)
{
    
    
    switch(matrix[0])
    {
        case '+':
        {
            for (int igroup = 0; igroup < Chunk.groupkmernums.size() ; ++igroup)
            {
                if (ifingroup[igroup]) Chunk.groupkmernums[igroup] += ifingroup_counter;
            }
            
            std::fill(ifingroup.begin(), ifingroup.end(), 0);
            
            for (int igroup = 0; igroup <matrix[1] ; ++igroup)
            {
                ifingroup[   Chunk.genegroups[ matrix[FIXCOL + igroup] ]  ] = 1;
            }
            
            ifingroup_counter = 1;
        
            break;
            
        }
            
        case '-':
        {
            for (int igroup = 0; igroup < Chunk.groupkmernums.size() ; ++igroup)
            {
                if (ifingroup[igroup]) Chunk.groupkmernums[igroup] += ifingroup_counter;
            }
            
            std::fill(ifingroup.begin(), ifingroup.end(), 0);
            
            vector<bool> row_reverse (gnum, 1);
            for (int igroup = 0; igroup <matrix[1]  ; ++igroup)
            {
                row_reverse[ matrix[FIXCOL + igroup] ] = 0;
            }
            
            for (int igroup = 0; igroup < gnum; ++igroup)
            {
                if (row_reverse[igroup]) ifingroup[   Chunk.genegroups[ igroup]  ] = 1;
            }
            
            ifingroup_counter = 1;
            break;
        }
            
        default:
        {
            ifingroup_counter ++;
            break;
        }
    }
}

void PriorData::LoadMatrix(PriorChunk &Chunk, size_t new_kmer_matrix_allocsize)
{
    string StrLine;
    StrLine.resize(MAX_LINE);
    
    uint16*& kmer_matrix = Chunk.kmer_matrix;
    const size_t kmernum = Chunk.kmervec_size;
    size_t& kmer_matrix_size = Chunk.kmer_matrix_allocsize;
    
    if (new_kmer_matrix_allocsize != kmer_matrix_size )
    {
        kmer_matrix = (uint16 *) realloc(kmer_matrix, sizeof(uint16) * new_kmer_matrix_allocsize + 10 );
        
        kmer_matrix_size = new_kmer_matrix_allocsize;
    }
    
    Chunk.pathsizes.resize(0);
    
    
    vector<bool> ifingroup (Chunk.numgroups,0);
    uint ifingroup_counter = 0;
        
    uint16* matrix = kmer_matrix;
    size_t rindex =0;
    for (rindex =0; rindex < kmernum; ++ rindex)
    {
        uint16 rsize = LoadRow(matrix , rindex, StrLine, Chunk.pathsizes);
        
        GetGroupKmerNum(Chunk, matrix, ifingroup_counter, ifingroup, Chunk.genenum);

        matrix = &matrix[rsize];
    }
    
    matrix[0] = '+';
    matrix[1] = 0;
    GetGroupKmerNum(Chunk, matrix,  ifingroup_counter, ifingroup, Chunk.genenum);
    
    /*
    size_t total2 = 0;
    for (size_t rindex =0; rindex < kmernum; ++ rindex)
    {
        total2 += kmer_matrix[total2 + 1] + 2;
    }
    */

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

	static float pow10[7] = {1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001};
    
    size_t count = Chunk.nodenum;
    count = std::count_if( StrLine.begin(), StrLine.begin()+len, []( char c ){return c ==':';}) + 1;
    
    size_t& phylo_tree_allocsize = Chunk.phylo_tree_allocsize;
    node*& phylo_tree = Chunk.phylo_tree;
    
    
    if (count > phylo_tree_allocsize || 2*count < phylo_tree_allocsize)
    {
        phylo_tree = (node*) realloc(phylo_tree, sizeof(node) * count);
        
        phylo_tree_allocsize = count;
        
    }
        
    for (size_t index =0; index < count; ++index)
    {
        phylo_tree[index].clear();
    }
    
    node *current_node = &phylo_tree[0];
    float current_num = 0.0;
    uint16 current_index = 1;
    float ifdeci = 0;
    
    int notuselast = (StrLine[len-2] == ';') ;
    
	float sign = 1;
    bool ifsci_e = 0;
    int sci_e = 0;
    char c;
    for (int pos=1; pos < len - notuselast; ++pos)
    {
        c = StrLine[pos];
        switch(c)
        {
            case ' ': case '\n':
                break;
            case '(':
                current_node = current_node->add(&phylo_tree[current_index++]);
                break;
            case ')': case ';':
                if (sci_e>5) current_num = 0;
                else if (sci_e > 0) current_num *= pow10[sci_e];
                current_node->dist = sign * current_num  ;
                ifdeci = 0;
                current_num = 0;
                sign = 1;
                sci_e = 0;
                ifsci_e = 0;
                current_node = current_node->parent;
                break;
            case ',':
                if (sci_e>5) current_num = 0;
                else if (sci_e > 0) current_num *= pow10[sci_e];
                current_node->dist = sign *  current_num ;
                ifdeci = 0;
                current_num = 0;
                sign = 1;
                sci_e = 0;
                ifsci_e = 0;
                current_node = current_node->parent->add(&phylo_tree[current_index++]);
                break;
            case ':':
                ifdeci = 0;
                ifsci_e = 0;
                sci_e = 0;
                break;
            case '.':
                ifdeci *= 0.1;
                break;
            case '-':
                if (ifdeci==1) sign = -1;
                break;
            case 'e':
                ifdeci = 0;
                ifsci_e = 1;
                sci_e = 0;
                break;
            default:
                if (ifdeci>0)
                {
                    if (ifdeci==1) current_num *= 10;
                    current_num += ifdeci * (c - '0');
                    if (ifdeci<1) ifdeci *= 0.1;
                }
                else if (ifsci_e)
                {
                    sci_e *=10;
                    sci_e += (c - '0');
                }
                
                break;
        }
    }
    
    const uint* gene_kmercounts = Chunk.gene_kmercounts;
 
    int leaveindex = 0;
    int nonleaveindex = -1;
    for (size_t index =0; index < count; ++index)
    {
        if (phylo_tree[index].numchildren == 0)
        {
            phylo_tree[index].size = gene_kmercounts[leaveindex];
            phylo_tree[index].index = leaveindex ++ ;
        

        }
        else
        {
            phylo_tree[index].index = nonleaveindex --;
        }
    }

}

PriorChunk* PriorData::getFreeBuffer(size_t Chunkindex)
{
    
    size_t i = 0;
	
    for (; i < buffer_size; ++i)
    {
        if (Buffer_indexes[i] == INT_MAX ) break;  //uninitialized block
    }
	
    for (i = 0; i < buffer_size; ++i)
    {
        if (Buffer_working_counts[i] == 0 ) break;  //not in using block
    }
    
    Buffer_working_counts[i]++;
    Buffer_indexes[i] = Chunkindex;
   
    return &Buffers[i];
}


PriorChunk* PriorData::getChunkData(size_t Chunkindex)
{
        
    PriorChunk &Chunk = *getFreeBuffer(Chunkindex);
    Chunk.index = Chunkindex;
    
    auto chunk_region = file_pos[Chunkindex];
    size_t chunk_start = chunk_region.first;
    file.Seek(chunk_start);
    
    auto &kmervec_range = kmervec_pos[Chunkindex];
    Chunk.kmervec_start = kmervec_range.first;
    Chunk.kmervec_size = kmervec_range.second;
        
    LoadHeader(Chunk);
    
    LoadSizes(Chunk);
 
    LoadTree(Chunk);
    
    LoadAlleles(Chunk);
    
    LoadNorm(Chunk);
    
    LoadGroups(Chunk);
    
    LoadMatrix(Chunk, indexed_matrix_sizes[Chunkindex] + 10);
   

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
    
    for (size_t i = 0 ; i < Buffer_indexes.size(); ++i)
    {
        auto buffer_index = Buffer_indexes[i];
                
        if (buffer_index != INT_MAX && finished[buffer_index] == 0)
        {
            Buffer_working_counts[i] ++;
            return &Buffers[i];
        }
    }
    
    size_t i = 0;
    for (; i < finished.size(); ++i)
    {
        if (finished[i] == 0) break;
    }
  

    if (i >= finished.size())
    {
        std::cerr << "Buffer Error" <<std::endl;
        std::_Exit(EXIT_FAILURE);
    }
    
    return getChunkData(i);
}
