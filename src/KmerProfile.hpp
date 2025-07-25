//
//  KmerProfile.hpp
//  Ctyper1
//
//  Created by walfred on 6/29/25.
//

#ifndef KmerProfile_hpp
#define KmerProfile_hpp

#include "KmerCounter.hpp"

constexpr int hotspot_freq_cutoff = 1000;
constexpr int hotspot_window_cutoff = 1000;
constexpr int extendsize = 150;

static inline uint findinsegments(uint index, const vector<size_t> &segments)
{
    auto it = std::upper_bound(segments.begin(), segments.end(), index);
    
    index = std::distance(segments.begin(), it);
    
    return index;
}

template <typename T1, typename T2, typename T3>
static void SearchGroup(T2 &kmer_hash, T3 &kmer_multi_hash, T1 &larger_kmer , vector<uint> &groupindex , uint &group_len, vector<size_t> &segments)
{
    uint hash_;
    auto map_find = kmer_hash.find(larger_kmer, hash_);
    
    group_len = 0;
    if (map_find != NULL)
    {
        uint index = *map_find;
        
        if (index < UINT_MAX )
        {
            auto gindex = findinsegments(index,segments);
            groupindex[0] = gindex;
            group_len = 1;
            
        }
        else
        {
            uint* const& data = kmer_multi_hash.find(larger_kmer)->second;
            group_len= data[0];
            
            for (uint j = 1; j <= group_len; ++j)
            {
                uint idx = data[j];
                auto gindex = findinsegments(idx, segments);
                groupindex[j -1] = gindex;
            }
            
        }
    }
    
    
}


template <int dictsize>
class KmerProfile
{
public:
    using kmer_int = typename std::conditional<(dictsize>32), u128, ull>::type;
    using kmer_hash_type = typename std::conditional<(dictsize>32), Kmer64_hash, Kmer32_hash>::type;
    using kmer_hash_type_mul = typename std::conditional<(dictsize>32), kmer64_dict_mul, kmer32_dict_mul>::type;

    KmerProfile(kmer_hash_type &hash1, kmer_hash_type_mul &hash2, vector<string>& groups, const char* outputfile)
        : kmer_hash(hash1), kmer_multi_hash(hash2), group_names(groups), outFile (outputfile)
    {
        if (!outFile.is_open())
        {
            std::cerr << "Error: Could not open file " << outputfile << " for writing." << std::endl;
            return;
        }
    }
    ~KmerProfile(){outFile.close();};
    
    void LoadRegion(std::vector<char *> &r, int hla, int unmap);
    
    template <class typefile>
    void count_kmer_(typefile &file, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg, const int nthreads);
    
    template <bool EnableLock>
    void count_kmer(FastaReader &file, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg);
   
    template <bool EnableLock>
    void count_kmer(FastqReader &file, uint16* samplevecs , ull_atom &nBases, ull_atom &nReads ,ull_atom &nBg);
    
    template <bool EnableLock>
    void count_kmer(CramReader &file, uint16* samplevecs , ull_atom &nBases, ull_atom &nReads ,ull_atom &nBg);
    
    void Call(const char* infile, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg, const int nthreads);
    
    template <bool EnableLock>
    void WriteHotspots(string &header,vector<vector<hotspot>>& allhotspots);
    
    void UpdateLocation(vector<vector<hotspot>> &allhotspots,vector<array<int, hotspot_freq_cutoff>> &kmer_lastposis, std::vector<int> &kmer_lastindex,  const vector<uint> &groups, const uint groups_len, const int posi);
    
    template <bool EnableLock>
    void PrintReaders(vector<uint> &groups, vector<string> &group_names,  string &Header);
    
    kmer_hash_type &kmer_hash;
    kmer_hash_type_mul &kmer_multi_hash;
    std::vector<char *>* regions=NULL;
    vector<size_t> kmer_ranges;
    vector<string>& group_names;
    unordered_map<string,uint16> group_index;
    int ifhla = 0;
    int ifunmap = 0;
    std::ofstream outFile;
    
    std::mutex lock;
};

template <int dictsize>
void KmerProfile<dictsize>::LoadRegion(std::vector<char *> &r, int hla,  int unmap)
{
    regions = &r;
    ifhla = hla;
    ifunmap = unmap;
}

template <int dictsize>
template <bool EnableLock>
void KmerProfile<dictsize>::PrintReaders(vector<uint> &groups, vector<string> &group_names, string &Header)
{
    if constexpr (EnableLock) std::lock_guard<std::mutex> guard(lock);
    
    for ( int i = 0 ; i < groups.size() ; ++i)
    {
        outFile << group_names[groups[i]] << ";";
    }
    outFile << "\t";
    
    outFile << Header << "\n";
    
}


template <int dictsize>
template <bool EnableLock>
void KmerProfile<dictsize>::WriteHotspots(string &header,vector<vector<hotspot>>& allhotspots)
{
    
    
    if constexpr (EnableLock)
    {
        std::lock_guard<std::mutex> guard(lock);
        // Loop through each group in group_names and allhotspots
        for (size_t i = 0; i < group_names.size(); ++i)
        {
            const std::string& group_name = group_names[i];
            
            std::vector<hotspot>& hotspots = allhotspots[i];  // List of hotspots
            
            // Loop through each hotspot
            for (const auto& h : hotspots)
            {
                int start = h.first;  // Start position of the hotspot
                int end = h.second;   // End position of the hotspot
                
                
                // Write the chromosome, start, end, and group name to the file
                outFile << header << '\t' << start << '\t' << end << '\t' << group_name << std::endl;
            }
            
            hotspots.clear();
        }
    }
    else
    {
        for (size_t i = 0; i < group_names.size(); ++i)
        {
            const std::string& group_name = group_names[i];
            
            std::vector<hotspot>& hotspots = allhotspots[i];  // List of hotspots
            
            // Loop through each hotspot
            for (const auto& h : hotspots)
            {
                int start = h.first;  // Start position of the hotspot
                int end = h.second;   // End position of the hotspot
                
                
                // Write the chromosome, start, end, and group name to the file
                outFile << header << '\t' << start << '\t' << end << '\t' << group_name << std::endl;
            }
            
            hotspots.clear();
        }
    }
    
}

template <int dictsize>
void KmerProfile<dictsize>::UpdateLocation(vector<vector<hotspot>> &allhotspots,vector<array<int, hotspot_freq_cutoff>> &kmer_lastposis, std::vector<int> &kmer_lastindex,  const vector<uint> &groups, const uint group_len,const int posi)
{
    for (uint index =0 ; index < group_len; ++index)
    {
        uint i = groups[index];
        vector<hotspot>& hotspots = allhotspots[i];
        auto& positions = kmer_lastposis[i];
        auto& index_start = kmer_lastindex[i];
        
        auto posi_start = positions[index_start];
        
        positions[index_start] = posi;
        index_start = (index_start + 1) % hotspot_freq_cutoff;
                    
        if (abs(posi - posi_start) < hotspot_window_cutoff)
        {
            if (hotspots.size() && posi_start - hotspots[hotspots.size()-1].second < hotspot_window_cutoff )
            {
                hotspots[hotspots.size()-1].second = posi + extendsize;
            }
            
            else
            {
                if (hotspots.size()%1000 == 0)
                {
                    hotspots.reserve(hotspots.size() + 1000);
                }
                
                hotspots.push_back(std::make_pair( MAX (0, posi_start -extendsize ), posi + extendsize));
                
            }
        }
    }
}


template <int dictsize>
template <bool EnableLock>
void KmerProfile<dictsize>::count_kmer(FastaReader &file, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg)
{
    std::size_t current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    vector<uint> groups(group_names.size()+1,0);
    uint groups_len = 0;

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
                Header = StrLine;
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
            
            SearchGroup(kmer_hash, kmer_multi_hash, larger_kmer, groups, groups_len, kmer_ranges);
            
            if (groups.size())
            {
                PrintReaders<EnableLock>(groups, group_names, Header);
                break;
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
void KmerProfile<dictsize>::count_kmer(FastqReader &file, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg)
{
    std::size_t current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    vector<uint> groups(group_names.size()+1,0);
    uint groups_len = 0;
    uint hash_ = 0;
    size_t rstart =0 , rlen = 0;
    char *StrLine = NULL;
    char *Header = NULL;
    while (file.nextLine(StrLine,rlen,Header))
    {
        current_size = 0;
        for (int pos = 0; pos < rlen; ++pos)
        {
            char base = StrLine[pos];
            
            //if (base=='\n') break;
            
            kmer_read_31(base, ++current_size, current_kmer, reverse_kmer);
            
            if (current_size < klen) continue;
            
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
            
            SearchGroup(kmer_hash, kmer_multi_hash, larger_kmer, groups, groups_len, kmer_ranges);
            
            if (groups.size())
            {
                string header_str(Header, strlen(Header));
                PrintReaders<EnableLock>(groups, group_names, header_str);
                break;
            }

        }
        
    nReads+=1;
    if (nReads % 10000000 == 0) {
      cerr << "processed " << nReads / 1000000 << "M reads." << endl;
    }
    }
        
    return ;
};

template <int dictsize>
template <bool EnableLock>
void KmerProfile<dictsize>::count_kmer(CramReader &file, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg)
{
    auto num_targets = MAX ( 1, group_names.size() );
    vector<vector<hotspot>> allhotspots(num_targets, vector<hotspot>());
    vector<array<int, hotspot_freq_cutoff>> kmer_lastposis(num_targets, array<int, hotspot_freq_cutoff>{});
    vector<int> kmer_lastindex(num_targets, 0);
    
    std::string Header = "";
    std::string lastHeader = "";
    
    std::size_t current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    vector<uint> groups(num_targets+1,0);
    uint groups_len = 0;
    uint hash_ ;
    size_t rlen = 0;
    uint8_t* StrLine = NULL;
    char base;
    auto SRread = bam_init1();
    while (file.nextLine(StrLine, rlen, SRread))
    {
        
        if ( ( (SRread->core.flag & BAM_FUNMAP) == 0 && SRread->core.tid >= 0 && SRread->core.tid < file.header->n_targets)  == false )
        {
            SRread->core.pos = 1000;
            Header = "Unmap";
        }
        else
        {
            Header = std::string(file.header->target_name[SRread->core.tid]);
        }
        
        if (lastHeader != Header)
        {
            WriteHotspots<EnableLock>(lastHeader, allhotspots);
            std::fill(kmer_lastindex.begin(), kmer_lastindex.end(), 0);
            for (auto& arr : kmer_lastposis)
            {
                arr.fill(-1001);
            }
            
            for (auto& arr: allhotspots)
            {
                arr.clear();
            }
        }
        lastHeader = Header;
              
        current_size = 0;
        for (int pos = 0; pos < rlen; ++pos)
        {
            base = seq_nt16_str[bam_seqi(StrLine,pos)];
            //if (base=='\n') break;
            kmer_read_31(base, ++current_size, current_kmer, reverse_kmer);
            
            if (current_size < klen) continue;
            
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
            
            SearchGroup(kmer_hash, kmer_multi_hash, larger_kmer, groups, groups_len, kmer_ranges);
            
            if (groups.size())
            {
                UpdateLocation(allhotspots,kmer_lastposis, kmer_lastindex, groups, groups_len , (int)SRread->core.pos);
            }
        }
        nReads+=1;
        if (nReads % 10000000 == 0)
        {
            cerr << "processed " << nReads / 1000000 << "M reads." << endl;
        }
    }
    
    
    
    WriteHotspots<EnableLock>(lastHeader, allhotspots);
    bam_destroy1(SRread);
    
    return ;
};


template <int dictsize>
template <class typefile>
void KmerProfile<dictsize>::count_kmer_(typefile &file, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg, const int nthreads)
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
void KmerProfile<dictsize>::Call(const char* inputfile, uint16* samplevecs, ull_atom &nBases, ull_atom &nReads, ull_atom &nBg, const int nthreads)
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
        
        if ( pathlen >= 3 && strcmp(file+(pathlen-3),".gz") == 0 )
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
            CramReader readsfile(file,"",nthreads);
            readsfile.LoadRegion(*regions, ifhla, ifunmap);
            readsfile.filtertag = 0;
            count_kmer_(readsfile, samplevecs, nBases, nReads, nBg, nthreads);
            readsfile.Close();
        }
    }
    
}

#endif /* KmerProfile_hpp */
