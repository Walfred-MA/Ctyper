//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License.
//  If you use this code, please cite our work.
//

#include "CramReader.hpp"
#include <string>
#include <chrono>
#include <thread>
#include <unistd.h>
#include <fcntl.h>
#include <mutex>

using namespace std;

constexpr const char* NULLSTDERR = "/dev/null";
std::mutex stderr_mutex;
bool stderr_is_redirected = false;
int old_stderr = -1;

void suppress_stderr()
{
    std::lock_guard<std::mutex> lock(stderr_mutex);

    if (stderr_is_redirected)
        return;

    int old_stderr_ = dup(STDERR_FILENO);  // Save current stderr
    if (old_stderr_ == -1)
        return;

    int dev_null = open(NULLSTDERR, O_WRONLY);
    if (dev_null == -1) {
        close(old_stderr_);
        return;
    }

    if (dup2(dev_null, STDERR_FILENO) == -1) {
        close(dev_null);
        close(old_stderr_);
        return;
    }

    close(dev_null);
    old_stderr = old_stderr_;
    stderr_is_redirected = true;
}

void restore_stderr()
{
    std::lock_guard<std::mutex> lock(stderr_mutex);

    if (!stderr_is_redirected)
        return;

    dup2(old_stderr, STDERR_FILENO);  // Restore original stderr
    close(old_stderr);
    old_stderr = -1;
    stderr_is_redirected = false;
}

void LoadHLARegion(hts_idx_t* index, bam_hdr_t* header, htsFile* samfile, std::vector<char *>& regions)
{
    unordered_set<string> HLA_chrs;
    
    regions.erase(
        std::remove_if(regions.begin(), regions.end(),
                       [](char* contig)
                       {
                           std::string s(contig);
                           if (s.find("HLA") != std::string::npos ||
                               s.find("hla") != std::string::npos ||
                               s.find("MHC") != std::string::npos ||
                               s.find("HSCHR6") != std::string::npos ||
                               s.find("GL000") != std::string::npos)
                           {
                               return true;
                           }
                           return false;
                       }),
        regions.end()
    );
    
    std::vector<char*> hla_regions;
    
    for (int i = 0; i < header->n_targets; ++i)
    {
        std::string contig(header->target_name[i]);

        if (contig.find("HLA") != std::string::npos ||
            contig.find("hla") != std::string::npos ||
            contig.find("MHC") != std::string::npos ||
            contig.find("HSCHR6") != std::string::npos ||
            contig.find("GL000") != std::string::npos)
        {
            if (HLA_chrs.find(contig) == HLA_chrs.end())
                regions.push_back(strdup(header->target_name[i]));
        }
    }
}


uint64_t CramReader::TotalReads()
{
    samFile* fp = sam_open(filepath, "r");

    if (refpath.length())
    {
        auto ret = hts_set_fai_filename(fp, refpath.c_str());
        if (ret != 0) {
                std::cerr << "Failed to set reference: " << refpath << std::endl;
            }
        
    }
    sam_hdr_t* header = sam_hdr_read(fp);
    hts_idx_t* idx = sam_index_load(fp, filepath);
    
    uint64_t mapped = 0, unmapped = 0;
    for (int tid = 0; tid < header->n_targets; ++tid)
    {
        
        uint64_t mapped_ = 0, unmapped_ = 0;
        hts_idx_get_stat(idx, tid, &mapped_, &unmapped_);
        mapped += mapped_;
        unmapped += unmapped_;
    }

    hts_idx_destroy(idx);
    sam_hdr_destroy(header);
    sam_close(fp);
    
    return mapped + unmapped;
}
/*
double CramReader::QuickCoverage()
{
    if (!filepath) return 0.0;
    
    CramReader tmp_fp(filepath, refpath);
    
    std::vector<char*> depth_region_dummy = DEPTH_REGION;
    
    tmp_fp.LoadRegion(depth_region_dummy);

    unsigned long long total_length = 0;
    uint8_t * StrLine;
    size_t rlen = 0;
    tmp_fp.useunmap = 0;
    tmp_fp.ifindex = 1;
    while (tmp_fp.nextLine_prt(StrLine,rlen))
    {
        total_length += rlen - 30;
    }

    tmp_fp.Close();
    
    return (double)total_length;
}
 
 bool hasHLARegion(std::vector<char*>& workregions)
 {
     for (auto it = workregions.begin(); it != workregions.end(); ++it)
     {
         if (string(*it) == "HLA")
         {
             workregions.erase(it);
             return true;
         }
     }
     return false;
 }


 bool hasUnmapRegion(std::vector<char*>& workregions)
 {
     for (auto it = workregions.begin(); it != workregions.end(); ++it)
     {
         if (*it && strncmp(*it, "Unmap", 5) == 0)
         {
             workregions.erase(it);
             return true;
         }
     }
     return false;
 }

*/
void CramReader::Load()
{
    if (strlen(filepath)<2) return;
    
    std::string indexpath;
    if (strlen(filepath) >= 5 && std::string(filepath).substr(strlen(filepath) - 5) == ".cram")
    {
        indexpath = std::string(filepath) + ".crai";
    }
    else if (strlen(filepath) >= 4 && std::string(filepath).substr(strlen(filepath) - 4) == ".bam")
    {
        indexpath = std::string(filepath) + ".bai";
    }
    else if (strlen(filepath) >= 4 && std::string(filepath).substr(strlen(filepath) - 4) == ".sam")
    {
        indexpath = std::string(filepath) + ".sai";  // placeholder
    }
    else
    {
        std::cerr << "Unsupported file type: " << filepath << std::endl;
        std::_Exit(EXIT_FAILURE);
    }
            
    int attemp = 200;
    do
    {
        
        if (samfile)
        {
            hts_close(samfile); // free previous handle
            samfile = NULL;
        }
        if (indexdata)
        {
            hts_idx_destroy(indexdata);
            indexdata = NULL;
        }
        
        try
        {
            
            samfile = hts_open(filepath, "r");
            hts_set_threads(samfile, nthread);
            
            if (!refpath.empty())
            {
        
                if (hts_set_fai_filename(samfile, refpath.c_str()) != 0)
                {
                    throw std::runtime_error("Failed to set reference FASTA: " + refpath);
                }
            }
            
            indexdata = sam_index_load2(samfile, filepath, indexpath.c_str());
        }
        catch (const std::bad_alloc&)
        {
             samfile = NULL;
             indexdata = NULL;
        }
            
            
        if (!samfile || !indexdata)
        {
            std::cerr << "ERROR: Could not open " << filepath << " for reading:  attemp: "<< attemp <<" \n" << std::endl;
            std::this_thread::sleep_for(std::chrono::seconds(10));
        }
        
    } while ( (!samfile || !indexdata) && attemp-- > 0) ;
        
    header = sam_hdr_read(samfile);
}


void CramReader::LoadRegion(std::vector<char *>& bedregions, int hla = 0, int unmap = 0)
{
    workregions = bedregions;
    
    if (workregions.size() )
    {
        ifindex = 1;
        if (hla) LoadHLARegion(indexdata, header, samfile, workregions);
        useunmap = unmap;
        //suppress_stderr();
        iter = sam_itr_regarray(indexdata, header , workregions.data(), (int)workregions.size());
        //restore_stderr();
        
    }
        
    //strcpy(  regions[0], "chr1:1209512-1409512");
    
    //hts_reglist_t * reglist = hts_reglist_create_cram(regions, 1, &rcount, indexdata);
    
}

bool CramReader::nextLine_prt(uint8_t * &StrLine, size_t &rlen)
{
    //lock_guard<mutex> IO(IO_lock);

    while (true)
    {
        if (ifindex)
        {
            if (startumap)
            {
                iter = sam_itr_queryi(indexdata, HTS_IDX_NOCOOR, 0, 0);
                if (!iter || sam_itr_next(samfile, iter, SRread) < 0)
                {
                    return false;
                }
            }
            else
            {
                auto get = sam_itr_multi_next(samfile, iter, SRread);
                if ( get < 0)
                {
                    
                    if (useunmap == 1)
                    {
                        hts_itr_destroy(iter);
                        iter = sam_itr_queryi(indexdata, HTS_IDX_NOCOOR, 0, 0);
                        useunmap = 2;
                        if (!iter || sam_itr_next(samfile, iter, SRread) < 0)
                        {
                            return false;
                        }
                        
                    }
                    else
                    {
                        return false;
                    }
                }
            }
        }
        else
        {
            if (sam_read1(samfile,header ,SRread)<0)
            {
                return false;
            }
        }
        
        if ((SRread->core.flag & filtertag) == 0)
        {
            break;
        }
    }

    rlen = SRread->core.l_qseq;


    //if (StrLine.size() <= readLength ) StrLine.resize(readLength+1);
    //uint8_t *q = bam_get_seq(SRread);

    //for (int i=0; i < readLength; i++) {StrLine[i]=seq_nt16_str[bam_seqi(q,i)];    }

    //StrLine[readLength] = '\0';
    StrLine = bam_get_seq(SRread);

    return true;
}
 
bool CramReader::nextLine(uint8_t * &StrLine, size_t &rlen, bam1_t* &SRread)
{
    std::lock_guard<std::mutex> IO(IO_lock);
    
    while (true)
    {
        if (ifindex)
        {
            if (startumap)
            {
                iter = sam_itr_queryi(indexdata, HTS_IDX_NOCOOR, 0, 0);
                if (!iter || sam_itr_next(samfile, iter, SRread) < 0)
                {
                    return false;
                }
            }
            else
            {
                auto get = sam_itr_multi_next(samfile, iter, SRread);
                
                if (get < 0)
                {
                    if (useunmap)
                    {
                        hts_itr_destroy(iter);
                        iter = sam_itr_queryi(indexdata, HTS_IDX_NOCOOR, 0, 0);
                        useunmap = 0;
                        if (!iter || sam_itr_next(samfile, iter, SRread) < 0)
                        {
                            return false;
                        }
                        
                    }
                    else
                    {
                        return false;
                    }
                }
            }
        }
        else
        {
            if (sam_read1(samfile,header ,SRread)<0)
            {
                return false;
            }
        }
        
        if ((SRread->core.flag & filtertag) == 0)
        {
            
            break;
        }
    }
        
    rlen = SRread->core.l_qseq;
    StrLine = bam_get_seq(SRread);

    return true;
}


bool CramReader::nextLine(std::string &StrLine)
{
    //lock_guard<mutex> IO(IO_lock);

    while (true)
    {
        if (ifindex)
        {
            if (startumap)
            {
                iter = sam_itr_queryi(indexdata, HTS_IDX_NOCOOR, 0, 0);
                if (!iter || sam_itr_next(samfile, iter, SRread) < 0)
                {
                    return false;
                }
            }
            else
            {
                auto get = sam_itr_multi_next(samfile, iter, SRread);
                
                if ( get < 0)
                {
                    
                    if (useunmap)
                    {
                        hts_itr_destroy(iter);
                        iter = sam_itr_queryi(indexdata, HTS_IDX_NOCOOR, 0, 0);
                        useunmap = 0;
                        if (!iter || sam_itr_next(samfile, iter, SRread) < 0)
                        {
                            return false;
                        }
                        
                    }
                    else
                    {
                        return false;
                    }
                }
            }
        }
        else
        {
            if (sam_read1(samfile,header ,SRread)<0)
            {
                return false;
            }
        }
        
        if ((SRread->core.flag & filtertag) == 0)
        {
            
            break;
        }
    }
    
    size_t readLength = SRread->core.l_qseq;

    if (StrLine.size() <= readLength ) StrLine.resize(readLength+1);
    uint8_t *q = bam_get_seq(SRread);

    for (int i=0; i < readLength; i++) {StrLine[i]=seq_nt16_str[bam_seqi(q,i)];    }

    StrLine[readLength] = '\0';

    return true;
}

