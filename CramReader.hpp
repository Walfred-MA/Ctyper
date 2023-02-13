//
//  CramReader.hpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 1/29/23.
//  Copyright © 2023 USC_Mark. All rights reserved.
//

/*
#ifndef CramReader_hpp
#define CramReader_hpp

#include "readsfile.hpp"
#include "include/samtools/htslib-1.16/htslib/sam.h"
#include "include/samtools/htslib-1.16/cram/cram_encode.h"
#include <zlib.h>
#include <vector>
#include <string>


class CramReader: public readsfile
{
    
public:
    
    CramReader(const char* inputfile):readsfile(inputfile)
    {
        kstring = new kstring_t();        
        SRread = bam_init1();
    };
    ~CramReader()
    {
        ks_free(kstring);
        bam_destroy1(SRread);
        hts_idx_destroy(indexdata);
        sam_hdr_destroy(header);
        hts_itr_multi_destroy(iter);
    }
    
    bool nextLine(std::string &StrLine);
    
    void Load(){};
    void Load(std::vector<char *>& regions);
    //void TotalReads();
    
    void Close(){ sam_close(samfile); };
       
private:
    
    htsFile *samfile;
    
    bam1_t * SRread;
    
    kstring_t* kstring;
    
    hts_itr_t * iter;
    
    sam_hdr_t *header;
    
    hts_idx_t *indexdata;
    
    std::vector<char*>* workregions;
    
    bool ifindex = 0;
         
    unsigned long long nreads = 0;
};

*/

#endif /* CramReader_hpp */
