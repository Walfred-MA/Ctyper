//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License. 
//  If you use this code, please cite our work.
//   

#ifndef CramReader_hpp
#define CramReader_hpp

#include "FileReader.hpp"

#include "htslib/sam.h"
#include <zlib.h>
#include <vector>
#include <string>
#include <mutex>
#include <string.h>

class CramReader: public FileReader {
   
public:
    
    CramReader(const char* inputfile): FileReader(inputfile)
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
    
    bool nextLine(std::string &StrLine){return 0;};
    bool nextLine_prt(uint8_t * &Str, size_t &rlen);
    
    void LoadRegion(std::vector<char *>& regions);
    void Load(){};
    void TotalReads(){};
    
    void Close(){ sam_close(samfile); };

private:
    
    htsFile *samfile = NULL;
    
    bam1_t * SRread = NULL;
    
    kstring_t* kstring = NULL;
    
    hts_itr_t * iter = NULL;
    
    sam_hdr_t *header = NULL;
    
    hts_idx_t *indexdata = NULL;
    
    std::vector<char*>* workregions;

    std::mutex IO_lock;
    

  bool ifindex = 0;
         
    unsigned long long nreads = 0;
};



#endif
