//
//  CramReader.cpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 1/29/23.
//  Copyright © 2023 USC_Mark. All rights reserved.
//

#include "CramReader.hpp"
#include <string>

using namespace std;
void CramReader::LoadRegion(std::vector<char *>& bedregions)
{
        
    workregions = &bedregions;
    
    if (strlen(filepath)<2) return;
    
    std::string indexpath = std::string(filepath)+".crai";
            
    samfile = hts_open(filepath, "r");
                
    if (!samfile) {
        
        std::cerr << "ERROR: Could not open " << filepath << " for reading.\n" << std::endl;
        std::_Exit(EXIT_FAILURE);
        
        return;
    }

        
    indexdata = sam_index_load2(samfile, filepath, indexpath.c_str());
    
    header = sam_hdr_read(samfile);
    
    if (workregions->size() ) 

    {
	ifindex = 1;

	iter = sam_itr_regarray(indexdata, header , workregions->data(), (int)workregions->size());
    }
        
    //strcpy(  regions[0], "chr1:1209512-1409512");
    
    //hts_reglist_t * reglist = hts_reglist_create_cram(regions, 1, &rcount, indexdata);
    
    
}
 
bool CramReader::nextLine(std::string &StrLine)
{
    
    if (ifindex)
    {
        if (sam_itr_multi_next(samfile , iter, SRread)<0) return false;
    }
    else
    {
        if (sam_read1(samfile,header ,SRread)<0) 
	{
		return false;
	}


    }
    int readLength= SRread->core.l_qseq;
    StrLine.resize(readLength+1);
    uint8_t *q = bam_get_seq(SRread);

    for (int i=0; i < readLength; i++) {StrLine[i]=seq_nt16_str[bam_seqi(q,i)];	}

    return true;
}


