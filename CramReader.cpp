//
//  CramReader.cpp
//   _haplotyping
//
//  Created by Wangfei MA on 1/29/23.
//  Copyright © 2023 USC_Mark. All rights reserved.
//

/*
#include "CramReader.hpp"
#include <string>

void CramReader::Load(std::vector<char *>& bedregions)
{
    ifindex = 1;
        
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
    
    iter = sam_itr_regarray(indexdata, header , workregions->data(), (int)workregions->size());
        
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
        if (cram_get_seq(samfile)<0) return false;

    }
    
    if (sam_format1(header, SRread, kstring)<0) return false;
        
    int string_count = 0;
    int col_counter = 0 ;
    

    for (char StrLine_c:std::string(ks_str(kstring)))
    {
        if (StrLine_c =='\t' )
        {
            col_counter++;
            continue;
        }
        //std::cout<<StrLine_c<<":"<<col_counter <<",";
        switch (col_counter)
        {
            case 9:
                StrLine[string_count++] = StrLine_c;
                break;
            case 10:
                StrLine[string_count++] = '\0';
                goto exit_loop;
            default:
                break;
        }
        
    }
    exit_loop: ;
        
    if (col_counter == 10)  return true;
    
    return false;
}
*/

/*
void CramReader::TotalReads()
{
    Load();
    
    nreads = 0;
    
    while(cram_get_seq(samfile))
    {
        nreads++;
    }
    
    Close();
    
}
 */
