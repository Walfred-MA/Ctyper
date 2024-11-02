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

using namespace std;
void CramReader::LoadRegion(std::vector<char *>& bedregions)
{
        
    workregions = &bedregions;
    
    if (strlen(filepath)<2) return;
    
    std::string indexpath = std::string(filepath)+".crai";
            
    int attemp = 200;
    do
    {
	try
	{
	     samfile = hts_open(filepath, "r");
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
 
bool CramReader::nextLine_prt(uint8_t * &StrLine, size_t &rlen)
{
    //lock_guard<mutex> IO(IO_lock);

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

    rlen = SRread->core.l_qseq;


    //if (StrLine.size() <= readLength ) StrLine.resize(readLength+1);
    //uint8_t *q = bam_get_seq(SRread);

    //for (int i=0; i < readLength; i++) {StrLine[i]=seq_nt16_str[bam_seqi(q,i)];	}

    //StrLine[readLength] = '\0';
    StrLine = bam_get_seq(SRread);

    return true;
}


