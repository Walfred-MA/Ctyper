//
//  fastq.hpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 12/1/21.
//  Copyright © 2021 USC_Mark. All rights reserved.
//

#ifndef fastq_hpp
#define fastq_hpp

#include "readsfile.hpp"

using namespace std;

#include <zlib.h>

class fastq: public readsfile
{
    
public:
    
    fastq(const char* inputfile):readsfile(inputfile){Load();};
    
    bool nextLine(std::string &StrLine);
    
    void Load();
    
    void Close();
    
private:
    
    gzFile fafile;
        
};
#endif /* fastq_hpp */

