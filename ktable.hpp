//
//  table.hpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 1/30/23.
//  Copyright © 2023 USC_Mark. All rights reserved.
//

#ifndef ktable_hpp
#define ktable_hpp


#include "readsfile.hpp"

using namespace std;

#include <zlib.h>
#include <vector>

class ktable: public readsfile
{
    
public:
    
    ktable(const char* inputfile):readsfile(inputfile)
    {
        Load();
    };
    ~ktable()
    {
        Close();
    }
    
    bool nextLine(std::string &StrLine);
    
    void Load();
    
    void Close();
    
    void Indexposi();
            
private:
        
    std::fstream fafile;
    
};
#endif /* fastq_hpp */
