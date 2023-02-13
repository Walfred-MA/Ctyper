//
//  readsfile.h
//  kmer_haplotyping
//
//  Created by Wangfei MA on 12/1/21.
//  Copyright © 2021 USC_Mark. All rights reserved.
//

#ifndef readsfile_h
#define readsfile_h

#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <utility>
#include <cstring>

#define spair std::pair<std::string,std::string>

using namespace std;



class readsfile
{
    
public:
    
    const char *filepath;
    int num_seq = 0;
    
    std::vector<spair> seqs;
    
    readsfile(const char* inputfile, int index = -1):filepath(inputfile) {};
    
    virtual bool nextLine(std::string &StrLine)=0;
    
    virtual void Load()=0;
    
    virtual void Close()=0;
        
};


#endif /* readsfile_hpp */
