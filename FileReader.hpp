//
//  ReadsFile.hpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#ifndef FileReader_hpp
#define FileReader_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <utility>
#include <cstring>

#include "config.hpp"


using namespace std;

class FileReader
{
    
public:
    
    const char *filepath;
    int num_seq = 0;
    
    std::vector<spair> seqs;
    FileReader() { };    
    FileReader(const char* inputfile, int index = -1):filepath(inputfile) {};
    void Init(const char* inputfile, int index=-1) {
      filepath=inputfile;
    }
    
    virtual bool nextLine(std::string &StrLine)=0;
    
    virtual void Load()=0;
    
    virtual void Close()=0;
        
};


#endif /* FileReader_hpp */
