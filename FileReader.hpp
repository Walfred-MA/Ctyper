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
    
    FileReader(const char* inputfile, int index = -1):filepath(inputfile) {};
    
    virtual bool nextLine(std::string &StrLine)=0;
    
    virtual void Load()=0;
    
    virtual void Close()=0;

  int GetLine(FILE* fptr, string &line) { 
    char buf[1024];
    line="";
    int res;
    while ((res=fscanf(fptr, "%1024[^\n]", buf))) {
      if (res) {
	line+=buf;
      }
      else {
	break;
      }
    }
    return line.size();
  }
        
};


#endif /* FileReader_hpp */
