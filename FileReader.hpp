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
    
    int GetLine(FILE* fptr, std::string &line) {
        char buf[1024];
        line = "";
        
        while (true) {
            if (fgets(buf, sizeof(buf), fptr) == nullptr) {
                // fgets failed, check if it's an error or end-of-file
                if (feof(fptr)) {
                    // End-of-file reached
                    break;
                } else {
                    // I/O error occurred
                    perror("Error reading from file");
                    return -1;
                }
            }
            
            line += buf;
            
            if (line.back() == '\n') {
                // Remove the newline character from the end of the line
                line.pop_back();
                break;
            }
            
            // If there are no more characters in the file, break the loop
            if (feof(fptr)) {
                break;
            }
        }
        
        // Return 1 if the file has more lines, 0 if it's the last line, -1 if an error occurred
        return !feof(fptr);
    }

  
    /*int GetLine(FILE* fptr, string &line) {
    char buf[1025];
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
    char c = fgetc(fptr);
    return (!feof(fptr));
  } */
    
};


#endif /* FileReader_hpp */
