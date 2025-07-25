//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License.
//  If you use this code, please cite our work.
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
#include <htslib/bgzf.h>

#include "config.hpp"


using namespace std;

inline bool bgzf_fgets(char* line, size_t max_len, BGZF* bgzf)
{
    
    size_t count = 0;
    char c;
    
    while (count < max_len - 1)
    {
        int ret = bgzf_read(bgzf, &c, 1);
        if (ret <= 0)
        {
            if (count == 0)
            {
                return false;
            }
            break;
        }
        
        line[count++] = c;
        
        if (c == '\n') break;
    }
    
    line[count] = '\0';
    return true;
}

class FileReader
{
    
public:
    
    const char *filepath;
    int num_seq = 0;
    
    std::vector<spair> seqs;
    
    FileReader(const char* inputfile, int index = -1):filepath(inputfile) {};
    
    bool nextLine(std::string &StrLine);
    
    void Load();
    
    void Close();
    
    void Reset();
    
    int GetLine(FILE* fptr, std::string &line);
    
    int GetLine(BGZF* fptr, std::string &line);
    
private:

    std::fstream fafile;
    
};


#endif /* FileReader_hpp */

