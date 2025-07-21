//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License. 
//  If you use this code, please cite our work.
//   

#ifndef KtableReader_hpp
#define KtableReader_hpp

#include "config.hpp"
#include "FileReader.hpp"

using namespace std;

#include <zlib.h>
#include <vector>
#include <iostream>
#include <htslib/bgzf.h>
#include <unordered_set>

class KtableReader: public FileReader
{
    
public:
    
    KtableReader(const char* inputfile):FileReader(inputfile),
    ifbgzf(
        std::string(inputfile).size() >= 4 &&
        (std::string(inputfile).compare(std::string(inputfile).size() - 4, 4, ".bgz") == 0 ||
         std::string(inputfile).compare(std::string(inputfile).size() - 3, 3, ".bg") == 0)
    )
    {
        Load();
    };
    ~KtableReader()
    {
        //Close();
    }
    bool nextLine_bgzf(std::string &StrLine);
    bool nextLine(std::string &StrLine);

    bool nextLine_start(std::string &StrLine, const char startchar);
    bool nextLine_kmer(std::string &StrLine);
    
    bool Find(std::string &StrLine, const char startchar);
    bool Find(std::string &StrLine, const char startchar1, const char startchar2);
    
    void Load();
    
    void Close();
    
    int GetLine( std::string &line);
    
    void Seek(const size_t pos);
       
    const bool ifbgzf = 0;
    
    FILE* file= nullptr;
    BGZF* file_bgzf= nullptr;
};

#endif /* KtableReader_hpp */
