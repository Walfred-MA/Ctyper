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

class KtableReader: public FileReader
{
    
public:
    
    KtableReader(const char* inputfile):FileReader(inputfile)
    {
        Load();
    };
    ~KtableReader()
    {
        Close();
    }
    
    bool nextLine(std::string &StrLine);

    bool nextLine_start(std::string &StrLine, const char startchar);

    bool nextLine_kmer(std::string &StrLine);
    
    bool nextLine_norm(std::string &StrLine);
    
    bool nextLine_genename(std::string &StrLine);
    
    bool nextLine_genegroup(std::string &StrLine);
    
    void Load();
    
    void Close();
    
    void Seek(const size_t pos);
                    
private:
        
    FILE* file;
    
};

#endif /* KtableReader_hpp */
