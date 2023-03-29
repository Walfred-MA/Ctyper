//
//  KtableReader.hpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
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
    
    bool nextLine_kmer(std::string &StrLine);
    
    bool nextLine_norm(std::string &StrLine);
    
    bool nextLine_genename(std::string &StrLine);
    
    void Load();
    
    void Close();
    
    void Seek(const size_t pos);
                    
private:
        
    FILE* file;
    
};

#endif /* KtableReader_hpp */
