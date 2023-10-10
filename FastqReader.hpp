//
//  FastqReader.hpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#ifndef FastqReader_hpp
#define FastqReader_hpp

#include "config.hpp"
#include "FileReader.hpp"

using namespace std;

#include <zlib.h>

class FastqReader: public FileReader
{
    
public:
    
    FastqReader(const char* inputfile):FileReader(inputfile){Load();};
    
    bool nextLine(std::string &StrLine);
    
    void Load();
    
    void Close();
    
private:
    
    gzFile fafile;
        
};

#endif /* FastqReader_hpp */
