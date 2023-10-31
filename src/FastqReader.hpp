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
#include <mutex>

class FastqReader: public FileReader
{
    
public:
    
    FastqReader(const char* inputfile):FileReader(inputfile){Load();};
    
    bool nextLine(std::string &StrLine){return 0;};
    bool nextLine_prt(const char* &StrLine, size_t &rlen);
    
    void Load();
        
    void Close();
    
private:
    
    std::mutex IO_lock;
    gzFile fafile;
    string buffer;
    size_t bytes_read = 0;
    size_t pos_read = 0;
};

#endif /* FastqReader_hpp */
