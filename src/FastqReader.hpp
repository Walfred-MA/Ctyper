//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License. 
//  If you use this code, please cite our work.
//   

#ifndef FastqReader_hpp
#define FastqReader_hpp

#include "config.hpp"
#include "FileReader.hpp"

using namespace std;

#include <zlib.h>
#include <mutex>
#include <vector>




class FastqReader
{
    
public:
    
    FastqReader(const char* inputfile)
    {
        std::string all(inputfile);
        size_t start = 0;
        size_t pos;
        while ((pos = all.find(',', start)) != std::string::npos)
        {
            std::string part = all.substr(start, pos - start);
            if (!part.empty()) filepathes.push_back(part);
            start = pos + 1;
        }
        if (start < all.size()) {
            std::string part = all.substr(start);
            if (!part.empty()) filepathes.push_back(part);
        }
        filepath = filepathes[0].c_str();
        Load();
    };
    
    bool nextLine(std::string &StrLine){return 0;};
    bool nextLine(char*& strLine, size_t& rlen, char*& header);
    
    void Load();
    bool LoadBuffer();
    void Close();
    
    const char* filepath;
    uint fileindex = 0;
    std::vector<string> filepathes;
private:
    
    string buffer = string(10000000, '\0');
    size_t buffer_remain =0 ;
    size_t buffer_pos = 0;
    std::mutex IO_lock;
    gzFile fafile = nullptr;
    FILE* file= nullptr;
    bool ifgz ;
    
};

#endif /* FastqReader_hpp */
