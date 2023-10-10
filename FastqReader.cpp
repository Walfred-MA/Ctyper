//
//  FastqReader.cpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#include "FastqReader.hpp"

void FastqReader::Load()
{
    
    if (strlen(filepath)<2) return;
    
    fafile = gzopen(filepath, "rb");
    
    if (fafile == NULL)
    {
        std::cerr << "ERROR: Could not open " << filepath << " for reading.\n" << std::endl;
        std::_Exit(EXIT_FAILURE);
    }
    
}

bool FastqReader::nextLine(std::string &StrLine)
{
    if (StrLine.size() < 501) {StrLine.resize(501);}
        
    int bytes_read = gzread(fafile, &StrLine[0], 500);
    
    if (bytes_read <= 0) return -1;
        
    StrLine[bytes_read] = '\0';

    std::string line;
    
    size_t pos_read = 0;
    while (1)
    {
        if (pos_read == bytes_read)
        {
            int new_bytes_read = gzread(fafile, &StrLine[0], 500);
            
            if (new_bytes_read <= 0 ) return -1;
            
            bytes_read += new_bytes_read;
        }
        ++pos_read;
        if (StrLine[pos_read-1] == '\n') break;
    }
    
    while (1)
    {
        if (pos_read == bytes_read)
        {
            int new_bytes_read = gzread(fafile, &StrLine[0], 500);
            
            if (new_bytes_read <= 0 ) return -1;
            
            bytes_read += new_bytes_read;
        }
        ++pos_read;
        if (StrLine[pos_read-1] == '\n') break;
    }
    
    size_t start = pos_read;
    while (1)
    {
        if (pos_read == bytes_read)
        {
            int new_bytes_read = gzread(fafile, &StrLine[0], 500);
            
            if (new_bytes_read <= 0 ) return -1;
            
            bytes_read += new_bytes_read;
        }
        ++pos_read;
        if (StrLine[pos_read-1] == '\n') break;
    }
    
    StrLine = StrLine.substr(start, pos_read - start - 1);
    
    gzseek(fafile, - bytes_read + pos_read, SEEK_CUR);
    
    return 1;
}

void FastqReader::Close()
{
    gzclose(fafile);
}
