//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License. 
//  If you use this code, please cite our work.
//   

#include "FastqReader.hpp"

void FastqReader::Load()
{
    
    buffer.resize(MAX_LINE);
    
    if (strlen(filepath)<2) return;
    
    fafile = gzopen(filepath, "rb");
    
    if (fafile == NULL)
    {
        std::cerr << "ERROR: Could not open " << filepath << " for reading.\n" << std::endl;
        std::_Exit(EXIT_FAILURE);
    }
    
}

bool FastqReader:: nextLine_prt(const char* &StrLine, size_t &rlen)
{
    lock_guard<mutex> IO(IO_lock);
    
    //if (buffer.size() < 501) {buffer.resize(501);}
    
    if (buffer.size() - bytes_read < 2000000)
    {
        gzseek(fafile, - bytes_read + pos_read, SEEK_CUR);
        
        bytes_read = gzread(fafile, &buffer[0], 1000000);
        
        if (bytes_read <= 0) return 0;
        
        pos_read = 0;
    }
    
    
    while (1)
    {
        if (pos_read == bytes_read)
        {
            int new_bytes_read = gzread(fafile, &buffer[bytes_read], 1000000);
            
            if (new_bytes_read <= 0 ) return 0;
            
            bytes_read += new_bytes_read;
        }
        ++pos_read;
        if (buffer[pos_read-1] == '\n') break;
    }
    
    size_t start = pos_read;
    
    while (1)
    {
        if (pos_read == bytes_read)
        {
            int new_bytes_read = gzread(fafile, &buffer[bytes_read], 1000000);
            
            if (new_bytes_read <= 0 ) return 0;
            
            bytes_read += new_bytes_read;
        }
        ++pos_read;
        if (buffer[pos_read-1] == '\n') break;
    }
    
    rlen = pos_read - start;
    
    while (1)
    {
        if (pos_read == bytes_read)
        {
            int new_bytes_read = gzread(fafile, &buffer[bytes_read], 1000000);
            
            if (new_bytes_read <= 0 ) return 1;
            
            bytes_read += new_bytes_read;
        }
        ++pos_read;
        if (buffer[pos_read-1] == '\n') break;
    }
    
    StrLine =  &buffer[start];
    
    //gzseek(fafile, - bytes_read + pos_read, SEEK_CUR);
    
    return 1;
}

void FastqReader::Close()
{
    gzclose(fafile);
}
