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
    
    if (StrLine.size() < 10000) {StrLine.resize(10000);}
    
    char* gzget = gzgets(fafile,(char*)StrLine.c_str(),10000);
    
    bool ifgets = (gzget != NULL && !gzeof(fafile));
    
    if (ifgets && StrLine[0]!='@')
    {
        int c;
        while((c = gzgetc(fafile)) != '\n' && c != -1);
        while((c = gzgetc(fafile)) != '\n' && c != -1);
    }
    
    return ifgets;
}

void FastqReader::Close()
{
    
    gzclose(fafile);
    
}
