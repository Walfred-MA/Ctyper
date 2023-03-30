//
//  KtableReader.cpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#include "KtableReader.hpp"

void KtableReader::Load()
{
    if (strlen(filepath)<2) return;
    
    file = fopen(filepath, "r");
    
    if (!file)
    {
        
        std::cerr << "ERROR: Could not open " << filepath << " for reading.\n" << std::endl;
        std::_Exit(EXIT_FAILURE);
        
        return;
    }
}

bool KtableReader::nextLine(std::string &StrLine)
{
    return (bool)(fgets((char*)StrLine.c_str(), MAX_LINE, file));
}

bool KtableReader::nextLine_genename(std::string &StrLine)
{
    bool ifget = 0;
    
    while (fgets((char*)StrLine.c_str(), MAX_LINE, file) !=NULL )
    {
        if (StrLine[0] == '>')
        {
            ifget = 1;
            break;
        }
    }
    
    return ifget;
}

bool KtableReader::nextLine_kmer(std::string &StrLine)
{
    bool ifget = 0;
    
    while (fgets((char*)StrLine.c_str(), MAX_LINE, file) !=NULL )
    {
        if (StrLine[0] == '&')
        {
            ifget = 1;
            break;
        }
    }
    
    return ifget;
}

bool KtableReader::nextLine_norm(std::string &StrLine)
{
    bool ifget = 0;
    
    while (GetLine(file, StrLine)) 
    {
        if (StrLine[0] == '$')
        {
            ifget = 1;
            break;
        }
    }
    
    return ifget;
}

void KtableReader::Close()
{
  if (file != NULL) {
    fclose(file);
  }
  file=NULL;    
}

void KtableReader::Seek(const size_t pos)
{
    fseek(file, pos, SEEK_SET);
}
