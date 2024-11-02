//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License. 
//  If you use this code, please cite our work.
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

bool KtableReader::nextLine_start(std::string &StrLine, const char startchar)
{
    bool ifget = 0;
    auto position = ftell(file);
    
    if (fgets((char*)StrLine.c_str(), MAX_LINE, file) !=NULL )
    {
        if (StrLine[0] == startchar)
        {
            ifget = 1;
        }
        else
        {
            fseek(file, position, SEEK_SET);
        }
    }
    
    return ifget;
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

bool KtableReader::nextLine_genegroup(std::string &StrLine)
{
    bool ifget = 0;
    auto position = ftell(file);
    
    if (fgets((char*)StrLine.c_str(), MAX_LINE, file) !=NULL )
    {
        if (StrLine[0] == '@')
        {
            ifget = 1;            
        }
        else
        {
            fseek(file, position, SEEK_SET);
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
    
    StrLine += '\n';
    
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
