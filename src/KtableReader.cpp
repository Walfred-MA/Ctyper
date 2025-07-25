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
    
    if (ifbgzf)
    {
        file_bgzf = bgzf_open(filepath, "r");
        if (!file_bgzf)
        {
            std::cerr << "ERROR: Could not open " << filepath << " with bgzf_open.\n";
            std::_Exit(EXIT_FAILURE);
        }
    }
    else
    {
        file = fopen(filepath, "r");
        if (!file)
        {
            std::cerr << "ERROR: Could not open " << filepath << " with fopen.\n";
            std::_Exit(EXIT_FAILURE);
        }
    }
}

bool KtableReader::nextLine(std::string &StrLine)
{
    if (ifbgzf)
    {
        return (bool)(bgzf_fgets((char*)&StrLine[0], MAX_LINE, file_bgzf));
    }
    else
    {
        return (bool)(fgets((char*)&StrLine[0], MAX_LINE, file));
    }
    
}


bool KtableReader::nextLine_start(std::string &StrLine, const char startchar)
{
    bool ifget = 0;
    size_t position;
    if (ifbgzf)
    {
        position = bgzf_tell(file_bgzf);
    }
    else
    {
        position = ftell(file);
    }
    
    if ( nextLine(StrLine) )
    {
        if (StrLine[0] == startchar)
        {
            ifget = 1;
        }
        else
        {
            if (ifbgzf)
            {
                bgzf_seek(file_bgzf, position, SEEK_SET);
            }
            else
            {
                fseek(file, position, SEEK_SET);
            }
        }
    }
    
    return ifget;
}

bool KtableReader::Find(std::string &StrLine, const char startchar)
{
    bool ifget = 0;
    size_t position;
    if (ifbgzf)
    {
        position = bgzf_tell(file_bgzf);
    }
    else
    {
        position = ftell(file);
    }
    
    while ( nextLine(StrLine) )
    {
        if (StrLine[0] == startchar)
        {
            ifget = 1;
            
            if (ifbgzf)
            {
                bgzf_seek(file_bgzf, position, SEEK_SET);
            }
            else
            {
                fseek(file, position, SEEK_SET);
            }
            
            break;
        }
        
        if (ifbgzf)
        {
            position = bgzf_tell(file_bgzf);
        }
        else
        {
            position = ftell(file);
        }

    }
    
    return ifget;
}

bool KtableReader::Find(std::string &StrLine, const char startchar1, const char startchar2)
{
    bool ifget = 0;
    size_t position;
    if (ifbgzf)
    {
        position = bgzf_tell(file_bgzf);
    }
    else
    {
        position = ftell(file);
    }
    
    while ( nextLine(StrLine) )
    {
        if (StrLine[0] == startchar1 || StrLine[0] == startchar2)
        {
            ifget = 1;
            
            if (ifbgzf)
            {
                bgzf_seek(file_bgzf, position, SEEK_SET);
            }
            else
            {
                fseek(file, position, SEEK_SET);
            }
            
            break;
        }
        
        if (ifbgzf)
        {
            position = bgzf_tell(file_bgzf);
        }
        else
        {
            position = ftell(file);
        }

    }
    
    return ifget;
}

bool KtableReader::nextLine_kmer(std::string &StrLine)
{
    bool ifget = 0;
    
    while ( nextLine(StrLine) )
    {
        
        if (StrLine[0] == '&' || StrLine[0] == '^')
        {
            ifget = 1;
            break;
        }
    }
    
    return ifget;
}
void KtableReader::Close()
{
    if (ifbgzf )
    {
        if (file_bgzf != NULL)
        {
            bgzf_close(file_bgzf);
        }
    }
    else
    {
        if (file != NULL)
        {
            fclose(file);
        }
    }
}

void KtableReader::Seek(const size_t pos)
{
    if (ifbgzf)
    {
        bgzf_seek(file_bgzf, pos, SEEK_SET);
    }
    else
    {
        fseek(file, pos, SEEK_SET);
    }
}


int KtableReader::GetLine(std::string &line)
{
    char buf[1024];
    line.clear();
    
    if (ifbgzf)
    {
        return FileReader::GetLine(file_bgzf, line);
    }
    else
    {
        return FileReader::GetLine(file, line);
    }
    
    return 0;
}



