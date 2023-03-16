//
//  Fasta.cpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#include "FastaReader.hpp"


void FastaReader::Load()
{
    if (strlen(filepath)<2) return;
    
    fafile.open(filepath, ios_base::in);
    
    if (!fafile) {
        
        std::cerr << "ERROR: Could not open " << filepath << " for reading.\n" << std::endl;
        std::_Exit(EXIT_FAILURE);
        
        return;
    }
    
}

bool FastaReader::nextLine(std::string &StrLine)
{
    return (bool)( getline(fafile,StrLine));
}

void FastaReader::Close()
{
    
    fafile.close();
    
}
