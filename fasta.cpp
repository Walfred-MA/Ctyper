//
//  fasta.cpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 10/17/21.
//  Copyright © 2021 USC_Mark. All rights reserved.
//

#include "fasta.hpp"

//reading reference
void fasta::Load()
{
    if (strlen(filepath)<2) return;
    
    fafile.open(filepath, ios_base::in);
    
    if (!fafile) {
        
        std::cerr << "ERROR: Could not open " << filepath << " for reading.\n" << std::endl;
        std::_Exit(EXIT_FAILURE);
        
        return;
    }
    
}

bool fasta::nextLine(std::string &StrLine)
{
    return (bool)( getline(fafile,StrLine));
}

void fasta::Close()
{
    
    fafile.close();
    
}
