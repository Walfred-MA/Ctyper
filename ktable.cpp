//
//  table.cpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 1/30/23.
//  Copyright © 2023 USC_Mark. All rights reserved.
//

#include "ktable.hpp"

//reading reference
void ktable::Load()
{
    if (strlen(filepath)<2) return;
    
    fafile.open(filepath, ios_base::in);
    
    if (!fafile) {
        
        std::cerr << "ERROR: Could not open " << filepath << " for reading.\n" << std::endl;
        std::_Exit(EXIT_FAILURE);
        
        return;
    }
}

bool ktable::nextLine(std::string &StrLine)
{
    return (bool)( getline(fafile,StrLine));
}

void ktable::Close()
{
    fafile.close();
}
