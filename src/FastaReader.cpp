//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License. 
//  If you use this code, please cite our work.
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
