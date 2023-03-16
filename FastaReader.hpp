//
//  Fasta.hpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#ifndef Fasta_hpp
#define Fasta_hpp

#include "config.hpp"
#include "FileReader.hpp"

using namespace std;

class FastaReader: public FileReader
{
    
public:
    
    FastaReader(const char* inputfile, int index = -1):FileReader(inputfile){Load();};
    
    bool nextLine(std::string &StrLine);
    
    void Load();
        
    void Close();
    
private:
    
    std::fstream fafile;
            
};


#endif /* Fasta_hpp */
