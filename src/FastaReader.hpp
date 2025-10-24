//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License. 
//  If you use this code, please cite our work.
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
    bool nextLine_prt(const char* &StrLine, size_t &rlen){return 0;};
    
    void Load();
        
    void Close();
    
private:
    
    std::fstream fafile;
            
};


#endif /* Fasta_hpp */
