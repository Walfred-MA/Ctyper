//
//  fastq.cpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 12/1/21.
//  Copyright © 2021 USC_Mark. All rights reserved.
//

#include "fastq.hpp"

//reading reference
void fastq::Load()
{
    
    if (strlen(filepath)<2) return;
    
    fafile = gzopen(filepath, "rb");
    
    if (fafile == NULL)
    {
        std::cerr << "ERROR: Could not open " << filepath << " for reading.\n" << std::endl;
        std::_Exit(EXIT_FAILURE);
    }
    
}

bool fastq::nextLine(std::string &StrLine)
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

void fastq::Close()
{
    
    gzclose(fafile);
    
}
