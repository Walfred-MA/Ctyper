//
//  main.cpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 10/13/21.
//  Copyright © 2021 USC_MarkLab. All rights reserved.
//

#include <iostream>
#include <string>

#include "KmerCounter.hpp"
//#include "CramReader.hpp"

template <int dictsize>
void run(std::vector<std::string> inputfiles, std::vector<std::string> outputfiles, std::vector<std::string> prefixes, std::string targetfile, std::string matrixfile, std::vector<char *> regions, std::vector<float> depths, const int nthreads)
{
    
    const int kmer_size = 31;
    
    kmer_counter<dictsize> counter(kmer_size);
    
    if (targetfile.size())
    {
        counter.read_target(targetfile.c_str(), matrixfile.c_str());
    }
    
    counter.read_files(inputfiles, outputfiles, prefixes, depths, nthreads);
    
}


int main(int argc, const char * argv[]) {
    
    
    std::vector<std::string> inputfiles;
    
    std::vector<std::string> prefixes;
    
    std::vector<std::string> targetfiles;
    
    std::vector<std::string> outputfiles;
    
    std::vector<char *> regions;
    
    std::vector<float> depths;
        
    std::string matrixfile="";
    std::string targetfile="";
    
    const char* Argument;
    int nthreads = 1;
    
    
    for (int i = 1; i < argc ; i++)
    {
        
        if (argv[i][0] == '-')
        {
            Argument = argv[i];
        }
        else if (strcmp(Argument, "-i")==0 or strcmp(Argument, "--input")==0)
        {
            inputfiles.push_back(argv[i]);
        }
        else if (strcmp(Argument, "-I")==0 or strcmp(Argument, "--Inputs")==0)
        {
            std::ifstream pathfile(argv[i]);

            if(!pathfile)
            {
                std::cout<<"Error opening input file"<<std::endl;
                return -1;
            }
            std::string line;
            while (std::getline(pathfile, line))
            {
                inputfiles.push_back(line);
            }
        }
        
        else if (strcmp(Argument, "-t")==0 or strcmp(Argument, "--target")==0)
        {
            targetfile = string(argv[i]);
        }
        
        /*
        else if (strcmp(Argument, "-T")==0 or strcmp(Argument, "--Targets")==0)
        {
            
            
            std::ifstream pathfile(argv[i]);
            std::string line;
            if(!pathfile)
            {
                std::cout<<"Error opening target file"<<std::endl;
                return -1;
            }
            while (std::getline(pathfile, line))
            {
                targetfiles.push_back(line);
            }
        }
        */
        
        else if (strcmp(Argument, "-o")==0 or strcmp(Argument, "--output")==0)
        {
            outputfiles.push_back(argv[i]);
        }

        else if (strcmp(Argument, "-O")==0 or strcmp(Argument, "--Outputs")==0)
        {
            std::ifstream pathfile(argv[i]);

            if(!pathfile)
            {
                std::cout<<"Error opening output file"<<std::endl;
                return -1;
            }
            std::string line;
            while (std::getline(pathfile, line))
            {
                outputfiles.push_back(line);
            }
        }
        
        else if (strcmp(Argument, "-p")==0 or strcmp(Argument, "--pref")==0)
        {
            prefixes.push_back(argv[i]);
        }
        
        else if (strcmp(Argument, "-P")==0 or strcmp(Argument, "--Prefs")==0)
        {
            std::ifstream pathfile(argv[i]);
            std::string line;
            if(!pathfile)
            {
                std::cout<<"Error opening target file"<<std::endl;
                return -1;
            }
            while (std::getline(pathfile, line))
            {
                prefixes.push_back(line);
            }
        }

        else if (strcmp(Argument, "-b")==0 or strcmp(Argument, "--bed")==0)
        {
            regions.push_back( (char*) malloc( (strlen(argv[i])+1)*sizeof(char) ) );
            
            strcpy(regions[0], argv[i]);
            
        }
        else if (strcmp(Argument, "-B")==0 or strcmp(Argument, "--Bed")==0)
        {
            std::ifstream pathfile(argv[i]);

            if(!pathfile)
            {
                std::cout<<"Error opening bed file"<<std::endl;
                return -1;
            }
            std::string line;
            while (std::getline(pathfile, line))
            {
                regions.push_back((char*) malloc((line.length()+1)*sizeof(char))) ;
                strcpy(regions[regions.size()-1], line.c_str());
            }
        }
        
        else if (strcmp(Argument, "-d")==0 or strcmp(Argument, "--depth")==0)
        {
            depths.push_back((float)atof(argv[i]));
        }
        
        else if (strcmp(Argument, "-D")==0 or strcmp(Argument, "--depths")==0)
        {
            std::ifstream pathfile(argv[i]);
            std::string line;
            if(!pathfile)
            {
                std::cout<<"Error opening target file"<<std::endl;
                return -1;
            }
            while (std::getline(pathfile, line))
            {
                depths.push_back((float)atof(line.c_str()));
            }
        }
        
        else if (strcmp(Argument, "-m")==0 or strcmp(Argument, "--depth")==0)
        {
            matrixfile = string(argv[i]);
        }
        
        else if (strcmp(Argument, "-n")==0 or strcmp(Argument, "--nthreads")==0)
        {
            nthreads=(int)atoi(argv[i]);
        }
        
    }
    
    if (!inputfiles.size()) return 1;
    

    run<32>(inputfiles, outputfiles, prefixes, targetfile, matrixfile, regions, depths, nthreads);
    /*
    if (kmer_size<=32)
    {
        run<32>(inputfiles, targetfiles, outputfiles, regions, kmer_size);
    }
    
    else if (kmer_size<=64)
    {
        run<64>(inputfiles, targetfiles, outputfiles, regions, kmer_size);
    }
    */
    
    return 0;
}
