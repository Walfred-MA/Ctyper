//
//  main.cpp
//  CTyper
//
//  Created by Wangfei MA on 2/13/23.
//

#include <iostream>
#include <string>
#include <unordered_set>
#include <chrono>

#include "Processor.hpp"




int main(int argc, const char * argv[]) {
    
    
    std::vector<std::string> inputfiles;
        
    std::vector<std::string> outputfiles;
    
    std::vector<float> depths;
    
    std::unordered_set<std::string> genes;
    
    std::vector<char *> regions;
    
    std::string kmatrixfile="";
    
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
        
        else if (strcmp(Argument, "-m")==0 or strcmp(Argument, "--matrix")==0)
        {
            kmatrixfile = string(argv[i]);
        }
        
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
        
        else if (strcmp(Argument, "-g")==0 or strcmp(Argument, "--gene")==0)
        {
            genes.insert(argv[i]);
        }

        else if (strcmp(Argument, "-G")==0 or strcmp(Argument, "--Genes")==0)
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
                genes.insert(line);
            }
        }
        
        else if (strcmp(Argument, "-n")==0 or strcmp(Argument, "--nthreads")==0)
        {
            nthreads=(int)atoi(argv[i]);
        }
        
        else if (strcmp(Argument, "-d")==0 or strcmp(Argument, "--depth")==0)
        {
            depths.push_back(atof(argv[i]));
        }

        else if (strcmp(Argument, "-D")==0 or strcmp(Argument, "--Depth")==0)
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
                depths.push_back(atof(line.c_str()));
            }
        }
        
    }
    
    if (!inputfiles.size()) return 1;
    
    auto begin = std::chrono::high_resolution_clock::now();
    
    Processor<32> *processor = new Processor<32>(inputfiles, outputfiles, depths, kmatrixfile, genes, regions,  nthreads);

    processor->Run();
    
    auto end = std::chrono::high_resolution_clock::now();
    
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        
    cout<<"finished running at time: "<<elapsed.count()* 1e-9 <<endl;
    
    
    return 0;
}
