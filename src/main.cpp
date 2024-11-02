//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License. 
//  If you use this code, please cite our work.
//   

#include <iostream>
#include <string>
#include <unordered_set>
#include <chrono>
#include <filesystem>

#include "Processor.hpp"
#include "config.hpp"

extern bool optioncorr;
bool optioncorr = 0;


void printHelp() 
{
    std::cout << "Usage: CTyper [options]\n\n"
              << "Options:\n"
              << "  -i, --input <file>         Input NGS file, allows CRAM(indexed), BAM(indexed), SAM, FASTQ and FASTA formats\n"
              << "  -I, --Inputs <file>        File with a list of input files\n"
              << "  -m, --matrix <file>        Path to the matrix database, need to come with its index file \n"
              << "  -o, --output <file>        Output file\n"
              << "  -O, --Outputs <file>       File with a list of output files\n"
              << "  -n, --nthreads <number>    Number of threads to use (default: 1) \n"
              << "  -N, --Nsubthreads <number> Number of sub-threads (default: 1) \n"
              << "  -d, --depth <value>        Predetermined depth value, not compatible with -b options\n"
              << "  -D, --Depth <file>         File with a list of predetermined depth values, not compatible with -b options\n"
              << "  -b, --background <file>    Background k-mer file, used to determine NGS coverage, not compatible with -d/-D options \n"
              << "  -c, --corr <0/1>           Set correction for biased k-mer option for NGS data (default: 0)\n"
              << "  -h, --help                 Show this help message\n"
              << std::endl;
}

int main(int argc, const char * argv[])
{
    
    std::vector<std::string> inputfiles;
        
    std::vector<std::string> outputfiles;
    
    std::vector<float> depths;
    
    std::unordered_set<std::string> genes;
    
    std::vector<char *> regions;
    
    std::string kmatrixfile="";
    std::string backgroundfile="";
    
    const char* Argument;
    int nthreads = 1;
    int Nsubthreads = 1;
    int window = 30;

    if (argc == 1 and (strcmp(Argument, "-h")==0 or strcmp(Argument, "--help")==0) )
    {
        printHelp();
        return 0;
    }
    
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

        else if (strcmp(Argument, "-bed") ==0 )
        {
            regions.push_back( (char*) malloc( (strlen(argv[i])+1)*sizeof(char) ) );
            
            strcpy(regions[0], argv[i]);
            
        }
        else if (strcmp(Argument, "-BED") ==0 )
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
        
        else if (strcmp(Argument, "-N")==0 or strcmp(Argument, "--Nsubthreads")==0)
        {
            Nsubthreads=(int)atoi(argv[i]);
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
        
        else if (strcmp(Argument, "-b")==0 or strcmp(Argument, "--background")==0)
        {
            backgroundfile = string(argv[i]);
        }
        
        
        else if (strcmp(Argument, "-w")==0 or strcmp(Argument, "--window")==0)
        {
            window = (int)atoi(argv[i]);
        }
        
        else if (strcmp(Argument, "-c")==0 or strcmp(Argument, "--corr")==0)
        {
            optioncorr = (int)atoi(argv[i]);
        }
        
    }
    
    if (!inputfiles.size()) return 1;
    
    auto begin = std::chrono::high_resolution_clock::now();
    
    Processor<32> *processor = new Processor<32>(inputfiles, outputfiles, depths, kmatrixfile, backgroundfile,genes, regions, window, nthreads, Nsubthreads);

    processor->Run();
    
    auto end = std::chrono::high_resolution_clock::now();
    
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        
    cout<<"finished running at time: "<<elapsed.count()* 1e-9 <<endl;
    
    
    return 0;
}
