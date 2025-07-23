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
#include <sstream>
#include <vector>

#include "Processor.hpp"
#include "config.hpp"
#include "Interactive.hpp"

extern bool optioncorr;
bool optioncorr = 1;

extern int forceunmap;
int forceunmap = 0;

extern bool ifgeneregion;
bool ifgeneregion = 0;

extern int ifprofile;
int ifprofile = 0;

extern string vdatabase;
string vdatabase = "unknown";

extern string refpath;
string refpath = "";

extern string profilesummaryfile;
string profilesummaryfile = "";

extern int inputdep;
int inputdep = 0;


                                
void printHeader()
{
    std::cout << R"(
                            _____  _______  __     __  ______   ________  ______  
                           / ____| |__ __|  \ \   / /  |  __ \  |  ____|  |  __ \ 
                           | |       | |     \ \_/ /   | |__) | | |__     | |__) |
                           | |       | |      \   /    |  ___/  |  __|    |  _  / 
                           | |____   | |       | |     | |      | |____   | | \ \   
                           \_____|   |_|       |_|     |_|      |______|  |_|  \_\  
            )" << std::endl;
    
    std::cout << "Ctyper " + string(vctyper) + "\nPangenome allele database version: " + vdatabase + "\n"
        << "Author: Walfred (Wangfei) MA (walfredma0@gmail.com), Mark Chaisson Lab (mchaisso@usc.edu). Please contact me if you find any bug or usage question.\n"
        << "Please cite our work if you use this tool.\n"
        << "Source code: https://github.com/Walfred-MA/Ctyper\n"
        << "Using pangenome annotation, results can be converted into readable tables, FASTA files, VCF files, HLA/KIR/CYP2D6 nomenclatures, and mutation plots. "
        << "Tool scripts: https://github.com/Walfred-MA/Ctyper/tree/main/tools\n"
        << "Databases (genotyping matrices and pangenome allele annotations): https://zenodo.org/records/15779311\n\n";
}


void printHelp()
{
    printHeader();

    std::cout << "Usage: Ctyper " + string(vctyper) + " [options]\n\n";
    std::cout << "Required:\n";
    std::cout << "  -m, --matrix <file>          Path to the matrix database (requires <file>.index). If not provided, the tool runs in dry-run mode to only estimate NGS read depth.\n\n";

    std::cout << "Options:\n";
    std::cout << "  -i, --input <value>          Input NGS file; supports CRAM, BAM, SAM, FASTA, FASTQ, or Jellyfish formats.\n";
    std::cout << "  -I, --Inputs <file>          Path to a file listing multiple input files.\n\n";

    std::cout << "  -o, --output <value>         Output file (append if file exits)\n";
    std::cout << "  -O, --Outputs <file>         Path to a file listing output files corresponding to each input file.\n\n";

    std::cout << "  -n, --nthreads <number>      Number of threads to run different samples in parallel, parallel between sample is more efficient for large cohorts, especially when reading files on slower hard drives (default: 1).\n";
    std::cout << "  -N, --subthreads <number>    Number of threads to run each sample, due to file I/O bottleneck, suggest use only 1-4 on hdd drive, but parallele within sample is memory friendly (default: 1).\n\n";
    
    std::cout << "  -d, --depth <value>          Fixed 31-mer depth value (incompatible with -b/--background). IMPORTANT: 31-mer depth = (1 - 30/read_length) × sequencing_depth. For 150 bp reads, this = 0.8 × sequencing_depth.\n";
    std::cout << "  -D, --Depth <file>           File of depth values corresponding to each input (incompatible with -b/--background).\n";
    std::cout << "  -b, --background <file>      Background k-mer file to estimate NGS coverage (incompatible with -d/-D). In target runs, randomly generated 1M regions are used for coverage estimation. Defaults: <matrix>.bgd \n";

    std::cout << "  -T <ref fasta path>          Reference FASTA for reading CRAM files (default: use REF_CACHE and REF_PATH environment variables).\n";
    std::cout << "  -w, --window <number>        Window size for k-mer coverage report (default: 30).\n";
    std::cout << "  -c, --corr <0/1>             Enable NGS k-mer bias correction (default: 1).\n";

    std::cout << "\nTarget run options (for using reads aligned to specific regions or genes, will run all genes in database without specifying):\n\n";
    std::cout << "  -g, --gene <name>            Target gene name, prefix (ending with '*', remember to quote escape, e.g., 'HLA*' ), or matrix (starting with '#', e.g., #SMN_group1). Can be specified multiple times.\n";
    std::cout << "  -G, --Genes <file>           File listing target genes or matrices.\n\n";

    std::cout << "  -B <file>                    BED file to restrict region analysis. Make sure its medium name field matches your reference genome md5 values (md5sum $reference). One profiled on EBI/GRCh38_full_analysis_set_plus_decoy_hla.fa is included in https://github.com/Walfred-MA/Ctyper/tree/main/data. With -g/-G, only BED entries with names found in -g/-G are used.\n";
    std::cout << "  -r  <chr:start-end>          Add a specific region for analysis. Can be specified multiple times. Regions with be merged if multiple regions provided, can work with -B.\n";
    std::cout << "  -r  'gene'                   A special key for -r, add regions from matrix database (must be combined with -g/-G; incompatible with -BED). Less accurate than profiled regions. Not recommend if profile file avaiable or running global mode is possible. No need to add if using a frofiling bedfile \n\n";
    std::cout << "  -r  'Unmap'                  A special key for -r, to include all unmap reads. No need to add if using a frofiling bedfile.\n\n";
    std::cout << "  -r  'HLA'                    A special key for -r, to include reads on all HLA decoys. No need to add if using a frofiling bedfile. \n\n";
    std::cout << "  --unmap <0/1>                Force include:1 or exclude:0 unmapped reads (only valid in target run). More like a debug function and in most cases, no need to specified. No need to add if using a frofiling bedfile \n\n";

    std::cout << "Target run profiling options (locate mismapped reads and generate BED files for future use):\n";
    std::cout << "  -p, --profile <file>         Input aligned NGS file for profiling.\n";
    std::cout << "  -P, --Profile <file>         File listing multiple aligned NGS files for profiling. Can be used with both -O and -o; individual results go to -O paths, and a summary is saved to -o.\n\n\n";
    
    std::cout << "Example target mode genotyping command using 16*4 cores on 16 samples to genotype all HLAs:\n  ctyper -B TargetRegions.64b32de2fc934679c16e83a2bc072064.bed -m HprcCpcHgsvc_cmr_matrix.txt -I ./ALL_Input_Paths -O ./ALL_Output_Paths -g 'HLA*' -n 16 -N 4 -T GRCh38_full_analysis_set_plus_decoy_hla.fa 2> log.txt \n";
    std::cout << std::endl;
}



int main(int argc, char * argv_[])
{
    
    std::vector<std::string> inputfiles;
        
    std::vector<std::string> outputfiles;
    
    std::vector<float> depths;
    
    std::unordered_set<std::string> genes;
    
    std::vector<char *> regions;
    
    std::string kmatrixfile="";
    std::string backgroundfile="";
    std::string bedfile = "";
    
    
    bool checkp1 = 0, checkp2= 0;
    const char* Argument;
    int nthreads = 1;
    int Nsubthreads = 1;
    int window = 30;
    bool cramnoteprint = 0;
    
    if ( argc > 1 && (strcmp(argv_[1], "-h")==0 or strcmp(argv_[1], "--help")==0) )
    {
        printHelp();
        return 0;
    }
    
    char** argv = &argv_[0];
    if (argc == 1)
    {
        printHelp();
        argv = nullptr;
        Interactive interact(argc, argv);
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
            if (cramnoteprint == 0 && inputfiles.back().size() >= 5 && inputfiles.back().compare(inputfiles.back().size() - 5, 5, ".cram") == 0)
            {
                cramnoteprint = 1;
            }
        }
        else if (strcmp(Argument, "-I")==0 or strcmp(Argument, "--Inputs")==0)
        {
            std::ifstream pathfile(argv[i]);

            if(!pathfile)
            {
                std::cerr<<"Error opening input file"<<std::endl;
                return -1;
            }
            std::string line;
            while (std::getline(pathfile, line))
            {
                inputfiles.push_back(line);
                if (cramnoteprint == 0 && inputfiles.back().size() >= 5 && inputfiles.back().compare(inputfiles.back().size() - 5, 5, ".cram") == 0)
                {
                    cramnoteprint = 1;
                }
            }
        }
        
        else if (strcmp(Argument, "-m")==0 or strcmp(Argument, "--matrix")==0)
        {
            kmatrixfile = string(argv[i]);
            if (backgroundfile.size() == 0 && inputdep == 0 && ifprofile == 0 && std::filesystem::exists(kmatrixfile + ".bgd"))
            {
                backgroundfile = kmatrixfile + ".bgd";
            }
        }
        
        else if (strcmp(Argument, "-o")==0 or strcmp(Argument, "--output")==0)
        {
            outputfiles.push_back(argv[i]);
            if (checkp2) profilesummaryfile = argv[i];
            checkp1 = 1;
        }
        
        else if (strcmp(Argument, "-p")==0 or strcmp(Argument, "--profile")==0)
        {
            ifprofile = 1;
            backgroundfile = "";
            inputfiles.push_back(argv[i]);
        }
        
        else if (strcmp(Argument, "-P")==0 or strcmp(Argument, "--Profile")==0)
        {
            ifprofile = 1;
            backgroundfile = "";
            std::ifstream pathfile(argv[i]);

            if(!pathfile)
            {
                std::cerr<<"Error opening input file"<<std::endl;
                return -1;
            }
            std::string line;
            while (std::getline(pathfile, line))
            {
                inputfiles.push_back(line);
            }
            if (checkp1) profilesummaryfile = argv[i];
            checkp2 = 1;
        }

        else if (strcmp(Argument, "-O")==0 or strcmp(Argument, "--Outputs")==0)
        {
            std::ifstream pathfile(argv[i]);

            if(!pathfile)
            {
                std::cerr<<"Error opening output file"<<std::endl;
                return -1;
            }
            std::string line;
            while (std::getline(pathfile, line))
            {
                outputfiles.push_back(line);
            }
        }
        else if (strcmp(Argument, "--unmap") ==0 )
        {
            forceunmap = 2 - (bool)atoi(argv[i]);
        }
        
        else if (strcmp(Argument, "-r") ==0 )
        {
            
            if (strcmp(argv[i], "gene") == 0)
            {
                ifgeneregion = 1;
                continue;
            }
            
            char* copied = new char[strlen(argv[i]) + 1];
            std::strcpy(copied, argv[i]);
            regions.push_back(copied);
        }
        
        else if (strcmp(Argument, "-B") == 0)
        {
            bedfile = argv[i];
            std::ifstream bedfile_(argv[i]);

            if(!bedfile_)
            {
                std::cerr<<"Error opening bed file"<<std::endl;
                return -1;
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
                std::cerr<<"Error opening output file"<<std::endl;
                return -1;
            }
            std::string line;
            while (std::getline(pathfile, line))
            {
                genes.insert(line);
            }
        }
        
        else if (strcmp(Argument, "-T")==0 or strcmp(Argument, "--ref")==0)
        {
            refpath = argv[i];
        }
        
        else if (strcmp(Argument, "-n")==0 or strcmp(Argument, "--nthreads")==0)
        {
            nthreads=(int)atoi(argv[i]);
        }
        
        else if (strcmp(Argument, "-N")==0 or strcmp(Argument, "--subthreads")==0)
        {
            Nsubthreads=(int)atoi(argv[i]);
        }
        
        else if (strcmp(Argument, "-d")==0 or strcmp(Argument, "--depth")==0)
        {
            inputdep = 1;
            backgroundfile = "";
            depths.push_back(atof(argv[i]));
        }
        
        else if (strcmp(Argument, "-D")==0 or strcmp(Argument, "--Depth")==0)
        {
            inputdep = 1;
            backgroundfile = "";
            std::ifstream pathfile(argv[i]);

            if(!pathfile)
            {
                std::cerr<<"Error opening output file"<<std::endl;
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
    
    std::ifstream versionfile(kmatrixfile + ".index");
    if (versionfile)
    {
        std::string versionline;
        if (std::getline(versionfile, versionline) && versionline[0] == '@')
        {
            vdatabase = versionline.substr(1);

            if (vdatabase.find("support:"+string(vctyper)) == std::string::npos &&
                vdatabase.find(","+string(vctyper)) == std::string::npos)
            {
                std::cerr << "Error: Database version mismatch (expected "+string(vctyper)+" or support:"+string(vctyper)+")." << std::endl;
                return 1;
            }
        }
        else
        {
            std::cerr << "Warning: Missing database version number in index file. Error might happen if database mismatches." << std::endl;
        }
    }
    else
    {
        std::cerr << "Warning: Cannot open matrix index file to check version. Setting as dry-run to determine sequencing depth" << std::endl;
    }
    
    if (backgroundfile.size() == 0 && inputdep == 0)
    {
        std::cerr << "Warnning: providing backgroud kmers as -b, background kmers can be download at ctyper's github repository inside data folder or provide yourself as -d/D" << std::endl;
        
        return 1;
    }
    
    printHeader();
    
    if (cramnoteprint) cout << "\nNote: Using refernce for CRAM file: " << GET_CRAM_REF_PATH() << endl<< endl;
    
    if (!inputfiles.size()) return 1;
    
    auto begin = std::chrono::high_resolution_clock::now();
    
    Processor<32> *processor = new Processor<32>(inputfiles, outputfiles, depths, kmatrixfile, backgroundfile, bedfile,genes, regions, window, nthreads, Nsubthreads);

    processor->Run();
    
    auto end = std::chrono::high_resolution_clock::now();
    
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        
    cerr<<"finished running at time: "<<elapsed.count()* 1e-9 <<endl;
    
    
    return 0;
}
