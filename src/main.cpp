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
#include <algorithm>
#include <cstring>

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

std::string matrix_version = "";
MatrixLayout matrix_layout = MatrixLayout::UNKNOWN;

extern string refpath;
string refpath = "";

extern string profilesummaryfile;
string profilesummaryfile = "";

extern int inputdep;
int inputdep = 0;

static vector<string> split_matrix_row(const string& row)
{
    vector<string> fields;
    const size_t length = strlen(row.c_str());
    size_t start = 0;
    while (start < length)
    {
        size_t end = row.find('\t', start);
        if (end == string::npos || end > length) end = length;
        fields.push_back(row.substr(start, end - start));
        start = end + 1;
    }
    return fields;
}

static MatrixLayout detect_matrix_layout(const string& matrix_path)
{
    KtableReader reader(matrix_path.c_str());
    string row(MAX_LINE, '\0');
    bool all_v1_fixed = true;
    bool all_v2_fixed = true;
    bool saw_extended_delimiters = false;
    size_t rows_checked = 0;

    while (rows_checked < 32 && reader.nextLine_kmer(row))
    {
        const vector<string> fields = split_matrix_row(row);
        if (fields.size() < 5)
        {
            reader.Close();
            return MatrixLayout::UNKNOWN;
        }

        all_v2_fixed = all_v2_fixed &&
            fields[1].size() == MATRIX_ENCODING_V2.tag_width &&
            fields[2].size() ==
                MATRIX_ENCODING_V2.path_width + MATRIX_ENCODING_V2.query_index_width &&
            fields[3].size() == MATRIX_ENCODING_V2.size_width +
                MATRIX_ENCODING_V2.query_position_width +
                MATRIX_ENCODING_V2.reference_position_width &&
            fields[4].size() == MATRIX_ENCODING_V2.kmer_width;

        all_v1_fixed = all_v1_fixed &&
            fields[1].size() == MATRIX_ENCODING_V1.tag_width &&
            fields[2].size() == MATRIX_ENCODING_V1.path_width &&
            fields[3].size() == MATRIX_ENCODING_V1.reference_position_width &&
            fields[4].size() == MATRIX_ENCODING_V1.kmer_width;

        saw_extended_delimiters = saw_extended_delimiters ||
            fields[2].find_first_of("+-") != string::npos ||
            std::count(fields[3].begin(), fields[3].end(), '|') == 2;
        ++rows_checked;
    }
    reader.Close();

    if (rows_checked == 0) return MatrixLayout::UNKNOWN;
    if (all_v2_fixed) return MatrixLayout::V2_FIXED;
    if (all_v1_fixed) return MatrixLayout::V1_FIXED;
    if (saw_extended_delimiters) return MatrixLayout::V2_DELIMITED;
    return MatrixLayout::V1_DELIMITED;
}

static const char* matrix_layout_name(MatrixLayout layout)
{
    switch (layout)
    {
        case MatrixLayout::V1_FIXED: return "v1.0.0 fixed-width";
        case MatrixLayout::V1_DELIMITED: return "v1.0.0 legacy delimited";
        case MatrixLayout::V2_FIXED: return "v2.0.1 fixed-width";
        case MatrixLayout::V2_DELIMITED: return "legacy extended delimited";
        default: return "unknown";
    }
}


                                
void printHeader()
{
    std::cerr << R"(
                            _____  _______  __     __  ______   ________  ______  
                           / ____| |__ __|  \ \   / /  |  __ \  |  ____|  |  __ \ 
                           | |       | |     \ \_/ /   | |__) | | |__     | |__) |
                           | |       | |      \   /    |  ___/  |  __|    |  _  / 
                           | |____   | |       | |     | |      | |____   | | \ \   
                           \_____|   |_|       |_|     |_|      |______|  |_|  \_\  
            )" << std::endl;
    
    std::cerr << "Ctyper " + string(vctyper) + "\nPangenome allele database version: " + vdatabase + "\n"
        << "Author: Walfred (Wangfei) MA (walfredma0@gmail.com), Mark Chaisson Lab (mchaisso@usc.edu). Please contact me if you find any bug or usage question.\n"
        << "Please cite our work if you use this tool.\n"
        << "Source code: https://github.com/Walfred-MA/Ctyper\n"
        << "Using pangenome annotation, results can be converted into readable tables, FASTA files, VCF files, HLA/KIR/CYP2D6 nomenclatures, and mutation plots. "
        << "Tool scripts: https://github.com/Walfred-MA/Ctyper/tree/main/tools\n"
        << "Databases (genotyping matrices and pangenome allele annotations): https://zenodo.org/records/16340156\n\n";
}


void printHelp()
{
    printHeader();

    std::cerr << "Usage: Ctyper " + string(vctyper) + " [options]\n\n";
    std::cerr << "Required:\n";
    std::cerr << "  -m, --matrix <file>          Path to the matrix database (requires <file>.index). If not provided, the tool runs in dry-run mode to only estimate NGS read depth.\n\n";

    std::cerr << "Options:\n";
    std::cerr << "  -i, --input <value>          Input NGS file; supports CRAM, BAM, SAM, FASTA, FASTQ, or Jellyfish formats.\n";
    std::cerr << "  -I, --Inputs <file>          Path to a file listing multiple input files.\n\n";

    std::cerr << "  -o, --output <value>         Output file (append if file exits, default: stdout)\n";
    std::cerr << "  -O, --Outputs <file>         Path to a file listing output files corresponding to each input file.\n\n";

    std::cerr << "  -n, --nthreads <number>      Number of threads to run different samples in parallel, parallel between sample is more efficient for large cohorts, especially when reading files on slower hard drives (default: 1).\n";
    std::cerr << "  -N, --subthreads <number>    Number of threads to run each sample, due to file I/O bottleneck, suggest use only 1-4 on hdd drive, but parallele within sample is memory friendly (default: 1).\n\n";
    
    std::cerr << "  -d, --depth <value>          Fixed 31-mer depth value (incompatible with -b/--background). IMPORTANT: 31-mer depth = (1 - 30/read_length) × sequencing_depth. For 150 bp reads, this = 0.4 × sequencing_depth.\n";
    std::cerr << "  -D, --Depth <file>           File of depth values corresponding to each input (incompatible with -b/--background).\n";
    std::cerr << "  -b, --background <file>      Background k-mer file to estimate NGS coverage (incompatible with -d/-D). In target runs, randomly generated 1M regions are used for coverage estimation. Defaults: <matrix>.bgd \n";

    std::cerr << "  -T <ref fasta path>          Reference FASTA for reading CRAM files (default: use REF_CACHE and REF_PATH environment variables).\n";
    std::cerr << "  -w, --window <number>        Window size for k-mer coverage report (default: 30).\n";
    std::cerr << "  -c, --corr <0/1>             Enable NGS k-mer bias correction (default: 1).\n";

    std::cerr << "\nTarget run options (for using reads aligned to specific regions or genes, will run all genes in database without specifying):\n\n";
    std::cerr << "  -g, --gene <name>            Target gene name, prefix (ending with '*', remember to quote escape, e.g., 'HLA*' ), or matrix (starting with '#', e.g., #SMN_group1). Can be specified multiple times.\n";
    std::cerr << "  -G, --Genes <file>           File listing target genes or matrices.\n\n";

    std::cerr << "  -B <file>                    BED file to restrict region analysis. Make sure its medium name field matches your reference genome md5 values (md5sum $reference). One made from profiling run on EBI/GRCh38_full_analysis_set_plus_decoy_hla.fa is included in https://github.com/Walfred-MA/Ctyper/tree/main/data. With -g/-G, only BED entries with names found in -g/-G are used.\n";
    std::cerr << "  -r  <chr:start-end>          Add a specific region for analysis. Can be specified multiple times. Regions with be merged if multiple regions provided, can work with -B.\n";
    std::cerr << "  -r  'gene'                   A special key for -r, add regions from matrix database (must be combined with -g/-G; incompatible with -BED). Less accurate than profiled regions. Not recommend if profile file avaiable or running global mode is possible. No need to add if using a frofiling bedfile \n\n";
    std::cerr << "  -r  'Unmap'                  A special key for -r, to include all unmap reads. No need to add if using a frofiling bedfile.\n\n";
    std::cerr << "  -r  'HLA'                    A special key for -r, to include reads on all HLA decoys. No need to add if using a frofiling bedfile. \n\n";
    std::cerr << "  --unmap <0/1>                Force include:1 or exclude:0 unmapped reads (only valid in target run). More like a debug function and in most cases, no need to specified. No need to add if using a frofiling bedfile \n\n";

    std::cerr << "Target run profiling options (locate mismapped reads and generate BED files for future use):\n";
    std::cerr << "  -p, --profile <file>         Input aligned NGS file for profiling.\n";
    std::cerr << "  -P, --Profile <file>         File listing multiple aligned NGS files for profiling. Can be used with both -O and -o; individual results go to -O paths, and a summary is saved to -o.\n\n\n";
    
    std::cerr << "Example target mode genotyping command using 16*4 cores on 16 samples to genotype all HLAs:\n  ctyper -B TargetRegions.64b32de2fc934679c16e83a2bc072064.bed -m HprcCpcHgsvc_cmr_matrix.txt -I ./ALL_Input_Paths -O ./ALL_Output_Paths -g 'HLA*' -n 16 -N 4 -T GRCh38_full_analysis_set_plus_decoy_hla.fa 2> log.txt \n";
    std::cerr  << std::endl;
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
        if (std::getline(versionfile, versionline) &&
            !versionline.empty() && versionline[0] == '@')
        {
            vdatabase = versionline.substr(1);

            if (vdatabase.find("support:"+string(vctyper)) == std::string::npos &&
                vdatabase.find(","+string(vctyper)) == std::string::npos)
            {
                std::cerr << "Error: Database version mismatch (expected "+string(vctyper)+" or support:"+string(vctyper)+")." << std::endl;
                return 1;
            }

            const size_t version_end = vdatabase.find_first_of(",\t ");
            const string advertised_version = vdatabase.substr(0, version_end);
            if (advertised_version == MATRIX_VERSION_V1)
            {
                matrix_version = MATRIX_VERSION_V1;
                matrix_layout = MatrixLayout::V1_FIXED;
            }
            else if (advertised_version == MATRIX_VERSION_V2)
            {
                matrix_version = MATRIX_VERSION_V2;
                matrix_layout = MatrixLayout::V2_FIXED;
            }

            if (!kmatrixfile.empty() && std::filesystem::exists(kmatrixfile))
            {
                const MatrixLayout detected_layout = detect_matrix_layout(kmatrixfile);
                if (matrix_layout != MatrixLayout::UNKNOWN &&
                    detected_layout != MatrixLayout::UNKNOWN &&
                    detected_layout != matrix_layout)
                {
                    std::cerr << "Error: Matrix index advertises " << advertised_version
                              << ", but the top k-mer rows have "
                              << matrix_layout_name(detected_layout) << " encoding."
                              << std::endl;
                    return 1;
                }
                if (matrix_layout == MatrixLayout::UNKNOWN)
                {
                    matrix_layout = detected_layout;
                    if (detected_layout == MatrixLayout::V2_FIXED)
                        matrix_version = MATRIX_VERSION_V2;
                    else if (detected_layout == MatrixLayout::V1_FIXED ||
                             detected_layout == MatrixLayout::V1_DELIMITED)
                        matrix_version = MATRIX_VERSION_V1;
                    std::cerr << "Warning: Unrecognized advertised matrix version '"
                              << advertised_version << "'; using the auto-detected "
                              << matrix_layout_name(matrix_layout) << " layout."
                              << std::endl;
                }
            }
        }
        else
        {
            if (!kmatrixfile.empty() && std::filesystem::exists(kmatrixfile))
                matrix_layout = detect_matrix_layout(kmatrixfile);

            if (matrix_layout == MatrixLayout::V2_FIXED)
                matrix_version = MATRIX_VERSION_V2;
            else if (matrix_layout == MatrixLayout::V1_FIXED ||
                     matrix_layout == MatrixLayout::V1_DELIMITED)
                matrix_version = MATRIX_VERSION_V1;

            if (matrix_layout == MatrixLayout::UNKNOWN)
            {
                std::cerr << "Warning: Missing matrix version in index and unable "
                          << "to recognize the k-mer row layout." << std::endl;
            }
            else
            {
                vdatabase = matrix_version.empty()
                    ? "legacy-extended,auto-detected"
                    : matrix_version + ",auto-detected";
                std::cerr << "Warning: Missing matrix version in index; auto-detected "
                          << matrix_layout_name(matrix_layout) << "." << std::endl;
            }
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
    
    if (cramnoteprint) cerr << "\nNote: Using refernce for CRAM file: " << GET_CRAM_REF_PATH() << endl<< endl;
    
    if (!inputfiles.size()) return 1;
    if (!outputfiles.size() && nthreads == 1) outputfiles.push_back("stdout");
    
    auto begin = std::chrono::high_resolution_clock::now();
    
    Processor<32> *processor = new Processor<32>(inputfiles, outputfiles, depths, kmatrixfile, backgroundfile, bedfile,genes, regions, window, nthreads, Nsubthreads);

    processor->Run();
    
    auto end = std::chrono::high_resolution_clock::now();
    
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        
    cerr<<"finished running at time: "<<elapsed.count()* 1e-9 <<endl;
    
    
    return 0;
}
