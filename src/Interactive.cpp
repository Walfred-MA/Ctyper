//
//  Interactive.cpp
//  Ctyper1
//
//  Created by walfred on 7/14/25.
//

#include "Interactive.hpp"


// Helper function to check file existence
bool file_exists(const std::string& filename)
{
    return std::filesystem::exists(filename);
}




bool InputEither(std::vector<std::string> &args,
                 const std::string &arg1, const std::string &print1,
                 const std::string &arg2, const std::string &print2)
{
    std::string input;

    // First field
    std::cout << print1;
    std::getline(std::cin, input);
    if (!input.empty())
    {
        // Must keep trying until a valid file is provided
        while (!file_exists(input))
        {
            std::cout << "File does not exist. Please try again:" << endl;
            std::getline(std::cin, input);
        }
        args.push_back(arg1); args.push_back(input);
        return true;
    }

    // Second field (if first skipped)
    std::cout << print2;
    std::getline(std::cin, input);
    if (!input.empty()) {
        while (!file_exists(input))
        {
            std::cout << "File does not exist. Please try again: " << endl;
            std::getline(std::cin, input);
        }
        args.push_back(arg2); args.push_back(input);
        return true;
    }
    
    return false;
    // Both skipped
    // (If you want to force at least one, add a loop here. If skipping is OK, then just return.)
}


// Main interactive function
void Interactive::RUN(int& argc, char**& argv)
{
    std::vector<std::string> args = {"Ctyper"};
    std::string matrix;
    std::string input;
    std::string output;
    std::string depthfile;
    std::string gene_file, bedfile, regionfile;
    std::string val;

    // 1. Ask for database matrix, require .index file as well
    while (true)
    {
        std::cout << "Path to matrix database (-m): \n";
        std::getline(std::cin, matrix);
        if (file_exists(matrix) && file_exists(matrix + ".index")) break;
        std::cout << "Matrix or matrix.index does not exist. Please try again.\n";
    }
    args.push_back("-m"); args.push_back(matrix);

    // check if .bgd exists
    bool bgd_exists = file_exists(matrix + ".bgd");

    // 2. Ask user for input (-i or -I)
    
    string argsi = "-i";
    string Inputi = "Path to input NGS file (-i): ,empty to switch use -I \n";
    string argsI = "-I";
    string InputI = "Not using -i, please input a file with list of input files (-I): \n";
    
    while (InputEither(args,argsi,Inputi,argsI,InputI) == false){};
    
    string argso = "-o";
    string Inputo = "Path to appended output file (-o): , empty to switch use -O \n";
    string argsO = "-O";
    string InputO = "Not using -o, please input a file with list of output files (-O), corrspond to each of input files: \n";
    
    while (InputEither(args,argso,Inputo,argsO,InputO) == false){};

    if (!bgd_exists)
    {
        std::cout << "No background file (" << matrix << ".bgd) found.\n";
        std::cout << "You MUST provide -d (fixed depth) or -D (file of depths).\n";
        while (true) {
            std::cout << "Fixed depth (-d) or file of depths (-D) corrspond to each of input files \n";
            std::getline(std::cin, depthfile);
            if (depthfile.empty()) continue;
            if (depthfile.find_first_not_of("0123456789.") == std::string::npos)
            {
                args.push_back("-d"); args.push_back(depthfile);
                break;
            }
            else if (file_exists(depthfile))
            {
                args.push_back("-D"); args.push_back(depthfile);
                break;
            }
            std::cout << "File does not exist or not a valid number. Please try again.\n";
        }
    }

    // 5. Threads
    std::cout << "Number of samples to run in parallel (-n), default=1: \n";
    std::getline(std::cin, val);
    if (!val.empty()) { args.push_back("-n"); args.push_back(val); }
    std::cout << "Number of threads per sample (-N), default=1: \n";
    std::getline(std::cin, val);
    if (!val.empty()) { args.push_back("-N"); args.push_back(val); }

    // 6. Ask for mode
    std::string mode;
    while (true) {
        std::cout << "Target mode (Y/N): \n";
        std::getline(std::cin, mode);
        if (mode == "y" || mode == "Y" || mode == "n" || mode == "N") break;
        std::cout << "Type 'y' or 'n'.\n";
    }

    // 7. Target mode options
    if (mode == "Y" || mode == "y")
    {
        // -g or -G gene info
        while (true)
        {
            std::cout << "Gene name/prefix/group (-g), or existing file with list of genes (-G): \n";
            std::getline(std::cin, gene_file);
            if (gene_file.empty()) break;
            if (file_exists(gene_file)) { args.push_back("-G"); args.push_back(gene_file); break; }
            else { args.push_back("-g"); args.push_back(gene_file); break; }
        }
        // BED file
        std::cout << "Profiled BED file (-B), leave blank to skip: \n";
        std::getline(std::cin, bedfile);
        if (!bedfile.empty())
        {
            while (!file_exists(bedfile))
            {
                std::cout << "BED file does not exist. Please try again: \n";
                std::getline(std::cin, bedfile);
                if (bedfile.empty()) break;
            }
            if (!bedfile.empty()) { args.push_back("-B"); args.push_back(bedfile); }
        }
        // Regions
        std::cout << "Input region (-r), leave blank to skip: \n";
        std::getline(std::cin, regionfile);
        if (!regionfile.empty()) { args.push_back("-r"); args.push_back(regionfile); }
        // Warn if both BED and -r are missing
        if (bedfile.empty() && regionfile.empty())
        {
            args.push_back("-r"); args.push_back("gene");
            std::cout << "Warning: Both -BED and -r are skipped. Ctyper will use leftover coordinates in database (not recommended).\n";
        }
    }

    // Now build argc/argv
    argc = args.size();
    argv = new char*[argc];
    for (int i = 0; i < argc; ++i) {
        argv[i] = new char[args[i].size() + 1];
        strcpy(argv[i], args[i].c_str());
    }
}
