//
//  BedReader.hpp
//  Ctyper1
//
//  Created by walfred on 7/1/25.
//

#ifndef BedReader_hpp
#define BedReader_hpp

#include <stdio.h>
#include <iostream>  // Required for std::cerr, std::cout, etc.
#include <fstream>   // Required for std::ofstream, std::ifstream
#include <unordered_map>
#include <vector>
#include <string>
#include <unordered_set>
#include <sstream>  // for std::istringstream
#include <algorithm> // for std::sort
#include "config.hpp"
#include <sstream>
#include <zlib.h>

namespace std {
    template <>
    struct hash<pair<string, string>> {
        size_t operator()(const pair<string, string>& p) const {
            return hash<string>()(p.first) ^ hash<string>()(p.second);
        }
    };
}

using namespace std;
class BedReader
{
public:
    BedReader() {};
    
    void mergeHotspotsForKey(std::vector<hotspot>& hotspots, const int mergedis = 1);
    void MergeFiles(const vector<string> & inputfiles, const string &outputfile);
    void ReaderFile(const string & file, const std::unordered_set<std::string>& Genes, vector<char*> &regions);
    void MergeRegions(vector<char*> &regions);
};


#endif /* BedReader_hpp */
