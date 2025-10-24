//
//  BedReader.cpp
//  Ctyper1
//
//  Created by walfred on 7/1/25.
//

#include "BedReader.hpp"

constexpr int mergeextendsize = 150;
// Function to merge hotspots for each key (name, chr)
void BedReader::mergeHotspotsForKey(std::vector<hotspot>& hotspots, const int mergedis)
{
    // Step 1: Sort the hotspots by their start position
    sort(hotspots.begin(), hotspots.end(), [](const hotspot& a, const hotspot& b) {
        return a.first < b.first;  // Sort by start position
    });

    // Step 2: Merge overlapping or adjacent hotspots
    std::vector<hotspot> merged;
    
    for (const auto& h : hotspots)
    {
        if (merged.empty())
        {
            // If merged is empty, just add the first hotspot
            merged.push_back(h);
        }
        else
        {
            // Get the last merged hotspot
            hotspot& last = merged.back();
            
            // If the current hotspot starts before the last one ends or is adjacent, merge them
            if (h.first <= last.second + mergedis)
            {
                last.second = MAX(last.second, h.second);  // Extend the end if necessary
            }
            else
            {
                merged.push_back(h);  // Otherwise, start a new range
            }
        }
    }

    // Step 3: Replace the original hotspots with the merged ones
    hotspots = std::move(merged);  // Move merged hotspots into the original std::vector
}

void BedReader::ReaderFile(const string & file, const unordered_set<string>& Genes, vector<char*> &regions)
{
    unordered_map<string, vector<hotspot> > mergedHotspots;
    
    if (file.size() > 3 && file.substr(file.size() - 3) == ".gz")
    {
        // Open gzipped file using zlib
        gzFile infile = gzopen(file.c_str(), "r");
        if (!infile)
        {
            std::cerr << "Error: Could not open gzipped file " << file << std::endl;
            return;
        }

        char buffer[1000];  // Buffer to read lines
        while (gzgets(infile, buffer, sizeof(buffer)) != NULL)
        {
            string line(buffer);
            // Skip empty lines
            if (line.empty()) continue;

            std::istringstream stream(line);
            std::string chr, name;
            int start, end;

            if ((stream >> chr >> start >> end))
            {
                if (!(stream >> name)) {
                    name = "";  // Accept lines with only 3 fields
                }

                // If Genes is provided, filter only when name is non-empty
                if (!name.empty() && !Genes.empty() && Genes.find(name) == Genes.end())
                    continue;

                mergedHotspots[chr].emplace_back(start, end);
            }
            else
            {
                std::cerr << "Error: Invalid line format in file " << file << ": " << line << std::endl;
            }
        }
        gzclose(infile);
    }
    else
    {
        // Open as a regular text file
        std::ifstream infile(file);
        if (!infile.is_open())
        {
            std::cerr << "Error: Could not open file " << file << std::endl;
            return;
        }

        std::string line;
        while (std::getline(infile, line))
        {
            // Skip empty lines
            if (line.empty()) continue;

            std::istringstream stream(line);
            std::string chr, name;
            int start, end;

            if ((stream >> chr >> start >> end))
            {
                if (!(stream >> name)) {
                    name = "";  // Accept lines with only 3 fields
                }

                // If Genes is provided, filter only when name is non-empty
                if (!name.empty() && !Genes.empty() && Genes.find(name) == Genes.end())
                    continue;

                mergedHotspots[chr].emplace_back(start, end);
            }
            else
            {
                std::cerr << "Error: Invalid line format in file " << file << ": " << line << std::endl;
            }
        }
        
        infile.close();
    }
    
    
    for (auto& entry : mergedHotspots)
    {
        mergeHotspotsForKey(entry.second);  // Merge the intervals for this key
        for (const auto& m : entry.second)
        {
            std::string region_str = entry.first + ":" +  to_string(m.first)  + "-" + to_string(m.second);
            
            char* region_cstr = (char*) malloc((region_str.size() + 1) * sizeof(char));
            std::strcpy(region_cstr, region_str.c_str());
            regions.push_back(region_cstr);
        }
    }
    
    
    
}

void BedReader::MergeRegions(vector<char*>& regions)
{
    unordered_map<string, vector<hotspot>> mergedHotspots;
    vector<string> passthroughRegions;
    // Step 1: Load regions into mergedHotspots from the char* regions vector
    for (char* region_cstr : regions)
    {
        std::string region_str(region_cstr);  // Convert char* to std::string
        free(region_cstr);  // Free the memory since we no longer need the original char*

        // Parse the region string (format: chr:start-end)
        size_t colon_pos = region_str.find_last_of(":");
        size_t dash_pos = region_str.find_last_of("-");
        
        if (colon_pos == std::string::npos || dash_pos == std::string::npos || colon_pos >= dash_pos)
        {
            passthroughRegions.push_back(region_str);  // Store for final output
            continue;
        }
        
        std::string chr = region_str.substr(0, colon_pos);  // Extract chr
        int start = std::stoi(region_str.substr(colon_pos + 1, dash_pos - colon_pos - 1));  // Extract start
        int end = std::stoi(region_str.substr(dash_pos + 1));  // Extract end
        start += 1;
        // Insert {start, end} into the vector corresponding to this chr
        mergedHotspots[chr].emplace_back(start, end);
    }

    regions.clear();
    // Step 2: Merge the hotspots for each chromosome
    for (auto& entry : mergedHotspots)
    {
        mergeHotspotsForKey(entry.second,mergeextendsize);  // Merge the intervals for this chr

        // Step 3: Rebuild the merged regions and store them back in the regions vector
        for (const auto& m : entry.second)
        {
            std::string region_str = entry.first + ":" + std::to_string(m.first) + "-" + std::to_string(m.second);
            char* region_cstr = (char*) malloc((region_str.size() + 1) * sizeof(char));  // Allocate memory for C-string
            std::strcpy(region_cstr, region_str.c_str());  // Copy the string into the allocated space
            regions.push_back(region_cstr);  // Add the region to the regions vector
        }
    }
    for (const std::string& region_str : passthroughRegions)
    {
        char* region_cstr = (char*) malloc(region_str.size() + 1);
        std::strcpy(region_cstr, region_str.c_str());
        regions.push_back(region_cstr);
    }
}

void BedReader::MergeFiles(const vector<string> & inputfiles, const string &outputfile)
{
    
    unordered_map<pair<string,string>, vector<hotspot> > mergedHotspots;
    
    for (const auto& file : inputfiles)
    {
        std::ifstream infile(file);
        if (!infile.is_open())
        {
            std::cerr << "Error: Could not open file " << file << std::endl;
            continue;
        }

        std::string line;
        while (std::getline(infile, line))
        {
            // Skip empty lines
            if (line.empty()) continue;

            // Parse each line into chr, start, end, and name
            std::istringstream stream(line);
            std::string chr, name;
            int start, end;

            // Assuming the format: chr start end name
            if (stream >> chr >> start >> end >> name)
            {
                // Create the key as {name, chr}
                pair<string, string> key = {name, chr};

                // Insert {start, end} into the vector corresponding to this key
                mergedHotspots[key].emplace_back(start, end);
            }
            else
            {
                std::cerr << "Error: Invalid line format in file " << file << ": " << line << std::endl;
            }
        }

        infile.close();  // Close the file
    }
    
    ofstream outFile(outputfile);
    if (!outFile.is_open())
    {
        cerr << "Error: Could not open file " << outputfile << " for writing." << endl;
        return;
    }
    for (auto& entry : mergedHotspots)
    {
        mergeHotspotsForKey(entry.second);  // Merge the intervals for this key
        for (const auto& m : entry.second)
        {
            outFile << entry.first.second << '\t'    // chr
                        << m.first << '\t' // start
                        << m.second << '\t' // end
                        << entry.first.first << '\n';    // name
        }
    }
    
    outFile.close();
}


