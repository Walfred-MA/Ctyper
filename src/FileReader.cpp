//
//  FileReader.cpp
//  Ctyper2
//
//  Created by walfred on 6/24/25.
//

#include "FileReader.hpp"


void FileReader::Load()
{
    if (strlen(filepath)<2) return;
    
    fafile.open(filepath, ios_base::in);
    
    if (!fafile) {
        
        std::cerr << "ERROR: Could not open " << filepath << " for reading.\n" << std::endl;
        std::_Exit(EXIT_FAILURE);
        
        return;
    }
    
}

bool FileReader::nextLine(std::string &StrLine)
{
    return (bool)( getline(fafile,StrLine));
}

void FileReader::Close()
{
    
    fafile.close();
    
}

void FileReader::Reset()
{
    fafile.seekg(0, std::ios::beg);
}

int FileReader::GetLine(FILE* fptr, std::string &line)
{
    char buf[1024];
    line = "";
    
    while (true)
    {
        if (fgets(buf, sizeof(buf)-1, fptr) == nullptr)
        {
            // fgets failed, check if it's an error or end-of-file
            if (feof(fptr))
            {
                // End-of-file reached
                break;
            }
            else
            {
                // I/O error occurred
                perror("Error reading from file");
                return -1;
            }
        }
        
        line += buf;
        
        if (line.back() == '\n')
        {
            // Remove the newline character from the end of the line
            line.pop_back();
            break;
        }
        
        // If there are no more characters in the file, break the loop
        if (feof(fptr))
        {
            break;
        }
    }
    
    // Return 1 if the file has more lines, 0 if it's the last line, -1 if an error occurred
    return !feof(fptr);
}
int FileReader::GetLine(BGZF* fptr, std::string &line)
{
    char buf[1024];
    line.clear();
    
    while (true)
    {
        if (!bgzf_fgets(buf, sizeof(buf), fptr))
        {
            // If line is empty, we reached EOF
            return line.empty() ? 0 : 1;
        }
        
        line += buf;
        
        // Check if the line contains the full line ending
        if (!line.empty() && line.back() == '\n')
        {
            line.pop_back();  // remove newline
            return 1;         // successfully read a line
        }
    }
    return false;
}
