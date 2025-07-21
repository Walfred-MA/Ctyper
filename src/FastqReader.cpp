//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License. 
//  If you use this code, please cite our work.
//   

#include "FastqReader.hpp"

void FastqReader::Load()
{
    if (strlen(filepath)<2) return;
    
    ifgz = (std::string(filepath).size() >= 3 &&
            std::string(filepath).compare(std::string(filepath).size() - 3, 3, ".gz") == 0);
    if (ifgz)
    {
        fafile = gzopen(filepath, "rb");
        if (fafile == NULL)
        {
            std::cerr << "ERROR: Could not open " << filepath << " for reading.\n" << std::endl;
            std::_Exit(EXIT_FAILURE);
        }
    }
    else
    {
        file = fopen(filepath, "r");
        if (file == NULL)
        {
            std::cerr << "ERROR: Could not open " << filepath << " for reading.\n" << std::endl;
            std::_Exit(EXIT_FAILURE);
        }
    }
    
    
}

bool FastqReader::LoadBuffer()
{
    size_t leftover = 0;

    // Carry over unfinished data
    if (buffer_pos < buffer_remain)
    {
        leftover = buffer_remain - buffer_pos;
        std::memmove(buffer.data(), buffer.data() + buffer_pos, leftover);
    }

    buffer_pos = 0;
    buffer_remain = leftover;

    while (true)
    {
        size_t bytes_to_read = buffer.size() - buffer_remain;
        size_t bytes_read = 0;

        if (ifgz)
        {
            bytes_read = gzread(fafile, buffer.data() + buffer_remain, bytes_to_read);
            if (bytes_read <= 0)
            {
                if (gzeof(fafile))
                {
                    // try next file
                    fileindex++;
                    if (fileindex < filepathes.size())
                    {
                        filepath = filepathes[fileindex].c_str();
                        Close();
                        Load();
                        continue;  // try reading again
                    }
                    return (buffer_remain > 0);  // return true if leftover still to process
                }
                else
                {
                    int errnum;
                    std::cerr << "gzread error: " << gzerror(fafile, &errnum) << std::endl;
                    return false;
                }
            }
        }
        else
        {
            bytes_read = fread(buffer.data() + buffer_remain, 1, bytes_to_read, file);
            if (bytes_read == 0)
            {
                if (feof(file))
                {
                    fileindex++;
                    if (fileindex < filepathes.size())
                    {
                        filepath = filepathes[fileindex].c_str();
                        Close();
                        Load();
                        continue;
                    }
                    return (buffer_remain > 0);  // process leftover
                }
                else if (ferror(file))
                {
                    perror("fread error");
                    return false;
                }
            }
        }

        buffer_remain += bytes_read;
        return true;
    }
}



bool FastqReader::nextLine(char*& strLine, size_t& rlen, char*& header)
{
    std::lock_guard<std::mutex> IO(IO_lock);

    while (true)
    {
        if (buffer_pos >= buffer_remain)
        {
            if (!LoadBuffer()) return false;
        }

        size_t line_start = buffer_pos;
        int line_stage = 0;  // 0 = @header, 1 = sequence, 2 = +, 3 = quality
        rlen = 0;
        char* header_line = nullptr;

        for (size_t i = buffer_pos; i < buffer_remain; ++i)
        {
            if (buffer[i] == '\n')
            {
                switch (line_stage)
                {
                    case 0:
                        if (buffer[line_start] == '@')
                        {
                            header_line = buffer.data() + line_start;
                            line_stage = 1;
                        }
                        break;

                    case 1:
                        strLine = buffer.data() + line_start;
                        rlen = i - line_start;
                        header = header_line;
                        line_stage = 2;
                        break;

                    case 2:
                        line_stage = 3;
                        break;

                    case 3:
                        // End of record — return after reading full 4 lines
                        buffer_pos = i + 1;
                        return true;
                }

                buffer_pos = i + 1;
                line_start = buffer_pos;
            }
        }
        
        auto ifload = LoadBuffer();
        // Need to load more data
        if (!ifload) return false;
    }
}



void FastqReader::Close()
{
    if (fafile)
    {
        gzclose(fafile);
        fafile = nullptr;
    }
    if (file)
    {
        fclose(file);
        file = nullptr;
    }
}
