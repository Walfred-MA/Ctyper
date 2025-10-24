//
//  Interactive.hpp
//  Ctyper1
//
//  Created by walfred on 7/14/25.
//

#ifndef Interactive_hpp
#define Interactive_hpp

#include <stdio.h>
#include <iostream>
#include <string>
#include <unordered_set>
#include <chrono>
#include <filesystem>
#include <sstream>
#include <vector>
#include <cstring>
using namespace std;

class Interactive
{
public:
    Interactive(int& argc, char**& argv)
    {
        RUN(argc, argv);
    }
    ~Interactive(){};
    void RUN(int& argc, char**& argv);
    
};

#endif /* Interactive_hpp */
