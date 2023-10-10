//
//  KmerWindow.hpp
//  CTyper
//
//  Created by Wangfei MA on 6/30/23.
//

#ifndef WindowCover_hpp
#define WindowCover_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <tuple>
#include <mutex>
#include <memory>
#include <numeric>
using namespace std;


#include "config.hpp"


class KmerWindow
{
    
public:
    
    KmerWindow(const int w = 30): window(w)
    {
    };
    ~KmerWindow()
    {};
    
    void resize(const vector<uint>& sizes)
    {
        windowcovers.resize(sizes.size());
        for (int i =0; i < sizes.size(); ++i)
        {
            windowcovers[i].resize(sizes[i]/window + 1);
            
            std::fill(windowcovers[i].begin(), windowcovers[i].end(), make_tuple(0,0,0));
            
        }
    }
    
    void WindowCovers(const uint16* kmervec, const uint16* kmermatrix, const FLOAT_T depth, const uint16 gnum, const uint knum, const int genenum, const int* results, ull &total_obs, ull &total_exp);
    
    vector<vector<tuple<int,int,int>>> windowcovers;
    
private:
        
    const int window;
    
    
    
};

#endif /* KmerWindow_hpp */
