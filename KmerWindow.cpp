//
//  KmerWindow.cpp
//  CTyper
//
//  Created by Wangfei MA on 6/30/23.
//

#include "KmerWindow.hpp"


void KmerWindow::WindowCovers(const uint16* kmervec, const uint16* kmermatrix, const FLOAT_T depth, const uint16 gnum, const uint knum, const int genenum, const int* results, ull &total_obs, ull &total_exp)
{
    
    const uint16* rowdata = kmermatrix;
    
    int copynum = 0;
    for (int i = 0; i < genenum; ++i)
    {
        if (results[i] > 0) copynum += results[i];
    }
    
    uint loc;
    int genecounter = 0;
    for (size_t i = 0; i < knum; ++i)
    {
        
        loc= ( rowdata[3] << 16 ) + rowdata[4];
        loc /= window;

        auto& thiswindow = windowcovers[rowdata[2]][loc];

        get<0>(thiswindow) += (int)kmervec[i];
        
        switch (rowdata[0])
        {
            case '_': case '=':
                
                get<1>(thiswindow) += genecounter;
                get<2>(thiswindow) = MAX(get<2>(thiswindow), genecounter);
                break;

            case '-':
                genecounter = copynum;
                for (int i = 0 ; i < rowdata[1] ; i++)
                {
                    genecounter -= results[rowdata[ FIXCOL + i]];
                }
                
                get<1>(thiswindow) += genecounter;
                
                get<2>(thiswindow) = MAX(get<2>(thiswindow), genecounter);
                
                break;
            case '+':
                genecounter = 0;
                for (int i = 0 ; i < rowdata[1] ; i++)
                {
                    genecounter += results[rowdata[ FIXCOL + i]];
                }
                
                get<1>(thiswindow) += genecounter;
                get<2>(thiswindow) = MAX(get<2>(thiswindow), genecounter);
                break;
            default:
                break;
        }
        
        total_exp += genecounter;
        total_obs += kmervec[i];
        
        rowdata = &rowdata[rowdata[1] + FIXCOL];
    }
    
}
