//
//  KmerWindow.cpp
//  CTyper
//
//  Created by Wangfei MA on 6/30/23.
//

#include "KmerWindow.hpp"


void KmerWindow::WindowCovers(const uint16* kmervec, const uint16* kmermatrix, const FLOAT_T depth, const uint16 gnum, const uint knum, const int genenum, const int* results)
{
    
    const uint16* rowdata = kmermatrix;
    
    int copynum = 0;
    for (int i = 0; i < genenum; ++i)
    {
        copynum += results[i];
    }
    
    uint loc;
    int kmernum = 0;
    for (size_t i = 0; i < knum; ++i)
    {
        
        
        loc= ( rowdata[3] << 16 ) + rowdata[4];
        loc /= window;

        auto& thiswindow = windowcovers[rowdata[2]][loc];

        get<0>(thiswindow) += (int)kmervec[i];
        
        switch (rowdata[0])
        {
            case '_': case '=':
                
                get<1>(thiswindow) += kmernum;
                get<2>(thiswindow) = MAX(get<2>(thiswindow), kmernum);
                break;

            case '-':
                kmernum = copynum;
                for (int i = 0 ; i < rowdata[1] ; i++)
                {
                    kmernum -= results[rowdata[ FIXCOL + i]];
                }
                
                get<1>(thiswindow) += kmernum;
                get<2>(thiswindow) = MAX(get<2>(thiswindow), kmernum);
                
                break;
            case '+':
                kmernum = 0;
                for (int i = 0 ; i < rowdata[1] ; i++)
                {
                    kmernum += results[rowdata[ FIXCOL + i]];
                }
                
                get<1>(thiswindow) += kmernum;
                get<2>(thiswindow) = MAX(get<2>(thiswindow), kmernum);
                break;
            default:
                break;
        }
        
        rowdata = &rowdata[rowdata[1] + FIXCOL];
    }
    
}
