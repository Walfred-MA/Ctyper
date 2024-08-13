//
//  KmerWindow.cpp
//  CTyper
//
//  Created by Wangfei MA on 6/30/23.
//

#include "KmerWindow.hpp"

void distract_segments(const vector<float>& covers_, vector<tuple<int,int,float>> &results)
{
    vector<int> covers;
    for (auto cover_: covers_)
    {
        covers.push_back(floor(cover_+0.5));
    }
    covers.push_back(0);
    
    vector<int> missingpoints;
    missingpoints.reserve(covers.size());
    
    int lastvalid = 0;
    for (int index = 0; index < covers.size(); ++index)
    {
        auto cover = covers[index];
        
        if (cover == MAX_UINT16)
        {
            missingpoints.push_back(index);
            continue;
        }
        
        if (cover == lastvalid)
        {
            for (int misspoint: missingpoints)
            {
                covers[misspoint] = cover;
            }
        }
        else
        {
            for (int misspoint: missingpoints)
            {
                covers[misspoint] = 0;
            }
        }
        
        missingpoints.clear();
        lastvalid = cover;
    }
    
    
    unordered_map<int, int> startposes;
    int lastcover = 0;
    for (int index = 0; index < covers.size(); ++index)
    {
        
        auto cover = covers[index];
        
        if (cover == lastcover)
        {
            continue;
        }
        
        if (cover * lastcover >= 0 && abs(cover) < abs(lastcover))
        {
            for (int j = min(cover, lastcover) ; j < max(cover, lastcover); ++j )
            {
                if (startposes[j] > 0)
                {
                    results.push_back(make_tuple(startposes[j]-1,index,j));
                    startposes[j] = 0;
                }
            }
            
            if (startposes[lastcover] > 0)
            {
                results.push_back(make_tuple(startposes[lastcover]-1,index,lastcover));
                startposes[lastcover] = 0;
            }
        }
        else if (cover * lastcover < 0)
        {
            for (int j = min( lastcover, 0) ; j < max(0, lastcover ); ++j )
            {
                if (startposes[j] > 0)
                {
                    results.push_back(make_tuple(startposes[j]-1,index,j));
                    startposes[j] = 0 ;
                }
            }
            
            if (startposes[lastcover] > 0)
            {
                results.push_back(make_tuple(startposes[lastcover]-1,index,lastcover));
                startposes[lastcover] = 0;
            }
        }
        
        if (cover != 0 )
        {
            startposes[cover] = index+1;
        }
        
        lastcover = cover;
    }
}


void refine_edges(const vector<tuple<int,int,int>> &windowcover, const vector<float>& extras, vector<tuple<int,int,float>> &results, const float depth)
{
    for (auto &result: results)
    {
        int startindex = get<0>(result);
        int endindex = get<1>(result);
        int copynum = floor(get<2>(result) + 0.5);
        int sign = (int)(copynum >= 0);
        
        int startwindow = startindex * windowmerge + 1;
        int endwindow = endindex * windowmerge;
        
        auto startrange = make_pair(max(1, startwindow - windowmerge),min(startwindow + windowmerge, (int)windowcover.size()));
        auto endrange = make_pair(max(1, endwindow - windowmerge),min(endwindow + windowmerge, (int)windowcover.size()));
        
        float score_sum = 0.0;
        pair<float, int> leftscore_max = make_pair(0.0, startwindow);
        for (int index = startrange.second; index >= startrange.first; --index)
        {
            auto &thewindow = windowcover[index];
            
            int query = get<0>(thewindow);
            int expect = get<1>(thewindow);
            int cnum = get<2>(thewindow);
            
            if (cnum == 0) continue;
            
            score_sum +=   sign * ( cnum * query / depth / expect - cnum  - copynum );
            if (score_sum > leftscore_max.first)
            {
                leftscore_max = make_pair(score_sum, index);
            }
        }
        
        pair<float, int> rightscore_max = make_pair(0.0, endwindow);
        for (int index = endrange.first; index < endrange.second; ++index)
        {
            auto &thewindow = windowcover[index];
            
            int query = get<0>(thewindow);
            int expect = get<1>(thewindow);
            int cnum = get<2>(thewindow);
            
            if (cnum == 0) continue;
            
            score_sum +=   sign * ( cnum * query / depth / expect - cnum  - copynum );
            if (score_sum > rightscore_max.first)
            {
                rightscore_max = make_pair(score_sum, index);
            }
        }
        
        result = make_tuple(leftscore_max.second, rightscore_max.second + 1, copynum);

    }
    
}

void KmerWindow::PartialCopy(vector<vector<tuple<int, int, float>>>& results, const float depth)
{
    for (int path = 1; path < windowcovers.size(); ++path )
    {
        
        auto &windowcover = windowcovers[path];
        vector<tuple<int, int, float>> &result = results[path];
        
        ull totalwindow = 0;
        for (int index = 1; index < windowcover.size(); ++index)
        {
            auto &thepair = windowcover[index];
            totalwindow += get<1>(thepair);
        }
        if (totalwindow < 100) continue;
        
        auto windowsize = windowcover.size() - 1;
        
        auto mergenum = (windowsize-1)/windowmerge + 2;
        
        vector<float> covers(mergenum);
        vector<int> covers_int(mergenum);
        vector<float> extras(mergenum);
        vector<float> expts(mergenum);
        vector<float> cpnums(mergenum);
        vector<int> elenums(mergenum);
        
        for (int index = 0; index < windowsize; ++index)
        {
            auto &thepair = windowcover[index + 1];
            
            int query = get<0>(thepair);
            int expect = get<1>(thepair);
            int cnum = get<2>(thepair);
            
            if (cnum == 0) continue;
            
            auto &elenum = elenums[index/windowmerge];
            auto &elenum2 = elenums[index/windowmerge+1];
            
            extras[index/windowmerge]  += (query / depth);
            expts[index/windowmerge] +=  expect   ;
            cpnums[index/windowmerge] += cnum;
            
            extras[index/windowmerge+1] += (query / depth);
            expts[index/windowmerge+1] += expect   ;
            cpnums[index/windowmerge+1] += cnum;
            
            elenum ++;
            elenum2 ++;
        }
        
        for (int index = 0; index < mergenum; ++index)
        {
            auto &extra = extras[index];
            auto &expt = expts[index];
            auto &cpnum = cpnums[index];
            
            int elenum = elenums[index];
                        
            float cpnum_mean = cpnum/elenum;
            
            if (expt <= 30)
            {
                covers[index] = MAX_UINT16;
                continue;
            }
            
            if (abs( extra/expt - 1.0) > 0.25 )
            {
                covers[index] = cpnum_mean  * extra/expt - cpnum_mean;
            }
            else
            {
                covers[index] = 0;
            }
                                    
        }
        
        distract_segments(covers, result);
        
        refine_edges(windowcover, extras, result, depth);
    }
}

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

	float count_f = (int) kmervec[i];
        
        const uint16 flag = rowdata[5];
        if (optioncorr && (flag & 0x3F) >= errorcutoff1)
        {
            float corr = 0.01 * (flag & 0xFFC0)/64;
            count_f *= corr;
        }
        
        get<0>(thiswindow) += (int)(count_f+0.5);

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
