//
//  Created by Walfred (Wangfei) MA at the University of Southern California,
//  Mark Chaisson Lab on 2/13/23.
//
//  Licensed under the MIT License. 
//  If you use this code, please cite our work.
//   

#include "KmerWindow.hpp"

vector<string> namesplit(string& str, char deli)
{
    vector<string> eles;
    
    string::size_type start = 0;
    size_t end = str.find(deli);
    size_t len = str.length();
    
    while (end != std::string::npos && end < len)
    {
        eles.push_back(str.substr(start, end - start));
        start = end + 1;
        end = str.find(deli, start);
    }
    eles.push_back(str.substr(start, len));
    
    return eles;
}

void distract_segments(const vector<float>& covers_, vector<tuple<int,int,float,string>> &results)
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
                if (startposes[j] > 0 && index - startposes[j] + 1 >=  minwindowcutoff)
                {
                    results.push_back(make_tuple(startposes[j]-1,index,j,"NA"));
                    startposes[j] = 0;
                }
            }
            
            if (startposes[lastcover] > 0 && index - startposes[lastcover]+1 >=  minwindowcutoff)
            {
                results.push_back(make_tuple(startposes[lastcover]-1,index,lastcover,"NA"));
                startposes[lastcover] = 0;
            }
        }
        else if (cover * lastcover < 0)
        {
            for (int j = min( lastcover, 0) ; j < max(0, lastcover ); ++j )
            {
                if (startposes[j] > 0 && index - startposes[j] + 1 >= minwindowcutoff)
                {
                    results.push_back(make_tuple(startposes[j]-1,index,j,"NA"));
                    startposes[j] = 0 ;
                }
            }
            
            if (startposes[lastcover] > 0 && index - startposes[lastcover] + 1 >= minwindowcutoff)
            {
                results.push_back(make_tuple(startposes[lastcover]-1,index,lastcover,"NA"));
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


void refine_edges(const vector<tuple<int,int,int>> &windowcover, const vector<float>& extras, vector<tuple<int,int,float,string>> &results, const float depth)
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
        
        result = make_tuple(leftscore_max.second, rightscore_max.second + 1, copynum, -1);

    }
    
}

inline int getoverlap(vector<tuple<int,int,int,int>> &overlaps, int start,int end, string& coordi)
{
    auto segs = namesplit(coordi, '~');
    
    int total_overlap_size = 0 ;
    for (auto seg: segs)
    {
        auto coordis_ = namesplit(seg, '_');
        
        vector<int> coordis(coordis_.size());
        for (int i = 0; i < coordis_.size(); ++i)
        {
            coordis[i] = atoi(coordis_[i].c_str());
        }
        
        int overlap_size =  (end - start) + (coordis[1] - coordis[0]) - ( (MAX(end, coordis[1])) - (MIN(start, coordis[0])) );
        
        if (overlap_size > 100)
        {
            total_overlap_size += overlap_size;
            int arr[4] = {start, end, coordis[0], coordis[1]};
            std::sort(arr, arr + 4);
            
            int overlap_ps = arr[1];
            int overlap_pe = arr[2];
            
            int overlap_qs = coordis[2] + overlap_ps - coordis[0];
            int overlap_qe = coordis[3] + overlap_pe - coordis[1];
            
            overlaps.push_back(make_tuple(overlap_ps, overlap_pe, overlap_qs, overlap_qe));
            
        }
        
    }
    
    return total_overlap_size;

}

inline string getbestoverlap(const float cpchange, const FLOAT_T* reminders, const vector<string>& genenames, const string& pathname, const int start, const int end,const int size)
{
    
    int useindex = -1;
    FLOAT_T usereminder = 0.0;
    
    string result_str = "";
    result_str.reserve(100);
    for (int i = 0; i < genenames.size(); ++i)
    {
        if ( cpchange * reminders[i] < 0.1) continue;
        
        
        auto genename = genenames[i];
        
        auto elements = namesplit(genename, '\t');
        if ( elements.size() < 3) continue;
        auto element = elements[2];
        
        if (element == "NA:NA" || element == "NA" || element.find("_") == std::string::npos) continue;
        
        auto segments = namesplit(element, '&');
        
        int overlap_size = 0;
        for (auto &segment:segments )
        {
            auto cutpos = segment.find(':', 0);
            auto contigname = segment.substr(0, cutpos);
            auto coordi = segment.substr(cutpos+1, segment.length());
                        
            if (strncmp(pathname.c_str(), contigname.c_str(),contigname.length()) != 0) continue;
            
            vector<tuple<int,int,int,int>> overlaps ;
            
            auto total_overlap_size = getoverlap(overlaps, start,end, coordi);
            
            if (total_overlap_size > 0.5*size && cpchange * reminders[i] > usereminder)
            {
                usereminder = cpchange * reminders[i];
                useindex = i;
                
                result_str = genename.substr(0,genename.find('\t', 0))+":";
                for (size_t i = 0; i < overlaps.size(); ++i) 
                {
                    int val2 = std::get<2>(overlaps[i]);
                    int val3 = std::get<3>(overlaps[i]);
                    
                    result_str += std::to_string(val2) + "_" + std::to_string(val3);
                    
                    if (i < overlaps.size() - 1)  result_str += "~";
                }
                
            }
            
        }

    }
    
    return result_str;
    
}

void locate_partial(vector<tuple<int,int,float,string>> &results, const FLOAT_T* reminders, const vector<string>& genenames, const string pathname)
{
    
    for (auto &patial: results)
    {
        auto start =  30*(get<0>(patial)-1);
        auto end = 30*(get<1>(patial)-1);
        auto cpchange = (int) get<2>(patial);
        
        auto size = end - start;
        if (size < 3000) continue;
        
        auto result_str = getbestoverlap(cpchange, reminders, genenames, pathname, start, end, size);
        
        std::get<3>(patial) = result_str;
        
    }
    
}

void KmerWindow::PartialCopy(vector<vector<tuple<int, int, float, string>>>& results, const FLOAT_T* reminders, const vector<string>& genenames, const vector<string>& pathnames, const float depth)
{
    
    
    for (int path = 1; path < windowcovers.size(); ++path )
    {
        
        auto &windowcover = windowcovers[path];
        auto &result = results[path];
        result.clear();
        
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
            
            if (expt <= 200)
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
        
        locate_partial(result, reminders, genenames, pathnames[path]);
        
    }
}

void KmerWindow::WindowCovers(const uint16* kmervec, const uint16* kmermatrix, const FLOAT_T depth, const uint16 gnum, const uint knum, const int genenum, const int* results, ull &total_obs, ull &total_exp)
{
    
    const uint16* rowdata = kmermatrix;
    vector<uint64_t> kmer_samepos;
    kmer_samepos.reserve(knum);
    
    int copynum = 0;
    for (int i = 0; i < genenum; ++i)
    {
        if (results[i] > 0) copynum += results[i];
    }
    
    uint loc;
    int genecounter = 0;
    for (size_t i = 0; i < knum; ++i)
    {
        
        loc= ( ((uint)rowdata[3]) << 16 ) + rowdata[4];
        uint64_t kmer_pos = ( ((uint64_t)rowdata[2]) << 32 ) + loc;
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
                kmer_samepos.push_back(kmer_pos);
                
                //get<2>(thiswindow) = MAX(get<2>(thiswindow), genecounter);
                break;

            case '-':
                genecounter = copynum;
                for (int i = 0 ; i < rowdata[1] ; i++)
                {
                    genecounter -= results[rowdata[ FIXCOL + i]];
                }
                
                get<1>(thiswindow) += genecounter;
                kmer_samepos.push_back(kmer_pos);
                
                //get<2>(thiswindow) = MAX(get<2>(thiswindow), genecounter);
                
                break;
            case '+':
                genecounter = 0;
                for (int i = 0 ; i < rowdata[1] ; i++)
                {
                    genecounter += results[rowdata[ FIXCOL + i]];
                }
                
                get<1>(thiswindow) += genecounter;
                kmer_samepos.push_back(kmer_pos);
                //get<2>(thiswindow) = MAX(get<2>(thiswindow), genecounter);
                break;
            default:
                break;
        }
        
        total_exp += genecounter;
        total_obs += kmervec[i];
        
        rowdata = &rowdata[rowdata[1] + FIXCOL];
    }
    
    std::sort(kmer_samepos.begin(), kmer_samepos.end());
    
    ull lastpos = 0;
    for (auto &kmer_pos: kmer_samepos)
    {
        if (kmer_pos != lastpos)
        {
            uint16 first = (uint16) ( kmer_pos >> 32 );
            uint second = (uint)( kmer_pos & 0xFFFFFFFF );
            get<2>(windowcovers[first][second/window])++;
            lastpos = kmer_pos;
        }
    }
    
    for (auto &path: windowcovers)
    {
        for (auto &windowcover: path)
        {
            get<2>(windowcover) = (int)(1.0 * (float)get<1>(windowcover)/ (MAX(1,get<2>(windowcover))) + 0.5);
        }
    }
    
        
}
