#!/usr/bin/env python3

import os
import argparse
import collections as cl
import re
import gzip
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
from numba import njit
from multiprocessing import Pool


@njit
def merge_numba_flat(coordinates, mergedistance=150):
    n = len(coordinates)
    sorted_index = np.argsort(coordinates)
    merged = np.empty(n, dtype=np.uint32)
    out_i = -2
    current_start = -mergedistance - 1
    current_end = -mergedistance - 1
    
    for i in range(n):
        idx = sorted_index[i]
        if idx % 2 == 0:
            start = coordinates[idx]
            end = coordinates[idx + 1]
            if current_end + mergedistance >= start:
                current_end = max(current_end, end)
            else:
                if out_i >= 0:
                    merged[out_i] = current_start
                    merged[out_i + 1] = current_end
                out_i += 2
                current_start = start
                current_end = end
                
                
    merged[out_i] = current_start
    merged[out_i + 1] = current_end
    out_i += 2
    
    return merged[:out_i]

def merge_one(item):
    key, all_coords = item
    flat_coords = np.concatenate(all_coords)
    merged = merge_numba_flat(flat_coords, mergedistance=400)
    return key, [merged]

def MergeAllregions(merged_regions, nthreads=4):
    """
    Merges all coordinate lists in parallel using the specified number of threads.
    
    Args:
        merged_regions (dict): A dictionary of {key: list of numpy arrays of coordinates}.
        nthreads (int): Number of parallel processes to use.
    
    Returns:
        dict: A dictionary of {key: [merged numpy array]}.
    """
    
    items = list(merged_regions.items())
    
    with Pool(nthreads) as pool:
        results = pool.map(merge_one, items)
        
    return cl.defaultdict(list, results)

def process_file(file_path):
    # Get the line count from the file using wc -l
    
    counts = cl.defaultdict(int)
    with open(file_path, mode='r') as f:
        for line in f:
            line = line.strip().split()
            counts[(line[0],line[3])] += 1
            
    # Allocate numpy array with enough space for both start and end coordinates
    coordinates = { key :np.empty(count * 2, dtype=np.uint32) for key, count in counts.items()}  # *2 for both start and end
    
    counts = cl.defaultdict(int)
    idx = 0  # To keep track of where to place data in the array
    
    # Process each file to extract and store regions
    with open(file_path, mode='r') as f:
        for line in f:
            line = line.strip().split()
            start, end = int(line[1]), int(line[2])
            key = (line[0],line[3])
            idx = counts [key]
            coordinates[key][idx] = start
            coordinates[key][idx+1] = end
            counts[key] += 2
            
    return coordinates


def process_batch( merged_regions, files, nthreads):
    # Use Pool to read and process all files in parallel
    
    with Pool(nthreads) as pool:
        local_results = pool.map(process_file, files)
        
    # Merge the results from all threads
    for result in local_results:
        for key, regions in result.items():
            merged_regions[key].append(regions)
            
    # After processing the batch, merge the regions
    return MergeAllregions(merged_regions, nthreads)


    
def summary(inputfolder, outputfile, nthreads):
    findex = 0
    regions = cl.defaultdict(list)
    
    files = os.listdir(inputfolder)
    # Divide files into batches of 16
    batch_size = 16
    for i in range(0, len(files), batch_size):
        batch_files = files[i:i+batch_size]
        print(f"Processing batch {i // batch_size + 1} with files: {batch_files}")
        
        # Process the batch of files in parallel
        regions = process_batch(regions, [os.path.join(inputfolder, file) for file in batch_files], nthreads)
        
    keys_sort = sorted(regions.keys())
    
    with open(outputfile, mode='w') as f:
        for key in sorted(regions.keys()):
            region_chunks = regions[key]  # This is a list of numpy arrays
            region_array = np.concatenate(region_chunks)  # Flatten to one array
            
            for i in range(0, len(region_array), 2):
                f.write(f"{key[0]}\t{region_array[i]}\t{region_array[i+1]}\t{key[1]}\n")
               
def main(args):
    summary(args.input, args.output, args.nthreads)

def run():
    parser = argparse.ArgumentParser(description="Program to merge ctyper profiling bed files")
    parser.add_argument("-i", "--input", help="Path to input data folder", dest="input", type=str, required=True)
    parser.add_argument("-o", "--output", help="Path to output file", dest="output", type=str, required=True)
    parser.add_argument("-n", "--nthreads", help="Number of threads to use", dest="nthreads", type=int, default=8)
    
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    dummy = np.array([0, 1], dtype=np.uint32)
    _ = merge_numba_flat(dummy)
    run()
