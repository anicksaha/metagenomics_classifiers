## This file contains all the necessary functions 
## to generate regions from variabilities

import numpy
def get_regions_info(variabilities):
    increase_windows = get_windows(variabilities, 20, True)
    decrease_windows = get_windows(variabilities, 20, False)
    regions_info = merge_windows(increase_windows, decrease_windows, 10)
    return regions_info

## For a window, get net-increasing/decreasing score
def get_score(seq):
    res = 0
    for i in range(1,len(seq)):
        res += seq[i] - seq[i - 1]
    return res

def get_windows(variabilities, num_regions, increasing = True):
    n = len(variabilities)
    windows = []
    vis = [0 for i in range(n)]
    X = []
    Y = []
    processed = 0
    
    for window_size in range(20, 31):
        for i in range(n - window_size + 1):
            end = i + window_size - 1
            curr = []
            for j in range(i, end + 1):
                curr.append(variabilities[j])
            
            score = get_score(curr)
            windows.append([score, [i, window_size]])

    windows = sorted(windows, reverse = increasing)
    
    for segment in windows:
        start = segment[1][0]
        end = segment[1][0] + segment[1][1] - 1
        valid = 1
        
        for i in range(start, end + 1):
            if vis[i] == 1:
                valid = 0
                break

        if valid == 1:
            x = []
            y = []
            for i in range(start, end + 1):
                vis[i] = 1
                x.append(i)
                y.append(variabilities[i])
            processed += 1
            X.append(x)
            Y.append(y)
            
        if processed == num_regions:
            return [X, Y]

def merge_windows(increase_windows, decrease_windows, num_regions):
    regions_info = []
    vis = [0 for i in range(len(increase_windows[0]))]
    processed = 0
    
    for sequence_1 in range(len(decrease_windows[0])):
        start_val1 = decrease_windows[0][sequence_1][0]
        minVal = 1e9
        idx1 = -1
        idx2 = -1
        L1 = L2 = []
        for sequence_2 in range(len(increase_windows[0])):
            end_val1 = increase_windows[0][sequence_2][-1]
            if vis[sequence_2] == 0 and end_val1 - start_val1 + 1 < minVal and end_val1 >= start_val1:
                minVal = end_val1 - start_val1 + 1
                idx1 = sequence_1
                idx2 = sequence_2
                start = start_val1
                end = end_val1
        if (idx2 == -1):
            continue

        regions_info.append([start, end])
        vis[idx2] = 1
        processed += 1

        if(processed == num_regions):
            break

    regions_info = sorted(regions_info)
    return regions_info