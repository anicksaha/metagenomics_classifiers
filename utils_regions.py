## This file contains all the necessary functions 
## to generate regions from variabilities

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
    
def merge_windows(increase_windows, decrease_windows, num_regions):
    regions_info = []
    vis = [0 for i in range(len(increase_windows[0]))]
    done = 0
    
    for seq1 in range(len(decrease_windows[0])):
        start_val1 = decrease_windows[0][seq1][0]
        minVal = 1e9
        idx1 = -1
        idx2 = -1
        L1 = []
        L2 = []
        for seq2 in range(len(increase_windows[0])):
            end_val1 = increase_windows[0][seq2][-1]
            if vis[seq2] == 0 and end_val1 - start_val1 + 1 < minVal and end_val1 >= start_val1:
                minVal = end_val1 - start_val1 + 1
                idx1 = seq1
                idx2 = seq2
                start = start_val1
                end = end_val1
        if (idx2 == -1):
            continue

        regions_info.append([start, end])
        vis[idx2] = 1
        done += 1

        if(done == num_regions):
            break

    regions_info = sorted(regions_info)
    return regions_info

def get_windows(Y, num_regions, increasing = True):
    n = len(Y)
    windows = []
    vis = [0 for i in range(n)]
    retX = []
    retY = []
    ret = 0
    
    # Limits the window size of a decreasing region to be between 20 and 30.
    for window_size in range(20, 31):
        for i in range(n - window_size + 1):
            end = i + window_size - 1
            prev = []
            for j in range(i, end + 1):
                val_at_j = 0
                prev.append(Y[j])
            
            score = get_score(prev)
            windows.append([score, [i, window_size]])

    windows = sorted(windows, reverse=increasing)
    
    for L in windows:
        start = L[1][0]
        end = L[1][0] + L[1][1] - 1
        flag = 0
        for i in range(start, end + 1):
            if vis[i] == 1:
                flag = 1
                break
        if flag == 0:
            x = []
            y = []
            for i in range(start, end + 1):
                vis[i] = 1
                x.append(i)
                y.append(Y[i])
            ret += 1
            retX.append(x)
            retY.append(y)
            
        if ret == num_regions:
            return [retX, retY]