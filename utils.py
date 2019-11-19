import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
import collections

##################################################################################

## Read in the Sequence from input file
def read_fasta_file(filepath):
    input_sequence = {} # Mapping from {seq_id to sequence}
    with open(filepath) as f:
        input_lines = f.read().splitlines()
    seq_id = None
    seq_length = None
    for i, line in enumerate(input_lines):
        if i % 2 == 0:
            seq_id = line[1:]
        else:
            input_sequence[seq_id] = line
            seq_length = len(line)
    return input_sequence, seq_length

##################################################################################
  
def create_variabilities_file(sequences, seq_length):
    variabilities = []
    for i in range(seq_length):
        M = collections.defaultdict(int) # Frequency map for BASES
        # count frequency skipping the gaps
        for key, val in sequences.items():
            if val[i] == '-':
                continue
            M[val[i]] += 1

        if len(M) > 0:
            # Divide Max frequency with total frequency
            total_freq = len(sequences)
            max_base = max(M, key=M.get)
            # for val in M.values():
            #     total_freq += val
            variabilities.append(float(M[max_base]) / total_freq)
        else:
            variabilities.append(float(0))

    ## store in a file
    with open('variability.txt', 'w') as f:
        for variability in variabilities:
            f.write(str(variability) + '\n')

    return variabilities

##################################################################################

import statistics
def get_smoothing(variabilities):
    window_len = 15
    variabilities_2 = []
    N = len(variabilities)
    for i in range(window_len,N):
        variabilities_2.append(statistics.mean(variabilities[i-window_len:i+1]))
    return variabilities_2

##################################################################################

def create_regions_file():
    regions_info = get_regions_info()
    with open('regions.txt', 'w') as f:
        for start, end in regions_info:
            f.write(str(start) + '\t' + str(end) + '\n')
    return regions_info

def get_regions_info():
    regions_info = [ (0, 20), (150, 210), (400, 450), (530, 600), (670, 710),
    (780, 830), (950, 1000), (1070, 1120), (1200, 1240), (1380, 1420)]
    return regions_info

##################################################################################

def write_fasta_file(seq, filepath):
    with open(filepath, 'w') as f:
        for key, val in seq.items():
            f.write('>' + key + '\n')
            f.write(val + '\n')


def create_fasta_files(seq, seq_length, regions_info):
    seq_subset = {key: seq[key] for key in list(seq)[:100]} # random 100 sequences
    whole_16s = collections.defaultdict(str)

    # Fill whole_16s seq by ignoring gaps
    for i in range(seq_length):
        num_gaps = 0
        for key, val in seq_subset.items():
            if val[i] == '-':
                num_gaps += 1
        if num_gaps == 100:
            continue
        for key, val in seq_subset.items():
            whole_16s[key] += val[i]

    write_fasta_file(whole_16s, 'whole_16s.fna')

    # Fill region_1 and region_4 seq
    region_2_seq = {}
    region_4_seq = {}
    region_2 = regions_info[1]
    region_4 = regions_info[3]

    for key, val in whole_16s.items():
        region_2_seq[key] = val[region_2[0]:region_2[1]+1]
        region_4_seq[key] = val[region_4[0]:region_4[1]+1]

    write_fasta_file(region_2_seq, 'region2.fna')
    write_fasta_file(region_4_seq, 'region4.fna')