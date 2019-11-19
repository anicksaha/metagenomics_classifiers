import sys
import utils
import utils_vis

def main(filepath):
    # First step is to read the file and get the sequences
    sequences, seq_length = utils.read_fasta_file(filepath)

    variabilities = utils.create_variabilities_file(sequences, seq_length)

    # smoothing of variabilities
    variabilities_2 = utils.get_smoothing(variabilities)

    utils_vis.generate_variability_plot(variabilities_2)

    # Eye balling | TO_DO    
    regions_info = utils.create_regions_file()
    
    # Create fasta_files for phylo tree
    utils.create_fasta_files(sequences, seq_length, regions_info)

if __name__ == '__main__':
    main(sys.argv[1])