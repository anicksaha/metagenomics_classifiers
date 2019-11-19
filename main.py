import sys
import utils

def main(filepath):
    # First step is to read the file and get the sequences
    sequences, seq_length = utils.read_fasta_file(filepath)

    # Get variabilities
    variabilities = utils.create_variabilities_file(sequences, seq_length)

    # Smoothening of variabilities
    variabilities_2 = utils.get_smoothing(variabilities)

    # Genarate the variability plot using smoothened variabilities
    utils.generate_variability_plot(variabilities_2)

    # Create regions using the smoothened variabilities list
    regions_info = utils.create_regions_file(variabilities_2)
    
    # Create fasta_files for phylo tree [Regions: 2 and 4]
    utils.create_fasta_files(sequences, seq_length, regions_info)

if __name__ == '__main__':
    main(sys.argv[1])