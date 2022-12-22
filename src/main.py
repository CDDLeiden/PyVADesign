import os
import sys
import argparse
import design_gene_blocks as dgb
import design_IVA_primers as dip

def read_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input_gene", help="FASTA file containing the gene of interest")
    parser.add_argument("-m", "--mutations", help="TXT file containing the mutations to make")
    parser.add_argument("-o", "--output_location", help="Location where to store the output of the script")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    
    # First design gene blocks
    args = read_arguments()
    dgb.main(args)

    # Next; design IVA primers
    result_file = os.path.join(args.output_location, "gene_blocks.npy")
    dip.main(result_file, args.input_gene, args.output_location)

    # Finally; write results to files that SnapGene can open

    # TODO Remove unnecessary files when done