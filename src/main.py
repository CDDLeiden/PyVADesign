import os
import sys
import argparse
# import snapgene_output as sno
from design_gene_blocks import DesignEblocks
# from design_IVA_primers import DesignPrimers

def read_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input_gene", required=True, help="FASTA file containing the gene of interest")
    parser.add_argument("-sp", "--species", type=str, default="Mycobacterium Smegmatis", help="Species to make calculations for")
    parser.add_argument("-m", "--mutations", required=True, help="TXT file containing the mutations to make")
    parser.add_argument("-o", "--output_location", required=True, help="Location where to store the output of the script")
    parser.add_argument("-s", "--snapgene_file", required=False, help="Snapgene DNA file of the vector for which mutations will be made")
    args = parser.parse_args()
    return args

def cleanup(args):
    to_remove = [os.path.join(args.output_location, "mut_gene_blocks.npy"),
                 os.path.join(args.output_location, "wt_gene_blocks.npy")]
    for fp in to_remove:
        if os.path.exists(fp):
            os.remove(fp)
        else:
            print(f"{fp} does not exist")
    

if __name__ == "__main__":
    
    # First design gene blocks
    args = read_arguments()
    instance = DesignEblocks(sequence_fp=args.input_gene,
                             mutations_fp=args.mutations,
                             output_fp=args.output_location,
                             species=args.species)
    instance.run()


    sys.exit()

    # Next; design IVA primers
    mut_gene_blocks_fp = os.path.join(args.output_location, "mut_gene_blocks.npy")
    wt_gene_blocks_fp = os.path.join(args.output_location, "wt_gene_blocks.npy")
    dip.main(mut_gene_blocks_fp, args.input_gene, args.output_location, args.snapgene_file)

    # Also write results to files that SnapGene can open
    primers_fp = os.path.join(args.output_location, "IVA_primers.csv")
    gene_blocks_mutation_info_fp = os.path.join(args.output_location, "gene_blocks.txt")
    if args.snapgene_file:
        sno.main(wt_gene_blocks_fp,
                 mut_gene_blocks_fp,
                 primers_fp,
                 args.output_location,
                 args.snapgene_file,
                 gene_blocks_mutation_info_fp)

    # Finally; clean-up files when done
    cleanup(args)
