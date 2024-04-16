#!/usr/bin/env python

import os
import sys
import argparse
from .eblocks import EblockDesign, Eblocks
from .primer import DesignPrimers
from .mutation import Mutation
from .sequence import Plasmid
from .utils import SnapGene



def eBlocksArgParser():
    parser = argparse.ArgumentParser()
    """
    Define and read command line arguments.
    """
    parser.add_argument("-g", 
                        "--gene",
                        type=str,
                        required=True,
                        default=None,
                        help="path to sequence file (.fasta) containing the gene of interest")
    
    parser.add_argument("-m",
                        "--mutations",
                        type=str,
                        required=True,
                        default=None,
                        help="path to file containing mutations to make")
    
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        required=True,
                        default=None,
                        help="path to output directory")
    
    parser.add_argument("-s",
                        "--species",
                        type=str,
                        required=False,
                        default="Mycobacterium Smegmatis",
                        help="species to make calculations for")
    
    parser.add_argument("-v",
                        "--vector",
                        type=str,
                        required=False,
                        default=None,
                        help="path to SnapGene file of the vector for which mutations will be made (.dna)")
    
    args = parser.parse_args()
    return args



if __name__ == "__main__":
    
    args = eBlocksArgParser()

    # Check if output file exists
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Parse mutations
    mutation_instance = Mutation()
    mutation_instance.parse_mutations(args.mutations)

    # Parse sequence
    sequence_instance = Plasmid()
    sequence_instance.parse_sequence(args.gene)
    sequence_instance.parse_vector(args.vector)

    snapgene_instance = SnapGene()
    snapgene_instance.parse_snapgene(args.vector)







    design_eblocks = DesignEblocks(sequence_fp=args.input_gene,
                                   mutations_fp=args.mutations,
                                   output_fp=args.output_location,
                                   species=args.species)
    design_eblocks.run()

    # Next; design IVA primers
    mut_gene_blocks_fp = os.path.join(args.output_location, "mut_gene_blocks.npy")
    wt_gene_blocks_fp = os.path.join(args.output_location, "wt_gene_blocks.npy")

    design_primers = DesignPrimers(wt_gene_blocks_fp, 
                                   mut_gene_blocks_fp, 
                                   args.output_location,
                                   args.input_gene,
                                   args.snapgene_file)
    design_primers.run()

    # Also write results to files that SnapGene can open
    primers_fp = os.path.join(args.output_location, "IVA_primers.csv")
    gene_blocks_mutation_info_fp = os.path.join(args.output_location, "gene_blocks.txt")
    
    if args.snapgene_file:
        snapgene_output = SnapGeneOutput(wt_gene_blocks_fp = wt_gene_blocks_fp,
                                         mut_gene_blocks_fp = mut_gene_blocks_fp,
                                         primers_fp = primers_fp,
                                         output_location = args.output_location,
                                         snapgene_file = args.snapgene_file,
                                         gene_blocks_info_fp = gene_blocks_mutation_info_fp)
        snapgene_output.run()

    # Finally; clean-up files when done
    # cleanup(args)

    print("Done!")
