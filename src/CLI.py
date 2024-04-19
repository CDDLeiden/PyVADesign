#!/usr/bin/env python

import os
import sys
import argparse

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from src.eblocks import EblockDesign
from src.primer import DesignPrimers
from src.mutation import Mutation
from src.sequence import Plasmid
from src.utils import SnapGene
from src.plot import Plot

# run command with tutorial data (from root directory):
# python src/CLI.py -g tutorial/files/A0QX55.fasta -m tutorial/files/mutations.txt -o tutorial/output -v tutorial/files/vector.dna -s data/codon_usage/Mycobacterium_smegmatis.csv
# python src/CLI.py -g tutorial/files/A0QX55.fasta -o tutorial/output -v tutorial/files/vector.dna -s data/codon_usage/Mycobacterium_smegmatis.csv -m tutorial/files/mutations_begin_end_of_gene.txt



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
                        default="./data/codon_usage/Escherichia_coli.csv",
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
    sequence_instance.parse_vector(args.vector)
    sequence_instance.parse_sequence(args.gene)


    snapgene_instance = SnapGene(sequence_instance=sequence_instance,
                                     output_dir=args.output)
    
    design_instance = EblockDesign(sequence_instance=sequence_instance,
                                   mutation_instance=mutation_instance,
                                   output_dir=args.output,
                                   codon_usage=args.species)
    
    plot_instance = Plot(mutation_instance=mutation_instance,
                         sequence_instance=sequence_instance,
                         eblocks_design_instance=design_instance,
                         output_dir=args.output)
    
    
    # Check the input vector
    plot_instance.plot_vector(figsize=(7, 7))

    # Next; design eblocks
    design_instance.run_design_eblocks()

    plot_instance.plot_eblocks_mutations(figure_length=20,
                                         figure_width=5)
    plot_instance.plot_histogram_mutations()

    # Next; design IVA primers
    primers_instance = DesignPrimers(mutation_instance=mutation_instance,
                                     eblocks_design_instance=design_instance,
                                     sequence_instance=sequence_instance, 
                                     output_dir=args.output,
                                     snapgene_instance=snapgene_instance)
    
    primers_instance.run_design()