#!/usr/bin/env python3

import os
import argparse

from src.DNABlocks import DNABlockDesign
from src.primer import DesignPrimers
from src.mutation import Mutation
from src.sequence import Gene, Vector
from src.plot import Plot

default_settings_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), 'settings')
)

def PyVADesignParser():
    parser = argparse.ArgumentParser(
        description="PyVADesign: Design of dsDNA fragments and primers for the introduction of mutations in a gene of interest")

    parser.add_argument("-g", 
                        "--gene",
                        type=str,
                        required=True,
                        default=None,
                        help="path to gene file (.fasta) containing the gene of interest")
    
    parser.add_argument("-m",
                        "--mutations",
                        type=str,
                        required=True,
                        default=None,
                        help="path to file containing mutations (.txt) to be made")
    
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        required=True,
                        default=None,
                        help="path to output directory")
    
    parser.add_argument("-c",
                        "--codon_usage",
                        type=str,
                        required=False,
                        default="U00096",
                        help="Codon usage to use for the design of the mutations")
    
    parser.add_argument("-v",
                        "--vector",
                        type=str,
                        required=True,
                        default=None,
                        help="path to the Plasmid file (.gb/.dna) containing the gene of interest")
    
    parser.add_argument("-ds",
                        "--design_settings",
                        type=str,
                        required=False,
                        default=None,
                        help="path to the settings file (.txt) containing the settings for the design of dsDNA fragments")

    parser.add_argument("-sps",
                        "--sequence_primer_settings",
                        type=str,
                        required=False,
                        default=f"{default_settings_path}/primer3-seq-settings.txt",
                        help="path to the settings file (.txt) containing the settings for the design of sequencing primers")
    
    parser.add_argument("-ips",
                        "--iva_primer_settings",
                        type=str,
                        required=False,
                        default=f"{default_settings_path}/primer3-settings.txt",
                        help="path to the settings file (.txt) containing the settings for the design of IVA primers")
    
    args = parser.parse_args()
    return args



if __name__ == "__main__":
    
    args = PyVADesignParser()

    # Check if output file exists, otherwise create it
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Parse mutations
    mutation_instance = Mutation()
    mutation_instance.parse_mutations(args.mutations)

    # Parse sequence files
    gene_instance = Gene(stopcodon=False)
    gene_instance.parse_sequence(args.gene)

    vector_instance = Vector(gene=gene_instance)
    vector_instance.parse_vector(args.vector)
    
    # Create design instance and run the design of the dsDNA fragments
    design_instance = DNABlockDesign(mutation_instance=mutation_instance,
                                    vector_instance=vector_instance,
                                    gene_instance=gene_instance,
                                    codon_usage=args.codon_usage,
                                    output_dir=args.output,
                                    settings_file=args.design_settings)
    
    design_instance.run_design_DNABlocks()

    # Create plot instance
    plot_instance = Plot(DNABlocks_design_instance=design_instance,
                         mutation_instance=mutation_instance,
                         vector_instance=vector_instance,
                         gene_instance=gene_instance,
                         show=False,
                         save=True,
                         output_dir=args.output)
    # Visualize the dsDNA fragments
    plot_instance.plot_DNABlocks_mutations(figure_length=20, figure_width=5)

    # Design the primers
    primers_instance = DesignPrimers(mutation_instance=mutation_instance,
                                     DNABlocks_design_instance=design_instance,
                                     vector_instance=vector_instance,
                                     primers_settingsfile=args.iva_primer_settings,
                                     seqprimers_settingsfile=args.sequence_primer_settings,
                                     output_dir=args.output)
    
    primers_instance.run_design()