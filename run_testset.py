import sys
import os
import random
from datetime import datetime

from src.sequence import Vector, Gene
from tests.CreateRandomMutations import RandomMutations
from src.mutation import Mutation
from src.eblocks import EblockDesign, Eblock

settingsfile = 'settings/primer3-settings.txt'
seq_settingsfile = 'settings/primer3-seq-settings.txt'

if __name__ == '__main__':

    n_repeats = 10
    end_range = 500
    optimization = ["cost", "amount"]

    inpath="/zfsdata/data/rosan/eBlocks/testset_input"
    outpath="/zfsdata/data/rosan/eBlocks/testset_output"

    folders = ["Mtb_pINIT_DnaE1", "Mav_pACE_MmpL3", "Mabs_pCPF3_05_RNAP-core"]
    sequencefiles = ["A0QX55.fasta", "A0A0H2ZYQ2.fasta", "rpoB-rpoC_complex.fasta"]
    vectorfiles = ["dnae1-pinit.gb", "pACE_mmpL3.dna", "pCPF3_05_RNAP-core.dna"]

    for f,seq,vec in zip(folders, sequencefiles, vectorfiles):

        sequencefile=f"{inpath}/{f}/{seq}"
        vectorfile=f"{inpath}/{f}/{vec}"

        # Create instancesfor Gene and Vector
        gene_instance = Gene(stopcodon=False)  # One of the genes contains histag, so stopcodon is set to False
        gene_instance.parse_sequence(sequencefile)

        vector_instance = Vector(gene=gene_instance)
        vector_instance.parse_vector(vectorfile)

        # Loop over the files and try to make the eblocks
        for i in range(10, end_range, 10):

            curdir = f"{inpath}/{f}/datasets/N{i}"

            # Create output directories
            os.makedirs(f"{outpath}/{f}/N{i}", exist_ok=True)
            # Create output file
            fo = open(f"{outpath}/{f}/N{i}/results.txt", "w")
            fo.write(f"{f}\tN{i}\tfile\toptimization\teblock_finished\tnum_eblocks\tlength_eblocks\tcost\tprimer_finished\n")

            # Loop over the files in the directory
            for file in os.listdir(curdir):

                # Calculate eblocks for each file
                print(f, i, file)

                # Read mutations
                mutation_instance = Mutation()
                mutation_instance.parse_mutations(f"{curdir}/{file}")

                for opt in optimization:

                    if opt == "cost":
                        cost_optimization = True
                        amount_optimization = False
                    elif opt == "amount":
                        cost_optimization = False
                        amount_optimization = True

                    # Create Eblocks
                    design_instance = EblockDesign(mutation_instance=mutation_instance,
                                                gene_instance=gene_instance,
                                                vector_instance=vector_instance,
                                                output_dir=f"{outpath}/{f}/N{i}",
                                                verbose=True,
                                                clone_files=False,
                                                cost_optimization=cost_optimization,
                                                amount_optimization=amount_optimization)
                    try:
                        design_instance.run_design_eblocks()
                        eblock_finished = True
                        num_eblocks = len(design_instance.wt_eblocks)
                        lentgh_eblocks = [len(eblock) for eblock in design_instance.wt_eblocks]
                        cost = design_instance.cost
                        # run primers


primers_instance = DesignPrimers(mutation_instance=mutation_instance,
                                 eblocks_design_instance=design_instance,
                                 primers_settingsfile=settingsfile,
                                 seqprimers_settingsfile=seq_settingsfile,
                                 vector_instance=vector_instance,
                                 output_dir=output_dir)

primers_instance.run_design()




                    except:
                        print(f"Failed {f} {i} {file} {opt}")
                    # TODO Delete files and clones (do not store everything will become too big)