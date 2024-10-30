import sys
import os
import random
from datetime import datetime

from src.sequence import Vector, Gene
from tests.CreateRandomMutations import RandomMutations
from src.mutation import Mutation
from src.eblocks import EblockDesign, Eblock
from src.primer import DesignPrimers

# Primer3 settings
settingsfile = '/home/rosan/git/design_gene_blocks/settings/primer3-settings-mild.txt'
seq_settingsfile = '/home/rosan/git/design_gene_blocks/settings/primer3-seq-settings.txt'

# Paths
inpath="/zfsdata/data/rosan/eBlocks/testset_input"
outpath="/zfsdata/data/rosan/eBlocks/testset_output"

# Files
folders = ["Mtb_pINIT_DnaE1", "Mav_pACE_MmpL3", "Mabs_pCPF3_05_RNAP-core"]
sequencefiles = ["A0QX55.fasta", "A0A0H2ZYQ2.fasta", "rpoB-rpoC_complex.fasta"]
vectorfiles = ["dnae1-pinit.gb", "pACE_mmpL3.dna", "pCPF3_05_RNAP-core.dna"]

# Settings
n_repeats = 10
end_range = 500
optimization = ["amount", "cost"]

if __name__ == '__main__':

    # Loop over testdata and do eblock calculations for each of the files
    for f,seq,vec in zip(folders, sequencefiles, vectorfiles):

        sequencefile=f"{inpath}/{f}/{seq}"
        vectorfile=f"{inpath}/{f}/{vec}"

        # Create instance sfor Gene and Vector
        gene_instance = Gene(stopcodon=False)  # One of the genes contains histag, so stopcodon is set to False to avoid throwing an no-stopcodon error
        gene_instance.parse_sequence(sequencefile)

        vector_instance = Vector(gene=gene_instance)
        vector_instance.parse_vector(vectorfile)

        # Loop over the files and try to make the eblocks
        for i in range(10, end_range, 10):

            curdir = f"{inpath}/{f}/datasets/N{i}"
            print(f"Processing {curdir}")

            # Create output directories and output files
            os.makedirs(f"{outpath}/{f}/N{i}", exist_ok=True)
            with open(f"{outpath}/{f}/N{i}/results.txt", "w") as fo:
                fo.write(f"gene\tN\tfile\toptimization\teblock_finished\tnum_eblocks\tlength_eblocks\tcost\tprimer_finished\n")

                # Loop over the files in the directory
                for file in os.listdir(curdir):

                    # Calculate eblocks for each file
                    print(f, i, file)

                    # Create directory for output
                    outputdata = f"{outpath}/{f}/N{i}/{file}"
                    os.makedirs(outputdata, exist_ok=True)

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

                        os.makedirs(f"{outputdata}/{opt}", exist_ok=True)

                        # Create Eblocks
                        design_instance = EblockDesign(mutation_instance=mutation_instance,
                                                    gene_instance=gene_instance,
                                                    vector_instance=vector_instance,
                                                    output_dir=f"{outputdata}/{opt}",
                                                    verbose=True,
                                                    clone_files=True,
                                                    cost_optimization=cost_optimization,
                                                    amount_optimization=amount_optimization)
                        # Store results
                        eblock_finished = False
                        num_eblocks = -1
                        lentgh_eblocks = -1
                        cost = -1
                        primer_finished = False


                        try:
                            design_instance.run_design_eblocks()
                            eblock_finished = True
                            num_eblocks = len(design_instance.wt_eblocks)
                            cost = design_instance.cost
                            length_eblocks = ','.join([str(len(eblock.sequence)) for eblock in design_instance.wt_eblocks])
                            
                            # run primer design
                            primers_instance = DesignPrimers(mutation_instance=mutation_instance,
                                                            eblocks_design_instance=design_instance,
                                                            primers_settingsfile=settingsfile,
                                                            seqprimers_settingsfile=seq_settingsfile,
                                                            vector_instance=vector_instance,
                                                            output_dir=f"{outputdata}/{opt}")
                            try:
                                primers_instance.run_design()
                                primer_finished = True
                            except:
                                primer_finished = False

                            # Write results to file
                            fo.write(f"{f}\tN{i}\t{file}\t{opt}\t{eblock_finished}\t{str(num_eblocks)}\t{length_eblocks}\t{cost}\t{primer_finished}\n")
                            fo.flush()

                        except:
                            eblock_finished = False
                            num_eblocks = -1
                            lentgh_eblocks = -1
                            cost = -1
                            primer_finished = False
                            fo.write(f"{f}\tN{i}\t{file}\t{opt}\t{eblock_finished}\t{num_eblocks}\t{length_eblocks}\t{cost}\t{primer_finished}\n")
                            fo.flush()

            sys.exit()