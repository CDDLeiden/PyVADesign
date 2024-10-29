import sys
import os
import random
from datetime import datetime

from src.sequence import Vector, Gene
from tests.CreateRandomMutations import RandomMutations


if __name__ == '__main__':

    n_repeats = 10
    end_range = 500

    path="/zfsdata/data/rosan/eBlocks/testset_input"

    folders = ["Mtb_pINIT_DnaE1", "Mav_pACE_MmpL3", "Mabs_pCPF3_05_RNAP-core"]
    sequencefiles = ["A0QX55.fasta", "A0A0H2ZYQ2.fasta", "rpoB-rpoC_complex.fasta"]
    vectorfiles = ["dnae1-pinit.gb", "pACE_mmpL3.dna", "pCPF3_05_RNAP-core.dna"]

    for f,seq,vec in zip(folders, sequencefiles, vectorfiles):

        print(f)

        sequencefile=f"{path}/{f}/{seq}"
        vectorfile=f"{path}/{f}/{vec}"

        # Output directory
        output_dir=f"{path}/{f}"
        os.makedirs(f"{output_dir}/datasets", exist_ok=True)

        # Create instances
        gene_instance = Gene(stopcodon=False)
        gene_instance.parse_sequence(sequencefile)

        vector_instance = Vector(gene=gene_instance)
        vector_instance.parse_vector(vectorfile)

        # Make random mutations
        for i in range(10, end_range, 10):
            
            print(i)

            # Create output directories
            directory=f"N{i}"
            os.makedirs(f"{output_dir}/datasets/{directory}", exist_ok=True)

            for j in range(n_repeats): # for each number of mutations, create 10 files
                random_mutations = RandomMutations()
                random_mutations.make_random_mutations(nucleotide_sequence=gene_instance.sequence,
                                                    n_single_mutations=i, 
                                                    n_multiple_mutations=i,
                                                    n_deletions=i,
                                                    n_insertions=i, 
                                                    output_dir=f"{output_dir}/datasets/N{i}",
                                                    filename=f"N{i}_{j}_mutations.txt")