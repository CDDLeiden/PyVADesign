import os
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from utils import read_codon_usage, DNA_Codons
from design_gene_blocks import DesignEblocks

# Read DNA sequence
# Read mutations
# Calculate number of mutations needed to get from DNA sequence to mutated DNA sequence using codon usage Mtb (and how often do these codons occur in Mtb)
    # For each mutation:
        # Find codon in DNA sequence (index)
        # Find codon in mutated DNA sequence (index)
        # Find codon usage in Mtb (index)
        # Calculate number of mutations needed to get from codon to codon usage in Mtb
        # Calculate how often codon usage in Mtb occurs in Mtb

def read_fasta(fp):
    count = 0
    for record in SeqIO.parse(fp, "fasta"):
        seq = record.seq
        count += 1
    if count == 1:
        return seq
    else:
        print("More than one sequence in fasta file")
        sys.exit()

def main(dna_sequence_fp, mutations_list_fp, prot_sequence_fp):

    eblocks = DesignEblocks(sequence_fp=dna_sequence_fp, 
                            mutations_fp=mutations_list_fp, # Not needed
                            output_fp=None, # Not needed
                            species="Mycobacterium Smegmatis")  # Not needed
    
    # Read DNA/protein sequence
    dna_sequence = eblocks.read_seq(dna_sequence_fp)
    prot_sequence = read_fasta(prot_sequence_fp)
    num_mutations = 0

    # Find number of mutations needed to get from DNA sequence to mutated DNA
    # Find number of mutations
    count = 1
    i = 3  # Skip start codon
    mutation_idexes = {}  # Mutation_index = {G1: [TAA], ...} (key:mutation, value:WT codon)
    while True:
        codon = dna_sequence[i:i+3]
        # Break if stop codon or end of sequence
        if (codon == "TAA" or codon == "TAG" or codon == "TGA") or (count == len(prot_sequence)):
            break
        if DNA_Codons[codon.upper()] != prot_sequence[count]:
            mutation_idexes[prot_sequence[count] + str(count)] = [str(codon.upper())]
            num_mutations += 1
        i += 3
        count += 1
    
    print(num_mutations)
    print(mutation_idexes)

    # Find number of mutations needed to get from DNA sequence to mutated DNA
    for mutation, wt_codon in mutation_idexes.items():
        
        # Find all codons that encode the same amino acid as the mutation
        mut_codons = []
        for k, v in DNA_Codons.items():
            if v == mutation[0]:
                mut_codons.append(k)
        
        # Find out how many mutations are needed to get from wt_codon to mut_codon
        mut_codons_operations_needed = []
        for mut_codon in mut_codons:
            num_mutations = 0
            for i in range(3):
                if mut_codon[i] != wt_codon[0][i]:
                    num_mutations += 1
            mut_codons_operations_needed.append(num_mutations)

        # Minimum number of point mutations in the DNA needed to change the amino acid
        # TODO USE THIS TO INSPECT RESISTANCE MUTATIONS IN MMPL3    
        min_number_of_mutations = min(mut_codons_operations_needed)

        print(mutation, wt_codon, mut_codons, min_number_of_mutations)



        

if __name__ == '__main__':

    dna_sequence_fp = r"C:\Users\Rosan\Documents\projects\two-entopy-analysis\data\Mtb_DnaE1_DNA.txt"
    mutations_list_fp = r"C:\Users\Rosan\Documents\git\my_repositories\design_gene_blocks\tests\T1\T1_in\mutations_ed.txt"
    prot_sequence_fp = r"C:\Users\Rosan\Documents\projects\two-entopy-analysis\data\Mtb_DnaE1_prot_mut.txt"
    main(dna_sequence_fp, mutations_list_fp, prot_sequence_fp)