import os
import re
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from utils import DNA_Codons
import matplotlib.pyplot as plt
from utils import read_codon_usage

# TODO FIND OPTIMAL NUMBER OF BLOCKS (AS FEW AS POSSIBLE WHILE COSTS ARE OKAY)

def check_type_input_sequence(sequence):
    """
    Check that the type of sequence is DNA

    Args:
        sequence (str): sequence

    Returns:
        Bool: _description_
    """    
    valid = 'actg'
    if all(i.lower() in valid for i in sequence):
        return True
    else:
        print("Please provide a DNA sequence")
        sys.exit()

def check_for_start_stop_codon(sequence):
    # Make sure that the sequence starts with a start codon and ends with a stop codon
    if (sequence[0:3] == "ATG".lower()) and ((sequence[-3:] == "TAA".lower()) | (sequence[-3:] == "TAG".lower()) | (sequence[-3:] =="TGA".lower())):
        return True
    else:
        print("Sequence does not start with a start codon or end with an stop codon")
        sys.exit()

def read_seq(fp):
    """
    Read DNA sequence in FASTA format using biopython

    Args:
        fp (str):   Location of the FASTA file

    Returns:
        sequence (str): DNA sequence
    """
    count = 0
    for record in SeqIO.parse(fp, "fasta"):
        sequence = record.seq
        count += 1
    result_type = check_type_input_sequence(sequence)  # Check that the sequence if of type DNA
    result_start_stop = check_for_start_stop_codon(sequence)
    if (count == 1) and (result_type) and (result_start_stop):      
        return sequence
    else:
        print("Check input sequence")


def index_mutations(mut_list):
    """_summary_

    Args:
        mut_list (list): _description_

    Returns:
        _type_: _description_
    """    
    idx_dna, idx_aa = [], []
    for mut in mut_list:
        idx_aa.append(int(mut[1:-1]))
        idx_dna.append(int(mut[1:-1]) * 3)  # A residue consists of 3 nucleotides
    return idx_dna, idx_aa

def check_input_mutations(mutations):
    valid = 'acdefghiklmnpqrstvwy'
    if all(i[0].lower() in valid for i in mutations):
        if all(i[-1].lower() in valid for i in mutations):
            return True
    else:
        print("Mutations contain non-natural amino acids")
        sys.exit()

def read_mutations(fp):
    mutations = []
    with open(fp, 'r') as f:
        content = f.readlines()
        for line in content:
            line = line.split()
            mutations.append(line[0])
    if check_input_mutations:  # Make sure there are NO non-natural amino acids in the mutations
        return mutations

def translate_sequence(dna_seq):    
    return dna_seq.translate()

def make_bins(data, binwidth=280):
    bins=np.arange(min(data) -20, max(data) + binwidth, binwidth)
    return bins

def make_histogram(data, outpath, fname="hist.png"):
    # TODO: ADD NUMBER OF INSTANCES PER BIN
    outname = os.path.join(outpath, fname)
    plt.hist(data, bins=make_bins(data))
    plt.savefig(outname)

def name_block(num, bins):
    return f'Block_{num}_pos_{bins[num]}_{bins[num+1]}'

def make_gene_block(bins, dna_sequence):
    # TODO: Make sure that a block starts in the correct reading frame
    gene_blocks = {}
    for num in range(len(bins) - 1):
        name = name_block(num, bins)
        block = dna_sequence[bins[num]:bins[num+1]]
        gene_blocks[name] = str(block)
    return gene_blocks

def map_codons_aas(protein_sequence, dna_sequence):
    residues_codons = {}
    r = 0
    ib = 0
    it = 3  # Amino acid consists of 3 nucleotides
    while r < len(protein_sequence):
        residues_codons[protein_sequence[r] + f'{r + 1}'] = str(dna_sequence[ib:it])
        r += 1
        ib += 3
        it += 3
    return residues_codons

def extract_wt_codon(mutation: str, mapped_residues: dict):
    original_codon = mapped_residues[mutation[:-1]]
    return original_codon

def select_mut_codon(codon_list):
    """
    Choose codon from list of codons based on occurance of codon in nature.
    """
    codon_dict = read_codon_usage()
    highest_freq = 0
    most_occuring_codon = 'xxx'
    for codon in codon_list:
        codon_freq = codon_dict[codon]
        if codon_freq > highest_freq:
            highest_freq = codon_freq
            most_occuring_codon = codon
    return most_occuring_codon

def extract_mut_codons(mutation: str):
    mut_codons = []
    mut_residue = mutation[-1]
    for key, value in DNA_Codons.items():  # d[codon] = amino acid
        if value == mut_residue:
            mut_codons.append(key.lower())
    return mut_codons

def gene_block_range(gene_block_name):
    begin_range = int(gene_block_name.split('_')[3])
    end_range = int(gene_block_name.split('_')[4])
    return begin_range, end_range

def find_gene_block(gene_blocks, mutation_idx):
    for key, value in gene_blocks.items():
        begin_range, end_range = gene_block_range(key)
        if begin_range < int(mutation_idx) < end_range:
            result = (key, value)
            return result
    
def check_position_gene_block(gene_block, idx_mutation, min_nonmutated=15):
    # TODO: Think of way to fix this
    begin_range, end_range = gene_block_range(gene_block)
    if (idx_mutation - begin_range) > min_nonmutated:
        print("Mutation is too close to beginning of gene block")
        sys.exit()
    elif (end_range - idx_mutation) < min_nonmutated:
        print("Mutation is too close to final part of gene block")
        sys.exit()
    
def find_mutation_index_in_gene_block(gene_block, idx_mutation):
    begin_range, _ = gene_block_range(gene_block)
    # Find index of mutation within geneblock
    index = idx_mutation - begin_range
    # Check that mutation is not in the first or final 15 residues of the gene block (this would make it very difficult to design IVA primers)
    check_position_gene_block(gene_block, index)
    return index

def mutate_gene_block(mut_codon, mut_index, gene_block_seq):
    # Change this codon in the gene_block
    mut_block = gene_block_seq[:mut_index -3] + mut_codon + gene_block_seq[mut_index:]
    return mut_block

def write_gene_blocks_to_txt(gene_block_dict, 
                             outpath, 
                             fname="code_blocks.txt"):
    header = ['mutation', 'gene block name', 'length gene block', 'gene block sequence', 'index mutation', 'mut codon']
    outfile = os.path.join(outpath, fname)
    with open(outfile, 'w+') as out:
        out.write('\t'.join(header) + '\n')
        for key, value in gene_block_dict.items():
            len_gene_block = length_gene_block(value[1])
            out.write(key + '\t' + value[0] + '\t' + str(len_gene_block) + '\t' + value[1] + '\t' + str(value[2]) + '\t' + value[3] + '\n')

def length_gene_block(gene_block):
    return len(gene_block)


if __name__ == "__main__":
    
    gene = r"example_data\mtb_DnaE1_seq.txt"
    mutations_fp = r"example_data\mutations.txt"
    outpath = r"C:\Users\Rosan\Documents\git\my_repositories\design_gene_blocks\example_output"

    # Read input DNA sequence
    dna_seq = read_seq(gene)
    print(dna_seq)

    # Convert DNA sequence to protein sequence
    prot_seq = translate_sequence(dna_seq)
    print(prot_seq)

    # Read mutations
    mutations = read_mutations(mutations_fp)
    print(mutations)

    # Find indexes in sequence where mutation occures
    idx_dna, idx_res = index_mutations(mutations)
    print(idx_dna)

    # Make histogram with bins
    make_histogram(idx_dna, outpath)

    bins = make_bins(idx_dna)
    print(bins)

    gene_blocks = make_gene_block(bins, dna_seq)
    print(gene_blocks)
    
    residues_codons_mapped = map_codons_aas(prot_seq, dna_seq)
    print(residues_codons_mapped)

    # Make mutations
    # - loop over mutations to make
        # - find which code block it should be in
        # - mutate in block
        # - save output in dictionary

    results = {}
    for num, mut in enumerate(mutations):
        print("*********")
        print("processing mutation:", mut)
        print("*********")
        mut_res = mut[-1]
        
        # Find codon of WT residue
        wt_codon = extract_wt_codon(mut, residues_codons_mapped)
        print("original_codon:", wt_codon)

        mut_codons = extract_mut_codons(mut)
        print(mut_codons)

        # Find most occuring mutant codon
        mut_codon = select_mut_codon(mut_codons)
        print(mut_codon)

        # Find gene block
        mut_idx = idx_dna[num]
        mut_gene_block_name, mut_gene_block_value = find_gene_block(gene_blocks, mut_idx)
        print(mut_gene_block_name, mut_gene_block_value)

        # Mutate gene block
        idx = find_mutation_index_in_gene_block(mut_gene_block_name, mut_idx)
        mut_gene_block = mutate_gene_block(mut_codon, mut_idx, mut_gene_block_value)
        results[mut] = [mut_gene_block_name, mut_gene_block, idx, mut_codon]
        

    for key, value in results.items():
        print(key, value)

    write_gene_blocks_to_txt(results, outpath)