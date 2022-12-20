import os
import sys
import random
import argparse
import numpy as np
from Bio import SeqIO
from utils import DNA_Codons
from random import randrange
from operator import add, sub
import matplotlib.pyplot as plt
from utils import read_codon_usage
from scipy.signal import argrelextrema
from sklearn.cluster import MeanShift, estimate_bandwidth


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
    bins=np.arange(min(data) - min_bin_overlap(), max(data) + min_bin_overlap(), binwidth)
    print(bins)
    return bins

def make_bins_from_cluster(clusters):
    bins = []
    for key, value in clusters.items():
        bins.append(min(value) - min_bin_overlap())
        bins.append(max(value) + min_bin_overlap())
    bins.sort()
    return bins

def make_histogram(data, outpath, bins, fname="hist.png"):
    # TODO: ADD NUMBER OF INSTANCES PER BIN
    # TODO: ADD NAME OF GENE BLOCK
    outname = os.path.join(outpath, fname)
    plt.hist(data, bins=bins)
    plt.savefig(outname)

def calculate_max_number_of_bases(max_price: int, base_price=0.055):
    num_bases = (max_price) // (base_price)
    return num_bases

def calculate_max_length_fragment(num_bases, num_mutations):
    max_fragment_length = (num_bases) // (num_mutations)
    if (idt_min_length_fragment() < max_fragment_length) and (idt_max_length_fragment() > max_fragment_length):
        return int(max_fragment_length)
    elif (idt_min_length_fragment() > max_fragment_length):
        return idt_min_length_fragment()
    elif (idt_max_length_fragment() < max_fragment_length):
        return idt_max_length_fragment()

def calculate_cost(clusters):
    total_cost = 0
    for key, value in clusters.items():
        min_val = min(value)
        max_val = max(value)
        len_gene_block = (max_val - min_val) + 2 * min_bin_overlap()
        cost = len_gene_block * 0.05 * len(value)
        total_cost += cost
    return round(total_cost, 2)

def check_fragment_sizes(clusters, bandwidth):
    for _, value in clusters.items():
        
        min_val = min(value)
        max_val = max(value)
        len_gene_block = (max_val - min_val) + 2 * min_bin_overlap()
        
        # size of gene block is too small > increasing bandwidth
        if len_gene_block < idt_min_length_fragment():
            newbandwidth = bandwidth + 1
            return newbandwidth
        # size of gene block is too large > decreasing bandwidth
        elif len_gene_block > idt_max_length_fragment():
            newbandwidth = bandwidth - 1
            return newbandwidth 
        else:
            continue
    
    return bandwidth
    
def optimize_bins(x):

    lowest_cost = 10000
    num_iterations = 200
    optimal_bandwidth = 0
    bandwidth = 200

    # TODO CHANGE SO THAT ALGORITHM STOPS WHEN NO IMPROVEMENT AFTER X ITERATIONS
    for i in range(num_iterations):

        # Start with default value
        clusters = meanshift(x, bandwidth)

        # Calculate costs and check size of fragments
        cost = calculate_cost(clusters)
        new_bandwidth = check_fragment_sizes(clusters, bandwidth)
        
        print(bandwidth, new_bandwidth)

        if bandwidth == new_bandwidth:
            if lowest_cost > cost:
                lowest_cost = cost
                optimal_bandwidth = new_bandwidth

                print(lowest_cost, optimal_bandwidth)
            
            ops = (add, sub)
            operation = random.choice(ops)
            random_int = random.randint(1, 15)
            bandwidth = operation(bandwidth, random_int)
        else:
            bandwidth = new_bandwidth

    return optimal_bandwidth, lowest_cost


def meanshift(x: list, bandwidth: int):
    
    # https://stackoverflow.com/questions/18364026/clustering-values-by-their-proximity-in-python-machine-learning
    X = np.array(list(zip(x, np.zeros(len(x)))), dtype=np.int64)

    bandwidth = bandwidth # estimate_bandwidth(X, quantile=0.3)
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    ms.fit(X)
    labels = ms.labels_
    cluster_centers = ms.cluster_centers_

    labels_unique = np.unique(labels)
    n_clusters_ = len(labels_unique)

    clusters = {}
    for label in labels_unique:
        clusters[f'cluster {label}'] = []
    for num, i in enumerate(labels):
        clusters[f'cluster {i}'].append(x[num])

    return clusters

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

def min_bin_overlap():
    return 15

def idt_max_length_fragment():
    return 1500

def idt_min_length_fragment():
    return 300

def read_arguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input_gene", help="FASTA file containing the gene of interest")
    parser.add_argument("-m", "--mutations", help="TXT file containing the mutations to make")
    parser.add_argument("-o", "--output_location", help="Location where to store the output of the script")
    args = parser.parse_args()
    return args

def main(args):

    # Read input DNA sequence
    dna_seq = read_seq(args.input_gene)

    # Convert DNA sequence to protein sequence
    prot_seq = translate_sequence(dna_seq)

    # Read mutations
    mutations = read_mutations(args.mutations)

    # Find indexes in sequence where mutation occures
    idx_dna, _ = index_mutations(mutations)

    # Optimize bins
    print("Optimizing bin sizes ..")
    optimal_bandwidth, lowest_cost = optimize_bins(idx_dna)
    print(optimal_bandwidth, lowest_cost)

    clusters = meanshift(idx_dna, optimal_bandwidth)
    print(clusters)

    # Make bins (TODO: OPTIMIZE THIS)
    # bins = make_bins(idx_dna, binwidth=350)
    # print(bins)

    bins = make_bins_from_cluster(clusters)
    print(bins)

    # Make histogram with bins
    make_histogram(idx_dna, args.output_location, bins)

    # Make gene blocks
    gene_blocks = make_gene_block(bins, dna_seq)
    residues_codons_mapped = map_codons_aas(prot_seq, dna_seq)

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
        # mut_res = mut[-1]
        
        # Find codon of WT residue
        # wt_codon = extract_wt_codon(mut, residues_codons_mapped)
        # print("original_codon:", wt_codon)

        mut_codons = extract_mut_codons(mut)
        # print(mut_codons)

        # Find most occuring mutant codon
        mut_codon = select_mut_codon(mut_codons)
        # print(mut_codon)

        # Find gene block
        mut_idx = idx_dna[num]
        mut_gene_block_name, mut_gene_block_value = find_gene_block(gene_blocks, mut_idx)
        # print(mut_gene_block_name, mut_gene_block_value)

        # Mutate gene block
        idx = find_mutation_index_in_gene_block(mut_gene_block_name, mut_idx)
        mut_gene_block = mutate_gene_block(mut_codon, mut_idx, mut_gene_block_value)

        # Store output in dictionary
        results[mut] = [mut_gene_block_name, mut_gene_block, idx, mut_codon]
        
    # for key, value in results.items():
    #     print(key, value)

    write_gene_blocks_to_txt(results, args.output_location)
    print(f"Estimated costs are {lowest_cost} euros plus the costs of the primers")

    
if __name__ == "__main__":

    arguments = read_arguments()
    main(arguments)
    print("Finished without any problems")
    
    # input_gene = r"example_data\mtb_DnaE1_seq.txt"
    # mutations = r"example_data\mutations.txt"
    # output_location = r"C:\Users\Rosan\Documents\git\my_repositories\design_gene_blocks\example_output"