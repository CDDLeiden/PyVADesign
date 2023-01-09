import os
import sys
import random
import pickle
import openpyxl
import numpy as np
import pandas as pd
from Bio import SeqIO
from operator import add, sub
import matplotlib.pyplot as plt
from sklearn.cluster import MeanShift
from utils import read_codon_usage, DNA_Codons


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
        print("Sequence does not start with a start codon or end with an stop codon, \
               this is very likely to result in a shifted reading frame. Please correct your input sequence.")
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
    """
    Make sure that none of the input mutations contains a unusual amino acid

    Args:
        mutations (list): List of mutations in format [WT residue][index][mutant residue], such as G432W

    Returns:
        Bool: Returns true if all amino acids are valid
    """    
    valid = 'acdefghiklmnpqrstvwy'
    if all(i[0].lower() in valid for i in mutations):
        if all(i[-1].lower() in valid for i in mutations):
            return True
    else:
        print("Mutations contain non-natural amino acids")
        sys.exit()

def check_number_input_mutations(mutations):
    if len(mutations) < idt_min_order():
        print(f"Minimum number of mutations {len(mutations)} is lower than the minimum amount of {idt_min_order()}. \
                Please make sure you have enough mutations in your input.")
        sys.exit()
    elif len(mutations) >= idt_min_order():
        return True

def read_mutations(fp: str):
    """
    Read mutations file in TXT format. 
    Each line of the file should contain one amino acid in the format [WT residue][index][mutant residue], such as G432W
    
    Args:
        fp (string): filepath of input mutations TXT

    Returns:
        list: list of mutations that were extracted from the input file
    """    
    mutations = []
    with open(fp, 'r') as f:
        content = f.readlines()
        for line in content:
            line = line.split()
            mutations.append(line[0])
    # (1) Check there are NO non-natural amino acids in the mutations
    # (2) Check that there are enough mutations to process
    if (check_input_mutations) and (check_number_input_mutations):  
        return mutations

def translate_sequence(dna_seq):    
    """
    Translate DNA sequence to protein sequence
    """
    return dna_seq.translate()

def make_bins(clusters: dict):
    """
    Create bins, based on optimal mutation clusters

    Args:
        clusters (dict): dictionary of the format d['cluster X'] = [list of mutation indexes belonging to cluster]

    Returns:
        list: list of bins
    """    
    bins = []
    for _, value in clusters.items():
        bins.append(min(value) - min_bin_overlap())
        bins.append(max(value) + min_bin_overlap())
    bins.sort()
    return bins

def make_histogram(data, outpath, bins, labels, fname="hist.png"):
    # TODO ADD LABELS TO HISTOGRAM
    # TODO ADD COUNTS TO HISTOGRAM
    outname = os.path.join(outpath, fname)
    plt.hist(data, bins=bins)
    # plt.xticks(labels)
    plt.savefig(outname)

def calculate_cost(clusters: dict) -> float:
    """
    Calculate the total cost of all fragments, based on clusters

    Args:
        clusters (dict): dictionary of the format d['cluster X'] = [list of mutation indexes belonging to cluster]

    Returns:
        float: cost in euros
    """    
    total_cost = 0
    for _, value in clusters.items():
        min_val = min(value)
        max_val = max(value)
        len_gene_block = (max_val - min_val) + 2 * min_bin_overlap()  # on both size of the gene block there should be a number of non-mutated basepairs for IVA primer design
        cost = len_gene_block * 0.05 * len(value)  # 0.05 euros per base pair
        total_cost += cost
    return round(total_cost, 2)

def check_fragment_sizes(clusters, bandwidth):
    """
    Check that the size of the fragments in the clusters are within bounds

    Args:
        clusters (dict): clusters (dict): dictionary of the format d['cluster X'] = [list of mutation indexes belonging to cluster]
        bandwidth (int): sklearn.cluster.estimate_bandwidth

    Returns:
        bandwidth (int): new bandwith value based on size of gene blocks
    """    
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
    
def optimize_bins(x, bandwidth=200, num_iterations=1000):
    """
    Optimize the bins using a meanshift algorithm

    Args:
        x (list): indexes of mutations

    Returns:
        optimal_bandwidth (int): 
        lowest_cost (float): estimated costs of all gene blocks together
    """    

    lowest_cost = np.inf
    optimal_bandwidth = bandwidth
    
    for i in range(num_iterations):

        # Start with default value
        clusters = meanshift(x, bandwidth)

        # Calculate costs and check size of fragments
        cost = calculate_cost(clusters)
        new_bandwidth = check_fragment_sizes(clusters, bandwidth)
        
        if bandwidth == new_bandwidth:
            if lowest_cost > cost:
                lowest_cost = cost
                optimal_bandwidth = new_bandwidth
            
            ops = (add, sub)
            operation = random.choice(ops)
            random_int = random.randint(1, 50)
            bandwidth = operation(bandwidth, random_int)
        else:
            bandwidth = new_bandwidth

    return optimal_bandwidth, lowest_cost


def meanshift(x: list, bandwidth: int):
    """
    Meanshift algorithm for finding clusters of mutations that fit in a gene block

    Args:
        x (list): _description_
        bandwidth (int): _description_

    Returns:
        clusters (dict): _description_
    """    
    # https://stackoverflow.com/questions/18364026/clustering-values-by-their-proximity-in-python-machine-learning
    X = np.array(list(zip(x, np.zeros(len(x)))), dtype=np.int64)
    bandwidth = bandwidth
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    ms.fit(X)
    labels = ms.labels_
    labels_unique = np.unique(labels)
    clusters = {}
    for label in labels_unique:
        clusters[f'cluster {label}'] = []
    for num, i in enumerate(labels):
        clusters[f'cluster {i}'].append(x[num])
    return clusters

def name_block(num, bins):
    return f'Block_{num}_pos_{bins[num]}_{bins[num+1]}'

def short_name(name):
    short_name = '_'.join(name.split('_')[0:2])
    return short_name

def make_gene_block(bins, dna_sequence):
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
    
def check_position_gene_block(gene_block, idx_mutation):
    begin_range, end_range = gene_block_range(gene_block)
    if (idx_mutation - begin_range) > idt_min_length_fragment():
        print("Mutation is too close to beginning of gene block")
        sys.exit()
    elif (end_range - idx_mutation) < idt_min_length_fragment():
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
                             fname="gene_blocks.txt"):
    header = ['mutation', 'gene block name', 'length gene block', 'gene block sequence', 'index mutation', 'mut codon']
    outfile = os.path.join(outpath, fname)
    with open(outfile, 'w+') as out:
        out.write('\t'.join(header) + '\n')
        for key, value in gene_block_dict.items():
            len_gene_block = length_gene_block(value[1])
            out.write(key + '\t' + value[0] + '\t' + str(len_gene_block) + '\t' + value[1] + '\t' + str(value[2]) + '\t' + value[3] + '\n')

def write_gene_blocks_to_template(gene_block_dict, outpath, fname="eblocks-plate-upload-template-96-filled.xlsx", template='data\eblocks-plate-upload-template-96.xlsx'):
    outfile = os.path.join(outpath, fname)
    df = pd.read_excel(template)
    names = []
    seqs = []
    for key, value in gene_block_dict.items():
        mutation = key
        block = value[0].split('_')[0:1]
        name = mutation + '_' + block
        names.append(name)
        seqs.append(value[1])
    df['Name'] = names
    df['Sequence'] = seqs
    df.to_excel(outfile, index=False)

def write_pickle(obj,
                 outpath,
                 fname="mut_gene_blocks.npy"):
    with open(os.path.join(outpath, fname), 'wb') as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
def length_gene_block(gene_block):
    return len(gene_block)

def min_bin_overlap():
    return 25

def idt_max_length_fragment():
    return 1500

def idt_min_length_fragment():
    return 300

def idt_min_order():
    return 24

def main(args):

    # Read input DNA sequence
    dna_seq = read_seq(args.input_gene)

    # Read mutations
    mutations = read_mutations(args.mutations)

    # Find indexes in sequence where mutation occures
    idx_dna, _ = index_mutations(mutations)

    # Optimize bin size and bin position using meanshift
    print("Optimizing bin sizes ..")
    optimal_bandwidth, lowest_cost = optimize_bins(idx_dna)

    clusters = meanshift(idx_dna, optimal_bandwidth)
    print("Lowest cost: ", str(lowest_cost), f"with {len(clusters)} clusters")

    bins = make_bins(clusters)

    # Make gene blocks
    gene_blocks = make_gene_block(bins, dna_seq)
    write_pickle(gene_blocks, args.output_location, fname="wt_gene_blocks.npy")

    # Make histogram with bins
    labels = gene_blocks.keys()
    make_histogram(idx_dna, args.output_location, bins, labels)

    results = {}
    for num, mut in enumerate(mutations):
        
        mut_codons = extract_mut_codons(mut)

        # Find most occuring mutant codon
        mut_codon = select_mut_codon(mut_codons)

        # Find gene block
        mut_idx = idx_dna[num]
        mut_gene_block_name, mut_gene_block_value = find_gene_block(gene_blocks, mut_idx)

        # Mutate gene block
        idx = find_mutation_index_in_gene_block(mut_gene_block_name, mut_idx)
        mut_gene_block = mutate_gene_block(mut_codon, idx, mut_gene_block_value)

        # Store output in dictionary
        results[mut] = [mut_gene_block_name, mut_gene_block, idx, mut_codon]
        
    # Store output
    write_gene_blocks_to_txt(results, args.output_location)
    # write_gene_blocks_to_template(results, args.output_location)
    write_pickle(results, args.output_location)
