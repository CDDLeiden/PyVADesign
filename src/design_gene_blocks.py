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


class DesignEblocks:
    """
    """

    def __init__(self, sequence_fp: str, mutations_fp: str, output_fp: str):
        
        self.dna_seq = self.read_seq(sequence_fp)
        self.mutations = self.read_mutations(mutations_fp)
        self.out_fp = output_fp

        # Default settings based on IDT requirements
        self.min_bin_overlap = 25
        self.idt_max_length_fragment = 1500
        self.idt_min_length_fragment = 300
        self.idt_min_order = 24

    def run(self):

        # Find indexes in sequence where mutation occures
        idx_dna, _ = self.index_mutations()

        # Optimize bin size and bin position using meanshift
        print("Optimizing bin sizes ..")
        optimal_bandwidth, lowest_cost = self.optimize_bins(idx_dna)

        clusters = self.meanshift(idx_dna, optimal_bandwidth)
        print("Lowest cost: ", str(lowest_cost), f"with {len(clusters)} clusters")

        bins = self.make_bins(clusters)

        # Make gene blocks
        gene_blocks = self.make_gene_block(bins, self.dna_seq)
        self.write_pickle(gene_blocks, self.output_fp, fname="wt_gene_blocks.npy")

        # Make histogram with bins
        labels = gene_blocks.keys()
        self.make_histogram(idx_dna, self.output_fp, bins, labels)

        results = {}
        for num, mut in enumerate(self.mutations):
            
            mut_codons = self.extract_mut_codons(mut)

            # Find most occuring mutant codon
            mut_codon = self.select_mut_codon(mut_codons)

            # Find gene block
            mut_idx = idx_dna[num]
            mut_gene_block_name, mut_gene_block_value = self.find_gene_block(gene_blocks, mut_idx)

            # Mutate gene block
            idx = self.find_mutation_index_in_gene_block(mut_gene_block_name, mut_idx)
            mut_gene_block = self.mutate_gene_block(mut_codon, idx, mut_gene_block_value)

            # Store output in dictionary
            results[mut] = [mut_gene_block_name, mut_gene_block, idx, mut_codon]
            
        # Store output
        self.write_gene_blocks_to_txt(results, self.out_fp)
        self.write_gene_blocks_to_template(results, self.out_fp)
        self.write_pickle(results, self.out_fp)

    def check_type_input_sequence(self, sequence):
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

    def check_for_start_stop_codon(self, sequence):
        # Make sure that the sequence starts with a start codon and ends with a stop codon
        if (sequence[0:3] == "ATG".lower()) and ((sequence[-3:] == "TAA".lower()) | (sequence[-3:] == "TAG".lower()) | (sequence[-3:] =="TGA".lower())):
            return True
        else:
            print("Sequence does not start with a start codon or end with an stop codon, \
                this is very likely to result in a shifted reading frame. Please correct your input sequence.")
            sys.exit()

    def read_seq(self, fp):
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
        result_type = self.check_type_input_sequence(sequence)  # Check that the sequence is of type DNA
        result_start_stop = self.check_for_start_stop_codon(sequence)
        if (count == 1) and (result_type) and (result_start_stop):      
            return sequence
        else:
            print("Check input sequence")

    def index_mutations(self, mut_list):
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

    def check_input_mutations(self, mutations):
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

    def check_number_input_mutations(self, mutations):
        if len(mutations) < self.idt_min_order:
            print(f"Minimum number of mutations {len(mutations)} is lower than the minimum amount of {self.idt_min_order}. \
                    Please make sure you have enough mutations in your input.")
            sys.exit()
        elif len(mutations) >= self.idt_min_order:
            return True

    def read_mutations(self, fp: str):
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
        if (self.check_input_mutations) and (self.check_number_input_mutations):  
            return mutations

    def translate_sequence(self, dna_seq):    
        """
        Translate DNA sequence to protein sequence
        """
        return dna_seq.translate()

    def make_bins(self, clusters: dict):
        """
        Create bins, based on optimal mutation clusters

        Args:
            clusters (dict): dictionary of the format d['cluster X'] = [list of mutation indexes belonging to cluster]

        Returns:
            list: list of bins
        """    
        bins = []
        for _, value in clusters.items():
            bins.append(min(value) - self.min_bin_overlap)
            bins.append(max(value) + self.min_bin_overlap)
        bins.sort()
        return bins

    def make_histogram(self, data, outpath, bins, labels, fname="hist.png"):
        # TODO ADD LABELS TO HISTOGRAM
        # TODO ADD COUNTS TO HISTOGRAM
        outname = os.path.join(outpath, fname)
        plt.hist(data, bins=bins)
        # plt.xticks(labels)
        plt.savefig(outname)

    def calculate_cost(self, clusters: dict) -> float:
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
            len_gene_block = (max_val - min_val) + 2 * self.min_bin_overlap  # on both size of the gene block there should be a number of non-mutated basepairs for IVA primer design
            cost = len_gene_block * 0.05 * len(value)  # 0.05 euros per base pair
            total_cost += cost
        return round(total_cost, 2)

    def check_fragment_sizes(self, clusters, bandwidth):
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
            len_gene_block = (max_val - min_val) + 2 * self.min_bin_overlap
            
            # size of gene block is too small > increasing bandwidth
            if len_gene_block < self.idt_min_length_fragment:
                newbandwidth = bandwidth + 1
                return newbandwidth

            # size of gene block is too large > decreasing bandwidth
            elif len_gene_block > self.idt_max_length_fragment:
                newbandwidth = bandwidth - 1
                return newbandwidth 
            else:
                continue

        return bandwidth
        
    def optimize_bins(self, x, bandwidth=200, num_iterations=1000):
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
            clusters = self.meanshift(x, bandwidth)

            # Calculate costs and check size of fragments
            cost = self.calculate_cost(clusters)
            new_bandwidth = self.check_fragment_sizes(clusters, bandwidth)
            
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


    def meanshift(self, x: list, bandwidth: int):
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

    def name_block(self, num, bins):
        return f'Block_{num}_pos_{bins[num]}_{bins[num+1]}'

    def short_name(self, name):
        short_name = '_'.join(name.split('_')[0:2])
        return short_name

    def make_gene_block(self, bins, dna_sequence):
        gene_blocks = {}
        for num in range(len(bins) - 1):
            name = self.name_block(num, bins)
            block = dna_sequence[bins[num]:bins[num+1]]
            gene_blocks[name] = str(block)
        return gene_blocks

    def map_codons_aas(self, protein_sequence, dna_sequence):
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

    def extract_wt_codon(self, mutation: str, mapped_residues: dict):
        original_codon = mapped_residues[mutation[:-1]]
        return original_codon

    def select_mut_codon(self, codon_list):
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

    def extract_mut_codons(self, mutation: str):
        mut_codons = []
        mut_residue = mutation[-1]
        for key, value in DNA_Codons.items():  # d[codon] = amino acid
            if value == mut_residue:
                mut_codons.append(key.lower())
        return mut_codons

    def gene_block_range(self, gene_block_name):
        begin_range = int(gene_block_name.split('_')[3])
        end_range = int(gene_block_name.split('_')[4])
        return begin_range, end_range

    def find_gene_block(self, gene_blocks, mutation_idx):
        for key, value in gene_blocks.items():
            begin_range, end_range = self.gene_block_range(key)
            if begin_range < int(mutation_idx) < end_range:
                result = (key, value)
                return result
        
    def check_position_gene_block(self, gene_block, idx_mutation):
        begin_range, end_range = self.gene_block_range(gene_block)
        if (idx_mutation - begin_range) > self.idt_min_length_fragment:
            print("Mutation is too close to beginning of gene block")
            sys.exit()
        elif (end_range - idx_mutation) < self.idt_min_length_fragment:
            print("Mutation is too close to final part of gene block")
            sys.exit()
        
    def find_mutation_index_in_gene_block(self, gene_block, idx_mutation):
        begin_range, _ = self.gene_block_range(gene_block)
        # Find index of mutation within geneblock
        index = idx_mutation - begin_range
        # Check that mutation is not in the first or final 15 residues of the gene block (this would make it very difficult to design IVA primers)
        self.check_position_gene_block(gene_block, index)
        return index

    def mutate_gene_block(self, mut_codon, mut_index, gene_block_seq):
        # Change this codon in the gene_block
        mut_block = gene_block_seq[:mut_index -3] + mut_codon + gene_block_seq[mut_index:]
        return mut_block

    def write_gene_blocks_to_txt(self, gene_block_dict, 
                                outpath, 
                                fname="gene_blocks.txt"):
        header = ['mutation', 'gene block name', 'length gene block', 'gene block sequence', 'index mutation', 'mut codon']
        outfile = os.path.join(outpath, fname)
        with open(outfile, 'w+') as out:
            out.write('\t'.join(header) + '\n')
            for key, value in gene_block_dict.items():
                len_gene_block = self.length_gene_block(value[1])
                out.write(key + '\t' + value[0] + '\t' + str(len_gene_block) + '\t' + value[1] + '\t' + str(value[2]) + '\t' + value[3] + '\n')

    def extract_wells_from_template(self, template='data\eblocks-plate-upload-template-96.xlsx'):
        df = pd.read_excel(template)
        wells = df['Well Position'].tolist()
        return wells

    def write_gene_blocks_to_template(self, gene_block_dict, outpath, fname="eblocks-plate-upload-template-96-filled.xlsx"):
        outfile = os.path.join(outpath, fname)
        names = []
        seqs = []
        for key, value in gene_block_dict.items():
            mutation = key
            block = value[0].split('_')[0:2]
            block = block[0] + '-' + block[1]
            name = mutation + '_' + block
            names.append(name)
            seqs.append(value[1])
        wells = self.extract_wells_from_template()
        wells = wells[:len(names)]
        df = pd.DataFrame()
        df['Well Position'] = wells
        df['Name'] = names
        df['Sequence'] = seqs
        df.to_excel(outfile, index=False)

    def write_pickle(self, obj,
                    outpath,
                    fname="mut_gene_blocks.npy"):
        with open(os.path.join(outpath, fname), 'wb') as handle:
            pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
    def length_gene_block(self, gene_block):
        return len(gene_block)
