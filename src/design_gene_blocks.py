import os
import sys
import random
# import openpyxl
import numpy as np
import pandas as pd
from Bio import SeqIO
from operator import add, sub
import matplotlib.pyplot as plt
from sklearn.cluster import MeanShift
from utils import read_codon_usage, DNA_Codons, write_pickle


class DesignEblocks:
    """
    Class to design eBlocks
    """

    def __init__(self, 
                 sequence_fp: str, 
                 mutations_fp: str,
                 output_fp: str,  # Path to store output files
                 species: str,
                 min_bin_overlap = 25,
                 idt_max_length_fragment = 1500,
                 idt_min_length_fragment = 300,
                 idt_min_order = 24):
        
        # Mutation types supported
        self.type_mutation = "Mutation"
        self.type_insert = "Insert"
        self.type_deletion = "Deletion"
        self.type_combined = "Combined"

        self.output_fp = output_fp
        self.species = species
        self.min_bin_overlap = min_bin_overlap
        self.idt_max_length_fragment = idt_max_length_fragment
        self.idt_min_length_fragment = idt_min_length_fragment
        self.idt_min_order = idt_min_order
        self.n_iterations = 100  # Number of iterations to run mean shift clustering
        self.gene_blocks = None
        
        self.dna_seq = self.read_seq(sequence_fp)
        self.mutations, self.mutation_types = self.read_mutations(mutations_fp)  # Types > Mutation, Insert, Deletion
        self.codon_usage = self.check_existance_codon_usage_table()  # Check if codon usage table is present for selected species

    def run(self):

        # TODO Describe somewhere how this should be formatted in input and describe which codon usage tables are present at the moment
        # TODO Add an option to include own codon usage table in CLI
        # TODO Also add option for default codon usage table?

        # Find indexes in sequence where mutation occures
        idx_dna, idx_dna_tups = self.index_mutations(self.mutations, self.mutation_types)

        # Optimize bin size and bin position using meanshift
        print("Optimizing bin sizes ...")
        optimal_bandwidth, lowest_cost = self.optimize_bins(idx_dna)

        clusters = self.meanshift(idx_dna, optimal_bandwidth)
        print(clusters)
        print("Lowest cost: ", str(lowest_cost), f"with {len(clusters)} clusters")

        # Write expected cost to file
        with open(os.path.join(self.output_fp, "expected_cost.txt"), "w") as f:
            f.write(f"expected costs: {str(lowest_cost)}" + '\n')
            f.write(f"number of cluters: {str(len(clusters))}" + '\n')

        bins = self.make_bins(clusters)
        
        # Make gene blocks (WT DNA sequences cut to the correct size)
        gene_blocks = self.make_gene_block(bins, self.dna_seq)
        self.gene_blocks = gene_blocks

        # Make histogram with bins
        labels = gene_blocks.keys()
        self.make_histogram(idx_dna, self.output_fp, bins, labels)

        # TODO Maybe first check for all gene blocks whether the mutation is in a correct place, e.g. not the beginning or end of a gene block
        # TODO Restart option is error was found?

        results = {}
        should_restart = True
        while should_restart:
            should_restart = False
            for num, mut in enumerate(self.mutations):

                mut_type = self.mutation_types[num]
                
                print(mut, mut_type)

                # Change gene block in selected position for mutation
                if mut_type == self.type_mutation:

                    # Find gene block and index of insert/deletion/mutation
                    mut_idx = idx_dna_tups[num][1]
                    mut_gene_block_name, mut_gene_block_value = self.find_gene_block(gene_blocks, mut_idx)
                    idx = self.find_mutation_index_in_gene_block(mut_gene_block_name, mut_idx)

                    # Check if WT codon at index is same residue as mutation
                    self.check_wt_codon(mut_gene_block_value, idx, mut[0])
                    
                    mut_codons = self.extract_mut_codons(mut[-1])
                    # Find most occuring mutant codon based on codon usage for species
                    mut_codon = self.select_mut_codon(mut_codons)

                    # Mutate gene block
                    mut_gene_block = self.mutate_gene_block(mut_codon, idx, mut_gene_block_value, mut_type)

                    # Store output in dictionary
                    results[mut] = [mut_gene_block_name, mut_gene_block, idx, mut_codon, mut_type]

                elif mut_type == self.type_insert:

                    res_insert = mut.split('-')[1]
                    codon_insert = ''  # Sequence to insert in gene block
                    for res in res_insert:
                        codons = self.extract_mut_codons(res)
                        codon = self.select_mut_codon(codons)
                        codon_insert += codon

                    # Find gene block and index of insert/deletion/mutation
                    mut_idx = idx_dna_tups[num][1]
                    mut_gene_block_name, mut_gene_block_value = self.find_gene_block(gene_blocks, mut_idx)
                    idx = self.find_mutation_index_in_gene_block(mut_gene_block_name, mut_idx)

                    # Check if WT codon at index is same residue as mutation
                    self.check_wt_codon(mut_gene_block_value, idx, mut[0])

                    # Mutate gene block
                    mut_gene_block = self.mutate_gene_block(codon_insert, idx, mut_gene_block_value, mut_type)

                    # Check if eBlock is too long / too short
                    if self.check_eblock_length(mut_gene_block):
                        # Store output in dictionary
                        results[mut] = [mut_gene_block_name, mut_gene_block, idx, codon_insert, mut_type]

                elif mut_type == self.type_deletion:

                    idx_del_start = idx_dna_tups[num][1]
                    idx_del_end = int(mut.split('-')[1][1:]) * 3
                    print(idx_del_start, idx_del_end)
                    
                    mut_gene_block_name, mut_gene_block_value = self.find_gene_block(gene_blocks, idx_del_start)
                    idx = self.find_mutation_index_in_gene_block(mut_gene_block_name, idx_del_start)
                    idx_end = self.find_mutation_index_in_gene_block(mut_gene_block_name, idx_del_end)

                    # Check if WT codon at index is same residue as mutation
                    self.check_wt_codon(mut_gene_block_value, idx, mut[0])

                    # Mutate gene block
                    mut_gene_block = self.mutate_gene_block('', idx, mut_gene_block_value, mut_type, idx_end)

                    # Check if eBlock is too long / too short
                    if self.check_eblock_length(mut_gene_block):
                        # Store output in dictionary
                        results[mut] = [mut_gene_block_name, mut_gene_block, idx, '', mut_type]  # TODO Maybe do it based on a string to remove 

                elif mut_type == self.type_combined:

                    # TODO Check if the different mutations are in the same eblock, otherwise throw error

                    # Find gene block and index of insert/deletion/mutation
                    mut_idx = idx_dna_tups[num][0][1]
                    print(mut_idx)

                    mut_gene_block_name, mut_gene_block_value = self.find_gene_block(gene_blocks, mut_idx)
                    print(mut_gene_block_name, mut_gene_block_value)

                    idxs = []
                    codons = []
                    
                    for mut_i in idx_dna_tups[num]:

                        idx = self.find_mutation_index_in_gene_block(mut_gene_block_name, mut_i[1])

                        idxs.append(idx)
                        mut_codons = self.extract_mut_codons(mut_i[0][-1])
                        # Find most occuring mutant codon based on codon usage for species
                        mut_codon = self.select_mut_codon(mut_codons)
                        codons.append(mut_codon)

                        # Check if WT codon at index is same residue as mutation
                        self.check_wt_codon(mut_gene_block_value, idx, mut_i[0])

                        # Mutate gene block
                        mut_gene_block_value = self.mutate_gene_block(mut_codon, idx, mut_gene_block_value, mut_type)

                    # Store output in dictionary
                    results[mut] = [mut_gene_block_name, mut_gene_block_value, idx, mut_codon, mut_type]
            
        # Store output
        self.write_gene_blocks_to_txt(results, self.output_fp)
        self.write_gene_blocks_to_template(results, self.output_fp)
        write_pickle(results, self.output_fp)
        write_pickle(gene_blocks, self.output_fp, fname="wt_gene_blocks.npy")

    def check_wt_codon(self, gene_block_value, idx, mut):
        for key, value in DNA_Codons.items():
                if key.lower() == gene_block_value[idx-3:idx]:
                    result = value
        if not result == mut[0]:
            print(f"WT codon does not match residue {mut}, but is {result}, the codon is {gene_block_value[idx-3:idx]}")
            sys.exit()

    def check_existance_codon_usage_table(self):
        codon_usage_present = [file for file in os.listdir(r"data/codon_usage")]
        organisms_present = [file.split('.')[0:-1] for file in codon_usage_present]
        organisms_present = [i[0] for i in organisms_present]
        organisms_present_format = []
        for i in organisms_present:
            inew = i.split('_')
            inew = ' '.join([i for i in inew])
            organisms_present_format.append(inew.lower())
        if self.species.lower() in organisms_present_format:
            i = organisms_present_format.index(self.species.lower())
            return os.path.join(r"data/codon_usage", codon_usage_present[i])
        else:
            print("It looks like the codon usage table for the specified organism is not present.")
            sys.exit()

    def check_eblock_length(self, mut_gene_block):
            # Make sure that the codon insert is not too long or too short for eBlock
            if len(mut_gene_block) > self.idt_max_length_fragment:
                print("Codon insert is too long for eBlock")
                sys.exit()
            elif len(mut_gene_block) < self.idt_min_length_fragment:
                print(f"Codon insert is too short for eBlock, length is {len(mut_gene_block)}")
                sys.exit()
            else:
                return True

    def index_mutations(self, mut_list, mut_types):
        """_summary_

        Args:
            mut_list (list): _description_

        Returns:
            _type_: _description_
        """    
        idx_dna, idx_dna_tups = [], []
        for mut, type in zip(mut_list, mut_types):
            if type == self.type_mutation:
                idx_dna.append(int(mut[1:-1]) * 3)  # A residue consists of 3 nucleotides
                idx_dna_tups.append([mut, int(mut[1:-1]) * 3])
            elif (type == self.type_insert) or (type == self.type_deletion):
                mut_i = mut.split('-')[0]
                print(mut_i)
                idx_dna.append(int(mut_i[1:]) * 3)  # A residue consists of 3 nucleotides
                idx_dna_tups.append([mut_i, int(mut_i[1:]) * 3])
            elif type == self.type_combined:
                muts = mut.split('-')
                tmp = []
                for i in muts:
                    idx_dna.append(int(i[1:-1]) * 3)  # A residue consists of 3 nucleotides
                    tmp.append([i, int(i[1:-1]) * 3])
                idx_dna_tups.append(tmp)
        return idx_dna, idx_dna_tups
    
    def check_input_mutations(self, mutations, mutation_types):
        """
        Make sure that none of the input mutations contains a unusual amino acid

        Args:
            mutations (list): List of mutations in format [WT residue][index][mutant residue], such as G432W
            mutation_types (list): List with type of mutations, either: Mutation, Insert or Deletion

        Returns:
            Bool: Returns true if all amino acids are valid
        """
        valid = 'acdefghiklmnpqrstvwy'
        for mut, type in zip(mutations, mutation_types):
            print(mut, type)
            if type == self.type_mutation:
                if not self.check_mut_format(mut, self.type_mutation):
                    print(f"Input {mut} contain non-natural amino acids or incorrect formatting")
                    sys.exit()
            elif type == self.type_insert:
                mut_format = mut.split('-')[0]
                added_residues = mut.split('-')[1]
                if not self.check_mut_format(mut_format, self.type_insert):
                        print(f"Input {mut} contain non-natural amino acids or incorrect formatting")
                        sys.exit()
                for i in added_residues:
                    if not i.lower() in valid:
                        print(f"Input {i} contain non-natural amino acids or incorrect formatting")
                        sys.exit()
            elif type == self.type_deletion:
                mut_start = mut.split('-')[0]
                mut_end = mut.split('-')[0]
                if not self.check_mut_format(mut_start, self.type_deletion):
                    print(f"Input {mut_start} contain non-natural amino acids or incorrect formatting")
                    sys.exit()
                if not self.check_mut_format(mut_end, self.type_deletion):
                    print(f"Input {mut_end} contain non-natural amino acids or incorrect formatting")
                    sys.exit()
            elif type == self.type_combined:
                mut_list = mut.split('-')
                for i in mut_list:
                    if not self.check_mut_format(i, self.type_combined):
                        print(f"Input {i} contain non-natural amino acids or incorrect formatting")
                        sys.exit()
            else:
                print("Input contains non standard mutation")
                sys.exit()
        return True

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
            list: list of mutations that were extracted from the input file as well as the types of mutation (mutation, insert, deletion)
        """
        mutations = []
        mutation_types = []
        with open(fp, 'r') as f:
            content = f.readlines()
            for line in content:
                line = line.split()
                if len(line) == 1:
                    mutations.append(line[0])
                    mutation_types.append(self.type_mutation)
                elif (len(line) == 2) and (line[0] == self.type_deletion):
                    mutation_types.append(self.type_deletion)
                    mutations.append(line[1])
                elif (len(line) == 2) and (line[0] == self.type_insert):
                    mutation_types.append(self.type_insert)
                    mutations.append(line[1])
                elif (len(line) == 2) and (line[0] == self.type_combined):
                    mutation_types.append(self.type_combined)
                    mutations.append(line[1])
                else:
                    print(f"Please check format of mutation {line}, one mutation should be written per line and for inserts and deletions check the requirements.")
                    sys.exit() 
        # (1) Check there are NO non-natural amino acids in the mutations
        # (2) Check that there are enough mutations to process
        # (3) Check formatting of mutations
        # (4) check that there are no duplicate mutations
        if len(mutations) != len(set(mutations)):
            # TODO print duplicate mutations
            print("Duplicate mutations detected. Please remove and rerun.")
            sys.exit()
        if (self.check_input_mutations(mutations, mutation_types)) and (self.check_number_input_mutations(mutations)):  
            return mutations, mutation_types


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
        plt.hist(data, bins=bins, align=('mid'), label=labels)
        # plt.set_xticks(labels)
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
        
    def optimize_bins(self, x, bandwidth=200):
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
        
        for i in range(self.n_iterations):

            # Start with default value
            clusters = self.meanshift(x, bandwidth)

            # Calculate costs and check size of fragments
            cost = self.calculate_cost(clusters)
            new_bandwidth = self.check_fragment_sizes(clusters, bandwidth)
            
            print(new_bandwidth, cost)
            
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
        codon_dict = read_codon_usage(fp = self.codon_usage)
        highest_freq = 0
        most_occuring_codon = 'xxx'
        for codon in codon_list:
            codon_freq = codon_dict[codon]
            if codon_freq > highest_freq:
                highest_freq = codon_freq
                most_occuring_codon = codon
        return most_occuring_codon

    def extract_mut_codons(self, res: str):
        mut_codons = []
        for key, value in DNA_Codons.items():  # d[codon] = amino acid
            if value == res:
                mut_codons.append(key.lower())
        return mut_codons
    
    def find_gene_block(self, gene_blocks, mutation_idx):
        for key, value in gene_blocks.items():
            begin_range, end_range = self.gene_block_range(key)
            if begin_range < int(mutation_idx) < end_range:
                result = (key, value)
                return result
        
    def check_position_gene_block(self, gene_block, idx_mutation):
        begin_range, end_range = self.gene_block_range(gene_block)
        # SHOULD NOT BE NECCESSARY

        # print(begin_range, end_range, idx_mutation)
        # # TODO CHANGE
        # if (idx_mutation - begin_range) > self.idt_min_length_fragment:
        #     # self.alter_gene_block_range(gene_block, 'extend_begin')
        #     print("Mutation is too close to beginning of gene block")
        #     sys.exit()
        # elif (end_range - idx_mutation) < self.idt_min_length_fragment:
        #     print(end_range - idx_mutation)
        #     print("Mutation is too close to final part of gene block")
        #     sys.exit()

    def alter_gene_block_range(self, gene_block, alteration):
        if alteration == 'extend_begin':
            pass
            # TODO Add 20 nucleotides to beginning of gene block
            # TODO: Make this more general
            # TODO: Then check again if size is ok
            # TODO: If correct, then return new range

        
    def find_mutation_index_in_gene_block(self, gene_block, idx_mutation):
        begin_range, _ = self.gene_block_range(gene_block)
        # Find index of mutation within geneblock
        index = idx_mutation - begin_range
        # Check that mutation is not in the first or final X residues of the gene block (this would make it very difficult to design IVA primers)
        self.check_position_gene_block(gene_block, index)
        return index

    def mutate_gene_block(self, mut_codon, mut_index, gene_block_seq, change_type, end_idx = None):
        # end_idx neccessary for deletion
        if (change_type == self.type_mutation) or (change_type == self.type_combined):
            # Change this codon in the gene_block
            mut_block = gene_block_seq[:mut_index -3] + mut_codon + gene_block_seq[mut_index:]
        elif change_type == self.type_insert:
            mut_block = gene_block_seq[:mut_index] + mut_codon + gene_block_seq[mut_index:]
        elif change_type == self.type_deletion:
            mut_block = gene_block_seq[:mut_index -3] + gene_block_seq[end_idx -3:]
        return mut_block

    def write_gene_blocks_to_txt(self, gene_block_dict, 
                                outpath, 
                                fname="gene_blocks.txt"):
        header = ['mutation', 'gene block name', 'length gene block', 'gene block sequence', 'index mutation', 'mut codon', 'type']
        outfile = os.path.join(outpath, fname)
        with open(outfile, 'w+') as out:
            out.write('\t'.join(header) + '\n')
            for key, value in gene_block_dict.items():
                len_gene_block = len(value[1])
                out.write(key + '\t' + value[0] + '\t' + str(len_gene_block) + '\t' + value[1] + '\t' + str(value[2]) + '\t' + value[3] + '\t' + value[4] + '\n')

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
        
    @staticmethod
    def short_name(name):
        short_name = '_'.join(name.split('_')[0:2])
        return short_name
    
    @staticmethod
    def check_mut_format(mut, mut_type):
        valid = 'acdefghiklmnpqrstvwy'
        if (mut_type == "Mutation") or (mut_type == "Combined"):
            if (mut[0].lower() in valid) and (mut[-1].lower() in valid) and (isinstance(int(mut[1:-1]), int)):
                return True
        elif (mut_type == 'Insert') or (mut_type == 'Deletion'):
            if (mut[0].lower() in valid) and (isinstance(int(mut[1:]), int)):
                return True
        else:
            return False
    
    @staticmethod
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
        
    @staticmethod
    def check_for_start_stop_codon(sequence):
        # Make sure that the sequence starts with a start codon and ends with a stop codon
        if (sequence[0:3] == "ATG".lower()) and ((sequence[-3:] == "TAA".lower()) | (sequence[-3:] == "TAG".lower()) | (sequence[-3:] =="TGA".lower())):
            return True
        else:
            print("Sequence does not start with a start codon or end with an stop codon, \
                this is very likely to result in a shifted reading frame. Please correct your input sequence.")
            sys.exit()

    @staticmethod  
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
        result_type = DesignEblocks.check_type_input_sequence(sequence)  # Check that the sequence is of type DNA
        result_start_stop = DesignEblocks.check_for_start_stop_codon(sequence)
        if (count == 1) and (result_type) and (result_start_stop):      
            return sequence
        else:
            print("Check input sequence")

    @staticmethod
    def gene_block_range(gene_block_name):
        begin_range = int(gene_block_name.split('_')[3])
        end_range = int(gene_block_name.split('_')[4])
        return begin_range, end_range
