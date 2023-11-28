import os
import sys
import math
import random
# import openpyxl
import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from dna_features_viewer import GraphicFeature, GraphicRecord
from utils import read_codon_usage, DNA_Codons, write_pickle, natural_amino_acids

script_dir = os.path.dirname(os.path.abspath(__file__))
template_path_96 = os.path.join(script_dir, 'data/eblocks-plate-upload-template-96.xlsx')
template_path_384 = os.path.join(script_dir, 'data/eblocks-plate-upload-template-384.xlsx')

# TODO Add @classmethod function to some of the functions
# TODO Make a separate class for plotting maybe?

class DesignEblocks:
    """
    Class to design eBlocks

    Args:
        sequence_fp (str): Path to FASTA file with DNA sequence
        mutations_fp (str): Path to TXT file with mutations
        output_fp (str): Path to store output files
        species (str): Species for which codon usage table is available
        codon_usage_fp (str): Path to folder with codon usage tables
        min_bin_overlap (int): Minimum overlap between bins
        idt_max_length_fragment (int): Maximum length of gene block
        idt_min_length_fragment (int): Minimum length of gene block
        idt_min_order (int): Minimum number of mutations to process
        optimize (str): 'cost' or 'amount' for optimization

    Returns:
        _type_: _description_
    """
    # Mutation types supported
    type_mutation = "Mutation"
    type_insert = "Insert"
    type_deletion = "Deletion"
    type_combined = "Combined"

    valid_aas = 'acdefghiklmnpqrstvwy'
    valid_nucleotides = 'atcg'

    mutation_type_colors = {'Mutation': 'black', 'Insert': 'red', 'Deletion': 'blue', 'Combined': 'green'}

    def __init__(self, 
                 sequence_fp: str, 
                 mutations_fp: str,
                 output_fp: str,
                 codon_usage_fp: str,
                 optimize: str, 
                 species = "Escherichia coli",
                 bp_price = 0.05,
                 gene_name = None,
                 min_bin_overlap = 25,
                 eblock_colors = None,
                 idt_max_length_fragment = 1500,
                 idt_min_length_fragment = 300,
                 idt_min_order = 24):
        
        self.output_fp = output_fp
        self.species = species
        self.optimization = optimize
        self.min_bin_overlap = min_bin_overlap
        self.gene_name = gene_name
        self.bp_price = bp_price
        self.idt_max_length_fragment = idt_max_length_fragment
        self.idt_min_length_fragment = idt_min_length_fragment
        self.idt_min_order = idt_min_order
        self.codon_usage_fp = codon_usage_fp
        self.gene_blocks = None
        self.num_mutations = None
        self.eblock_colors = eblock_colors
        self.dna_seq = self.read_single_seq(sequence_fp)
        self.mutations, self.mutation_types = self.read_mutations(mutations_fp)  # Types = {Mutation, Insert, Deletion, Combined}
        self.codon_usage = self.check_existance_codon_usage_table()  # Check if codon usage table is present for selected species
        self.counts = None  # Number of mutations per eblock

    def __str__(self):
        return (
            f"DesignEblocks("
            # f"sequence_fp='{self.sequence_fp}', "
            # f"mutations_fp='{self.mutations_fp}', "
            f"output_fp='{self.output_fp}', "
            f"codon_usage_fp='{self.codon_usage_fp}', "
            f"optimize='{self.optimization}', "
            f"species='{self.species}', "
            f"bp_price={self.bp_price}, "
            f"gene_name={self.gene_name}, "
            f"min_bin_overlap={self.min_bin_overlap}, "
            f"eblock_colors={self.eblock_colors}, "
            f"idt_max_length_fragment={self.idt_max_length_fragment}, "
            f"idt_min_length_fragment={self.idt_min_length_fragment}, "
            f"idt_min_order={self.idt_min_order})"
        )


    def run(self, show=False):
        """
        Run the design of eBlocks
        """
        # Find indexes in sequence where mutation occures
        idx_dna_tups, idx_test, paired = self.index_mutations(self.mutations, self.mutation_types)
        self.num_mutations = len(idx_dna_tups)
     
        # length of ebLock is checked here and should be within bounds
        clusters = self.make_clusters(idx_dna_tups, idx_test, paired)
        print(f"expected cost, counting {self.bp_price} cent per bp: {self.calculate_cost(clusters)} euros")

        # Define the beginning and end of each gene block, based on the clusters and include the minimum overlap
        bins = self.make_bins(clusters)

        # Make gene blocks (WT DNA sequences cut to the correct size, according to the bins) and renumber them starting from 1
        self.gene_blocks = self.make_gene_block(bins, self.dna_seq)

        # Count number of mutations per eblock, for making barplot
        self.counts = self.count_mutations_per_eblock(idx_dna_tups)

        # Create legend, barplot and plot of eBlocks and mutations for visualization
        # Legend
        self.plot_legend(show=show)
    
        if not self.eblock_colors:  # If no colors are provided, randomly generate them
            eblock_hex_colors = self.generate_eblock_colors()
            self.eblock_colors = list(eblock_hex_colors.values())[0:len(self.counts)]
        elif len(self.eblock_colors) < len(self.gene_blocks):  # If not enough colors are provided, randomly generate the remaining ones (amount optimization uses less colors)
            eblock_hex_colors = self.generate_eblock_colors()
            self.eblock_colors = list(eblock_hex_colors.values())[0:len(self.counts)]

        # Barplot showing the number of mutations that can be made with each eBlock
        self.make_barplot(show=show)
        
        # Plot showing the eBlocks and which mutations are in which eBlock
        self.plot_eblocks_mutations(idx_dna_tups=idx_dna_tups, eblocks=True, mutations=True, genename=self.gene_name, show=show)

        # Loop over all mutations and create the eBlocks
        results = {}
        for num, mut in enumerate(self.mutations):
            mut_type = self.mutation_types[num]
            # Create mutated eblock, based on mutation type
            results = self.create_mutated_eblock(num, mut, mut_type, idx_dna_tups, results)
                                        
        # Store output
        self.write_gene_blocks_to_txt(results, self.output_fp)
        self.write_gene_blocks_to_template(results, self.output_fp)
        write_pickle(results, self.output_fp)
        write_pickle(self.gene_blocks, self.output_fp, fname="wt_gene_blocks.npy")

        print("Designed eBlocks and stored output in ", self.output_fp)


    def create_mutated_eblock(self, num: int, mut: str, mut_type: str, idx_dna_tups: list, results: dict) -> dict:
        """
        Generates a mutated version of an eBlock based on the given mutation information.

        Parameters:
        - num (int): Index indicating the specific mutation to be processed.
        - mut (str): Description of the mutation, containing relevant details.
        - mut_type (str): Type of mutation, including 'mutation', 'insert', 'deletion', or 'combined'.
        - idx_dna_tups (list): List of tuples containing mutation information.
        - results (dict): Dictionary to store the results of the mutation analysis.

        Returns:
        - dict: Updated results dictionary containing information about the mutated eBlock.

        Note:
        - This function handles different mutation types ('mutation', 'insert', 'deletion', 'combined').
        - It applies the corresponding mutation to the eBlock and checks for validity based on various criteria.
        - The results are stored in the provided dictionary with details such as gene block name, mutated eBlock,
        mutation index, mutant codon, and mutation type.
        """
        # Change gene block in selected position for mutation
        if mut_type == self.type_mutation:

            # Find gene block and index of insert/deletion/mutation
            mut_idx = idx_dna_tups[num][1]
            mut_gene_block_name, mut_gene_block_value = self.map_mutation_to_eblock(mut_idx)
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
            return results

        elif mut_type == self.type_insert:
        
            # Design insert
            res_insert = mut.split('-')[1]
            codon_insert = self.design_insert(res_insert)
            # Find gene block and index of insert/deletion/mutation
            mut_idx = idx_dna_tups[num][1]
            mut_gene_block_name, mut_gene_block_value = self.map_mutation_to_eblock(mut_idx)
            idx = self.find_mutation_index_in_gene_block(mut_gene_block_name, mut_idx)
            # Check if WT codon at index is same residue as mutation
            self.check_wt_codon(mut_gene_block_value, idx, mut[0])
            # Mutate gene block
            mut_gene_block = self.mutate_gene_block(codon_insert, idx, mut_gene_block_value, mut_type)
            # Check if eBlock is too long / too short
            self.check_eblock_length(mut_gene_block)
            results[mut] = [mut_gene_block_name, mut_gene_block, idx, codon_insert, mut_type]
            return results

        elif mut_type == self.type_deletion:

            idx_del_start = idx_dna_tups[num][1]
            idx_del_end = int(mut.split('-')[1][1:]) * 3
            mut_gene_block_name, mut_gene_block_value = self.map_mutation_to_eblock(idx_del_start)
            
            idx = self.find_mutation_index_in_gene_block(mut_gene_block_name, idx_del_start)
            idx_end = self.find_mutation_index_in_gene_block(mut_gene_block_name, idx_del_end)

            # Check if WT codon at index is same residue as mutation
            self.check_wt_codon(mut_gene_block_value, idx, mut[0])
            # Mutate gene block
            mut_gene_block = self.mutate_gene_block('', idx, mut_gene_block_value, mut_type, idx_end)
            # Check if eBlock is too long / too short
            self.check_eblock_length(mut_gene_block)
            results[mut] = [mut_gene_block_name, mut_gene_block, idx, '', mut_type]
            return results
        
        elif mut_type == self.type_combined:

            mut_gene_block_name  = None
            mut_gene_block_value = None
            lowest_count = None 

            for mut_i in idx_dna_tups[num]:
                possible_gene_blocks, counts = self.find_gene_block(self.gene_blocks, mut_i[1])
                if (counts == 1) and (mut_gene_block_name is None):
                    mut_gene_block_name = list(possible_gene_blocks.keys())[0]
                    mut_gene_block_value = possible_gene_blocks[mut_gene_block_name]

            if mut_gene_block_name is None:
                all_counts = [counts for _, counts in (self.find_gene_block(self.gene_blocks, mut_i[1]) for mut_i in idx_dna_tups[num])]
                lowest_count = min(all_counts)

            for mut_i in idx_dna_tups[num]:
                possible_gene_blocks, counts = self.find_gene_block(self.gene_blocks, mut_i[1])
                if (counts == lowest_count) and (mut_gene_block_name is None):
                    mut_gene_block_name = list(possible_gene_blocks.keys())[0]
                    mut_gene_block_value = possible_gene_blocks[mut_gene_block_name]
                    
                    # Try to find indexes of mutations, based on eblock. Check if they are too close to beginning or end of eblock
                    try:
                        for mut_i in idx_dna_tups[num]:
                            # Check too close to beginning or end
                            idx = self.find_mutation_index_in_gene_block(mut_gene_block_name, mut_i[1])
                            if (idx < self.min_bin_overlap) or (idx > (len(mut_gene_block_value) - self.min_bin_overlap)):
                                raise Exception("Mutation too close to beginning or end of eBlock")
                    except Exception:
                        continue

            idxs, codons = [], []
            for mut_i in idx_dna_tups[num]:

                idx = self.find_mutation_index_in_gene_block(mut_gene_block_name, mut_i[1])
                idxs.append(idx)
                # Find most occuring mutant codon based on codon usage for species
                mut_codons = self.extract_mut_codons(mut_i[0][-1])
                mut_codon = self.select_mut_codon(mut_codons)
                codons.append(mut_codon)
                # Check if WT codon at index is same residue as mutation
                self.check_wt_codon(mut_gene_block_value, idx, mut_i[0])
                # Mutate gene block
                mut_gene_block_value = self.mutate_gene_block(mut_codon, idx, mut_gene_block_value, mut_type)

            results[mut] = [mut_gene_block_name, mut_gene_block_value, idx, mut_codon, mut_type]
            return results


    def count_mutations_per_eblock(self, idx_dna_tups: list) -> dict:
        """
        Count the number of mutations in each eblock

        Args:
            idx_dna_tups (list): list of tuples with [mutation, index]

        Returns:
            dict: dictionary with eblock name as key and number of mutations as value
        """
        counts = {}
        for i in self.gene_blocks.keys():
            begin_range, end_range = self.gene_block_range(i)
            counts[i] = 1
            for idx in idx_dna_tups:
                if type(idx[1]) == int:  
                    if begin_range < idx[1] < end_range:
                        if i in counts.keys():
                            counts[i] += 1
                elif type(idx[1]) == list:  # Take the first mutation of a combined mutation, as they will all be in the same eblock
                    if begin_range < idx[1][1] < end_range:
                        if i in counts.keys():
                            counts[i] += 1
        return counts


    def check_wt_codon(self, eblock_seq: str, idx: int, mut: str):
        """
        Check if WT codon at given index is matching the mutation

        Args:
            gene_block_value (str): WT gene block DNA sequence
            idx (int): index of mutation in gene block
            mut (str): mutation in format [WT residue][index][mutant residue], such as G432W

        Returns:
            None: None
        """
        codon = eblock_seq[idx-3:idx]
        result = next((value for key, value in DNA_Codons.items() if key.lower() == codon), None)
        if result is not None and result != mut[0]:
            print(f"WT codon does not match residue {mut}, but is {result}, the codon is {codon}")
            print("This is probably due to the fact that paired mutations are not in the same eBlock")
            sys.exit()


    def check_existance_codon_usage_table(self):
        """
        Check if codon usage table is present for selected species
        """
        codon_files = [file for file in os.listdir(self.codon_usage_fp)]
        organisms_present = [file.split('.')[0] for file in codon_files]
        organisms_present_format = [' '.join(i.split('_')).lower() for i in organisms_present]

        if self.species.lower() in organisms_present_format:
            index = organisms_present_format.index(self.species.lower())
            return os.path.join(self.codon_usage_fp, codon_files[index])
        else:
            print("It looks like the codon usage table for the specified organism is not present.")
            sys.exit()


    def check_eblock_length(self, eblock_seq: str) -> bool:
        """
        Check if the length of the gene block is within bounds

        Args:
            eblock_seq (str): gene block DNA sequence
        """
        length_eblock = len(eblock_seq)
        if not self.idt_min_length_fragment <= length_eblock <= self.idt_max_length_fragment:
            if length_eblock > self.idt_max_length_fragment:
                print("Codon insert is too long for eBlock")
            else:
                print(f"Codon insert is too short for mutation eBlock, length is {length_eblock}, minimum length is {self.idt_min_length_fragment}")
                sys.exit()


    def check_input_mutations(self, mutations, mutation_types):
        """
        Make sure that none of the input mutations contains a unusual amino acid

        Args:
            mutations (list): List of mutations in format [WT residue][index][mutant residue], such as G432W
            mutation_types (list): List with type of mutations, either: Mutation, Insert or Deletion

        Returns:
            Bool: Returns true if all amino acids are valid
        """
        for mut, type in zip(mutations, mutation_types):
            if type == self.type_mutation:
                if not self.check_mut_format(mut, self.type_mutation):
                    print( f"Input {mut} contain non-natural amino acids or incorrect formatting")
                    sys.exit()
            elif type == self.type_insert:
                mut_format = mut.split('-')[0]
                added_residues = mut.split('-')[1]
                if not self.check_mut_format(mut_format, self.type_insert):
                        print( f"Input {mut} contain non-natural amino acids or incorrect formatting")
                        sys.exit()
                for i in added_residues:
                    if not i.lower() in self.valid_aas:
                        print( f"Input {i} contain non-natural amino acids or incorrect formatting")
                        sys.exit()
            elif type == self.type_deletion:
                mut_start = mut.split('-')[0]
                mut_end = mut.split('-')[0]
                if not self.check_mut_format(mut_start, self.type_deletion):
                    print( f"Input {mut_start} contain non-natural amino acids or incorrect formatting")
                    sys.exit()
                if not self.check_mut_format(mut_end, self.type_deletion):
                    print( f"Input {mut_end} contain non-natural amino acids or incorrect formatting")
                    sys.exit()
            elif type == self.type_combined:
                mut_list = mut.split('-')
                for i in mut_list:
                    if not self.check_mut_format(i, self.type_combined):
                        print( f"Input {i} contain non-natural amino acids or incorrect formatting")
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
        

    def map_mutation_to_eblock(self, mutation_index: int):
        name_val, _ = self.find_gene_block(self.gene_blocks, mutation_index)
        eblock_name = list(name_val.keys())[0]
        eblock_value = name_val[eblock_name]
        return eblock_name, eblock_value
    

    def design_insert(self, aas):
        codon_insert = ''  # Sequence to insert in gene block
        for res in aas:
            codons = self.extract_mut_codons(res)
            codon = self.select_mut_codon(codons)
            codon_insert += codon
        return codon_insert


    def index_mutations(self, mut_list: list, mut_types: list):
        """
        Find the indexes in the DNA sequence where mutations occur

        Args:
            mut_list (list): list of mutations in format [WT residue][index][mutant residue], such as G432W
            mut_types (list): list with type of mutations, either: Mutation, Combined, Insert or Deletion

        Returns:
            list: list of indexes in DNA sequence where mutations occur
        """
        mutname_idx = []  # List of tuples with [mutation, index], double mutation are split in multiple sublists
        idx_dna_list = []  # Single list containting all indexes, nested lists for double mutations and deletions
        paired = []  # Paired = double mutation and insertion
        for mut, mut_type in zip(mut_list, mut_types):
            if (mut_type == self.type_mutation) or (mut_type == self.type_insert):
                if (mut_type) == self.type_mutation:
                    idx = int(mut[1:-1]) * 3
                elif (mut_type) == self.type_insert:
                    idx = int(mut.split('-')[0][1:]) * 3
                mutname_idx.append([mut, idx])
                idx_dna_list.append(idx)

            elif (mut_type == self.type_deletion):
                mut_begin = int(mut.split('-')[0][1:])
                mut_end = int(mut.split('-')[1][1:])
                mut_length = int(mut_end) - int(mut_begin)
                
                mutname_idx.append([f"{mut}", mut_begin * 3])
                idx_end = (mut_begin * 3) + (mut_length * 3)

                tmp_list = list(range(mut_begin * 3, idx_end, 3))  # For a deletion add all indexes in between the start and end of the deletion
                idx_dna_list.append(tmp_list)
                paired.append(tuple(tmp_list))

            elif (mut_type == self.type_combined):
                mutations = mut.split('-')
                mutname_idx.append([[i, int(i[1:-1]) * 3] for i in mutations])
                idx_dna_list.append([int(i[1:-1]) * 3 for i in mutations])
                paired.append(tuple(idx_dna_list[-1]))

        return mutname_idx, idx_dna_list, paired
    

    def read_mutations(self, fp: str):
        """
        Read mutations file in TXT format.
        
        Each line of the file should contain one amino acid in the format [WT residue][index][mutant residue], such as G432W.
        
        The function performs the following checks:
        - Ensure there are NO non-natural amino acids in the mutations.
        - Verify that there are enough mutations to process.
        - Check the formatting of mutations.
        - Ensure there are no duplicate mutations.
        
        Args:
            fp (str): Filepath of the input mutations TXT file.

        Returns:
            tuple: A tuple containing two lists - 
                1. List of mutations extracted from the input file.
                2. List of corresponding mutation types ('mutation', 'insert', 'deletion').
                Returns None if checks fail.
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
        if len(mutations) != len(set(mutations)):
            print("Duplicate mutations detected. Please remove and rerun.")
            sys.exit()
        if (self.check_input_mutations(mutations, mutation_types)) and (self.check_number_input_mutations(mutations)):  
            return mutations, mutation_types


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
            bins.append(int(min(value) - self.min_bin_overlap))
            bins.append(int(max(value) + self.min_bin_overlap))
        return bins
    

    def plot_legend(self, legend_alpha=0.2, font_size='x-large', marker_size=10, linestyle='None', marker='o', loc='center', bbox_to_anchor=(0.5, 0.5), show=False):
        """
        Plot legend for eBlocks plot
        """
        # Create an empty plot with no data points
        fig, ax = plt.subplots(figsize=(3, 2))
        # Add mutation type colors to the legend
        handles = []
        for k, v in self.mutation_type_colors.items():
            handle = plt.Line2D([0], [0], marker=marker, color=f'{v}', label=f'{k}', markersize=marker_size, linestyle=linestyle, alpha=legend_alpha)
            handles.append(handle) 
        legend = ax.legend(handles=handles, loc=loc, bbox_to_anchor=bbox_to_anchor, fontsize=font_size, framealpha=legend_alpha)
        # Hide the axes
        ax.axis('off')
        fig.savefig(os.path.join(self.output_fp, 'legend.png'), dpi=100)
        if show:
            plt.show()
        else:
            plt.close()

    
    def plot_eblocks_mutations(self, idx_dna_tups=None, eblocks=True, mutations=True, genename=None, genecolor="#d3d3d3", show=False, figure_width=20, figure_length=10):
        """
        Plot mutations and selected eBlocks
        """
        features = []
        if not idx_dna_tups:
            idx_dna_tups, idx_test, paired = self.index_mutations(self.mutations, self.mutation_types)

        # Add gene to plot
        if genename:
            features.append(GraphicFeature(start=0, 
                                        end=len(self.dna_seq), 
                                        strand=+1, 
                                        color=genecolor, 
                                        label=f"{genename}"))

        # Add mutations to plot
        if mutations:
            for num, mut in enumerate(idx_dna_tups):
                if type(mut[1]) == int:
                    features.append(GraphicFeature(start=int(mut[1]), 
                                                end=int(mut[1]) + 3,
                                                strand=+1, 
                                                color=self.mutation_type_colors[self.mutation_types[num]], 
                                                label=f"{mut[0]}"))
                elif type(mut[1]) == list:
                        for m in mut:
                            features.append(GraphicFeature(start=int(m[1]), 
                                                    end=int(m[1]) + 3,
                                                    strand=+1, 
                                                    color=self.mutation_type_colors['Combined'], 
                                                    label=f"{m[0]}"))

        # Add eBlocks to plot
        if eblocks:
            for num, key in enumerate(self.gene_blocks.keys()):
                features.append(GraphicFeature(start=int(key.split('_')[3]), 
                                            end=int(key.split('_')[4]), 
                                            strand=+1, 
                                            color=self.eblock_colors[num], 
                                            label=f"Block {key.split('_')[1]}"))
        

        record = GraphicRecord(sequence_length=len(self.dna_seq), features=features)
        # if self.num_mutations:
        #     figure_length = np.ceil(len(self.dna_seq) / 150)
        #     figure_width = np.ceil(self.num_mutations / 10)
        fig_size = (figure_length, figure_width)
        fig, ax = plt.subplots(figsize=fig_size) 
        record.plot(ax=ax, figure_width=20)
        if show:
            plt.show()
        else:
            plt.close()
        if eblocks:
            fig.savefig(os.path.join(self.output_fp, f'eblocks_{self.gene_name}_N{self.num_mutations}_{self.optimization}.png'), dpi=100)
        if not eblocks:
            fig.savefig(os.path.join(self.output_fp, f'{self.gene_name}_N{self.num_mutations}.png'), dpi=100)


    def make_barplot(self, show, figure_width=5, figure_length=5):
        """
        Make barplot of bins
        """
        fig, ax = plt.subplots(figsize=(figure_width, figure_length))
        labels = []
        for k, v in self.counts.items():
            kn = k.split('_')[0] + ' ' + k.split('_')[1]
            labels.append(kn)
        ax.bar(range(len(self.counts)), list(self.counts.values()), align='center', color=self.eblock_colors)
        ax.set_xticks(range(len(self.counts)), labels)
        ax.set_ylabel('Number of mutants per eBlock')
        ax.set_xlabel('eBlock')
        ax.set_title('Number of mutants per eBlock')
        ax.bar_label(ax.containers[0])
        if show:
            plt.show()
        else:
            plt.close()
        fig.savefig(os.path.join(self.output_fp, f'counts_{self.gene_name}_N{self.num_mutations}_{self.optimization}.png'), dpi=100)


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
            cost = len_gene_block * self.bp_price * len(value)
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
    

    def optimize_clusters(self, possibilities: str) -> dict:
        # Choose the best clustering based on the optimization parameter
        if self.optimization == 'cost':
            print("Optimizing based on price per bp ...")
            lowest_cost, best_clustering = min((self.calculate_cost(value), value) for value in possibilities.values())
            print(f"Lowest cost: {lowest_cost} with cluster {best_clustering}")
            return best_clustering

        elif self.optimization == 'amount':
            print("Optimizing based on amount of eBlocks ...")
            fewest_blocks, best_clustering = min((len(value), value) for value in possibilities.values())
            print(f"Fewest blocks: {fewest_blocks} with cluster {best_clustering}")
            return best_clustering
        

    def make_clusters(self, idxs_tuple, idx_test, paired_mutations):        
        possibilities = {} # Store all possible clusterings
        n = 1
        valid_clusters = True

        while valid_clusters:
            
            clusters = {}
            cluster_labels, idxs_reordered = self.kmeans_clustering(idx_test, paired_mutations, n)

            # Calculate the size of each cluster
            for i, j in zip(cluster_labels, idxs_reordered):
                if i not in clusters:
                    clusters[i] = []
                clusters[i].append(j)

            # Check if the size of the clusters is within bounds
            cluster_sizes = [max(v) - min(v) for v in clusters.values()]

            # Find in which cluster the insertion and deleted is located and check if the size of the cluster is within bounds when the mutation is added
            for num, mut in enumerate(idxs_tuple):
                if self.mutation_types[num] == self.type_insert:
                    for key, value in clusters.items():
                        if mut[1] in value:
                            length_insert = len(mut[0].split('-')[1]) * 3
                            cluster_sizes[key] += length_insert
                elif self.mutation_types[num] == self.type_deletion:
                    for key, value in clusters.items():
                        if mut[1] in value:
                            del_begin = int(mut[0].split('-')[0][1:]) * 3
                            del_end = int(mut[0].split('-')[1][1:]) * 3
                            length_del = (del_end - del_begin) + 3
                            cluster_sizes[key] -= length_del

            max_cluster_size = max(cluster_sizes)
            min_cluster_size = min(cluster_sizes)

            if max_cluster_size > (self.idt_max_length_fragment - 2 * self.min_bin_overlap): # Take into account the OH on both sides of the gene block, so decrease max size with 2 times min length of OH
                n += 1
            elif min_cluster_size < (self.idt_min_length_fragment - 2 * self.min_bin_overlap):
                valid_clusters = False
            else:
                possibilities[f'cluster N={n}'] = clusters
                n += 1
        
        # TODO Fix this, it is not working properly
        if len(possibilities) == 0:  # No valid clusters found
            print("No valid clusterings found, please check your input mutations and make sure that the multiple mutants are not too far apart.")
            sys.exit()

        optimal_clustering = self.optimize_clusters(possibilities)
        return optimal_clustering

        
    def kmeans_clustering(self, idx_test, paired_mutations, num_clusters, n_init='auto', random_state=42, OMP_NUM_THREADS=1):
        # Extract the first index of each mutation in idx_test
        idx_first = [np.mean(i) if isinstance(i, list) else i for i in idx_test]
        
        # Create a numpy array for clustering
        mutation_arr = np.array(idx_first).reshape(-1, 1)

        # Initialize KMeans with the number of clusters
        kmeans = KMeans(n_clusters=num_clusters, random_state=random_state, n_init=n_init)

        # Fit the model and obtain cluster labels for connected indices
        cluster_labels = kmeans.fit_predict(mutation_arr)

        # Add the remaining indices of the pairs to the assigned cluster
        for pair in paired_mutations:
            mean_values = [np.mean(pair) for _ in pair]
            indices = [list(mutation_arr).index(mean) for mean in mean_values]
            cluster_label = cluster_labels[indices[0]]
            
            cluster_labels = np.concatenate((cluster_labels, np.full(len(pair), cluster_label)))
            idx_first.extend(pair)

        # Remove the mean values from the list of indices and cluster labels
        to_remove = sorted(set(indices), reverse=True)
        for i in to_remove:
            del idx_first[i]
            cluster_labels = np.delete(cluster_labels, i)

        return cluster_labels, idx_first
    

    def renumber_gene_blocks(self, gene_blocks):
        new_gene_blocks = {}
        for num, (k, v) in enumerate(sorted(gene_blocks.items(), key=lambda item: int(item[0].split('_')[3])), 1):
            new_gene_blocks[f'Block_{num}_pos_{k.split("_")[3]}_{k.split("_")[4]}'] = v
        return new_gene_blocks


    def make_gene_block(self, bins, dna_sequence):
        gene_blocks = {}
        block_num = 1
        for num in range(0, len(bins), 2):
            name = f'Block_{block_num}_pos_{bins[num]}_{bins[num+1]}'
            block = dna_sequence[bins[num]:bins[num+1]]
            gene_blocks[name] = str(block)
            block_num += 1
        return self.renumber_gene_blocks(gene_blocks)


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
        Choose codon from list of codons based on occurrence of codon in nature.
        """
        codon_dict = read_codon_usage(fp=self.codon_usage)
        most_occuring_codon = max(codon_list, key=codon_dict.get, default='xxx')
        return most_occuring_codon


    def extract_mut_codons(self, res: str):
        return [key.lower() for key, value in DNA_Codons.items() if value == res]
    

    def find_gene_block(self, gene_blocks: dict, mutation_idx: int):
        results = {key: value for key, value in gene_blocks.items() if self.is_within_gene_block(key, mutation_idx)}
        count = len(results)
        return results, count


    def is_within_gene_block(self, gene_block_name: str, mutation_idx: int) -> bool:
        begin_range, end_range = self.gene_block_range(gene_block_name)
        return begin_range < int(mutation_idx) < end_range


    def mutate_gene_block(self, mut_codon, mut_index, gene_block_seq, change_type, end_idx = None):
        """
        Mutate gene block based on mutation type

        Args:
            mut_codon (str): mutant codon
            mut_index (int): index of mutation in gene block
            gene_block_seq (str): gene block DNA sequence
            change_type (str): type of mutation, either: Mutation, Insert or Deletion
            end_idx (int, optional): index of end of deletion. Defaults to None.

        end_idx is only neccessary for deletions, as the end of the deletion is not the same as the index of the mutation

        Returns:
            str: mutated gene block
        """
        if (change_type == self.type_mutation) or (change_type == self.type_combined):
            mut_block = gene_block_seq[:mut_index -3] + mut_codon + gene_block_seq[mut_index:]
        elif change_type == self.type_insert:
            mut_block = gene_block_seq[:mut_index] + mut_codon + gene_block_seq[mut_index:]
        elif change_type == self.type_deletion:
            mut_block = gene_block_seq[:mut_index -3] + gene_block_seq[end_idx -3:]
        return mut_block


    def write_gene_blocks_to_txt(self, eblocks: dict, outpath):
        fname= f"gene_blocks_{self.gene_name}_N{self.num_mutations}_{self.optimization}.txt"
        header = ['mutation', 'gene block name', 'length gene block', 'gene block sequence', 'index mutation', 'mut codon', 'type']
        outfile = os.path.join(outpath, fname)
        with open(outfile, 'w+') as out:
            out.write('\t'.join(header) + '\n')
            for key, value in eblocks.items():
                len_gene_block = len(value[1])
                out.write(f"{key}\t{value[0]}\t{len_gene_block}\t{value[1]}\t{value[2]}\t{value[3]}\t{value[4]}\n")


    def determine_template(self, n_samples: int):
        if n_samples <= 96:
            return template_path_96, f"eblocks-plate-upload-96-{self.gene_name}_N{self.num_mutations}_{self.optimization}.xlsx"
        elif n_samples > 96 and n_samples <= 384:
            return template_path_384, f"eblocks-plate-upload-template-384-{self.gene_name}_N{self.num_mutations}_{self.optimization}.xlsx"
        elif n_samples > 384:
            return template_path_384, f"eblocks-plate-upload-template-384-{self.gene_name}_N{self.num_mutations}_{self.optimization}.xlsx" 


    def write_gene_blocks_to_template(self, eblocks: dict, outpath):
        n_samples = len(eblocks)
        template, fname = self.determine_template(n_samples)
        n_plates = self.determine_number_plates(n_samples)

        if n_plates == 1:
            outfile = os.path.join(outpath, fname)
            names, seqs = [], []

            for key, value in eblocks.items():
                mutation = key
                block = '-'.join(value[0].split('_')[0:2])
                names.append(f"{mutation}_{block}")
                seqs.append(value[1])

            wells = self.extract_wells_from_template(template=template)[:len(names)]
            df = pd.DataFrame({'Well Position': wells, 'Name': names, 'Sequence': seqs})
            df.to_excel(outfile, index=False)

        elif n_plates > 1:
            max_samples = 384
            fnames = [f"eblocks-plate-upload-template-384-filled_{i}.xlsx" for i in range(1, n_plates + 1)]

            for num, outfile in enumerate([os.path.join(outpath, i) for i in fnames], 1):
                names, seqs = [], []

                for index, (key, value) in enumerate(eblocks.items(), 1):
                    if (max_samples * (num - 1) + 1) <= index <= (max_samples * num):
                        mutation = key
                        block = '-'.join(value[0].split('_')[0:2])
                        names.append(f"{mutation}_{block}")
                        seqs.append(value[1])

                wells = self.extract_wells_from_template(template=template)[:len(names)]
                df = pd.DataFrame({'Well Position': wells, 'Name': names, 'Sequence': seqs})
                df.to_excel(outfile, index=False)


    def find_mutation_index_in_gene_block(self, gene_block, idx_mutation: int) -> int:
        begin_range, _ = self.gene_block_range(gene_block)
        return idx_mutation - begin_range


    def extract_wells_from_template(self, template):
        df = pd.read_excel(template)
        wells = df['Well Position'].tolist()
        return wells


    @staticmethod
    def determine_number_plates(n_samples: int):
        if n_samples <= 96:
            return 1
        elif n_samples > 96 and n_samples <= 384:
            return 1
        elif n_samples > 384:
            return math.ceil(n_samples / 384)


    @staticmethod
    def check_mut_format(mut: str, mut_type: str) -> bool:
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
        """    
        valid = 'actg'
        if all(i.lower() in valid for i in sequence):
            return True
        else:
            print("Please provide a DNA sequence")
            sys.exit()


    @staticmethod
    def check_for_start_stop_codon(dnaseq):
        start_codon = 'atg'
        stop_codons = ['taa', 'tag', 'tga']
        if (dnaseq[:3] == start_codon) and (dnaseq[-3:].lower() in stop_codons):
            return True
        else:
            print("Sequence does not start with a start codon or end with a stop codon. This is very likely to result in a shifted reading frame. Please correct your input sequence.")
            sys.exit()


    @staticmethod  
    def read_single_seq(fp):
        """
        Read single DNA sequence in FASTA format using biopython

        Args:
            fp (str):   Location of the FASTA file

        Returns:
            sequence (str): DNA sequence
        """
        for record in SeqIO.parse(fp, "fasta"):
            sequence = record.seq

        if (DesignEblocks.check_type_input_sequence(sequence)) and (DesignEblocks.check_for_start_stop_codon(sequence)):
            return sequence
        else:
            print("Check input sequence. This should be a single sequence in FASTA format.")
            sys.exit()


    @staticmethod
    def translate_sequence(dna_seq):    
        """
        Translate DNA sequence to protein sequence
        """
        return dna_seq.translate()


    @staticmethod
    def gene_block_range(eblock_name: str):
        begin_range = int(eblock_name.split('_')[3])
        end_range = int(eblock_name.split('_')[4])
        return begin_range, end_range
    

    @staticmethod
    def generate_eblock_colors() -> dict:
        """
        Create dictionary with colors for plotting eBlocks
        """
        return {i: '#%06X' % random.randint(0, 0xFFFFFF) for i in range(100)}