import os
import sys
import math
import random
# import openpyxl
import numpy as np
import pandas as pd
from Bio import SeqIO
from operator import add, sub
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from sklearn.cluster import MeanShift
from sklearn.cluster import KMeans
from dna_features_viewer import GraphicFeature, GraphicRecord
from utils import read_codon_usage, DNA_Codons, write_pickle, log_to_file_and_console, create_or_clear_file

script_dir = os.path.dirname(os.path.abspath(__file__))

# TODO Change this
template_path_96 = os.path.join(script_dir, 'data/eblocks-plate-upload-template-96.xlsx')
template_path_384 = os.path.join(script_dir, 'data/eblocks-plate-upload-template-384.xlsx')

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

    Returns:
        _type_: _description_
    """

    # TODO Add all parameters here
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
        self.optimization = optimize
        self.min_bin_overlap = min_bin_overlap
        self.gene_name = gene_name
        self.bp_price = bp_price  # Based on IDT pricing
        self.idt_max_length_fragment = idt_max_length_fragment
        self.idt_min_length_fragment = idt_min_length_fragment
        self.idt_min_order = idt_min_order
        self.codon_usage_fp = codon_usage_fp
        self.gene_blocks = None

        self.colors = None
        
        self.dna_seq = self.read_single_seq(sequence_fp)
        self.mutations, self.mutation_types = self.read_mutations(mutations_fp)  # Types = {Mutation, Insert, Deletion}
        self.codon_usage = self.check_existance_codon_usage_table()  # Check if codon usage table is present for selected species

        # TODO ADD THE PLOT HERE WITHOUT THE EBLOCKS, JUST THE MUTATIONS
        # TODO Put all check functions in maybe a check class or just close to each other

    def run(self, show=False):
        """
        Run the design of eBlocks
        """

        # TODO Add option IDT that formats it to IDT format

        # Find indexes in sequence where mutation occures
        idx_dna_tups, idx_test, paired = self.index_mutations(self.mutations, self.mutation_types)
     
        # length of ebLock is already checked here and should be within bounds
        # TODO Check here upon addition of the insert whether the gene block is still within bounds
        clusters = self.make_clusters(idx_dna_tups, idx_test, paired) # or 'amount', 'cost'
        print(f"expected cost, counting {self.bp_price} cent per bp: {self.calculate_cost(clusters)} euros")
        # print(clusters)

        bins = self.make_bins(clusters) # OH is added to both sides of the gene bloc
        # print(bins)
        # print(len(bins))

        # Make gene blocks (WT DNA sequences cut to the correct size)
        # Gene blocks are renumbered here
        self.gene_blocks = self.make_gene_block(bins, self.dna_seq)
        # print(self.gene_blocks)

        # Plot gene blocks
        # TODO Save plot to file?
        # TODO Make legend and save to file
        # TODO SEPARATE PLOTTING AND SAVING TO FILE
        # TODO CHECK WHAT HAPPENS WHEN USING COMMANDLINE FOR RUNNING (SHOULD NOT SHOW)
        # TODO Save colors of eblock to instance

        # Plot legend
        legend = self.plot_legend(show=show)
        legend.savefig(os.path.join(self.output_fp, 'legend.png'), dpi=100)

        # Make histogram with bins
        counts = self.count_mutations_per_eblock(idx_dna_tups)

        if not self.colors:
            eblock_hex_colors = self.eblock_colors()
            self.colors = list(eblock_hex_colors.values())[0:len(counts)]
        eblock_counts = self.make_barplot(counts, show=show)
        eblock_counts.savefig(os.path.join(self.output_fp, 'barplot.png'), dpi=100)

        record = self.plot_eblocks_mutations(idx_dna_tups=idx_dna_tups, eblocks=True, mutations=True, genename=self.gene_name)
        record.plot(figure_width=20)[0].figure.savefig(os.path.join(self.output_fp, 'eblocks.png'), dpi=100)
        if show:
            record.plot(figure_width=20)

        # TODO Plot the lengths of the gene blocks
        # for k, v in self.gene_blocks.items():
        #     print(k, len(v))
        
        # Create dictionary that contains the gene blocks and the mutations that are in the gene blocks
        # mut_idxs = {}
        # for key, value in self.gene_blocks.items():
        #     mut_idxs[key] = []
        #     for idx in idx_dna_tups:
        #         if type(idx[1]) == int:
        #             if int(key.split('_')[3]) < idx[1] < int(key.split('_')[4]):
        #                 mut_idxs[key].append(idx[0])
        #         elif type(idx[1]) == list:
        #     pass

        # TODO SIZES ALREADY HAVE BEEN CHECKED, SO NO NEED TO CHECK AGAIN HERE

        results = {}
        for num, mut in enumerate(self.mutations):
            mut_type = self.mutation_types[num]
                            
            # Change gene block in selected position for mutation
            if mut_type == self.type_mutation:

                # Find gene block and index of insert/deletion/mutation
                mut_idx = idx_dna_tups[num][1]
                name_val, _ = self.find_gene_block(self.gene_blocks, mut_idx)
                mut_gene_block_name = list(name_val.keys())[0]
                mut_gene_block_value = name_val[mut_gene_block_name]
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
                name_val, _ = self.find_gene_block(self.gene_blocks, mut_idx)
                mut_gene_block_name = list(name_val.keys())[0]
                mut_gene_block_value = name_val[mut_gene_block_name]
                idx = self.find_mutation_index_in_gene_block(mut_gene_block_name, mut_idx)

                # Check if WT codon at index is same residue as mutation
                self.check_wt_codon(mut_gene_block_value, idx, mut[0])

                # Mutate gene block
                mut_gene_block = self.mutate_gene_block(codon_insert, idx, mut_gene_block_value, mut_type)

                # Check if eBlock is too long / too short
                self.check_eblock_length(mut_gene_block)
                results[mut] = [mut_gene_block_name, mut_gene_block, idx, codon_insert, mut_type]

            elif mut_type == self.type_deletion:

                idx_del_start = idx_dna_tups[num][1]
                idx_del_end = int(mut.split('-')[1][1:]) * 3
                
                name_val, _ = self.find_gene_block(self.gene_blocks, idx_del_start)
                mut_gene_block_name = list(name_val.keys())[0]
                mut_gene_block_value = name_val[mut_gene_block_name]

                idx = self.find_mutation_index_in_gene_block(mut_gene_block_name, idx_del_start)
                idx_end = self.find_mutation_index_in_gene_block(mut_gene_block_name, idx_del_end)

                # Check if WT codon at index is same residue as mutation
                self.check_wt_codon(mut_gene_block_value, idx, mut[0])

                # Mutate gene block
                mut_gene_block = self.mutate_gene_block('', idx, mut_gene_block_value, mut_type, idx_end)

                # Check if eBlock is too long / too short
                self.check_eblock_length(mut_gene_block)
                results[mut] = [mut_gene_block_name, mut_gene_block, idx, '', mut_type]

            # TODO Cleanup this part
            elif mut_type == self.type_combined:

                mut_gene_block_name = None
                mut_gene_block_value = None

                all_counts = []
                for mut_i in idx_dna_tups[num]:
                    possible_gene_blocks, counts = self.find_gene_block(self.gene_blocks, mut_i[1])
                    # print("possible_gene_blocks", possible_gene_blocks, counts)
                    all_counts.append(counts)
                    if counts == 1:
                        mut_gene_block_name = list(possible_gene_blocks.keys())[0]
                        mut_gene_block_value = possible_gene_blocks[mut_gene_block_name]

                if not mut_gene_block_name:
                    lowest_count = min(all_counts)
                    for mut_i in idx_dna_tups[num]:
                        possible_gene_blocks, counts = self.find_gene_block(self.gene_blocks, mut_i[1])
                        if counts == lowest_count:
                            mut_gene_block_name = list(possible_gene_blocks.keys())[0]
                            mut_gene_block_value = possible_gene_blocks[mut_gene_block_name]
                            try:
                                # Try to find indexes of mutations, based on eblock. Check if they are too close to beginning or end of eblock
                                idxs = []
                                codons = []
                                for mut_i in idx_dna_tups[num]:
                                    # Check too close to beginning or end
                                    idx = self.find_mutation_index_in_gene_block(mut_gene_block_name, mut_i[1])
                                    if (idx < self.min_bin_overlap) or (idx > (len(mut_gene_block_value) - self.min_bin_overlap)):
                                        raise Exception("Mutation too close to beginning or end of eBlock")
                                    mut_codons = self.extract_mut_codons(mut_i[0][-1])
                                    mut_codon = self.select_mut_codon(mut_codons)
                                    codons.append(mut_codon)
                                    self.check_wt_codon(mut_gene_block_value, idx, mut_i[0])
                            except:
                                continue

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
        write_pickle(self.gene_blocks, self.output_fp, fname="wt_gene_blocks.npy")

        print("Designed eBlocks and stored output in ", self.output_fp)



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
                elif type(idx[1]) == list:  # Just take the first mutation of a combined mutation, as they will all be in the same eblock
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
        valid = 'acdefghiklmnpqrstvwy'
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
                    if not i.lower() in valid:
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
        Each line of the file should contain one amino acid in the format [WT residue][index][mutant residue], such as G432W

        Following checks are performed:
        - Check there are NO non-natural amino acids in the mutations
        - Check that there are enough mutations to process
        - Check formatting of mutations
        - check that there are no duplicate mutations
        
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
            bins.append(min(value) - self.min_bin_overlap)
            bins.append(max(value) + self.min_bin_overlap)
        return bins
    
    def plot_legend(self, legend_alpha=0.2, font_size='x-large', marker_size=10, linestyle='None', marker='o', loc='center', bbox_to_anchor=(0.5, 0.5), show=False):
        """
        Plot legend for eBlocks plot
        """
        mutation_type_colors = self.mutation_type_colors()
        # Create an empty plot with no data points
        fig, ax = plt.subplots(figsize=(3, 2))
        # Add mutation type colors to the legend
        handles = []
        for k, v in mutation_type_colors.items():
            handle = plt.Line2D([0], [0], marker=marker, color=f'{v}', label=f'{k}', markersize=marker_size, linestyle=linestyle, alpha=legend_alpha)
            handles.append(handle) 
        legend = ax.legend(handles=handles, loc=loc, bbox_to_anchor=bbox_to_anchor, fontsize=font_size, framealpha=legend_alpha)
        # Hide the axes
        ax.axis('off')
        if not show:
            plt.close()
        return fig
    
    def plot_eblocks_mutations(self, idx_dna_tups=None, eblocks=True, mutations=True, genename=None, genecolor="#d3d3d3"):
        """
        Plot mutations and selected eBlocks
        """
        mutation_type_colors = self.mutation_type_colors()
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
                                                end=int(mut[1]) + 3, # TODO MAKE MORE SPECIFIC FOR INSERTS ETC
                                                strand=+1, 
                                                color=mutation_type_colors[self.mutation_types[num]], 
                                                label=f"{mut[0]}"))
                elif type(mut[1]) == list:
                        for m in mut:
                            features.append(GraphicFeature(start=int(m[1]), 
                                                    end=int(m[1]) + 3, # TODO MAKE MORE SPECIFIC FOR INSERTS ETC
                                                    strand=+1, 
                                                    color=mutation_type_colors['Combined'], 
                                                    label=f"{m[0]}"))

        # Add eBlocks to plot
        if eblocks:
            for num, key in enumerate(self.gene_blocks.keys()):
                features.append(GraphicFeature(start=int(key.split('_')[3]), 
                                            end=int(key.split('_')[4]), 
                                            strand=+1, 
                                            color=self.colors[num], 
                                            label=f"Block {key.split('_')[1]}"))
            
        record = GraphicRecord(sequence_length=len(self.dna_seq), features=features)
        return record

    def make_barplot(self, data, show):
        """
        Make barplot of bins
        """
        fig, ax = plt.subplots(figsize=(8, 4))
        labels = []
        for k, v in data.items():
            kn = k.split('_')[0] + ' ' + k.split('_')[1]
            labels.append(kn)
        ax.bar(range(len(data)), list(data.values()), align='center', color=self.colors)
        ax.set_xticks(range(len(data)), labels)
        ax.set_ylabel('Number of mutants per eBlock')
        ax.set_xlabel('eBlock')
        ax.set_title('Number of mutants per eBlock')
        ax.bar_label(ax.containers[0])
        if not show:
            plt.close()
        return fig

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
        

    def make_clusters(self, idxs_tuple, idx_test, paired_mutations):
        # OTHER OPTIMIZE Possibilty = 'amount', 'cost'
        
        possibilities = {} # Store all possible clusterings
        n = 1
        valid_clusters = True

        # TODO CHECK IF PAIRED MUTATIONS ARE CORRECTLY ADDED TO THE CLUSTER > CHECK THIS
        # TODO IF DESIGN OF EBLOCK IS NOT POSSIBLE > SUGGEST REMOVAL OF MUTATIONS

        while valid_clusters:
            
            # print(f"Clustering with {n} clusters ...")

            clusters = {}
            cluster_labels, idxs_reordered = self.kmeans_clustering(idx_test, paired_mutations, n)
            # print("cluster_labels", cluster_labels)
            # print("idxs_reordered", idxs_reordered)

            # Calculate the size of each cluster
            for i, j in zip(cluster_labels, idxs_reordered):
                if i not in clusters:
                    clusters[i] = []
                clusters[i].append(j)

            # Check if the size of the clusters is within bounds
            cluster_sizes = [max(v) - min(v) for v in clusters.values()]
            # print(f"Cluster sizes: {cluster_sizes}", type(cluster_sizes))

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
                # print(f"N={n}, Cluster size is still too large, increasing the number of clusters")
                # print(f"Max cluster size: {max_cluster_size}, max allowed size: {self.idt_max_length_fragment - 2 * self.min_bin_overlap}")
                n += 1
            elif min_cluster_size < (self.idt_min_length_fragment - 2 * self.min_bin_overlap):
                # print(f"N={n}, Cluster size is too small, stop increasing the number of clusters")
                # print(f"Min cluster size: {min_cluster_size}, min allowed size: {self.idt_min_length_fragment - 2 * self.min_bin_overlap}")
                valid_clusters = False
            else:
                # print(f"N={n}, Cluster size is within bounds, increasing the number of clusters")
                # print(f"Min cluster size: {min_cluster_size}, min allowed size: {self.idt_min_length_fragment - 2 * self.min_bin_overlap}")
                possibilities[f'cluster N={n}'] = clusters # TODO ADD CLUSTERING PARAMS HERE? Atleast store them somwhere
                n += 1
        
        if len(possibilities) == 0:  # No valid clusters found
            # The double mutants are too far apart to be clustered together
            # TODO Think about what to do here and what to suggest to the user
            # TODO Find out which mutations are too far apart and suggest to remove them
            print("No valid clusterings found, please check your input mutations and make sure that the multiple mutants are not too far apart.")
            sys.exit()

        # Choose the best clustering based on the optimization parameter
        if self.optimization == 'cost':
            print("Optimizing based on price per bp ...")
            # Find the clustering with the lowest cost
            lowest_cost = np.inf
            for key, value in possibilities.items():
                cost = self.calculate_cost(value)
                if cost < lowest_cost:
                    lowest_cost = cost
                    best_clustering = value
            print(f"Lowest cost: {lowest_cost} with cluster {key}")
            return best_clustering
        
        elif self.optimization == 'amount':
            print("Optimizing based on amount of eBlocks ...")
            # Find the clustering with the lowest number of eBlocks
            fewest_blocks = np.inf
            for key, value in possibilities.items():
                n_blocks = len(value)
                if n_blocks < fewest_blocks:
                    fewest_blocks = n_blocks
                    best_clustering = value
            print(f"Fewest blocks: {fewest_blocks} with cluster {key}")
            return best_clustering
        

    def kmeans_clustering(self, idx_test, paired_mutations, num_clusters, visualize=False, n_init='auto', random_state=42):

        idx_first = []
        for i in idx_test:
            if isinstance(i, list):
                mean = np.mean(i)
                idx_first.append(mean)
            else:
                idx_first.append(i)

        mutation_arr = np.array(idx_first).reshape(-1, 1)

        # Initialize KMeans with the number of clusters
        kmeans = KMeans(n_clusters=num_clusters, random_state=random_state, n_init=n_init)

        # Fit the model and obtain cluster labels for connected indices
        cluster_labels = kmeans.fit_predict(mutation_arr)
        cluster_labels = list(cluster_labels)
        # print("length cluster label:", len(cluster_labels))

        # Add the remaining indices of the pairs to the assigned cluster
        to_remove = []
        for pair in paired_mutations:
            mean = np.mean(pair)
            index = list(mutation_arr).index(mean)
            to_remove.append(index)
            cluster_label = cluster_labels[index]
            L = [cluster_label] * len(pair)
            for i in L:
                cluster_labels.append(i)
            for i in pair:
                idx_first.append(i)
            # print(pair, mean, index, cluster_label, L)
            # Remove mean value from labels and indices
        
        # Remove the mean values from the list of indices and cluster labels
        to_remove = sorted(to_remove, reverse=True)
        for i in to_remove:
            del idx_first[i]
            del cluster_labels[i]
             
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


    def write_gene_blocks_to_txt(self, eblocks: dict, outpath, fname="gene_blocks.txt"):
        header = ['mutation', 'gene block name', 'length gene block', 'gene block sequence', 'index mutation', 'mut codon', 'type']
        outfile = os.path.join(outpath, fname)
        with open(outfile, 'w+') as out:
            out.write('\t'.join(header) + '\n')
            for key, value in eblocks.items():
                len_gene_block = len(value[1])
                out.write(f"{key}\t{value[0]}\t{len_gene_block}\t{value[1]}\t{value[2]}\t{value[3]}\t{value[4]}\n")


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
    def determine_template(n_samples: int):
        if n_samples <= 96:
            return template_path_96, "eblocks-plate-upload-template-96-filled.xlsx"
        elif n_samples > 96 and n_samples <= 384:
            return template_path_384
        elif n_samples > 384:
            return template_path_384, "eblocks-plate-upload-template-384-filled.xlsx" 


    @staticmethod
    def check_mut_format(mut: str, mut_type: str) -> bool:
        valid = 'acdefghiklmnpqrstvwy'
        # TODO Mutation and Combined is hardcoded here
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

        if (len(sequence) == 1) and (DesignEblocks.check_type_input_sequence(sequence)) and (DesignEblocks.check_for_start_stop_codon(sequence)):
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
    def eblock_colors() -> dict:
        """
        Create dictionary with colors for plotting eBlocks
        """
        return {i: '#%06X' % random.randint(0, 0xFFFFFF) for i in range(100)}
    

    @staticmethod
    def mutation_type_colors() -> dict:
        """
        Create dictionary with colors for plotting mutation types

        Returns:
            dict: dictionary with colors for plotting mutation types
        """
        return {'Mutation': 'black', 'Insert': 'red', 'Deletion': 'blue', 'Combined': 'green'}

