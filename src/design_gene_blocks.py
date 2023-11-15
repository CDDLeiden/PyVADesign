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
        
        self.dna_seq = self.read_seq(sequence_fp)
        self.mutations, self.mutation_types = self.read_mutations(mutations_fp)  # Types = {Mutation, Insert, Deletion}
        self.codon_usage = self.check_existance_codon_usage_table()  # Check if codon usage table is present for selected species

        # TODO ADD THE PLOT HERE WITHOUT THE EBLOCKS, JUST THE MUTATIONS

    def run(self, show=False):
        """
        Run the design of eBlocks
        """

        # Find indexes in sequence where mutation occures
        # TODO Remove/rename idx test
        # TODO Change this function so that you obtain for deletions and insertions all indexes (so the full range of mutations)
        # Paired = double mutation and insertion
        idx_dna, idx_dna_tups, idx_test, paired, idx_all, double_mut_pairs = self.index_mutations(self.mutations, self.mutation_types)

        # print("idx_dna_tups", len(idx_dna_tups), idx_dna_tups)
        # print("idx_dna", len(idx_dna), idx_dna)
        # print("idx_test", len(idx_test), idx_test)
        # print("paired", len(paired), paired)
        # print("idx_all", len(idx_all), idx_all)
        # print("double_mut_pairs", len(double_mut_pairs), double_mut_pairs)
        
        # length of ebLock is already checked here and should be within bounds
        # TODO Check here upon addition of the insert whether the gene block is still within bounds
        clusters = self.make_clusters(idx_all, idx_dna_tups, idx_test, paired) # or 'amount', 'cost'
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

        # for k, v in self.gene_blocks.items():
        #     print(k, len(v))

        # Plot legend
        legend = self.plot_legend(show=show)
        legend.savefig(os.path.join(self.output_fp, 'legend.png'), dpi=100)

        # Make histogram with bins
        # TODO Count double mutations in single eblock as single mutations in the histogram
        counts = self.count_mutations_per_eblock(idx_dna_tups)

        # TODO See if the total number of mutations adds to the value of the histogram
        num_mutations = sum(counts.values())
        # print("num mutations:", num_mutations)

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
        for key, value in self.gene_blocks.items():
            pass


        # TODO SIZES ALREADY HAVE BEEN CHECKED, SO NO NEED TO CHECK AGAIN HERE

        results = {}
        should_restart = True
        while should_restart:
            should_restart = False
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
                    # TODO This part should be moved more to the beginning
                    if self.check_eblock_length(mut_gene_block):
                        # Store output in dictionary
                        # TODO Check here upon addition of the insert whether the gene block is still within bounds
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
                    # TODO Check here upon addition of the insert whether the gene block is still within bounds
                    if self.check_eblock_length(mut_gene_block):
                        # Store output in dictionary
                        results[mut] = [mut_gene_block_name, mut_gene_block, idx, '', mut_type]  # TODO Maybe do it based on a string to remove 

                elif mut_type == self.type_combined:

                    # print(num, mut, mut_type)
                    # print("idx_dna_tups", idx_dna_tups[num])

                    # Find gene block and index of insert/deletion/mutation
                    # TODO Assuming here that there are at least some mutations not overlapping with each other in different eblocks
                    mut_gene_block_name = None
                    mut_gene_block_value = None

                    for mut_i in idx_dna_tups[num]:
                        possible_gene_blocks, counts = self.find_gene_block(self.gene_blocks, mut_i[1])
                        # print("possible_gene_blocks", possible_gene_blocks, counts)
                        if counts == 1:
                            mut_gene_block_name = list(possible_gene_blocks.keys())[0]
                            mut_gene_block_value = possible_gene_blocks[mut_gene_block_name]

                    idxs = []
                    codons = []
                    
                    for mut_i in idx_dna_tups[num]:

                        # print(mut_i, idx_dna_tups[num])
                        # print("mut_gene_block_name:", mut_gene_block_name)
                        # print("mut_gene_block_value:", mut_gene_block_value)

                        idx = self.find_mutation_index_in_gene_block(mut_gene_block_name, mut_i[1])
                        # print(idx)

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

        # Restore stdout to the original value
        # TODO FIX THIS
        # sys.stdout = sys.__stdout__
        # self.logfile.close()

    def count_mutations_per_eblock(self, idx_dna_tups):
        # Count the number of mutations in each eblock
        counts = {}
        for i in self.gene_blocks.keys():
            bbegin = int(i.split('_')[3])
            bend = int(i.split('_')[4])
            for idx in idx_dna_tups:
                if type(idx[1]) == int:
                    if bbegin < idx[1] < bend:
                        if i in counts.keys():
                            counts[i] += 1
                        else:
                            counts[i] = 1
                elif type(idx[1]) == list:  # just take the first mutation of a combined mutation, as they will all be in the same eblock
                    if bbegin < idx[1][1] < bend:
                        if i in counts.keys():
                            counts[i] += 1
                        else:
                            counts[i] = 1
        return counts

    def clean_gene_blocks(self, idx_dna, bins, gene_blocks):
        """
        Delete gene blocks that are not used for mutations
        """
        counts = {}
        for num, bin in enumerate(bins):
            if num < len(bins) - 1:
                counts[f"Block_{num}_pos_{str(bin)}_{str(bins[num + 1])}"] = 0
        
        for idx in idx_dna:
            for key, value in counts.items():
                if int(key.split('_')[3]) < idx < int(key.split('_')[4]):
                    counts[key] += 1

        # Delete gene blocks that are not used for mutations
        for key, value in counts.items():
            if value == 0:
                del gene_blocks[f'{key}']

        # Rename the gene blocks (starting from 1)
        nkeys = []
        nkey_gene_blocks = {}
        count = 1
        for key in gene_blocks.keys():
            splitkey = key.split('_')
            nkey = f'Block_{count}_pos_{splitkey[3]}_{splitkey[4]}'
            nkeys.append(nkey)
            count += 1
        for key, nkey in zip(gene_blocks.keys(), nkeys):
            nkey_gene_blocks[nkey] = gene_blocks[key]
            
        return nkey_gene_blocks, counts

    def check_wt_codon(self, gene_block_value, idx, mut):
        result = None
        for key, value in DNA_Codons.items():
                if key.lower() == gene_block_value[idx-3:idx]:
                    result = value
        if not result == mut[0]:
            print( f"WT codon does not match residue {mut}, but is {result}, the codon is {gene_block_value[idx-3:idx]}")
            print(f"This is probably due to the fact that paired mutations are not in same eBlock")
            sys.exit()

    def check_existance_codon_usage_table(self):
        codon_usage_present = [file for file in os.listdir(self.codon_usage_fp)]
        organisms_present = [file.split('.')[0:-1] for file in codon_usage_present]
        organisms_present = [i[0] for i in organisms_present]
        organisms_present_format = []
        for i in organisms_present:
            inew = i.split('_')
            inew = ' '.join([i for i in inew])
            organisms_present_format.append(inew.lower())
        if self.species.lower() in organisms_present_format:
            i = organisms_present_format.index(self.species.lower())
            return os.path.join(self.codon_usage_fp, codon_usage_present[i])
        else:
            print( "It looks like the codon usage table for the specified organism is not present.")
            sys.exit()

    def check_eblock_length(self, mut_gene_block):
            # Make sure that the codon insert is not too long or too short for eBlock
            if len(mut_gene_block) > self.idt_max_length_fragment:
                print( "Codon insert is too long for eBlock")
                sys.exit()
            elif len(mut_gene_block) < self.idt_min_length_fragment:
                print( f"Codon insert is too short for mutation eBlock, length is {len(mut_gene_block)}, minimum length is {self.idt_min_length_fragment}")
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
        # TODO CLEAN THIS FUNCTION
        # TODO Add check for insertions and deletions (so that the full range of mutations is added)
        idx_dna, idx_dna_tups, idx_test, paired, idx_all, double_mut_pairs = [], [], [], [], [], []
        for mut, type in zip(mut_list, mut_types):
            if type == self.type_mutation:
                idx_dna.append(int(mut[1:-1]) * 3)  # A residue consists of 3 nucleotides
                idx_dna_tups.append([mut, int(mut[1:-1]) * 3])
                idx_test.append(int(mut[1:-1]) * 3)
            elif (type == self.type_insert):
                mut_i = mut.split('-')[0]
                idx_dna.append(int(mut_i[1:]) * 3)  # A residue consists of 3 nucleotides
                idx_dna_tups.append([f"{mut}", int(mut_i[1:]) * 3])
                idx_test.append(int(mut_i[1:]) * 3)
            elif (type == self.type_deletion):
                mut_b = mut.split('-')[0]
                mut_e = mut.split('-')[1]
                idx_dna.append(int(mut_b[1:]) * 3)
                idx_dna_tups.append([f"{mut}", int(mut_b[1:]) * 3])
                length = int(mut_e[1:]) - int(mut_b[1:])
                idx_begin = int(int(mut_b[1:])*3)
                idx_end = idx_begin + (int(length)*3)
                tmp = range(int(mut_b[1:]) * 3, idx_end, 3)
                tmp_list = list(tmp)
                tmp_tuple = tuple(tmp_list)
                idx_test.append(tmp_list)
                paired.append(tmp_tuple)
            elif type == self.type_combined:
                muts = mut.split('-')
                tmp = []
                for i in muts:
                    idx_dna.append(int(i[1:-1]) * 3)  # A residue consists of 3 nucleotides
                    tmp.append([i, int(i[1:-1]) * 3])
                idx_dna_tups.append(tmp)
                tmp = []
                for i in muts:
                    tmp.append(int(i[1:-1]) * 3)
                idx_test.append(tmp)
                tmp_tuple = tuple(tmp)
                paired.append(tmp_tuple)
                double_mut_pairs.append(tmp_tuple)

        for i in idx_test:
            if isinstance(i, list):
                idx_all.extend(i)
            else:
                idx_all.append(i)

        return idx_dna, idx_dna_tups, idx_test, paired, idx_all, double_mut_pairs
    
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
            _, idx_dna_tups, _, _, _, _ = self.index_mutations(self.mutations, self.mutation_types)

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
        
    # def optimize_bins(self, x, bandwidth=200):
    #     """
    #     Optimize the bins using a meanshift algorithm

    #     Args:
    #         x (list): indexes of mutations

    #     Returns:
    #         optimal_bandwidth (int): 
    #         lowest_cost (float): estimated costs of all gene blocks together
    #     """    

    #     lowest_cost = np.inf
    #     optimal_bandwidth = bandwidth
        
    #     for i in range(self.n_iterations):

    #         # Start with default value
    #         clusters = self.meanshift(x, bandwidth)

    #         # Calculate costs and check size of fragments
    #         cost = self.calculate_cost(clusters)
    #         new_bandwidth = self.check_fragment_sizes(clusters, bandwidth)
                   
    #         if bandwidth == new_bandwidth:
    #             if lowest_cost > cost:
    #                 lowest_cost = cost
    #                 optimal_bandwidth = new_bandwidth
                
    #             ops = (add, sub)
    #             operation = random.choice(ops)
    #             random_int = random.randint(1, 50)
    #             bandwidth = operation(bandwidth, random_int)
    #         else:
    #             bandwidth = new_bandwidth

    #     return optimal_bandwidth, lowest_cost
    
    # def dbscan_clustering(self, data, epsilon):
    #     # Flatten the list of ranges to individual numbers
    #     print(data)
    #     flattened_data = [item if isinstance(item, int) else item for sublist in data for item in (sublist if isinstance(sublist, list) else [sublist])]
    #     print(flattened_data)
                
    #     X = np.array(flattened_data).reshape(-1, 1)
    #     dbscan = DBSCAN(eps=epsilon, min_samples=2, metric='euclidean')
    #     labels = dbscan.fit_predict(X)

    #     clusters = {}
    #     for idx, label in enumerate(labels):
    #         if label not in clusters:
    #             clusters[label] = []
    #         clusters[label].append(flattened_data[idx])

    #     return clusters
    
    def make_clusters(self, idxs, idxs_tuple, idx_test, paired_mutations):
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
             
        # print("length cluster labels:", len(cluster_labels))

        # print("length idx_first:", len(idx_first))

        # Extract cluster centers
        # Visualize the clusters
        # TODO Check what to do with this. Can be removed I guess
        # if visualize:
        #     plt.scatter(mutation_matrix, np.zeros_like(mutation_matrix),
        #                 c=cluster_labels, cmap='viridis')
        #     plt.scatter(cluster_centers, np.zeros_like(cluster_centers), c='red', marker='X', s=100, label='Cluster Centers')
        #     plt.xlabel('Mutation Indices')
        #     plt.title('K-means Clustering with Connected Mutations')
        #     plt.legend()
        #     plt.show()

        return cluster_labels, idx_first



    # def meanshift(self, x: list, bandwidth: int):
    #     """
    #     Meanshift algorithm for finding clusters of mutations that fit in a gene block

    #     Args:
    #         x (list): _description_
    #         bandwidth (int): _description_

    #     Returns:
    #         clusters (dict): _description_
    #     """    
    #     # https://stackoverflow.com/questions/18364026/clustering-values-by-their-proximity-in-python-machine-learning
    #     X = np.array(list(zip(x, np.zeros(len(x)))), dtype=np.int64)
    #     bandwidth = bandwidth
    #     ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    #     ms.fit(X)
    #     labels = ms.labels_
    #     labels_unique = np.unique(labels)
    #     clusters = {}
    #     for label in labels_unique:
    #         clusters[f'cluster {label}'] = []
    #     for num, i in enumerate(labels):
    #         clusters[f'cluster {i}'].append(x[num])
    #     return clusters

    def renumber_gene_blocks(self, gene_blocks):
        """
        """
        new_gene_blocks = {}
        first_positions = []
        for k, v in gene_blocks.items():
            pos1 = k.split('_')[3]
            first_positions.append(int(pos1))
        # Sort the first positions
        first_positions = sorted(first_positions)


        for num, i in enumerate(first_positions, 1):
            for k, v in gene_blocks.items():
                if int(k.split('_')[3]) == i:
                    new_gene_blocks[f'Block_{num}_pos_{k.split("_")[3]}_{k.split("_")[4]}'] = v
        return new_gene_blocks

    def make_gene_block(self, bins, dna_sequence):
        gene_blocks = {}
        num = 0
        block_num = 1
        while num < len(bins):
            # print(f"Creating gene block {block_num} ...")
            name = f'Block_{str(block_num)}_pos_{bins[num]}_{bins[num+1]}'
            block = dna_sequence[bins[num]:bins[num+1]]
            gene_blocks[name] = str(block)
            num += 2
            block_num += 1
        gene_blocks = self.renumber_gene_blocks(gene_blocks)
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
        count = 0
        results = {}
        for key, value in gene_blocks.items():
            begin_range, end_range = self.gene_block_range(key)
            if begin_range < int(mutation_idx) < end_range:
                results[key] = value
                count += 1
        return results, count
                         
    def find_mutation_index_in_gene_block(self, gene_block, idx_mutation):
        begin_range, _ = self.gene_block_range(gene_block)
        # Find index of mutation within geneblock
        index = idx_mutation - begin_range
        # Check that mutation is not in the first or final X residues of the gene block (this would make it very difficult to design IVA primers)
        # self.check_position_gene_block(gene_block, index)
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

    def extract_wells_from_template(self, template):
        df = pd.read_excel(template)
        wells = df['Well Position'].tolist()
        return wells

    def write_gene_blocks_to_template(self, gene_block_dict, outpath):
        # TODO Change this function when more than 96 eBlocks are designed
        # 384-well plate
        # Count number of samples
        n_samples = len(gene_block_dict)
        if n_samples <= 96: # Use 96-well plate
            fname="eblocks-plate-upload-template-96-filled.xlsx"
            n_plates = 1
            template = template_path_96
        elif n_samples > 96 and n_samples <= 384: # Use 384-well plate
            fname="eblocks-plate-upload-template-384-filled.xlsx"
            n_plates = 1
            template = template_path_384
        elif n_samples > 384:
            n_plates = math.ceil(n_samples / 384)
            fnames = [f"eblocks-plate-upload-template-384-filled_{i}.xlsx" for i in range(1, n_plates + 1)]
            template = template_path_384

        if n_plates == 1:
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
            wells = self.extract_wells_from_template(template=template)
            wells = wells[:len(names)]
            df = pd.DataFrame()
            df['Well Position'] = wells
            df['Name'] = names
            df['Sequence'] = seqs
            df.to_excel(outfile, index=False)

        elif n_plates > 1:
            # count = 1
            max_samples = 384
            for num, i in enumerate(fnames, 1):
                # print(num, i, count)
                outfile = os.path.join(outpath, i)
                names = []
                seqs = []
                for index, (key, value) in enumerate(gene_block_dict.items(), 1):
                    if index <= (max_samples * num) and index >= (max_samples * (num - 1) + 1):
                        mutation = key
                        block = value[0].split('_')[0:2]
                        block = block[0] + '-' + block[1]
                        name = mutation + '_' + block
                        names.append(name)
                        seqs.append(value[1])
                        # count += 1
                wells = self.extract_wells_from_template(template=template)
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
    
    @staticmethod
    def eblock_colors():
        """
        Create dictionary with colors for plotting eBlocks
        """
        # TODO Update this to nice colors, maybe specific color palette from sns or so
        colors = {}
        for i in range(100):
            colors[i] = '#%06X' % random.randint(0, 0xFFFFFF)
        return colors
    
    @staticmethod
    def mutation_type_colors():
        """
        Create dictionary with colors for plotting mutation types

        Returns:
            dict: dictionary with colors for plotting mutation types
        """
        colors = {}
        colors['Mutation'] = 'black'
        colors['Insert'] = 'red'
        colors['Deletion'] = 'blue'
        colors['Combined'] = 'green'
        return colors

