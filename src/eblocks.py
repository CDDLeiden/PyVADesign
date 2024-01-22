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
# from dna_features_viewer import GraphicFeature, GraphicRecord
# from utils import read_codon_usage, DNA_Codons, write_pickle, natural_amino_acids

# from mutation import Mutation
# TODO Change filepath (str) to Path object

from mutation import Mutation
from sequence import Sequence
from utils import Utils

class Eblocks:
    def __init__(self):
        self.eblock_parameters = {"common_param": "shared_value"}
        self.block_sequences = []

    def set_parameters(self, parameters):
        # TODO Store the information for the designed eblocks here (mutation, eblock nun, eblock sequences, etc.)
        self.eblock_parameters = parameters

    def set_block_sequences(self, block_sequences):
        self.block_sequences = block_sequences

    def display_parameters(self):
        print("Eblocks Parameters:", self.eblock_parameters)
        print("Block Sequences:", self.block_sequences)


class EblockDesign:
    def __init__(self, 
                 eblocks_instance: Eblocks,
                 mutation_instance: Mutation,
                 sequence_instance: Sequence,
                 # File paths for input files
                 output_fp: str = None,
                 verbose = True,
                 
                # IDT parameters
                 bp_price: float = 0.05,
                 max_eblock_length: int = 1500,
                 min_eblock_length: int = 300,
                 min_overlap: int = 25,
                 min_order: int = 24,
                 optimization_method = "cost"
                ):
        
        self.eblocks_instance = eblocks_instance
        self.mutation_instance = mutation_instance
        self.sequence_instance = sequence_instance
        self.verbose = verbose

        # IDT parameters
        self.bp_price = bp_price
        self.max_eblock_length = max_eblock_length
        self.min_eblock_length = min_eblock_length
        self.min_overlap = min_overlap
        self.min_order = min_order
        self.optimization_method = optimization_method


        self.block_sequences = []

    def run_design_eblocks(self):
        """
        This function runs the design process for eblocks.
        """

        valid_clusters = self.find_possible_clusters()
        optimal_clustering = self.choose_clusters(valid_clusters)
        print(optimal_clustering)

        # Define the beginning and end of each gene block, based on the clusters and include the minimum overlap
        bins = self.make_bins(optimal_clustering)

        # Make gene blocks (WT DNA sequences cut to the correct size, according to the bins) and renumber them starting from 1
        self.wt_eblocks = self.make_wt_eblocks(bins)

        # TODO ADD THE PLOTTING FUNCTIONALITY SOMEWHERE

        # Loop over all mutations and create the eBlocks
        results = {}  # TODO STORE IN EBLOCK CLASS
        for num, mutation in enumerate(self.mutation_instance.mutations):
            # Create mutated eblock, based on mutation type
            results = self.make_mutant_eblock(num, mutation, results)
                                        
        # Store output
        # TODO 
        # self.write_gene_blocks_to_txt(results, self.output_fp)
        # self.write_gene_blocks_to_template(results, self.output_fp)
        # write_pickle(results, self.output_fp)
        # write_pickle(self.wt_eblocks, self.output_fp, fname="wt_gene_blocks.npy")

        # print("Designed eBlocks and stored output in ", self.output_fp)

        # Count number of mutations per eblock, for making barplot
        # TODO 
        # self.counts = self.count_mutations_per_eblock(idx_dna_tups)

        # # Your design logic to generate block_sequences
        # self.block_sequences = ["sequence1", "sequence2", "sequence3"]

        # # Set the block_sequences in the Eblocks instance
        # self.eblocks_instance.set_block_sequences(self.block_sequences)

    def find_possible_clusters(self):        
        possibilities = {} # Store all possible clusterings
        n = 1
        valid_clusters = True

        while valid_clusters:
            
            clusters = {}
            cluster_labels, idxs_reordered = self.kmeans_clustering(num_clusters=n)

            # Calculate the size of each cluster
            for i, j in zip(cluster_labels, idxs_reordered):
                if i not in clusters:
                    clusters[i] = []
                clusters[i].append(j)

            # Check if the size of the clusters is within bounds
            cluster_sizes = [max(v) - min(v) for v in clusters.values()]

            # Find in which cluster insertions and deletions are located and check if the size of the cluster is within bounds when the mutation is added
            for mutation in self.mutation_instance.mutations:
                if mutation.type == "Insert":
                    for key, value in clusters.items():
                        if mutation.idx_dna in value:
                            cluster_sizes[key] += mutation.length_insert
                elif mutation.type == "Deletion":
                    for key, value in clusters.items():
                        if mutation.idx_dna in value:
                            cluster_sizes[key] -= mutation.length_deletion

            max_cluster_size = max(cluster_sizes)
            min_cluster_size = min(cluster_sizes)

            if max_cluster_size > (self.max_eblock_length - 2 * self.min_overlap): # Take into account the OH on both sides of the gene block, so decrease max size with 2 times min length of OH
                n += 1
            elif min_cluster_size < (self.min_eblock_length - 2 * self.min_overlap):
                valid_clusters = False
            else:
                possibilities[f'cluster N={n}'] = clusters  
                n += 1
        
        if len(possibilities) == 0:  # No valid clusters found
            print("No valid clusterings found. Please check your input mutations and make sure that the multi-mutants are not too far apart.")
            sys.exit()
        return possibilities
    
    def kmeans_clustering(self, num_clusters: int, n_init=10, random_state=42):
        # TODO FIX THIS FUNCTION!!! WHAT IF TWO SAME MULTIP MUTATIONS POSISITONS E.G> (L305I, L807I) (L305G, L807G) !!
        # Extract the first index of each mutation in idx_test
        mean_idxs_dna = self.mutation_instance.mean_idxs_dna()
        # Create a numpy array for clustering
        mutation_arr = np.array(mean_idxs_dna, dtype=float).reshape(-1, 1)
        # Initialize KMeans with the number of clusters
        kmeans = KMeans(n_clusters=num_clusters, random_state=random_state, n_init=n_init)
        # Fit the model and obtain cluster labels for connected indices
        cluster_labels = kmeans.fit_predict(mutation_arr)
        # Add the remaining indices of the pairs to the assigned cluster
        indices = []
        for mutation in self.mutation_instance.mutations:
            if len(mutation.idx_dna) > 1:
                pair = mutation.idx_dna
                mean_pair_value = np.mean(pair)
                idx = list(mutation_arr).index(mean_pair_value)
                indices.append(idx)
                cluster_label = cluster_labels[idx]
                cluster_labels = np.concatenate((cluster_labels, np.full(len(pair), cluster_label)))
                mean_idxs_dna.extend(pair)
        # Remove the mean values from the list of indices and cluster labels
        to_remove = sorted(set(indices), reverse=True)
        for i in to_remove:
            del mean_idxs_dna[i]
            cluster_labels = np.delete(cluster_labels, i)
        return cluster_labels, mean_idxs_dna
    
    def choose_clusters(self, clusters: dict) -> dict:
        if self.optimization_method == "cost":
            print("Optimizing based on price per bp ...")
            lowest_cost, best_clustering = min((self.calculate_cost(value), value) for value in clusters.values())
            print(f"Lowest cost: {lowest_cost} with cluster {best_clustering}")
            return best_clustering
        elif self.optimization_method == "amount":
            print("Optimizing based on number of eBlocks ...")
            fewest_blocks, best_clustering = min((len(value), value) for value in clusters.values())
            print(f"Fewest blocks: {fewest_blocks} with cluster {best_clustering}")
            return best_clustering

    def calculate_cost(self, clusters: dict) -> float:
        total_cost = 0
        for _, value in clusters.items():
            min_val = min(value)
            max_val = max(value)
            len_gene_block = (max_val - min_val) + 2 * self.min_overlap  # on both size of the gene block there should be a number of non-mutated basepairs for IVA primer design
            cost = len_gene_block * self.bp_price * len(value)
            total_cost += cost
        return round(total_cost, 2)
    
    def make_bins(self, clusters: dict):
        bins = []
        for _, value in clusters.items():
            bins.append(int(min(value) - self.min_overlap))
            bins.append(int(max(value) + self.min_overlap))
        return bins
    
    def make_wt_eblocks(self, bins):
        gene_blocks = {}
        block_num = 1
        for num in range(0, len(bins), 2):
            name = f'Block_{block_num}_pos_{bins[num]}_{bins[num+1]}'
            block = self.sequence_instance.sequence[bins[num]:bins[num+1]]
            gene_blocks[name] = str(block)
            block_num += 1
        return self.renumber_gene_blocks(gene_blocks)

    def renumber_gene_blocks(self, gene_blocks):
        renumbered_gene_blocks = {}
        for num, (k, v) in enumerate(sorted(gene_blocks.items(), key=lambda item: int(item[0].split('_')[3])), 1):
            renumbered_gene_blocks[f'Block_{num}_pos_{k.split("_")[3]}_{k.split("_")[4]}'] = v
        return renumbered_gene_blocks
    
    def find_eblock_index(self, gene_block, idx_mutation: int) -> int:
        begin_range, _ = self.gene_block_range(gene_block)
        return idx_mutation - begin_range
    
    def map_mutation_to_eblock(self, dna_idx: int):
        name_val, _ = self.find_gene_block(self.gene_blocks, dna_idx)
        eblock_name = list(name_val.keys())[0]
        eblock_value = name_val[eblock_name]
        return eblock_name, eblock_value

    def extract_mut_codons(self, res: str):
        return [key.lower() for key, value in Utils.DNA_Codons.items() if value == res]
    
    def check_wt_codon(self, eblock_seq: str, idx: int, mut: str):
        """
        This function checks if the WT codon at the mutation index is the same as in the mutation.
        """
        codon = eblock_seq[idx-3:idx]
        result = next((value for key, value in DNA_Codons.items() if key.lower() == codon), None)
        if result is not None and result != mut[0]:
            print(f"WT codon does not match residue {mut}, but is {result}, the codon is {codon}")
            print("This is probably due to the fact that paired mutations are not in the same eBlock")
            sys.exit()

    def select_mut_codon(self, codon_list):
        """
        Choose codon from list of codons based on occurrence of codon in nature.
        """
        codon_dict = read_codon_usage(fp=self.codon_usage)
        most_occuring_codon = max(codon_list, key=codon_dict.get, default='xxx')
        return most_occuring_codon

    
    def make_mutant_eblock(self, num: int, mutation, mut_type: str, idx_dna_tups: list, results: dict) -> dict:
        """
        This function creates the mutated eBlock, based on the WT eBlock and the mutation.
        """
        # Change gene block in selected position for mutation
        if mutation.type == "Mutation":
            # Find gene block
            mut_gene_block_name, mut_gene_block_value = self.map_mutation_to_eblock(mutation.idx_dna)
            eblock_index = self.find_eblock_index(mut_gene_block_name, mutation.idx_dna)
            # Check if WT codon at index is same residue as mutation
            self.check_wt_codon(mut_gene_block_value, mutation.idx_dna, mutation.mutation[0])
            mut_codons = self.extract_mut_codons(mutation.mutation[-1])
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
            idx = self.find_eblock_index(mut_gene_block_name, mut_idx)
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
            
            idx = self.find_eblock_index(mut_gene_block_name, idx_del_start)
            idx_end = self.find_eblock_index(mut_gene_block_name, idx_del_end)

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
                            idx = self.find_eblock_index(mut_gene_block_name, mut_i[1])
                            if (idx < self.min_bin_overlap) or (idx > (len(mut_gene_block_value) - self.min_bin_overlap)):
                                raise Exception("Mutation too close to beginning or end of eBlock")
                    except Exception:
                        continue

            idxs, codons = [], []
            for mut_i in idx_dna_tups[num]:

                idx = self.find_eblock_index(mut_gene_block_name, mut_i[1])
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