import os
import sys
import copy
import numpy as np
from sklearn.cluster import KMeans

from .mutation import Mutation
from .sequence import Plasmid
from .utils import Utils, SnapGene


class Eblock:
    """
    This class contains information of 
    """
    def __init__(self):

        self.name: str = None
        self.sequence: str = None
        self.start_index: int = None
        self.end_index: int = None

        self.mutation: str = None
        self.mutation_start_index: int = None
        self.mutation_end_index: int = None
        self.wt_codon: str = None
        self.mutant_codon: str = None
        self.insert: str = None

        # Multiple mutations
        self.mutation_idxs = []
        self.mutant_codons = []

    def set_name(self, name: str):
        self.name = name

    def set_wt_sequence(self, sequence: str):
        self.wt_sequence = sequence

    def set_mutant_sequence(self, sequence: str):
        self.mutant_sequence = sequence

    def set_start_index(self, start_index: int):
        self.start_index = start_index

    def set_end_index(self, end_index: int):
        self.end_index = end_index



class EblockDesign:
    """
    This class contains functions to design eBlocks based on mutations in a gene sequence.
    """
    def __init__(self, 
                 mutation_instance: Mutation,
                 sequence_instance: Plasmid,
                 output_dir: str = None,

                 verbose = True,
                 
                # IDT parameters
                 bp_price: float = 0.05,
                 max_eblock_length: int = 1500,
                 min_eblock_length: int = 300,
                 min_overlap: int = 25,
                 min_order: int = 24,
                 optimization_method = "cost",

                # TODO Change this
                 codon_usage: str = r"C:\Users\Rosan\Documents\git\my_repositories\design_gene_blocks\src\data\codon_usage\Escherichia_coli.csv",
                ):
        
        self.mutation_instance = mutation_instance
        self.sequence_instance = sequence_instance
        # self.eblocks_instance = eblocks_instance
        self.verbose = verbose
        self.output_dir = output_dir

        # IDT parameters
        self.bp_price = bp_price
        self.max_eblock_length = max_eblock_length
        self.min_eblock_length = min_eblock_length
        self.min_overlap = min_overlap
        self.min_order = min_order
        self.optimization_method = optimization_method

        self.codon_usage = codon_usage
        self.block_sequences = []
        self.eblock_colors = self.generate_eblock_colors()
        self.snapgene_files = True

        self.wt_eblocks = {}
        self.eblocks = {}

    def run_design_eblocks(self):
        """
        This function runs the design process for eblocks.
        """

        print("Starting eBlock design ...")

        valid_clusters = self.find_possible_clusters()
        optimal_clustering = self.choose_clusters(valid_clusters)

        # Define the beginning and end of each gene block, based on the clusters and include the minimum overlap
        bins = self.make_bins(optimal_clustering)

        # Make gene blocks (WT DNA sequences cut to the correct size, according to the bins) and renumber them starting from 1
        wt_eblocks = self.make_wt_eblocks(bins)
        self.wt_eblocks = wt_eblocks

        # Loop over all mutations and create the eBlocks
        results = {}
        for mutation in self.mutation_instance.mutations:
            results = self.make_mutant_eblock(mutation, results)  # Create mutated eblock, based on mutation type

        # Sort the eblocks based on the index of the first mutation in the eblock and the number of the eblock
        sorted_dict = dict(sorted(results.items(), key=lambda x: (x[1].name, x[1].start_index)))
        self.eblocks = sorted_dict

        # Create a GFF3 file for easy visualization of eBlocks in SnapGene
        if self.snapgene_files:
            snapgene_dict = {}
            for i in self.wt_eblocks:
                snapgene_dict[i.name] = [i.start_index, i.end_index, self.eblock_colors[int(i.name[-1])-1]]

            snapgene_instance = SnapGene(output_dir=self.output_dir, 
                                        sequence_instance=self.sequence_instance)
            snapgene_instance.eblocks_to_gff3(eblocks=snapgene_dict)

            # TODO ADD GENE DnaE1 (or whatever to features in snapgene as full gene)
            # TODO SAVE INDUVIDUAL MUTATIONS AS WELL HERE
                                    
        print("Completed eBlock design.")

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
            print(f"Lowest estimated cost: €{lowest_cost} (given price per bp of €{self.bp_price})")
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
    
    def make_wt_eblocks(self, bins) -> list:
        """
        Create wild-type eBlocks
        """
        gene_blocks = []
        block_num = 1
        for num in range(0, len(bins), 2):
            eblock = Eblock()
            eblock.start_index = bins[num]
            eblock.end_index = bins[num+1]
            eblock.sequence = self.sequence_instance.sequence[bins[num]:bins[num+1]]
            gene_blocks.append(eblock)
            block_num += 1
        gene_blocks = self.renumber_eblock(gene_blocks)
        return gene_blocks
        
    def eblock_index(self, eblock: Eblock, idx_mutation: int) -> int:
        """
        Find the index of a mutation in an eblock.
        """
        return idx_mutation - eblock.start_index
        
    def eblocks_within_range(self, mutation_idx: int):
        """
        For a given mutation index, find the number of overlapping eBlocks
        """
        eblocks = []
        for eblock in self.wt_eblocks:
            if eblock.start_index  < int(mutation_idx) < eblock.end_index:
                eblocks.append(eblock)
        count = len(eblocks)
        return copy.deepcopy(eblocks[0]), count  # Copy to avoid editing original class
    

    def check_wt_codon(self, eblock: Eblock, mut: str):
        """
        This function checks if the WT codon at the mutation index is the same as in the mutation.
        """
        all_codons = Utils.DNA_codons()
        codon = eblock.sequence[eblock.mutation_start_index-3:eblock.mutation_start_index]
        result = next((value for key, value in all_codons.items() if key.lower() == codon), None)
        if result is not None and result != mut[0]:
            print(f"WT codon does not match residue {mut}, but is {result}, the codon is {codon}")
            print("This is probably due to the fact that paired mutations are not in the same eBlock")
            sys.exit()

    def select_mut_codon(self, res: str):
        all_codons = Utils.DNA_codons()
        possible_codons = [key.lower() for key, value in all_codons.items() if value == res]
        codon_dict = Utils.read_codon_usage(fp=self.codon_usage)
        selected_codon = max(possible_codons, key=codon_dict.get, default='xxx')
        return selected_codon
    
    def mutate_eblock(self, mutation, eblock: Eblock, end_idx=None):
        """
        Mutate gene block based on mutation type
        """
        if (mutation.is_singlemutation) or (mutation.is_multiplemutation):
            mutated_eblock = eblock.sequence[:eblock.mutation_start_index -3] + eblock.mutant_codon + eblock.sequence[eblock.mutation_start_index:]
        elif mutation.is_insert:
            mutated_eblock = eblock.sequence[:eblock.mutation_start_index] + eblock.insert + eblock.sequence[eblock.mutation_start_index:]
        elif mutation.is_deletion:
            mutated_eblock = eblock.sequence[:eblock.mutation_start_index -3] + eblock.sequence[end_idx -3:]
        return mutated_eblock
    
    def design_insert(self, aas):
        codon_insert = ''  # Sequence to insert in gene block
        for res in aas:
            codon = self.select_mut_codon(res)
            codon_insert += codon
        return codon_insert
    
    def check_eblock_length(self, eblock_seq: str) -> bool:
        """
        Check if the length of the gene block is within bounds
        """
        length_eblock = len(eblock_seq)
        if not self.min_eblock_length <= length_eblock <= self.max_eblock_length:
            if length_eblock > self.max_eblock_length:
                print("Codon insert is too long for eBlock")
            else:
                print(f"Codon insert is too short for mutation eBlock, length is {length_eblock}, minimum length is {self.min_eblock_length}")
                sys.exit()
    
    def count_mutations_per_eblock(self) -> dict:
        """
        Count the number of mutations in each eblock
        """
        counts = {}
        for i in self.wt_eblocks.keys():
            counts[i] = 0
        for _, val in self.eblocks.items():
            counts[val[0]] += 1
        return counts
    
    def make_mutant_eblock(self, mutation: Mutation, results: dict) -> dict:
        """
        This function creates the mutated eBlock, based on the WT eBlocks and the mutation.
        """

        if mutation.is_singlemutation:
            eblock, _ = self.eblocks_within_range(mutation.idx_dna[0])
            eblock.mutation_start_index = self.eblock_index(eblock, mutation.idx_dna[0])    
            self.check_wt_codon(eblock, mutation.mutation[0][0])  # Check if WT codon at index is same residue as mutation
            eblock.mutant_codon = self.select_mut_codon(mutation.mutation[0][-1])
            eblock.sequence = self.mutate_eblock(mutation, eblock)

        elif mutation.is_insert:
            eblock, _ = self.eblocks_within_range(mutation.idx_dna[0])  # Find gene block and index of insert
            eblock.mutation_start_index = self.eblock_index(eblock, mutation.idx_dna[0])
            
            self.check_wt_codon(eblock, mutation.mutation[0][0])
            eblock.insert = self.design_insert(mutation.insert)
            eblock.sequence = self.mutate_eblock(mutation, eblock)
            self.check_eblock_length(eblock.sequence)  # Check if eBlock is too long including the insert

        elif mutation.is_deletion:
            eblock, _ = self.eblocks_within_range(mutation.idx_dna_deletion_begin)
            eblock.mutation_start_index = self.eblock_index(eblock, mutation.idx_dna_deletion_begin)
            idx_end = self.eblock_index(eblock, mutation.idx_dna_deletion_end)
            
            self.check_wt_codon(eblock, mutation.mutation[0][0])  # Check if WT codon at index is same residue as mutation
            eblock.sequence = self.mutate_eblock(mutation, eblock, idx_end)
            self.check_eblock_length(eblock.sequence)  # Check if eBlock is too short
     
        elif mutation.is_multiplemutation:
            selected_eblock = None
            for mut_i in mutation.idx_dna:
                eblock, counts = self.eblocks_within_range(mut_i)
                if counts == 1:  # If only one eBlock for mutation > should be correct eBlock
                    selected_eblock = eblock

            all_counts = [counts for _, counts in (self.eblocks_within_range(mut_i) for mut_i in mutation.idx_dna)]
            lowest_count = min(all_counts)
            
            for mut_i in mutation.idx_dna:
                eblock, counts = self.eblocks_within_range(mut_i)
                if counts == lowest_count:
                    selected_eblock = eblock
                    try:  # Try to find indexes of mutations, based on eblock. Check if they are too close to beginning or end of eblock
                        for mut_i in mutation.idx_dna:
                            eblock.mutation_start_index = self.eblock_index(selected_eblock, mut_i)  # Check too close to beginning or end
                            if (selected_eblock.mutation_start_index < self.min_overlap) or (selected_eblock.mutation_start_index > (len(selected_eblock.sequence) - self.min_overlap)):
                                raise Exception("Mutation too close to beginning or end of eBlock")
                    except Exception:
                        continue

            for num_i, mut_i in enumerate(mutation.idx_dna):

                print(mutation.mutation[num_i])

                selected_eblock.mutation_start_index = self.eblock_index(selected_eblock, mut_i)
                selected_eblock.mutant_codon = self.select_mut_codon(mutation.mutation[num_i][-1])  # Find most occuring mutant codon based on codon usage for species
                print(selected_eblock.mutant_codon)

                selected_eblock.mutation_idxs.append(selected_eblock.start_index)
                selected_eblock.mutant_codons.append(selected_eblock.mutant_codon)
                
                self.check_wt_codon(selected_eblock, mutation.mutation[num_i][0])  # Check if WT codon at index is same residue as mutation
                selected_eblock.sequence = self.mutate_eblock(mutation, selected_eblock)

            eblock = selected_eblock

        print('Mutated eBlock', eblock.sequence)
        results[mutation] = eblock
        return results
    
    @staticmethod
    def renumber_eblock(eblocks: list):
        sorted_list = sorted(eblocks, key=lambda x: x.start_index)
        for i, obj in enumerate(sorted_list, start=1):
            obj.name = f"eBlock-{i}"
        return sorted_list
                
    @staticmethod
    def generate_eblock_colors() -> dict:
        """
        Create dictionary with colors for plotting eBlocks using the tab10 color scheme.
        """
        tab10_colors = ['#1f77b4','#ff7f0e','#2ca02c', '#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf',
                        '#aec7e8','#ffbb78','#98df8a','#ff9896','#c5b0d5','#c49c94','#f7b6d2','#c7c7c7','#dbdb8d','#9edae5',
                        '#393b79','#ff7f0e','#2ca02c','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf']
        return {i: tab10_colors[i] for i in range(len(tab10_colors))}