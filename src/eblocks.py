import os

os.environ["OMP_NUM_THREADS"] = '1'  # KMeans is not parallelized, so set to 1 thread

import sys
import copy
import numpy as np
from sklearn.cluster import KMeans

from .mutation import Mutation
from .sequence import Plasmid
from .utils import Utils, SnapGene

script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(script_dir, os.pardir))  # Navigate one directory upwards to reach the parent directory
data_dir = os.path.join(parent_dir, "data")  # Navigate to the data directory



class Eblock:
    """
    This class contains information about eBlocks and mutations in the eBlocks.
    """
    def __init__(self,
                 name: str = None,
                 block_number: int = None,
                 sequence: str = None,
                 start_index: int = None,
                 end_index: int = None,
                 mutation_start_index: int = None,
                 mutation_end_index: int = None,
                 wt_codon: str = None,
                 mutant_codon: str = None,
                 insert: str = None):
        
        self.name = name
        self.block_number = block_number
        self.sequence = sequence
        self.start_index = start_index
        self.end_index = end_index
        self.mutation_start_index = mutation_start_index
        self.mutation_end_index = mutation_end_index
        self.wt_codon = wt_codon
        self.mutant_codon = mutant_codon
        self.insert = insert

    def set_name(self, name: str):
        self.name = name

    def set_sequence(self, sequence: str):
        self.sequence = sequence

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

                 cost_optimization: bool = True,
                 amount_optimization: bool = False,
                 to_snapgene: bool = True,  # TODO Change name to more general
                 verbose: bool = True,
                 codon_usage: str = os.path.join(data_dir, "codon_usage", "Escherichia_coli.csv"),

                # IDT values
                 bp_price: float = 0.05,
                 max_eblock_length: int = 1500,
                 min_eblock_length: int = 300,
                 min_overlap: int = 25,
                 min_order: int = 24,
                ):
        
        self.mutation_instance = mutation_instance
        self.sequence_instance = sequence_instance
        self.output_dir = output_dir

        self.cost_optimization = cost_optimization
        self.amount_optimization = amount_optimization
        self.eblock_colors = self.generate_eblock_colors()
        self.to_snapgene = to_snapgene
        self.verbose = verbose
        self.codon_usage = codon_usage

        # IDT parameters
        self.bp_price = bp_price
        self.max_eblock_length = max_eblock_length
        self.min_eblock_length = min_eblock_length
        self.min_overlap = min_overlap
        self.min_order = min_order

        # Store WT and mutated gene blocks
        self.wt_eblocks: dict = {}  # Wild-type gene blocks
        self.eblocks: dict = {}  # Mutated gene blocks

        self.validate_optimization_parameters()

    def run_design_eblocks(self):
        """
        This function runs the design process for eBlocks.
        """

        self.print_line("Starting eBlock design ...")

        # Divide the target gene into clusters based on the mutations
        valid_clusters = self.find_possible_clusters()
        optimal_clustering = self.choose_clusters(valid_clusters)

        # Define the beginning and end of each gene block, based on the clusters and include the minimum overlap
        bins = self.make_bins(optimal_clustering)

        # Make gene blocks (WT DNA sequences sliced to the correct size, according to the bins) and renumber them starting from 1
        self.wt_eblocks = self.make_wt_eblocks(bins)

        # Loop over all mutations and create the eBlocks, based on the WT eBlocks
        results = {}
        for mutation in self.mutation_instance.mutations:
            results = self.make_mutant_eblock(mutation, results)  # Create mutated eBlock, based on mutation type
        sorted_dict = dict(sorted(results.items(), key=lambda x: (x[1].name, x[1].start_index)))  # Sort the eblocks based on the index of the first mutation in the eblock and the number of the eblock
        self.eblocks = sorted_dict

        # Create a GFF3 file for easy visualization of eBlocks in SnapGene
        # TODO paramter save_clones?
        if self.to_snapgene:
            self.output_to_snapgene()

        self.print_line("Completed eBlock design.")
                                             
    def find_possible_clusters(self, threshold_small_clusters=2):
        """
        This function finds all possible clusterings of the mutations based on the index of the mutations.
        """
        possibilities = {} # Store all possible clusterings
        n = 1
        valid_clusters = True
        while (valid_clusters) and (n <= len(self.mutation_instance.mutations)):
            
            clusters = {}
            cluster_labels, idxs_reordered = self.kmeans_clustering(num_clusters=n)
            clusters = self.cluster_size(clusters, cluster_labels, idxs_reordered)  # Add size of clusters
            cluster_sizes = [max(v) - min(v) for v in clusters.values()]  # Check if the size of the clusters is within bounds
            cluster_sizes = self.check_clusters(clusters, cluster_sizes) # Check if the size of the clusters is within bounds (inserts + deletions)

            clusters_copy = copy.deepcopy(clusters)
            cluster_too_small = 0
            cluster_correct = 0
            cluster_too_big = 0

            # Loop over cluster and increase size if necessary
            for k, v in clusters.items():

                size = max(v) - min(v)                    
                if size > (self.max_eblock_length - 2 * self.min_overlap): # Take into account the OH on both sides of the gene block, so decrease max size with 2 times min length of OH
                    cluster_too_big += 1

                elif size < (self.min_eblock_length - 2 * self.min_overlap):  # Cluster size too small, increasing the eBlock length to fit requirements    
                    min_required_length_toadd = (self.min_eblock_length - size) - 2 * self.min_overlap       
                    if max(v) + min_required_length_toadd <= len(self.sequence_instance.sequence):
                        clusters_copy[k].append(max(v) + min_required_length_toadd)
                        cluster_too_small += 1
                        cluster_correct += 1
                    else:
                        clusters_copy[k].append(min(v) - min_required_length_toadd)
                        cluster_too_small += 1
                        cluster_correct += 1

                else:
                    cluster_correct += 1  

            if (cluster_correct == len(clusters)) and (cluster_too_small < threshold_small_clusters):
                possibilities[f'cluster N={n}'] = clusters_copy
                n += 1
            elif cluster_too_big > 0:
                n += 1
            else:
                valid_clusters = False
        
        if len(possibilities) == 0:  # No valid clusters found
            print("No valid clusterings found for current mutations. Please check your input mutations and make sure \
                   that the multi-mutants are not too far apart.")
            sys.exit()
        return possibilities
    
    def kmeans_clustering(self, num_clusters: int, n_init=10, random_state=42):
        mean_idxs_dna = self.mutation_instance.mean_idxs_dna()  # Extract the first index of each mutation in idx_test
        mutation_arr = np.array(mean_idxs_dna, dtype=float).reshape(-1, 1)
        kmeans = KMeans(n_clusters=num_clusters, random_state=random_state, n_init=n_init)  # Initialize KMeans with the number of clusters
        cluster_labels = kmeans.fit_predict(mutation_arr)
        indices = []  # Add the remaining indices of the pairs to the assigned cluster
        for mutation in self.mutation_instance.mutations:
            if len(mutation.idx_dna) > 1: 
                pair = mutation.idx_dna
                mean_pair_value = np.mean(pair)
                idx = list(mutation_arr).index(mean_pair_value)
                indices.append(idx)
                cluster_label = cluster_labels[idx]
                cluster_labels = np.concatenate((cluster_labels, np.full(len(pair), cluster_label)))
                mean_idxs_dna.extend(pair)
        to_remove = sorted(set(indices), reverse=True)  # Remove the mean values from the list of indices and cluster labels
        for i in to_remove:
            del mean_idxs_dna[i]
            cluster_labels = np.delete(cluster_labels, i)
        return cluster_labels, mean_idxs_dna
    
    def choose_clusters(self, clusters: dict) -> dict:
        if self.cost_optimization:
            self.print_line("Optimizing based on price per bp ...")
            lowest_cost, best_clustering = min((self.calculate_cost(value), value) for value in clusters.values())
            self.print_line(f"Lowest estimated cost: €{lowest_cost} (given price per bp of €{self.bp_price})")
            return best_clustering
        elif self.amount_optimization:
            self.print_line("Optimizing based on number of eBlocks ...")
            fewest_blocks, best_clustering = min((len(value), value) for value in clusters.values())
            self.print_line(f"Fewest blocks: {fewest_blocks} with cluster {best_clustering}")
            return best_clustering

    def calculate_cost(self, clusters: dict) -> float:
        total_cost = 0
        for _, value in clusters.items():
            min_val = min(value)
            max_val = max(value)
            len_gene_block = (max_val - min_val) + 2 * self.min_overlap  # on both size of the gene block there should be an overhang for IVA primer design
            cost = len_gene_block * self.bp_price * len(value)
            total_cost += cost
        return round(total_cost, 2)
    
    def make_bins(self, clusters: dict):
        bins = []
        for _, value in clusters.items():
            bins.append(int(min(value) - self.min_overlap))
            bins.append(int(max(value) + self.min_overlap))
        return bins
    
    def make_wt_eblocks(self, bins: list) -> list:
        """
        Create wild-type eBlocks
        """
        gene_blocks = []
        for num in range(0, len(bins), 2):
            eblock = Eblock(
                start_index=bins[num],
                end_index=bins[num+1],
                sequence=Plasmid.slice_circular_sequence(self.sequence_instance.vector.seq, bins[num], bins[num+1]))
            gene_blocks.append(eblock)
        gene_blocks = self.renumber_eblock(gene_blocks)
        return gene_blocks
        
    def eblock_index(self, eblock: Eblock, idx_mutation: int) -> int:
        """
        Find the index of a mutation in an eblock.
        """
        if self.sequence_instance.gene_start_idx > eblock.start_index:
            mutation_index_in_eblock = self.sequence_instance.gene_start_idx - eblock.start_index + idx_mutation
        else:
            mutation_index_in_eblock = (idx_mutation - eblock.start_index) + self.sequence_instance.gene_start_idx
        return mutation_index_in_eblock
        
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
                print(f"eBlock is too long, length is {length_eblock}, maximum length is {self.max_eblock_length}")
                sys.exit()
            else:
                print(f"Codon insert is too short for mutation eBlock, length is {length_eblock}, minimum length is {self.min_eblock_length}")
                sys.exit()
        
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

                selected_eblock.mutation_start_index = self.eblock_index(selected_eblock, mut_i)
                selected_eblock.mutant_codon = self.select_mut_codon(mutation.mutation[num_i][-1])  # Find most occuring mutant codon based on codon usage for species
                
                self.check_wt_codon(selected_eblock, mutation.mutation[num_i][0])  # Check if WT codon at index is same residue as mutation
                selected_eblock.sequence = self.mutate_eblock(mutation, selected_eblock)

            eblock = selected_eblock
        results[mutation] = eblock
        return results
    
    def output_to_snapgene(self):
            
            snapgene_instance = SnapGene(output_dir=self.output_dir, sequence_instance=self.sequence_instance)
            snapgene_instance.make_dir()  # Make clones-dir
            original_dir = self.output_dir
            self.set_output_dir(snapgene_instance.output_dir)
            self.sequence_instance.output_dir = snapgene_instance.output_dir

            # For test purposes, save wt eblocks to SnapGene
            # snapgene_dict = {}
            # for i in self.wt_eblocks:
            #     snapgene_dict[i.name] = [Plasmid.circular_index(i.start_index, len(self.sequence_instance.vector.seq)), 
            #                              Plasmid.circular_index(i.end_index, len(self.sequence_instance.vector.seq)), self.eblock_colors[i.block_number]]
            #     snapgene_instance.eblocks_to_gff3(eblocks=snapgene_dict, output_dir=original_dir, filename=f"eblocks.gff3")

            # Loop over all mutations and create mutated vector and features that can be read by snapgene
            for mut, eblock in self.eblocks.items():
                snapgene_dict = {}
                snapgene_dict[self.sequence_instance.seqid] = [self.sequence_instance.gene_start_idx, self.sequence_instance.gene_end_idx, self.sequence_instance.color]
                snapgene_dict[eblock.name] = [eblock.start_index, eblock.end_index, self.eblock_colors[eblock.block_number]]
                
                if mut.is_singlemutation:
                    filename = mut.mutation[0]
                    snapgene_dict[mut.mutation[0]] = [self.sequence_instance.gene_start_idx -2 + eblock.mutation_start_index,
                                                      self.sequence_instance.gene_start_idx + eblock.mutation_start_index,
                                                      self.mutation_instance.colors[mut.type]]

                elif mut.is_insert:
                    # TODO Remove Insert and Deletion, because in mutated vector they will be different
                    filename = mut.mutation
                    snapgene_dict[mut.mutation] = [self.sequence_instance.gene_start_idx -2 + eblock.mutation_start_index,
                                                      self.sequence_instance.gene_start_idx + eblock.mutation_start_index,
                                                      self.mutation_instance.colors[mut.type]]

                elif mut.is_deletion:
                    filename = mut.mutation 
                    snapgene_dict[mut.mutation] = [self.sequence_instance.gene_start_idx -2 + eblock.mutation_start_index,
                                                      self.sequence_instance.gene_start_idx + eblock.mutation_start_index + mut.length_deletion,
                                                      self.mutation_instance.colors[mut.type]]

                elif mut.is_multiplemutation:
                    filename = '-'.join(mut.mutation)
                    for i, _ in enumerate(mut.idx_dna):
                        snapgene_dict[mut.mutation[i]] = [self.sequence_instance.gene_start_idx -2 + mut.idx_dna[i],
                                                      self.sequence_instance.gene_start_idx + mut.idx_dna[i],
                                                      self.mutation_instance.colors[mut.type]]
                self.make_dir(dirname=filename)
                snapgene_instance.eblocks_to_gff3(eblocks=snapgene_dict, output_dir=os.path.join(snapgene_instance.output_dir, filename), filename=f"{filename}.gff3")
                mutated_vector = self.sequence_instance.mutate_vector(eblock.start_index, eblock.end_index, eblock.sequence)
                self.sequence_instance.save_vector(vector=mutated_vector, output_dir=os.path.join(snapgene_instance.output_dir, filename), filename=f"{filename}.dna")
            
            self.output_dir = original_dir
            self.sequence_instance.output_dir = original_dir
        
    def count_mutations_per_eblock(self) -> dict:
        """
        Count the numer of mutations in each eBlock
        """
        counts = {}
        for i in self.wt_eblocks:
            counts[i.name] = 0
        for _, val in self.eblocks.items():
            counts[val.name] += 1
        return counts
    
    def check_clusters(self, clusters, cluster_sizes):
        """
        Find in which cluster insertions and deletions are located and check if the size of the cluster
        is within bounds when the mutation is added
        """
        for mutation in self.mutation_instance.mutations:
            if mutation.is_insert:
                for key, value in clusters.items():
                    if mutation.idx_dna in value:
                        cluster_sizes[key] += mutation.length_insert
            elif mutation.is_deletion:
                for key, value in clusters.items():
                    if mutation.idx_dna in value:
                        cluster_sizes[key] -= mutation.length_deletion
        return cluster_sizes
    
    def set_output_dir(self, output_dir: str):
        self.output_dir = output_dir

    def make_dir(self, dirname: str):
        try:
            os.makedirs(os.path.join(self.output_dir, f"{dirname}"))
        except FileExistsError:
            pass

    def validate_optimization_parameters(self):
        """
        Check if the optimization parameters are set correctly.
        """
        if not self.cost_optimization and not self.amount_optimization:
            print("Please set either cost_optimization or amount_optimization to True.")
            sys.exit()
        if self.cost_optimization and self.amount_optimization:
            print("Please set either cost_optimization or amount_optimization to True, not both.")
            sys.exit()
    
    def print_line(self, txt):
        if self.verbose:
            print(txt)

    @staticmethod
    def renumber_eblock(eblocks: list):
        sorted_list = sorted(eblocks, key=lambda x: x.start_index)
        for i, obj in enumerate(sorted_list, start=1):
            obj.name = f"eBlock-{i}"
            obj.block_number = i
        return sorted_list

    @staticmethod
    def cluster_size(clusters, labels, indexes):
        for i, j in zip(labels, indexes):
            if i not in clusters:
                clusters[i] = []
            clusters[i].append(j)
        return clusters
                
    @staticmethod
    def generate_eblock_colors() -> dict:
        """
        Create dictionary with colors for plotting eBlocks using the tab10 color scheme.
        """
        tab10_colors = ['#1f77b4','#ff7f0e','#2ca02c', '#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf',
                        '#aec7e8','#ffbb78','#98df8a','#ff9896','#c5b0d5','#c49c94','#f7b6d2','#c7c7c7','#dbdb8d','#9edae5',
                        '#393b79','#ff7f0e','#2ca02c','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf']
        return {i: tab10_colors[i] for i in range(len(tab10_colors))}