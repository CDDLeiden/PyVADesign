#!/usr/bin/env python3

import os
import sys
import copy
import math
import random
import warnings
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from datetime import datetime
import biotite.sequence as seq
from collections import Counter
from Bio.SeqRecord import SeqRecord
from sklearn_extra.cluster import KMedoids
from scipy.spatial.distance import pdist, squareform
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

from .utils import CodonUsage
from .mutation import Mutation
from .sequence import Vector, Gene

warnings.filterwarnings("ignore", category=UserWarning, module="sklearn_extra.cluster._k_medoids")  # Ignore warnings from sklearn_extra.cluster._k_medoids



class DNABlock:
    """
    A class to represent a single dsDNA fragment

    Attributes
    ----------
    name : str
        The name of the dsDNA fragment
    block_number : int
        The number of the dsDNA fragment
    sequence : str
        The DNA sequence of the dsDNA fragment
    start_index : int
        The start index of the dsDNA fragment in the vector
    end_index : int
        The end index of the dsDNA fragment in the vector
    bin_start : int
        The start index of the dsDNA fragment in the gene
    bin_end : int
        The end index of the dsDNA fragment in the gene
    mutation_start_index : int
        The start index of the mutation in the dsDNA fragment
    mutation_end_index : int
        The end index of the mutation in the dsDNA fragment
    wt_codon : str
        The wild-type codon at the mutation index
    mutant_codon : str
        The mutant codon at the mutation index
    insert : str
        The insert sequence for the dsDNA fragment
    """
    def __init__(self,
                 name: str = None,
                 block_number: int = None,
                 sequence: str = None,
                 start_index: int = None,
                 end_index: int = None,
                 bin_start: int = None,
                 bin_end: int = None,
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
        self.bin_start = bin_start
        self.bin_end = bin_end
        self.mutation_start_index = mutation_start_index
        self.mutation_end_index = mutation_end_index
        self.wt_codon = wt_codon
        self.mutant_codon = mutant_codon
        self.insert = insert
        self.silent_mutations = []



class DNABlockDesign:
    """
    A class to design dsDNA fragments

    Attributes
    ----------
    mutation_instance : Mutation
        The mutation instance containing all mutations of interest
    vector_instance : Vector
        The vector instance containing the gene of interest and the vector sequence
    gene_instance : Gene
        The gene instance containing the gene sequence
    output_dir : str
        The output directory to save the results
    settings_file : str
        The settings file to set the optimization parameters
    cost_optimization : bool
        Optimize the cost of the dsDNA fragments
    amount_optimization : bool
        Optimize the amount of the dsDNA fragments
    DNABlock_colors : dict
        A dictionary with colors for plotting the dsDNA fragments
    clone_files : bool
        Create clone files for each dsDNA fragment
    verbose : bool
        Print progress messages
    codon_usage : str
        The codon usage table to use for codon optimization
    bp_price : float
        The price per basepair
    max_DNABlock_length : int
        The maximum length of a dsDNA fragment
    min_DNABlock_length : int
        The minimum length of a dsDNA fragment
    min_overlap : int
        The minimum overlap between vector and dsDNA fragment in basepairs
    min_order : int
        The minimum amount of dsDNA fragments to order
    silent_mutations : bool
        Add silent mutations to the dsDNA fragments to simplify sequencing
    cost : float
        The cost of all dsDNA fragments combined
    silent_mutation_index : int
        The index of the silent mutation in the dsDNA fragment for sequencing purposes
    """
    def __init__(self, 
                 mutation_instance: Mutation,
                 vector_instance: Vector,
                 gene_instance: Gene,
                 output_dir: str = None,
                 settings_file: str = None,

                 cost_optimization: bool = False,
                 amount_optimization: bool = True,
                 clone_files: bool = True,
                 verbose: bool = True,
                 codon_usage: str = "U00096",  # Escherichia coli str. K-12 substr. MG1655

                 bp_price: float = 0.05,
                 max_DNABlock_length: int = 1500,
                 min_DNABlock_length: int = 300,
                 min_overlap: int = 25,
                 min_order: int = 24,
                 silent_mutation: bool = True,
                 silent_mutation_index: int = 30,
                ):
        
        self.mutation_instance = mutation_instance
        self.vector_instance = vector_instance
        self.gene_instance = gene_instance
        self.output_dir = output_dir
        self.settings_file = settings_file

        self.cost_optimization = cost_optimization
        self.amount_optimization = amount_optimization
        self.DNABlock_colors = self.generate_DNABlock_colors()
        self.clone_files = clone_files
        self.verbose = verbose
        self.codon_usage = codon_usage

        self.bp_price = bp_price
        self.max_DNABlock_length = max_DNABlock_length
        self.min_DNABlock_length = min_DNABlock_length
        self.min_overlap = min_overlap
        self.min_order = min_order
        self.cost = -1
        self.silent_mutations = silent_mutation
        self.silent_mutation_index = silent_mutation_index

        # Store WT and mutated gene blocks
        self.wt_DNABlocks: list = []
        self.DNABlocks: dict = {} # key: mutation, value: DNABlock

        self.most_abundant_codons: dict = {}  # Store most abundant codons for the selected genome
        self.codon_frequences = None
        self.validate_optimization_parameters()

    def run_design_DNABlocks(self):
        """Function to run the dsDNA fragment design process."""
        if self.settings_file:
            settings = self.parse_settings(self.settings_file)
            self.update_settings(settings)
            if self.verbose:
                self.display_settings()

        # Calculate the most abundant codon for each aa for the selected genome
        codonusage = CodonUsage(
            genome_id=self.codon_usage,
            output_dir=self.output_dir)
        # Check if codon-usage exists, otherwise calculate
        if not codonusage.codon_usage_exists():
            self.print_line(f"Calculating relative codon frequencies, based on the selected genome id {self.codon_usage}")
            self.most_abundant_codons, self.codon_frequences = codonusage.run()
        else:
            self.print_line(f"Loading relative codon frequencies from file")
            self.most_abundant_codons, self.codon_frequences = codonusage.load_relative_frequencies()

        self.print_line("Clustering mutations on gene of interest using K-medioids clustering")

        # Divide the target gene into fragment regions, based on the position of the mutations
        cluster_instance = Clustering(  
            mutation_instance=self.mutation_instance,
            vector_instance=self.vector_instance,
            gene_instance=self.gene_instance,
            max_DNABlock_length=self.max_DNABlock_length,
            min_overlap=self.min_overlap,
            min_DNABlock_length=self.min_DNABlock_length,
            cost_optimization=self.cost_optimization,
            amount_optimization=self.amount_optimization,
            bp_price=self.bp_price,
            verbose=self.verbose)
        optimal_clustering = cluster_instance.run_clustering()

        self.cost = cluster_instance.cost

        self.print_line("Starting dsDNA fragment design")  

        # Define the beginning and end of each dsDNA fragment, based on the clusters and include a minimum overlap that is needed for IVA cloning
        bins = self.make_bins(optimal_clustering)

        # Make fragment regions (WT DNA sequences sliced to the correct size, according to the bins) and renumber them starting from 1
        self.wt_DNABlocks = self.make_wt_DNABlocks(bins)
        
        # Randomly add more hex-colors if there are more fragment regions than colors in the default color scheme
        if len(self.wt_DNABlocks) > len(self.DNABlock_colors):  
            num_colors_to_make = len(self.wt_DNABlocks) - len(self.DNABlock_colors)
            new_colors = []
            for i in range(num_colors_to_make):
                new_colors.append(self.random_hex_color())
            max_num = max(self.DNABlock_colors.keys())
            for num, i in enumerate(new_colors, start=max_num+1):
                self.DNABlock_colors[num] = i

        # Create barcode for each fragment region by adding a silent mutation at the beginning and end of DNABlock to simplify sequencing
        if self.silent_mutations:
            for DNABlock in self.wt_DNABlocks:

                # Get starting and ending residue of the fragment region
                first_residue = math.ceil((DNABlock.start_index - self.vector_instance.gene_start_idx)/3) +1 # Add 1 to make sure the full codon is included
                final_residue = math.floor((DNABlock.end_index - self.vector_instance.gene_start_idx)/3) -1

                self.add_silent_mutations(first_residue, DNABlock, location='start')
                self.add_silent_mutations(final_residue, DNABlock, location='end')
           
        # Loop over all mutations in the mutations instance and create the dsDNA fragments, based on the fragments regions
        results = {}
        for mutation in self.mutation_instance.mutations:
            results = self.make_mutant_DNABlock(mutation, results)  # Create mutated DNABlock, based on mutation type
        
        # Check if all mutations could be mapped and remove mutations that could not be processed from mutation instance
        self.check_DNABlocks(results)
 
        sorted_dict = dict(sorted(results.items(), key=lambda x: (x[1].name, x[1].start_index)))  # Sort the dsDNA fragments based on the index of the first mutation in the dsDNA fragment and the fragment region
        self.DNABlocks = sorted_dict

        # Create a GFF3/gb file for each clone for easy visualization of dsDNA fragments in sequence manager
        if self.clone_files:
            self.make_clones()
            self.make_wt_clones()

        self.DNABlocks_to_csv()  # Save dsDNA fragments to a CSV file

        self.print_line("Completed dsDNA fragment design. \n")

        self.print_line(f"Clone files (.dna/.gb/.gff3) per input mutation stored in {os.path.join(self.output_dir, 'clones')}")
        self.print_line(f"CSV file with all dsDNA fragments stored in {os.path.join(self.output_dir, 'DNABlocks.csv')}")
        
        if len(self.mutation_instance.unprocessed_mutations) > 0:
            self.print_line(f"{' '.join([self.mutation_instance.unprocessed_mutations])} could not be processed and no dsDNA fragment was created for these mutations.")

    def choose_block(self, count_start: int, count_end: int, DNABlock_start: list, DNABlock_end: list):
        """Choose the dsDNA fragment region that is most suitable for the mutation."""
        if (count_start > 1) or (count_end > 1):
                for i in DNABlock_start:
                    for j in DNABlock_end:
                        if i.name == j.name:
                            DNABlock = i
        elif (count_start == 1) and (count_end == 1):
            DNABlock = DNABlock_start[0]
        return DNABlock

    def make_mutant_DNABlock(self, mutation: Mutation, results: dict) -> dict:
        """function to create mutated dsDNA fragment, based on the fragment region and the mutation."""
        DNABlock = None
        if mutation.is_singlemutation:
            DNABlock_start, count_start = self.DNABlocks_within_range(mutation.idx_dna[0])
            DNABlock_end, count_end = self.DNABlocks_within_range(mutation.idx_dna[0] + 3)
            DNABlock = self.choose_block(count_start, count_end, DNABlock_start, DNABlock_end)

            DNABlock.mutation_start_index = self.DNABlock_index(DNABlock, mutation.idx_dna[0])
            self.check_wt_codon(DNABlock, mutation.mutation[0][0])  # Check if WT codon at index is same residue as mutation
            DNABlock.mutant_codon = self.select_mut_codon(mutation.mutation[0][-1])
            DNABlock.sequence = self.mutate_DNABlock(mutation, DNABlock)

        elif mutation.is_insert:
            DNABlock_start, count_start = self.DNABlocks_within_range(mutation.idx_dna[0])
            DNABlock_end, count_end = self.DNABlocks_within_range(mutation.idx_dna[0] + mutation.length_insert)
            DNABlock = self.choose_block(count_start, count_end, DNABlock_start, DNABlock_end)
       
            DNABlock.mutation_start_index = self.DNABlock_index(DNABlock, mutation.idx_dna[0])
            self.check_wt_codon(DNABlock, mutation.mutation[0][0])
            DNABlock.insert = self.design_insert(mutation.insert)
            DNABlock.sequence = self.mutate_DNABlock(mutation, DNABlock)
            self.check_DNABlock_length(DNABlock.sequence)  # Check if DNABlock is too long including the insert

        elif mutation.is_deletion:
            DNABlock_start, count_start = self.DNABlocks_within_range(mutation.idx_dna_deletion_begin)
            DNABlock_end, count_end = self.DNABlocks_within_range(mutation.idx_dna_deletion_end)
            DNABlock = self.choose_block(count_start, count_end, DNABlock_start, DNABlock_end)

            DNABlock.mutation_start_index = self.DNABlock_index(DNABlock, mutation.idx_dna_deletion_begin)
            idx_end = self.DNABlock_index(DNABlock, mutation.idx_dna_deletion_end)
            self.check_wt_codon(DNABlock, mutation.mutation[0][0])  # Check if WT codon at index is same residue as mutation
            DNABlock.sequence = self.mutate_DNABlock(mutation, DNABlock, idx_end)
            self.check_DNABlock_length(DNABlock.sequence)  # Check if DNABlock is too short
     
        elif mutation.is_multiplemutation:
            selected_DNABlock = None
            all_DNABlocks = []
            for mut_i in mutation.idx_dna:
                DNABlocks, _ = self.DNABlocks_within_range(mut_i)
                all_DNABlocks.append([i.name for i in DNABlocks])
            counter = Counter()

            for lst in all_DNABlocks:
                counter.update(set(lst))

            common_DNABlock = [item for item, count in counter.items() if count == len(all_DNABlocks)]
            if len(common_DNABlock) == 1:
                selected_DNABlock = [e for e in DNABlocks if e.name == common_DNABlock[0]][0]
            elif len(common_DNABlock) > 1:

                # Check which DNABlock falls better within the range of the mutations (check not too close to beginning or end)
                possible_DNABlocks = [e for e in DNABlocks if e.name in common_DNABlock]
                for DNABlock in possible_DNABlocks:
                    for mut_i in mutation.idx_dna:
                        DNABlock.mutation_start_index = self.DNABlock_index(DNABlock, mut_i)
                        if not (DNABlock.mutation_start_index < self.min_overlap) or (DNABlock.mutation_start_index > (len(DNABlock.sequence) - self.min_overlap)):
                            selected_DNABlock = DNABlock           
            else:
                print(f"Multiple DNABlocks for multiple mutations, not all mutations in the same DNABlock. Skip mutation {mutation.name} {mutation.idx_dna}.")
                return results
       
            for num_i, mut_i in enumerate(mutation.idx_dna):
                selected_DNABlock.mutation_start_index = self.DNABlock_index(selected_DNABlock, mut_i)
                selected_DNABlock.mutant_codon = self.select_mut_codon(mutation.mutation[num_i][-1])  # Find most occuring mutant codon based on codon usage for species
                self.check_wt_codon(selected_DNABlock, mutation.mutation[num_i][0])  # Check if WT codon at index is same residue as mutation
                selected_DNABlock.sequence = self.mutate_DNABlock(mutation, selected_DNABlock)

            DNABlock = selected_DNABlock
        results[mutation] = DNABlock
        return results
    
    def make_wt_clones(self):
        """Create clone file for the fragment regions of the WT gene."""
        original_dir = self.output_dir
        self.set_output_dir(os.path.join(self.output_dir, 'clones'))
        filename = 'wild-type'
        results = {}
        for DNABlock in self.wt_DNABlocks:
            results[DNABlock.name] = [DNABlock.start_index, DNABlock.end_index, self.DNABlock_colors[DNABlock.block_number]]
            results[self.gene_instance.seqid] = [self.vector_instance.gene_start_idx, self.vector_instance.gene_end_idx, self.vector_instance.color]
            
            if len(DNABlock.silent_mutations) > 0:
                for silmut in DNABlock.silent_mutations:
                    start = self.vector_instance.gene_start_idx -3 + (int(silmut[1:-1])*3)
                    end = self.vector_instance.gene_start_idx + (int(silmut[1:-1])*3)
                    results[silmut] = [start, end, self.mutation_instance.colors['Silent']]

        self.make_dir(dirname=filename)
        self.check_directory(os.path.join(self.output_dir, filename))
        self.vector_instance.save_vector(vector=str(self.vector_instance.vector.seq), output_dir=os.path.join(self.output_dir, filename), filename=f"{filename}.dna")
        self.DNABlocks_to_gff3(DNABlocks=results, output_dir=os.path.join(self.output_dir, filename), filename=f"{filename}.gff3")
        self.DNABlocks_to_genbank(mutvector=str(self.vector_instance.vector.seq), DNABlocks=results, output_dir=os.path.join(self.output_dir, filename), filename=f"{filename}.gb")
        self.output_dir = original_dir
        
    def make_clones(self):
        """Create separate clone file for each dsDNA fragment.""" 
        self.make_dir(dirname='clones')  # Make clones-dir
        original_dir = self.output_dir
        self.set_output_dir(os.path.join(self.output_dir, 'clones'))

        # Loop over all mutations and create mutated vector and features that can be read by snapgene
        for mut, DNABlock in self.DNABlocks.items():
            
            results = {}
            if not (mut.is_deletion) and not (mut.is_insert):
                results[DNABlock.name] = [DNABlock.start_index, DNABlock.end_index, self.DNABlock_colors[DNABlock.block_number]]
                results[self.gene_instance.seqid] = [self.vector_instance.gene_start_idx, self.vector_instance.gene_end_idx, self.vector_instance.color]

            filename = mut.name
            
            # Add silent mutation information in clones
            for silmut in DNABlock.silent_mutations:
                start = self.vector_instance.gene_start_idx -3 + (int(silmut[1:-1])*3)
                end = self.vector_instance.gene_start_idx + (int(silmut[1:-1])*3)
                results[silmut] = [start, end, self.mutation_instance.colors['Silent']]
            
            if mut.is_singlemutation:
                start = self.vector_instance.gene_start_idx -3 + mut.idx_dna[0]
                end = self.vector_instance.gene_start_idx + mut.idx_dna[0]
                results[mut.name] = [start, end, self.mutation_instance.colors[mut.type]]

            elif mut.is_insert:
                start = self.vector_instance.gene_start_idx -3 + mut.idx_dna[0]
                end = self.vector_instance.gene_start_idx + mut.idx_dna[0] + mut.length_insert
                results[mut.name] = [start, end, self.mutation_instance.colors[mut.type]]
                results[DNABlock.name] = [DNABlock.start_index, DNABlock.end_index + mut.length_insert, self.DNABlock_colors[DNABlock.block_number]]
                results[self.gene_instance.seqid] = [self.vector_instance.gene_start_idx, self.vector_instance.gene_end_idx + mut.length_insert, self.vector_instance.color]
        
            elif mut.is_deletion:
                start = self.vector_instance.gene_start_idx -6 + mut.idx_dna_deletion_begin
                end = self.vector_instance.gene_start_idx -3 + mut.idx_dna_deletion_begin
                results[mut.name] = [start, end, self.mutation_instance.colors[mut.type]]
                results[self.gene_instance.seqid] = [self.vector_instance.gene_start_idx, self.vector_instance.gene_end_idx - mut.length_deletion, self.vector_instance.color]
                if DNABlock.start_index < DNABlock.end_index:
                    results[DNABlock.name] = [DNABlock.start_index, DNABlock.end_index - mut.length_deletion, self.DNABlock_colors[DNABlock.block_number]]
                else:
                    restoend = len(self.vector_instance.vector.seq) - DNABlock.start_index
                    newlength = len(self.vector_instance.vector.seq) - mut.length_deletion
                    newstart = newlength - restoend
                    results[DNABlock.name] = [newstart, DNABlock.end_index - mut.length_deletion, self.DNABlock_colors[DNABlock.block_number]]

            elif mut.is_multiplemutation:
                for i, _ in enumerate(mut.idx_dna):
                    start = self.vector_instance.gene_start_idx -3 + mut.idx_dna[i]
                    end = self.vector_instance.gene_start_idx + mut.idx_dna[i]
                    results[mut.mutation[i]] = [start, end, self.mutation_instance.colors[mut.type]]
                    
            self.make_dir(dirname=filename)
            self.check_directory(os.path.join(self.output_dir, filename))
            mutated_vector = self.vector_instance.mutate_vector(DNABlock.start_index, DNABlock.end_index, DNABlock.sequence, mutation=mut)
            self.vector_instance.save_vector(vector=mutated_vector, output_dir=os.path.join(self.output_dir, filename), filename=f"{filename}.dna")
            self.DNABlocks_to_gff3(DNABlocks=results, output_dir=os.path.join(self.output_dir, filename), filename=f"{filename}.gff3")
            self.DNABlocks_to_genbank(mutvector=mutated_vector, DNABlocks=results, output_dir=os.path.join(self.output_dir, filename), filename=f"{filename}.gb")
            
        self.output_dir = original_dir

    def add_silent_mutations(self, residue_index: int, DNABlock: DNABlock, location):
        """Add silent mutations to the WT DNABlocks to simplify sequencing."""
        n = self.silent_mutation_index

        while True:

            if location == 'start':  # Silent mutation beginning of DNABlock
                idx_start = residue_index + n
            elif location == 'end': # Silent mutation end of DNABlock
                idx_start = residue_index - n
            
            codon_start = self.gene_instance.residues[idx_start]
                
            if not str(codon_start[0]) in ['W', 'M']: # Only has one option so cannot make silent mutation
                all_codons = self.find_codon_for_residue(codon_start[0].upper())
                options = [c.upper() for c in all_codons if c != codon_start[1].upper()]
                options_freq = []
                if len(options) > 1:
                    codon_freqs = self.codon_frequences[codon_start[0]]
                    cs = [i[0] for i in codon_freqs]
                    fs = [i[1] for i in codon_freqs]
                    for i in options:
                        freq = fs[cs.index(i)]
                        options_freq.append(freq)
                    # Select most abundant codon
                    newcodon = options[np.argmax(options_freq)]
                elif len(options) == 1:
                    newcodon = options[0]
                else:
                    raise Exception(f"No codon found for residue {codon_start[0]}")
                DNABlock.mutation_start_index = self.DNABlock_index(DNABlock, idx_start*3)
                mutated_DNABlock = DNABlock.sequence[:DNABlock.mutation_start_index -3] + newcodon + DNABlock.sequence[DNABlock.mutation_start_index:]
                DNABlock.sequence = mutated_DNABlock
                DNABlock.silent_mutations.append(str(codon_start[0]) + str(idx_start) + str(codon_start[0]))
                DNABlock.mutation_start_index = None

                silent_mutation = Mutation(name=str(codon_start[0]) + str(idx_start) + str(codon_start[0]), 
                                            type='Silent',
                                            idx_dna=[idx_start * 3])
                if not silent_mutation.name in [i.name for i in self.mutation_instance.silent_mutations]:
                    self.mutation_instance.silent_mutations.append(silent_mutation)
                break
            else:
                n += 1

    def DNABlocks_to_genbank(self, mutvector, DNABlocks: dict, output_dir, type='gene', filename='DNABlocks.gb', header=True, max_filename_length=16):
        """This function saves a vector to a GenBank (gb) file."""
        sequence = Seq(mutvector)
        record = SeqRecord(sequence, id=self.gene_instance.seqid, name=self.gene_instance.seqid, description="")
        record.annotations["molecule_type"] = "DNA"
        record.annotations["organism"] = self.vector_instance.organism
        record.annotations["date"] = datetime.today().strftime('%d-%b-%Y').upper()
        record.name = record.name + "_" + filename
        # Limit filename length characters
        if len(record.name) > max_filename_length:
            record.name = record.name[:max_filename_length]
        features = []  # Add DNABlock and mutations as features
        for k, v in DNABlocks.items():
            if v[0] > v[1]: # Start index is larger than end index
                joint_location = CompoundLocation([FeatureLocation(v[0], len(mutvector)), FeatureLocation(0, v[1])])
                joint_feature = SeqFeature(joint_location, type="gene", qualifiers={"gene": k, "color": v[2]})
                features.append(joint_feature)
            else:
                feature = SeqFeature(FeatureLocation(v[0], v[1]), type=type, qualifiers={"gene": k, "color": v[2]})
                features.append(feature)
        record.features.extend(features)     
        outpath = os.path.join(output_dir, filename)
        SeqIO.write(record, outpath, "genbank")

    def make_bins(self, clusters: dict):
        """Create bins based on the clusters and add overlap."""
        bins = []
        for _, value in clusters.items():
            bins.append(int(min(value) - self.min_overlap))  # add overlap
            bins.append(int(max(value) + self.min_overlap))
        return bins
    
    def make_wt_DNABlocks(self, bins: list) -> list:
        """Create fragment regions based bins."""
        gene_blocks = []
        for num in range(0, len(bins), 2):
            start_index = self.vector_instance.gene_start_idx + bins[num]
            end_index = self.vector_instance.gene_start_idx + bins[num+1]
            wt_DNABlock = DNABlock(
                start_index=self.vector_instance.circular_index(start_index, len(self.vector_instance.vector.seq)),  # start index in the vector
                end_index=self.vector_instance.circular_index(end_index, len(self.vector_instance.vector.seq)), 
                bin_start=bins[num],  # start index in the gene
                bin_end=bins[num+1],
                sequence=Vector.slice_circular_sequence(self.vector_instance.vector.seq, start_index, end_index))
            gene_blocks.append(wt_DNABlock)
        gene_blocks = self.renumber_DNABlock(gene_blocks)
        return gene_blocks
        
    def DNABlock_index(self, DNABlock: DNABlock, idx_mutation: int) -> int:
        """Find the index of a mutation in an fragment region."""
        if (DNABlock.start_index > DNABlock.end_index) and (self.vector_instance.gene_start_idx < DNABlock.start_index):
            mutation_idx_in_gene  = self.vector_instance.gene_start_idx + idx_mutation
            residues_to_end = len(self.vector_instance.vector.seq) - DNABlock.start_index
            mutation_index_in_DNABlock = residues_to_end + mutation_idx_in_gene
        elif (DNABlock.start_index < DNABlock.end_index) and (self.vector_instance.gene_start_idx < DNABlock.start_index):
            mutation_index_in_DNABlock = (idx_mutation - DNABlock.start_index) + self.vector_instance.gene_start_idx
        elif (DNABlock.start_index < DNABlock.end_index) and (self.vector_instance.gene_start_idx > DNABlock.start_index):
            mutation_index_in_DNABlock = self.vector_instance.gene_start_idx - DNABlock.start_index + idx_mutation
        else:
            raise Exception("Error in mutation index calculation")
        return mutation_index_in_DNABlock
        
    def DNABlocks_within_range(self, mutation_idx: int):
        """For a given mutation index, find the number of overlapping dsDNA fragments."""
        within_range = []
        for DNABlock in self.wt_DNABlocks:
            if DNABlock.start_index < DNABlock.end_index:
                if DNABlock.bin_start < mutation_idx < DNABlock.bin_end:
                    within_range.append(DNABlock)
            else:
                if (DNABlock.bin_start < mutation_idx) or (mutation_idx < DNABlock.bin_end):
                    within_range.append(DNABlock)
        count = len(within_range)
        return copy.deepcopy(within_range), count

    def select_mut_codon(self, res: str):
        """Select the most abundant codon for a given residue, based on the relative frequencies of the codons in the selected genome."""
        return self.most_abundant_codons[res][0]
    
    def check_wt_codon(self, DNABlock: DNABlock, mut: str):
        """Check whether the WT codon at the mutation index is the same as in the proposed mutation."""
        codon = DNABlock.sequence[DNABlock.mutation_start_index-3:DNABlock.mutation_start_index].upper()
        try:
            result = seq.CodonTable.default_table()[str(codon)]
        except:
            result = None
        if result != mut[0]:
            raise Exception(f"WT codon does not match residue {mut}, but is {result}, the codon is {codon}")

    def update_settings(self, settings):
        """Update class attributes based on settings from a settings file."""
        for key, value in settings.items():
            if key == 'cost_optimization':
                self.cost_optimization = value == "True"
            elif key == 'amount_optimization':
                self.amount_optimization = value == "True"
            elif key == 'bp_price':
                self.bp_price = float(value)
            elif key == 'max_DNABlock_length':
                self.max_DNABlock_length = int(value)
            elif key == 'min_DNABlock_length':
                self.min_DNABlock_length = int(value)
            elif key == 'min_overlap':
                self.min_overlap = int(value)
            elif key == 'min_order':
                self.min_order = int(value)
            elif key == 'codon_usage':
                self.codon_usage = str(value)
            elif key == 'output_dir':
                self.output_dir = str(value)
            elif key == 'verbose':
                self.verbose = value == "True"
            elif key == 'clone_files':
                self.clone_files = value == "True"
            elif key == 'silent_mutations':
                self.silent_mutations = value == "True"
            elif key == 'silent_mutation_index':
                self.silent_mutation_index = int(value)
            else:
                raise Exception(f"Key {key} not found in class.")
            
    def count_mutations_per_DNABlock(self) -> dict:
        """Count the numer of mutations in each DNABlock."""
        counts = {}
        for i in self.wt_DNABlocks:
            counts[i.name] = 0
        for _, val in self.DNABlocks.items():
            counts[val.name] += 1
        return counts
        
    def set_output_dir(self, output_dir: str):
        self.output_dir = output_dir

    def make_dir(self, dirname: str):
        dir_path = os.path.join(self.output_dir, f"{dirname}")
        os.makedirs(dir_path, exist_ok=True)

    def make_file(self, directory, filename, header=False):
        """Check if the file exists; if not, create it and optionally add a header."""
        file_path = os.path.join(directory, filename)
        if not os.path.exists(file_path):
            with open(file_path, 'w') as f:
                if header:
                    f.write("\n".join(self.gff3_header(length_sequence=self.vector_instance.length)))
                    f.write("\n")
            
    def mutate_DNABlock(self, mutation, DNABlock: DNABlock, end_idx=None):
        """Mutate gene block based on mutation type."""
        if (mutation.is_singlemutation) or (mutation.is_multiplemutation):
            mutated_DNABlock = DNABlock.sequence[:DNABlock.mutation_start_index -3] + DNABlock.mutant_codon + DNABlock.sequence[DNABlock.mutation_start_index:]
        elif mutation.is_insert:
            mutated_DNABlock = DNABlock.sequence[:DNABlock.mutation_start_index] + DNABlock.insert + DNABlock.sequence[DNABlock.mutation_start_index:]
        elif mutation.is_deletion:
            mutated_DNABlock = DNABlock.sequence[:DNABlock.mutation_start_index -3] + DNABlock.sequence[end_idx -3:]
        return mutated_DNABlock
        
    def check_DNABlock_length(self, DNABlock_seq: str) -> bool:
        """Check if the length of the gene block is within bounds."""
        length_DNABlock = len(DNABlock_seq)
        if not self.min_DNABlock_length <= length_DNABlock <= self.max_DNABlock_length:
            if length_DNABlock > self.max_DNABlock_length:
                raise Exception(f"dsDNA fragment is too long, length is {length_DNABlock}, maximum length is {self.max_DNABlock_length}. Try increasing or decreasing the maximum fragment length.")
            else:
                raise Exception(f"dsDNA fragment is too short, length is {length_DNABlock}, minimum length is {self.min_DNABlock_length}")

    def DNABlocks_to_csv(self, filename="DNABlocks.csv"):
        """Save DNABlocks to a CSV file."""
        file_path = os.path.join(self.output_dir, filename)
        with open(file_path, 'w') as f:
            f.write("DNABlock,Mutation,Start,End,Sequence\n")
            for key, value in self.DNABlocks.items():
                f.write(f"{value.name},{key.name},{value.start_index},{value.end_index},{value.sequence}\n")

    def DNABlocks_to_gff3(self, DNABlocks: dict, output_dir, type='gene', filename='DNABlocks.gff3', header=True):
        """This function converts the DNABlocks to features that can be read by SnapGene."""
        self.make_file(output_dir, filename, header=header)
        with open(os.path.join(output_dir, filename), 'a') as f:
            for k, v in DNABlocks.items():
                line = self.gff3_line(v[0], v[1], k, v[2], type)
                f.write('\t'.join(line) + '\n')

    def validate_optimization_parameters(self):
        """Check if the optimization parameters are set correctly."""
        if not self.cost_optimization and not self.amount_optimization:
            raise Exception("Please set either cost_optimization or amount_optimization to True.")
        if self.cost_optimization and self.amount_optimization:
            raise Exception("Please set either cost_optimization or amount_optimization to True, not both.")
    
    def check_directory(self, directory):
        """Check if a directory exists and is empty."""
        if not os.path.exists(directory):  # Check if the directory exists
            raise FileNotFoundError(f"Directory {directory} does not exist.")
        if os.listdir(directory):  # Check if the directory is empty
            if self.verbose:
                print(f"Warning: Directory {directory} is not empty. Files might get overwritten or appended to.")

    def check_DNABlocks(self, results: dict):
        """Check if all mutations could be mapped to an DNABlock and remove mutations that could not be processed from the mutation instance."""
        failed_mutations = []
        for num, mutation in enumerate(self.mutation_instance.mutations):
            if mutation not in results.keys():
                failed_mutations.append(num)
        if len(failed_mutations) > 0:
            failed_mutations = sorted(failed_mutations, reverse=True)  # Sort failed mutations from high to low index
            for idx in failed_mutations:
                self.mutation_instance.remove_index(idx)

    def parse_settings(self, file):
        settings = {}
        with open(file, 'r') as f:
            for line in f:
                if '=' in line:
                    key, value = line.strip().split('=')
                    settings[key] = value
        return settings
    
    def design_insert(self, aas):
        codon_insert = ''  # Sequence to insert in gene block
        for res in aas:
            codon = self.select_mut_codon(res)
            codon_insert += codon
        return codon_insert
    
    def display_settings(self):
        """Print all class attributes using the __dict__ attribute."""
        max_key_length = max(len(key) for key in self.__dict__.keys())
        do_not_print = ['mutation_instance', 'vector_instance', 'gene_instance', 'DNABlock_colors', 'wt_DNABlocks', 'DNABlocks', 'most_abundant_codons', 'codon_frequences']
        print(f"Settings for {self.__class__.__name__}:")
        for key, value in self.__dict__.items():
            if key in do_not_print:
                continue
            print(f"\t\t{key.ljust(max_key_length)}: {value}")
        print("\n")

    def print_line(self, txt):
        if self.verbose:
            print(txt)
                
    @staticmethod
    def renumber_DNABlock(DNABlocks: list):
        sorted_list = sorted(DNABlocks, key=lambda x: x.bin_start)
        for i, obj in enumerate(sorted_list, start=1):
            obj.name = f"DNABlock-{i}"
            obj.block_number = i
        return sorted_list
    
    @staticmethod
    def find_codon_for_residue(aa: str):
        codon_table = CodonTable.standard_dna_table
        return [codon for codon, residue in codon_table.forward_table.items() if residue == aa]
                
    @staticmethod
    def generate_DNABlock_colors() -> dict:
        """Create dictionary with colors for plotting fragment regions using the tab10 color scheme."""
        tab10_colors = ['#1f77b4','#ff7f0e','#2ca02c', '#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf',
                        '#aec7e8','#ffbb78','#98df8a','#ff9896','#c5b0d5','#c49c94','#f7b6d2','#c7c7c7','#dbdb8d','#9edae5',
                        '#393b79','#ff7f0e','#2ca02c','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf']
        return {i: tab10_colors[i] for i in range(len(tab10_colors))}
    
    @staticmethod
    def random_hex_color():
        return "#{:06x}".format(random.randint(0, 0xFFFFFF))
    
    @staticmethod
    def gff3_header(length_sequence, version="3.2.1", sequence_name="myseq"):
        result = [f"##gff-version {version}", f"##sequence-region {sequence_name} 1 {length_sequence}"]
        return result

    @staticmethod
    def gff3_line(begin_pos, end_pos, name, hex_color, type):
        line = ['myseq', '.', f"{type}", str(begin_pos), str(end_pos), '.', '.', '.', f"Name={name};color={hex_color}"]
        return line
    
    @staticmethod
    def gff3_colors():
        colors = {
            'mutation': '#FF0000',
            'combination': '#D8FF00',
            'insert': '#0017FF',
            'deletion': '#FF5900',
            'primer': '#FF00D4'}
        return colors
    


class Clustering:
    """
    This class contains functions to cluster input mutations based on their index in the gene sequence

    Attributes
    ----------
    mutation_instance : Mutation
        Instance of the Mutation class
    vector_instance : Vector
        Instance of the Vector class
    gene_instance : Gene
        Instance of the Gene class
    max_DNABlock_length : int
        Maximum length of a dsDNA fragment
    min_overlap : int
        Minimum overlap between dsDNA fragment and gene sequence
    min_DNABlock_length : int
        Minimum length of a dsDNA fragment
    cost_optimization : bool
        Optimize the cost of the dsDNA fragments
    amount_optimization : bool
        Optimize the amount of dsDNA fragments
    bp_price : float
        Price per base pair
    verbose : bool
        Print additional information
    """
    def __init__(self,
                 mutation_instance: Mutation,
                 vector_instance: Vector,
                 gene_instance: Gene,
                 max_DNABlock_length: int = None,
                 min_overlap: int = None,
                 min_DNABlock_length: int = None,
                 cost_optimization: bool = None,
                 amount_optimization: bool = None,
                 bp_price: float = None,
                 verbose: bool = None):
        
        self.mutation_instance = mutation_instance
        self.vector_instance = vector_instance
        self.gene_instance = gene_instance
        self.max_DNABlock_length = max_DNABlock_length
        self.min_overlap = min_overlap
        self.min_DNABlock_length = min_DNABlock_length
        self.cost_optimization = cost_optimization
        self.amount_optimization = amount_optimization
        self.bp_price = bp_price
        self.verbose = verbose
        self.cost = -1

        self.idxs_constraints = []  # Store indices of constraints
        self.X = []  # Store all indices of the mutations

    def run_clustering(self):
        self.X, self.idxs_constraints = self.mutation_instance.extract_indices()
        valid_clusters = self.find_possible_clusters()
        optimal_clustering = self.choose_cluster(valid_clusters)
        return optimal_clustering
    
    def find_possible_clusters(self, threshold_small_clusters=3):
        """This function finds all possible clusterings of the mutations based on the index of the mutations."""
        possibilities = {} # Store all possible ways of clustering the mutations
        n = 1
        valid_clusters = True
        
        while (valid_clusters) and (n <= len(self.mutation_instance.mutations)):
      
            cluster_labels, invalid_constraints = self.kmedoids_clustering(num_clusters=n)
            clusters = self.cluster_to_dict(cluster_labels, self.X)
            cluster_sizes = [max(v) - min(v) for v in clusters.values()]
       
            # Calculate the size, based on the largest/deletion/insertion per cluster
            cluster_sizes = self.calculate_max_min_cluster_sizes(clusters)

            # Fix constraints if there are any invalid constraints
            if (invalid_constraints > 0) and not any(min_size < (self.min_DNABlock_length - 2 * self.min_overlap) for min_size, _ in cluster_sizes.values()):
                
                clusters_copy = self.fix_constraints(clusters, cluster_labels)
                new_labels = self.cluster_to_labels(clusters_copy)
                invalid_constraints = self.count_invalid_constraints(self.idxs_constraints, new_labels)
            else:
                clusters_copy = copy.deepcopy(clusters)

            cluster_too_small = 0
            cluster_correct = 0
            cluster_too_big = False

            # Check if the cluster sizes are within bounds
            if any(max_size > (self.max_DNABlock_length - 2 * self.min_overlap) for _, max_size in cluster_sizes.values()):  # Cluster too big
                cluster_too_big = True

            elif any(min_size < (self.min_DNABlock_length - 2 * self.min_overlap) for min_size, _ in cluster_sizes.values()):  # Cluster too small
                exceeding_items = [min_size for min_size, _ in cluster_sizes.values() if min_size < (self.min_DNABlock_length - 2 * self.min_overlap)]  # Length cluster
                cluster_keys = [k for k, v in cluster_sizes.items() if v[0] in exceeding_items]  # cluster label

                min_length_to_add = [(self.min_DNABlock_length - size) - 2 * self.min_overlap for size in exceeding_items]
                highest_indexes = []
                for i in cluster_keys:
                    indices = [k for k, x in enumerate(cluster_labels) if x == i]
                    highest_indexes.append(max(self.X[indices]))

                # See if adding some length to the cluster will make it fit the requirements
                for i, j, k, l in zip(exceeding_items, min_length_to_add, cluster_keys, highest_indexes):
                    vals = clusters[k]
                    if l + j <= len(self.gene_instance.sequence):
                        clusters_copy[k].append(max(vals) + j)
                        cluster_too_small += 1
                        cluster_correct += 1
                    else:
                        clusters_copy[k].append(min(vals) - j)
                        cluster_too_small += 1
                        cluster_correct += 1
                cluster_correct += (len(clusters) - len(exceeding_items))
            else:
                cluster_correct += len(cluster_sizes.keys())

            if (cluster_too_small > threshold_small_clusters):  # Too many small clusters > stop searching
                valid_clusters = False
            
            if (cluster_correct == len(clusters)) and (invalid_constraints == 0):
                possibilities[f'cluster N={n}'] = clusters_copy
                n += 1
            elif cluster_too_big:
                n += 1
            else:
                n += 1

        if len(possibilities) == 0:  # No valid clusters found
            raise Exception("No valid clusterings found for current mutations. \
                            Check your input mutations and make sure that (1) the multi-mutants are not too far apart. \
                            and (2) insertions and deletions are not too large or too small and (3) the DNABlocks settings are not too strict.")
        return possibilities
    
    def kmedoids_clustering(self, num_clusters: int, dist_metric='euclidean', random_state=42):
        distance_matrix = squareform(pdist(self.X[:, np.newaxis], metric=dist_metric))  # Precompute the distance matrix

        for idx1, idx2 in self.idxs_constraints:  # Modify the distance matrix based on the constraints so that the distance between the constraints is halved
            distance_matrix[idx1, idx2] *= 0.5
            distance_matrix[idx2, idx1] *= 0.5

        kmedoids = KMedoids(n_clusters=num_clusters, metric='precomputed', random_state=random_state, init='k-medoids++')
        kmedoids.fit(distance_matrix)

        labels = kmedoids.labels_
        labels = labels.tolist()

        invalid_constraints = self.count_invalid_constraints(self.idxs_constraints, labels)
        return labels, invalid_constraints
    
    def choose_cluster(self, clusters: dict) -> dict:
        if self.cost_optimization == True:
            self.print_line("Optimizing based on price per bp")
            get_score = lambda c: self.calculate_cost(c)
        elif self.amount_optimization == True:
            self.print_line("Optimizing based on number of fragment regions")
            get_score = lambda c: len(c)
        else:
            raise Exception("Please set either cost_optimization or amount_optimization to True, but not both.")

        # Find the cluster with the minimum score
        best_clustering_k = min(clusters, key=lambda k: get_score(clusters[k]))
        self.cost = self.calculate_cost(clusters[best_clustering_k])
        return clusters[best_clustering_k]
    
    def calculate_max_min_cluster_sizes(self, clusters: dict) -> list:
        min_max_sizes = {}
        for k, v in clusters.items():
            size = max(v) - min(v)
            max_insert = 0
            max_deletion = 0
            
            # Get the mutation based on the index
            for mut_idx in v:
                for mut in self.mutation_instance.mutations:
                    if mut_idx in mut.idx_dna:
                        matching_mutation = mut
                        if matching_mutation.is_insert:  # Calculate size based on insert
                            if matching_mutation.length_insert > max_insert:
                                max_insert = matching_mutation.length_insert
                        elif matching_mutation.is_deletion:  # Calculate size based on deletion
                            if matching_mutation.length_deletion > max_deletion:
                                max_deletion = matching_mutation.length_deletion
            
            # Calculate the new size based on the largest/deletion/insertion
            if max_insert > 0:
                max_size = size + max_insert
            else:
                max_size = size
            
            if max_deletion > 0:
                min_size = size - max_deletion
            else:
                min_size = size

            min_max_sizes[k] = (min_size, max_size)
        return min_max_sizes
                     
    def fix_constraints(self, groups, labels):
        """Fix invalid constraints by adding values to clusters if constraints allow it."""
        groups_copy = copy.deepcopy(groups)
        invalid_constraints = self._get_invalid_constraints(groups, labels)  # Identify invalid constraints
        invalid_constraints = self._remove_empty_keys(invalid_constraints)  # Remove empty constraint keys
        values_to_adjust = self._find_values_to_adjust(invalid_constraints)
        filtered_values = self._filter_constraints(values_to_adjust, labels)
        clean_values = self._clean_values_to_add(filtered_values)
        groups_copy = self._apply_adjustments_to_groups(clean_values, groups_copy)  # Add values to clusters if constraints allow it
        return groups_copy

    def count_invalid_constraints(self, constraints_indices, labels) -> int:
        """
        Count the number of invalid constraints.

        A constraint is considered invalid if the data points it references are 
        not assigned to the same cluster
        """
        invalid = 0
        for con in constraints_indices:
            con_labels = [labels[i] for i in con]
            same = all(x == con_labels[0] for x in con_labels) # Check if all pairs are in the same cluster
            if not same:
                invalid += 1
        return invalid
    
    def calculate_cost(self, clusters: dict) -> float:
        total_cost = 0
        for _, value in clusters.items():
            min_val = min(value)
            max_val = max(value)
            len_gene_block = (max_val - min_val) + 2 * self.min_overlap
            cost = len_gene_block * self.bp_price * len(value)
            total_cost += cost
        return round(total_cost, 2)

    def print_line(self, txt: str):
        if self.verbose:
            print(txt)
    
    def _get_invalid_constraints(self, groups, labels):
        """Identify invalid constraints based on label mismatches."""
        invalid_constraints = {key: [] for key in groups.keys()}
        for con in self.idxs_constraints:
            con_labels = [labels[i] for i in con]
            if not all(x == con_labels[0] for x in con_labels):
                difference = self.X[con[0]] - self.X[con[1]]
                if difference > 0:
                    invalid_constraints[con_labels[0]].append(difference)
                    invalid_constraints[con_labels[1]].append(-1 * difference)
                else:
                    invalid_constraints[con_labels[0]].append(-1 * difference)
                    invalid_constraints[con_labels[1]].append(difference)
        return invalid_constraints

    def _remove_empty_keys(self, invalid_constraints):
        """Remove keys with empty lists from invalid constraints."""
        return {k: v for k, v in invalid_constraints.items() if len(v) > 0}

    def _find_values_to_adjust(self, invalid_constraints):
        """Find values to be adjusted based on invalid constraints."""
        values_to_adjust = {}
        for key1, list1 in invalid_constraints.items():
            for key2, list2 in invalid_constraints.items():
                if key1 != key2:
                    abs_list1 = [abs(i) for i in list1]
                    abs_list2 = [abs(i) for i in list2]
                    matches = [i for i in abs_list1 if i in abs_list2]
                    if matches:
                        highest = max(matches)
                        # Track key where the highest value was found
                        if highest in list1:
                            values_to_adjust[key1, key2] = [highest, key1]
                        else:
                            values_to_adjust[key1, key2] = [highest, key2]      
        return values_to_adjust

    def _filter_constraints(self, values_to_adjust, labels):
        """Filter the values to adjust based on existing constraints."""
        filtered_values = {}
        for k, v in values_to_adjust.items():
            for con in self.idxs_constraints:
                con_labels = [labels[i] for i in con]
                if k[0] in con_labels and k[1] in con_labels:
                    filtered_values[k] = v
        return filtered_values

    def _clean_values_to_add(self, filtered_values):
        """Clean up the values to adjust by verifying correct label associations."""
        clean_values = {}
        for k, v in filtered_values.items():
            if k[0] == v[1]:
                clean_values[k] = v[0]
        return clean_values

    def _apply_adjustments_to_groups(self, clean_values, groups_copy):
        """Apply adjustments to groups if constraints are met."""
        for k, v in clean_values.items():
            max_value = max(groups_copy[k[0]])
            if (max_value - min(groups_copy[k[0]]) + v <= self.max_DNABlock_length - 2 * self.min_overlap):  # Check that the group size doesn't exceed the allowed limit
                groups_copy[k[0]].append(max_value + v)
        return groups_copy
        
    @staticmethod
    def cluster_to_dict(labels, indexes):
        clusters = {}
        for i, j in zip(labels, indexes):
            if i not in clusters:
                clusters[i] = []
            clusters[i].append(j)
        return clusters
    
    @staticmethod
    def cluster_to_labels(clusters: dict) -> list:
        labels = []
        for k, v in clusters.items():
            labels.extend([k] * len(v))
        return labels