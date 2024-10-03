import os
import sys
import copy
import itertools
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
from datetime import datetime
import biotite.sequence as seq
from Bio.SeqRecord import SeqRecord
from sklearn_extra.cluster import KMedoids
from scipy.spatial.distance import pdist, squareform
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

from .utils import CodonUsage
from .mutation import Mutation
from .sequence import Vector, Gene



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



class EblockDesign:
    """
    This class contains functions to design eBlocks based on mutations in a gene sequence.
    """
    def __init__(self, 
                 mutation_instance: Mutation,
                 vector_instance: Vector,
                 gene_instance: Gene,
                 output_dir: str = None,
                 settings_file: str = None,

                 cost_optimization: bool = True,
                 amount_optimization: bool = False,
                 clone_files: bool = True,
                 verbose: bool = True,
                 codon_usage: str = "U00096",  # Escherichia coli str. K-12 substr. MG1655

                 bp_price: float = 0.05,
                 max_eblock_length: int = 1500,
                 min_eblock_length: int = 300,
                 min_overlap: int = 25,
                 min_order: int = 24,
                ):
        
        self.mutation_instance = mutation_instance
        self.vector_instance = vector_instance
        self.gene_instance = gene_instance
        self.output_dir = output_dir
        self.settings_file = settings_file

        self.cost_optimization = cost_optimization
        self.amount_optimization = amount_optimization
        self.eblock_colors = self.generate_eblock_colors()
        self.clone_files = clone_files
        self.verbose = verbose
        self.codon_usage = codon_usage

        self.bp_price = bp_price
        self.max_eblock_length = max_eblock_length
        self.min_eblock_length = min_eblock_length
        self.min_overlap = min_overlap
        self.min_order = min_order

        # Store WT and mutated gene blocks
        self.wt_eblocks: list = []
        self.eblocks: list = []

        self.most_abundant_codons: dict = {}  # Store most abundant codons for the selected genome
        self.validate_optimization_parameters()

    def run_design_eblocks(self):
        """
        This function runs the design process for eBlocks.
        """
        if self.settings_file:
            settings = self.parse_settings(self.settings_file)
            self.update_settings(settings)
            if self.verbose:
                self.display_settings()

        self.print_line(f"Calculating relative codon frequencies, based on the selected genome id {self.codon_usage} ...")
        codonusage = CodonUsage(
            genome_id=self.codon_usage,
            output_dir=self.output_dir)
        self.most_abundant_codons = codonusage.run()

        self.print_line("Starting eBlock design ...")  # Divide the target gene into clusters based on the position of the mutations
        cluster_instance = Clustering(  
            mutation_instance=self.mutation_instance,
            vector_instance=self.vector_instance,
            gene_instance=self.gene_instance,
            max_eblock_length=self.max_eblock_length,
            min_overlap=self.min_overlap,
            min_eblock_length=self.min_eblock_length,
            cost_optimization=self.cost_optimization,
            amount_optimization=self.amount_optimization,
            bp_price=self.bp_price,
            verbose=self.verbose)
        optimal_clustering = cluster_instance.run_clustering()
        print("optimal_clustering", optimal_clustering)

        # Define the beginning and end of each gene block, based on the clusters and include the minimum overlap
        bins = self.make_bins(optimal_clustering)

        # Make gene blocks (WT DNA sequences sliced to the correct size, according to the bins) and renumber them starting from 1
        self.wt_eblocks = self.make_wt_eblocks(bins)
        print("WT EBLOCKS!!!!!")
        for i in self.wt_eblocks:
            print(i.name, i.start_index, i.end_index, i.bin_start, i.bin_end, i.sequence)

        # Loop over all mutations and create the eBlocks, based on the WT eBlocks
        results = {}
        for mutation in self.mutation_instance.mutations:
            print("MUTATION", mutation.name, mutation.idx_dna)
            results = self.make_mutant_eblock(mutation, results)  # Create mutated eBlock, based on mutation type
        # Check if all mutations could be mapped and remove mutations that could not be processed from mutation instance
        self.check_eblocks(results)
        for k, v in results.items():
            print(k.name, v.name, v.start_index, v.end_index, v.sequence)

        # print("UNPROCESSED MUTATIONS")
        # print(self.mutation_instance.unprocessed_mutations)
        # sys.exit()

        sorted_dict = dict(sorted(results.items(), key=lambda x: (x[1].name, x[1].start_index)))  # Sort the eblocks based on the index of the first mutation in the eblock and the number of the eblock
        self.eblocks = sorted_dict

        # Create a GFF3/gb file for each clone for easy visualization of eBlocks in sequence editor tools
        if self.clone_files:
            self.make_clones()

        self.eblocks_to_csv()  # Save eBlocks to a CSV file

        self.print_line("Completed eBlock design.")
                                                     
    def make_bins(self, clusters: dict):
        bins = []
        for _, value in clusters.items():
            bins.append(int(min(value) - self.min_overlap))  # add overlap
            bins.append(int(max(value) + self.min_overlap))
        return bins
    
    def make_wt_eblocks(self, bins: list) -> list:
        """
        Create wild-type eBlocks
        """
        gene_blocks = []
        for num in range(0, len(bins), 2):
            start_index = self.vector_instance.gene_start_idx + bins[num]
            end_index = self.vector_instance.gene_start_idx + bins[num+1]
            eblock = Eblock(
                start_index=self.vector_instance.circular_index(start_index, len(self.vector_instance.vector.seq)),  # start index in the vector
                end_index=self.vector_instance.circular_index(end_index, len(self.vector_instance.vector.seq)), 
                bin_start=bins[num],  # start index in the gene
                bin_end=bins[num+1],
                sequence=Vector.slice_circular_sequence(self.vector_instance.vector.seq, start_index, end_index))
            gene_blocks.append(eblock)
        gene_blocks = self.renumber_eblock(gene_blocks)
        return gene_blocks
        
    def eblock_index(self, eblock: Eblock, idx_mutation: int) -> int:
        """
        Find the index of a mutation in an eblock.
        """
        if (eblock.start_index > eblock.end_index) and (self.vector_instance.gene_start_idx < eblock.start_index):
            mutation_idx_in_gene  = self.vector_instance.gene_start_idx + idx_mutation
            residues_to_end = len(self.vector_instance.vector.seq) - eblock.start_index
            mutation_index_in_eblock = residues_to_end + mutation_idx_in_gene
        elif (eblock.start_index < eblock.end_index) and (self.vector_instance.gene_start_idx < eblock.start_index):
            mutation_index_in_eblock = (idx_mutation - eblock.start_index) + self.vector_instance.gene_start_idx
        elif (eblock.start_index < eblock.end_index) and (self.vector_instance.gene_start_idx > eblock.start_index):
            mutation_index_in_eblock = self.vector_instance.gene_start_idx - eblock.start_index + idx_mutation
        else:
            raise Exception("Error in mutation index calculation")
        return mutation_index_in_eblock
        
    def eblocks_within_range(self, mutation_idx: int):
        """
        For a given mutation index, find the number of overlapping eBlocks
        """
        eblocks = []
        for eblock in self.wt_eblocks:
            if eblock.start_index < eblock.end_index:
                if eblock.start_index < mutation_idx < eblock.end_index:
                    eblocks.append(eblock)
            else:
                if (eblock.start_index < mutation_idx) or (mutation_idx < eblock.end_index):
                    eblocks.append(eblock)
        count = len(eblocks)
        return copy.deepcopy(eblocks), count
    
    def check_wt_codon(self, eblock: Eblock, mut: str):
        """
        Check whether the WT codon at the mutation index is the same as in the proposed mutation
        """
        codon = eblock.sequence[eblock.mutation_start_index-3:eblock.mutation_start_index].upper()
        try:
            result = seq.CodonTable.default_table()[str(codon)]
        except:
            result = None
        if result != mut[0]:
            raise Exception(f"WT codon does not match residue {mut}, but is {result}, the codon is {codon}")
    
    def select_mut_codon(self, res: str):
        """
        Select the most abundant codon for a given residue, 
        based on the relative frequencies of the codons in the selected genome
        """
        return self.most_abundant_codons[res][0]
    
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
                raise Exception(f"eBlock is too long, length is {length_eblock}, maximum length is {self.max_eblock_length}")
            else:
                raise Exception(f"eBlock is too short, length is {length_eblock}, minimum length is {self.min_eblock_length}")
        
    def make_mutant_eblock(self, mutation: Mutation, results: dict) -> dict:
        """
        This function creates the mutated eBlock, based on the WT eBlocks and the mutation.
        """

        if mutation.is_singlemutation:
            eblocks, _ = self.eblocks_within_range(mutation.idx_dna[0])
            eblock = eblocks[0]
            eblock.mutation_start_index = self.eblock_index(eblock, mutation.idx_dna[0])
            self.check_wt_codon(eblock, mutation.mutation[0][0])  # Check if WT codon at index is same residue as mutation
            eblock.mutant_codon = self.select_mut_codon(mutation.mutation[0][-1])
            eblock.sequence = self.mutate_eblock(mutation, eblock)

        elif mutation.is_insert:
            # TODO Move to function (=same als deletion)
            eblock_start, count_start = self.eblocks_within_range(mutation.idx_dna[0])
            eblock_end, count_end = self.eblocks_within_range(mutation.idx_dna[0] + mutation.length_insert)
            eblock = None
            if (count_start > 1) or (count_end > 1):
                for i in eblock_start:
                    for j in eblock_end:
                        if i.name == j.name:
                            eblock = i
            elif (count_start == 1) and (count_end == 1):
                eblock = eblock_start[0]
            print("selected eblock", eblock.name)        
            eblock.mutation_start_index = self.eblock_index(eblock, mutation.idx_dna[0])
            self.check_wt_codon(eblock, mutation.mutation[0][0])
            eblock.insert = self.design_insert(mutation.insert)
            eblock.sequence = self.mutate_eblock(mutation, eblock)
            self.check_eblock_length(eblock.sequence)  # Check if eBlock is too long including the insert

        elif mutation.is_deletion:
            eblock_start, count_start = self.eblocks_within_range(mutation.idx_dna_deletion_begin)
            print("EBLOCK START", eblock_start, count_start)
            eblock_end, count_end = self.eblocks_within_range(mutation.idx_dna_deletion_end)
            print("EBLOCK END", eblock_end, count_end)
            eblock = None
            if (count_start > 1) or (count_end > 1):
                for i in eblock_start:
                    for j in eblock_end:
                        if i.name == j.name:
                            eblock = i
            else:
                eblock = eblock_start[0]
            print(mutation.name, mutation.idx_dna_deletion_begin, mutation.idx_dna_deletion_end)
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
                    selected_eblock = eblock[0]
                    
            all_counts = [counts for _, counts in (self.eblocks_within_range(mut_i) for mut_i in mutation.idx_dna)]
            all_eblocks = []
            for mut_i in mutation.idx_dna:
                eblocks, _ = self.eblocks_within_range(mut_i)
                all_eblocks.append([i.name for i in eblocks])
            counter = Counter()
            for lst in all_eblocks:
                counter.update(set(lst))
            common_eblock = [item for item, count in counter.items() if count == len(all_eblocks)]
            print("common_eblock", common_eblock)
            if len(common_eblock) > 0:
                selected_eblock = [e for e in eblocks if e.name == common_eblock[0]][0]
                print("selected_eblock", selected_eblock)
            else:
                print(f"Multiple eBlocks for multiple mutations, not all mutations in the same eBlock. Skip mutation {mutation.name} {mutation.idx_dna}.")
                return results

            # all_eblocks = [eblock.name for eblock, _ in (self.eblocks_within_range(mut_i) for mut_i in mutation.idx_dna)]
            # all_eblocks = [i.name for i, _ in (self.eblocks_within_range(mut_i) for mut_i in mutation.idx_dna)]
            # lowest_count = min(all_counts)            
            for mut_i in mutation.idx_dna:
                eblock, counts = self.eblocks_within_range(mut_i)
                try:  # Try to find indexes of mutations, based on eblock. Check if they are too close to beginning or end of eblock
                    # TODO Move this to a separate function
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
    
    def make_clones(self):
            
            self.make_dir(dirname='clones')  # Make clones-dir
            original_dir = self.output_dir
            self.set_output_dir(os.path.join(self.output_dir, 'clones'))

            # Loop over all mutations and create mutated vector and features that can be read by snapgene
            for mut, eblock in self.eblocks.items():
                results = {}
                if not (mut.is_deletion) and not (mut.is_insert):
                    results[eblock.name] = [eblock.start_index, eblock.end_index, self.eblock_colors[eblock.block_number]]
                    results[self.gene_instance.seqid] = [self.vector_instance.gene_start_idx, self.vector_instance.gene_end_idx, self.vector_instance.color]

                filename = mut.name
                
                if mut.is_singlemutation:
                    start = self.vector_instance.gene_start_idx -3 + mut.idx_dna[0]
                    end = self.vector_instance.gene_start_idx + mut.idx_dna[0]
                    results[mut.name] = [start, end, self.mutation_instance.colors[mut.type]]

                elif mut.is_insert:
                    start = self.vector_instance.gene_start_idx -3 + mut.idx_dna[0]
                    end = self.vector_instance.gene_start_idx + mut.idx_dna[0] + mut.length_insert
                    results[mut.name] = [start, end, self.mutation_instance.colors[mut.type]]
                    results[eblock.name] = [eblock.start_index, eblock.end_index + mut.length_insert, self.eblock_colors[eblock.block_number]]
                    results[self.gene_instance.seqid] = [self.vector_instance.gene_start_idx, self.vector_instance.gene_end_idx + mut.length_insert, self.vector_instance.color]
           
                elif mut.is_deletion:
                    print(mut.name, mut.idx_dna_deletion_begin, mut.idx_dna_deletion_end)
                    start = self.vector_instance.gene_start_idx -6 + mut.idx_dna_deletion_begin
                    end = self.vector_instance.gene_start_idx -3 + mut.idx_dna_deletion_begin
                    results[mut.name] = [start, end, self.mutation_instance.colors[mut.type]]
                    results[self.gene_instance.seqid] = [self.vector_instance.gene_start_idx, self.vector_instance.gene_end_idx - mut.length_deletion, self.vector_instance.color]
                    if eblock.start_index < eblock.end_index:
                        results[eblock.name] = [eblock.start_index, eblock.end_index - mut.length_deletion, self.eblock_colors[eblock.block_number]]
                    else:
                        restoend = len(self.vector_instance.vector.seq) - eblock.start_index
                        newlength = len(self.vector_instance.vector.seq) - mut.length_deletion
                        newstart = newlength - restoend
                        results[eblock.name] = [newstart, eblock.end_index - mut.length_deletion, self.eblock_colors[eblock.block_number]]

                elif mut.is_multiplemutation:
                    for i, _ in enumerate(mut.idx_dna):
                        start = self.vector_instance.gene_start_idx -3 + mut.idx_dna[i]
                        end = self.vector_instance.gene_start_idx + mut.idx_dna[i]
                        results[mut.mutation[i]] = [start, end, self.mutation_instance.colors[mut.type]]
                        
                self.make_dir(dirname=filename)
                self.check_directory(os.path.join(self.output_dir, filename))
                mutated_vector = self.vector_instance.mutate_vector(eblock.start_index, eblock.end_index, eblock.sequence, mutation=mut)
                self.vector_instance.save_vector(vector=mutated_vector, output_dir=os.path.join(self.output_dir, filename), filename=f"{filename}.dna")
                self.eblocks_to_gff3(eblocks=results, output_dir=os.path.join(self.output_dir, filename), filename=f"{filename}.gff3")
                self.eblocks_to_genbank(mutvector=mutated_vector, eblocks=results, output_dir=os.path.join(self.output_dir, filename), filename=f"{filename}.gb")
                
            self.output_dir = original_dir
        
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
        
    def set_output_dir(self, output_dir: str):
        self.output_dir = output_dir

    def make_dir(self, dirname: str):
        try:
            os.makedirs(os.path.join(self.output_dir, f"{dirname}"))
        except FileExistsError:
            pass

    def make_file(self, directory, filename, header=False):
        try:
            with open(os.path.join(directory, filename), 'r') as f:
                pass
        except FileNotFoundError:
            with open(os.path.join(directory, filename), 'w') as f:
                if header:
                    f.write("\n".join(self.gff3_header(length_sequence=self.vector_instance.length)))
                    f.write("\n")
                else:
                    pass

    def validate_optimization_parameters(self):
        """
        Check if the optimization parameters are set correctly.
        """
        if not self.cost_optimization and not self.amount_optimization:
            raise Exception("Please set either cost_optimization or amount_optimization to True.")
        if self.cost_optimization and self.amount_optimization:
            raise Exception("Please set either cost_optimization or amount_optimization to True, not both.")
    
    def print_line(self, txt):
        if self.verbose:
            print(txt)

    def check_directory(self, directory):
        """
        Check if a directory exists and is empty
        """
        if not os.path.exists(directory):  # Check if the directory exists
            raise FileNotFoundError(f"Directory {directory} does not exist.")
        if os.listdir(directory):  # Check if the directory is empty
            if self.verbose:
                print(f"Directory {directory} is not empty. Files might get overwritten or appended to.")

    def eblocks_to_csv(self, filename="eblocks.csv"):
        """
        Save eBlocks to a CSV file.
        """
        with open(os.path.join(self.output_dir, filename), 'w') as f:
            f.write("eBlock,Start,End,Sequence\n")
            for key, value in self.eblocks.items():
                f.write(f"{value.name},{value.start_index},{value.end_index},{value.sequence}\n")

    def eblocks_to_gff3(self, eblocks: dict, output_dir, type='gene', filename='eblocks.gff3', header=True):
        """
        This function converts the eBlocks to features that can be read by SnapGene.
        """
        self.make_file(output_dir, filename, header=header)
        with open(os.path.join(output_dir, filename), 'a') as f:
            for k, v in eblocks.items():
                line = self.gff3_line(v[0], v[1], k, v[2], type)
                f.write('\t'.join(line) + '\n')

    def eblocks_to_genbank(self, mutvector, eblocks: dict, output_dir, type='gene', filename='eblocks.gb', header=True, max_filename_length=16):
        """
        This function saves a vector to a GenBank (gb) file
        """
        sequence = Seq(mutvector)
        record = SeqRecord(sequence, id=self.gene_instance.seqid, name=self.gene_instance.seqid, description="")
        record.annotations["molecule_type"] = "DNA"
        record.annotations["organism"] = self.vector_instance.organism
        record.annotations["date"] = datetime.today().strftime('%d-%b-%Y').upper()
        record.name = record.name + "_" + filename
        # Limit filename length characters
        if len(record.name) > max_filename_length:
            record.name = record.name[:max_filename_length]
        features = []  # Add eBlock and mutations as features
        for k, v in eblocks.items():
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

    def check_eblocks(self, results: dict):
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
    
    def display_settings(self):
        # Print all class attributes using the __dict__ attribute
        for key, value in self.__dict__.items():
            print(f"{key}: {value}")
    
    def update_settings(self, settings):
        for key, value in settings.items():
            if key == 'cost_optimization':
                self.cost_optimization = bool(value)
            elif key == 'amount_optimization':
                self.amount_optimization = bool(value)
            elif key == 'bp_price':
                self.bp_price = float(value)
            elif key == 'max_eblock_length':
                self.max_eblock_length = int(value)
            elif key == 'min_eblock_length':
                self.min_eblock_length = int(value)
            elif key == 'min_overlap':
                self.min_overlap = int(value)
            elif key == 'min_order':
                self.min_order = int(value)
            elif key == 'codon_usage':
                self.codon_usage = str(value)
            elif key == 'output_dir':
                self.output_dir = str(value)
            elif key == 'verbose':
                self.verbose = bool(value)
            elif key == 'clone_files':
                self.clone_files = bool(value)
            else:
                raise Exception(f"Key {key} not found in class.")
            
    @staticmethod
    def renumber_eblock(eblocks: list):
        sorted_list = sorted(eblocks, key=lambda x: x.bin_start)
        for i, obj in enumerate(sorted_list, start=1):
            obj.name = f"eBlock-{i}"
            obj.block_number = i
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
    This class contains functions to cluster mutations based on the index of the mutations
    """
    def __init__(self,
                 mutation_instance: Mutation,
                 vector_instance: Vector,
                 gene_instance: Gene,
                 max_eblock_length: int = None,
                 min_overlap: int = None,
                 min_eblock_length: int = None,
                 cost_optimization: bool = None,
                 amount_optimization: bool = None,
                 bp_price: float = None,
                 verbose: bool = None):
        
        self.mutation_instance = mutation_instance
        self.vector_instance = vector_instance
        self.gene_instance = gene_instance
        self.max_eblock_length = max_eblock_length
        self.min_overlap = min_overlap
        self.min_eblock_length = min_eblock_length
        self.cost_optimization = cost_optimization
        self.amount_optimization = amount_optimization
        self.bp_price = bp_price
        self.verbose = verbose

        self.idxs_constraints = []  # Store indices of constraints
        self.X = []  # Store all indices of the mutations

    def run_clustering(self):
        self.X, self.idxs_constraints = self.extract_indices()
        valid_clusters = self.find_possible_clusters()
        optimal_clustering = self.choose_cluster(valid_clusters)
        return optimal_clustering
        
    def find_possible_clusters(self, threshold_small_clusters=2):
        """
        This function finds all possible clusterings of the mutations based on the index of the mutations.
        """
        possibilities = {} # Store all possible clusterings
        n = 1
        valid_clusters = True
        while (valid_clusters) and (n <= len(self.mutation_instance.mutations)):
            
            print("n", n)
            
            cluster_labels, invalid_constraints = self.kmedoids_clustering(num_clusters=n)
            print("cluster_labels", cluster_labels)
            clusters = self.cluster_to_dict(cluster_labels, self.X)
            
            # clusters = self.add_insert_deletion_sizes(clusters) # TODO Think about insertions and deletions and what to do with them
            cluster_sizes = [max(v) - min(v) for v in clusters.values()]
            clusters_copy = copy.deepcopy(clusters)

            # Fix constraints if there are any invalid constraints
            if (invalid_constraints > 0) and not any(size < (self.min_eblock_length - 2 * self.min_overlap) for size in cluster_sizes):
                clusters_copy = self.fix_constraints(clusters_copy, cluster_labels)
                invalid_constraints = 0
                # TODO Check constraints again

            print("clusters", clusters)
            print("cluster_sizes", cluster_sizes)
            print("invalid_constraints", invalid_constraints)

            cluster_too_small = 0
            cluster_correct = 0
            cluster_too_big = False

            # Check if the cluster sizes are within bounds
            if any(size > (self.max_eblock_length - 2 * self.min_overlap) for size in cluster_sizes):  # Cluster too big
                cluster_too_big = True

            elif any(size < (self.min_eblock_length - 2 * self.min_overlap) for size in cluster_sizes):  # Cluster too small
                exceeding_items = [i for i in cluster_sizes if i < (self.min_eblock_length - 2 * self.min_overlap)]
                cluster_keys = [k for k, v in clusters.items() if max(v) - min(v) < (self.min_eblock_length - 2 * self.min_overlap)]
                min_length_to_add = [(self.min_eblock_length - size) - 2 * self.min_overlap for size in exceeding_items]

                # See if adding some length to the cluster will make it fit the requirements
                for i, j, k in zip(exceeding_items, min_length_to_add, cluster_keys):
                    vals = clusters[k]
                    if i + j <= len(self.gene_instance.sequence):  # Add to the end of the eblock
                        clusters_copy[k].append(max(vals) + j)
                        cluster_too_small += 1
                        cluster_correct += 1
                    else:
                        clusters_copy[k].append(min(vals) - j)
                        cluster_too_small += 1
                        cluster_correct += 1
            else:
                cluster_correct += len(cluster_sizes)

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
                            and (2) insertions and deletions are not too large or too small and (3) the eblocks settings are not too strict.")
        return possibilities
    
    def extract_indices(self):
        X = []  # Store all indices of the mutations
        constraints = []  # Store constraints 
        for i in self.mutation_instance.mutations:
            if i.is_singlemutation:
                X.append(i.idx_dna[0])
            elif i.is_insert:
                start_insert = i.idx_dna[0]
                end_insert = i.idx_dna[0] + i.length_insert
                constraints.append((start_insert, end_insert))
                X.append(start_insert)
                X.append(end_insert)
            elif i.is_deletion:
                X.append(i.idx_dna_deletion_begin)
                X.append(i.idx_dna_deletion_end)
            elif i.is_multiplemutation:
                idxs = []
                for j in i.idx_dna:
                    idxs.append(j)
                    X.append(j)
                constraints.append(tuple(idxs))
        X = np.asarray(X).reshape(-1, 1).flatten()
        constraints_indices = []  # Store constraint indices
        for con in constraints:
            for pair in itertools.combinations(con, 2):
                idx_a = np.where(X == pair[0])[0]
                index_a = int(idx_a[0])  # Take the first match

                idx_b = np.where(X == pair[1])[0]
                index_b = int(idx_b[0])
                constraints_indices.append((index_a, index_b))
        return X, constraints_indices

    def kmedoids_clustering(self, num_clusters: int):
        distance_matrix = squareform(pdist(self.X[:, np.newaxis], metric='euclidean'))  # Precompute the distance matrix

        for idx1, idx2 in self.idxs_constraints:  # Modify the distance matrix based on the constraints so that the distance between the constraints is halved
            distance_matrix[idx1, idx2] *= 0.5
            distance_matrix[idx2, idx1] *= 0.5

        kmedoids = KMedoids(n_clusters=num_clusters, metric='precomputed', random_state=42)
        kmedoids.fit(distance_matrix)

        labels = kmedoids.labels_
        labels = labels.tolist()

        invalid_constraints = self.check_constraints(self.idxs_constraints, labels)
        return labels, invalid_constraints
    
    def choose_cluster(self, clusters: dict) -> dict:
        if self.cost_optimization:
            self.print_line("Optimizing based on price per bp ...")
            costs = {}
            for n, c in clusters.items():
                costs[n] = self.calculate_cost(c)
            lowest_cost = min(costs.values())
            best_clustering_k = [k for k, v in costs.items() if v == lowest_cost][0]
            print("best_clustering_k", best_clustering_k)
            self.print_line(f"Lowest estimated cost: €{lowest_cost} (given price per bp of €{self.bp_price})")
            return clusters[best_clustering_k]
        elif self.amount_optimization:
            self.print_line("Optimizing based on number of eBlocks ...")
            fewest_blocks = {}
            for n, c in clusters.items():
                fewest_blocks[n] = len(c)
            best_clustering_k = [k for k, v in fewest_blocks.items() if v == min(fewest_blocks.values())][0]
            print(f"Fewest number of eBlocks: {min(fewest_blocks.values())}")
            self.print_line(f"Lowest number of eBlocks: {best_clustering_k}")
            return clusters[best_clustering_k]

    # TODO Implement this 
    def add_insert_deletion_sizes(self, clusters):
        """
        Loop over mutation idxs and see whether the insertions and deletions are within bounds
        """
        clusters_copy = copy.deepcopy(clusters)
        for mutation in self.mutation_instance.mutations:
            if mutation.is_insert:
                for key, value in clusters.items():
                    if mutation.idx_dna[0] in value:
                        clusters_copy[key].append(mutation.idx_dna[0] + mutation.length_insert)
            elif mutation.is_deletion:
                for key, value in clusters.items():
                    if mutation.idx_dna[0] in value:
                        clusters_copy[key].append(mutation.idx_dna[0] -  mutation.length_deletion)
        return clusters_copy
            
    def calculate_cost(self, clusters: dict) -> float:
        total_cost = 0
        for _, value in clusters.items():
            min_val = min(value)
            max_val = max(value)
            len_gene_block = (max_val - min_val) + 2 * self.min_overlap  # on both size of the gene block there should be an overhang for primer design
            cost = len_gene_block * self.bp_price * len(value)
            total_cost += cost
        return round(total_cost, 2)
    
    def print_line(self, txt):
        if self.verbose:
            print(txt)

    def check_constraints(self, constraints_indices, labels) -> int:
        """
        Check if the constraints are valid.
        """
        invalid = 0
        for con in constraints_indices:
            con_labels = [labels[i] for i in con]
            same = all(x == con_labels[0] for x in con_labels) # Check if all pairs are in the same cluster
            if not same:
                invalid += 1
        return invalid
        
    def fix_constraints(self, clusters, labels):
        clusterscopy = copy.deepcopy(clusters)
        invalid_constraints = {}
        # Add keys to new dict
        for key in clusters.keys():
            invalid_constraints[key] = []

        for con in self.idxs_constraints:
            con_labels = [labels[i] for i in con]
            if not all(x == con_labels[0] for x in con_labels):
                print("X", self.X[con[0]], self.X[con[1]], "con", con, "con_labels", con_labels)
                distance = self.X[con[0]] - self.X[con[1]]
                if distance > 0:
                    invalid_constraints[con_labels[0]].append(distance)
                    invalid_constraints[con_labels[1]].append(-1 * distance)
                else:
                    invalid_constraints[con_labels[0]].append(-1 * distance)
                    invalid_constraints[con_labels[1]].append(distance)
    
        # Remove empty keys
        to_remove = []
        for k, v in invalid_constraints.items():
            if len(v) == 0:
                to_remove.append(k)
        for k in to_remove:
            del invalid_constraints[k]

        values_to_add = {}
        for key1, list1 in invalid_constraints.items():
            for key2, list2 in invalid_constraints.items():
                if key1 != key2:
                    abs_list1 = [abs(i) for i in list1]
                    abs_list2 = [abs(i) for i in list2]
                    matches = [i for i in abs_list1 if i in abs_list2]
                    if len(matches) > 0:
                        highest = max(matches)
                        # Find in which list the highest value is found
                        if highest in list1:
                            values_to_add[key1, key2] = [highest, key1]
                        else:
                            values_to_add[key1, key2] = [highest, key2]

        # Check if the constraints exists in the cluster
        values_to_add_2 = {}
        for k, v in values_to_add.items():
            for con in self.idxs_constraints:
                con_labels = [labels[i] for i in con]
                if k[0] in con_labels and k[1] in con_labels:
                    values_to_add_2[k] = v

        clean_dict = {}
        for k, v in values_to_add_2.items():
            if k[0] == v[1]:
                clean_dict[k] = v[0]
            
        # Add value to the eblock
        for k, v in clean_dict.items():
            max_value = max(clusters[k[0]])
            print(k, "max_value", max_value, "v", v)
            print("size:", max(clusters[k[0]]) - min(clusters[k[0]]))
            if (max(clusters[k[0]]) - min(clusters[k[0]]) + v <= self.max_eblock_length - 2 * self.min_overlap):  # Check that the eblock is not becoming too big # TODO Move to separate functino                
                print("add value to eblock")
                clusterscopy[k[0]].append(max_value + v)
        
        # TODO Check that if by adding the value to the eblock, the constraints are met
        # TODO Update constraints
        # TODO Update labels
        # self.update_labels(clusterscopy, labels)
        # print(labels)
        # print(clusters)
        # print(clusterscopy)
        return clusterscopy
        
    @staticmethod
    def cluster_to_dict(labels, indexes):
        clusters = {}
        for i, j in zip(labels, indexes):
            if i not in clusters:
                clusters[i] = []
            clusters[i].append(j)
        return clusters