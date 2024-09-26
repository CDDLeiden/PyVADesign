import os
import sys
import math
import primer3
import difflib
import pandas as pd
import biotite.sequence as seq
from Bio.SeqUtils import gc_fraction

from .mutation import Mutation
from .sequence import Vector, Gene
from .eblocks import EblockDesign
from .utils import OutputToFile, SnapGene


# TODO add primer data to genbankl file feature is  "primer_bind" 
# TODO Add examples of complementary and hairpin structures to the tests directory


class DesignPrimers:
    """
    This class designs IVA primers to open-up the expression plasmid and Sanger sequencing primers for validation. 
    """

    def __init__(self,
                 eblocks_design_instance: EblockDesign,
                 mutation_instance: Mutation,
                 vector_instance: Vector,
                 gene_instance: Gene,
                 snapgene_instance = None,
                 output_dir: str = None,
                 verbose: bool = True):

        self.eblocks_design_instance = eblocks_design_instance
        self.mutation_instance = mutation_instance
        self.vector_instance = vector_instance
        self.gene_instance = gene_instance
        self.output_dir = output_dir
        self.snapgene_instance = snapgene_instance
        self.verbose = verbose

    def run_design(self):
        """
        Run the design of the primers
        """

        primerinstance = Primer()
        snapgene_instance = SnapGene(sequence_instance=self.sequence_instance, output_dir=self.output_dir)
        primers = {}
        # with OutputToFile(os.path.join(self.output_dir, 'primer-warnings.txt')):  # Save warnings to file
        ivaprimers = self.design_iva_primer()
        for i in ivaprimers:
            primers[i.name] = i.sequence_5to3

        # Add primers to Genbank file
        for mut, eblock in self.eblocks_design_instance.eblocks.items():
            matching_fw_primer = [i for i in ivaprimers if i.name == f"IVA_Fw_eBlock_{eblock.block_number}"]
            matching_rv_primer = [i for i in ivaprimers if i.name == f"IVA_Rv_eBlock_{eblock.block_number}"]
            snapgene_instance.add_primers_to_genbank_file(genbank_file=os.path.join(self.output_dir, 'clones', f"{mut.name}", f"{mut.name}.gb"), primer=matching_fw_primer[0])
            snapgene_instance.add_primers_to_genbank_file(genbank_file=os.path.join(self.output_dir, 'clones', f"{mut.name}", f"{mut.name}.gb"), primer=matching_rv_primer[0])

        seqprimers, mapping = self.design_seq_primer()
        for i in seqprimers:
            primers[i.name] = str(i.sequence_5to3)
        self.snapgene_instance.primers_to_fasta(primers=primers, directory=self.output_dir, filename='primers.fasta')
        for k, v in mapping.items():
            for i in v:
                for s in seqprimers:
                    if s.name == k:
                        snapgene_instance.add_primers_to_genbank_file(genbank_file=os.path.join(self.output_dir, 'clones', f"{i}", f"{i}.gb"), primer=s.sequence_5to3)

        # for k, v in primers.items():  # Check primers for hairpin formation and multiple binding sites
        #     hairpin = primerinstance.calc_hairpin(v)
        #     n_binding_sites = primerinstance.check_multiple_binding_sites(vector=self.sequence_instance.vector.seq, sequence=v)

        # TODO Check fot 
        # hairpin = primer3.calc_hairpin(seq2)
        # het = primer3.calc_heterodimer(seq1, seq2)
        # homo = primer3.calc_homodimer(seq1)
        # add to ebl
                        
    def design_iva_primer(self):
        """
        Design IVA primers to open-up the expression plasmid
        """

        fw_sequence = str(self.sequence_instance.vector.seq.lower())
        rv_sequence = str(seq.NucleotideSequence(fw_sequence).complement())
        
        ivaprimerdesign = IVAprimer()
        ivaprimers = []  # Store all IVA primers in a list
        ivaprimerdesign.vector_length = len(self.sequence_instance.vector.seq)

        for eblock in self.eblocks_design_instance.wt_eblocks:  # Loop over gene blocks and design IVA primers (starting with initial sequences that are optimized later on)

            eblock.start_index = self.sequence_instance.circular_index(eblock.start_index, len(self.sequence_instance.vector.seq))
            vector_length = len(self.sequence_instance.vector.seq)

            init_fw_oh = ivaprimerdesign.Fw_overhang(eblock.end_index, fw_sequence, size=ivaprimerdesign.init_size)
            size = ivaprimerdesign.optimize_size(ivaprimerdesign.max_overhang_temp_IVA, init_fw_oh, eblock.end_index, ivaprimerdesign.init_size, fw_sequence, ivaprimerdesign.Fw_overhang)
            final_fw_oh = ivaprimerdesign.Fw_overhang(eblock.end_index, fw_sequence, size=size)

            init_fw_template = ivaprimerdesign.Fw_template(eblock.end_index, fw_sequence, size=ivaprimerdesign.init_size)
            size = ivaprimerdesign.optimize_size(ivaprimerdesign.max_template_temp_IVA, init_fw_template, eblock.end_index, ivaprimerdesign.init_size, fw_sequence, ivaprimerdesign.Fw_template)
            final_fw_template = ivaprimerdesign.Fw_template(eblock.end_index, fw_sequence, size)

            init_rv_oh = ivaprimerdesign.Rv_overhang(eblock.start_index, rv_sequence, size=ivaprimerdesign.init_size)
            size = ivaprimerdesign.optimize_size(ivaprimerdesign.max_overhang_temp_IVA, init_rv_oh, eblock.start_index, ivaprimerdesign.init_size, rv_sequence, ivaprimerdesign.Rv_overhang)
            final_rv_oh = ivaprimerdesign.Rv_overhang(eblock.start_index, rv_sequence, size)

            init_rv_template = ivaprimerdesign.Rv_template(eblock.start_index, rv_sequence, size=ivaprimerdesign.init_size)
            size = ivaprimerdesign.optimize_size(ivaprimerdesign.max_template_temp_IVA, init_rv_template, eblock.start_index, ivaprimerdesign.init_size, rv_sequence, ivaprimerdesign.Rv_template)
            final_rv_template = ivaprimerdesign.Rv_template(eblock.start_index, rv_sequence, size)

            # Combine template and overhang sequences
            fw_combined = ivaprimerdesign.combine_Fw_primer(final_fw_oh, final_fw_template)
            rv_combined = ivaprimerdesign.combine_Rv_primer(final_rv_template, final_rv_oh)

            iva_fw_primer = IVAprimer(name=ivaprimerdesign.Fw_name(eblock.block_number),
                                      sequence_5to3=''.join(fw_combined),
                                      is_forward=True,
                                      template=''.join(final_fw_template),
                                      overhang=''.join(final_fw_oh))
            
            idx_start, idx_end = IVAprimer.find_index_in_vector(fw_sequence, ''.join(fw_combined))
            print(f"Index start {idx_start} and index end {idx_end}")
            iva_fw_primer.idx_start = idx_start[0]
            iva_fw_primer.idx_end = idx_end[0]

            ivaprimers.append(iva_fw_primer)
            
            iva_rv_primer = IVAprimer(name=ivaprimerdesign.Rv_name(eblock.block_number),
                                      sequence_5to3=str(seq.NucleotideSequence(''.join(rv_combined)).reverse()), # sequence_5to3=self.sequence_instance.invert_sequence(''.join(rv_combined)),
                                      is_reverse=True,
                                      template=''.join(final_rv_template),
                                      overhang=''.join(final_rv_oh))
            
            idx_start, idx_end = IVAprimer.find_index_in_vector(rv_sequence, ''.join(rv_combined))
            print(f"Index start {idx_start} and index end {idx_end}")
            iva_rv_primer.idx_start = idx_start[0]
            iva_rv_primer.idx_end = idx_end[0]

            ivaprimers.append(iva_rv_primer)
                
            # Check primers and make sure they are within the desired parameters
            ivaprimerdesign.Tm_difference(iva_fw_primer, iva_rv_primer)
            overlap = ivaprimerdesign.check_complementarity(iva_fw_primer, iva_rv_primer)

        # Save primer information
        df = ivaprimerdesign.primers_to_dataframe(ivaprimers)
        df.to_csv(os.path.join(self.output_dir, 'IVAprimers.csv'), index=False)

        return ivaprimers

    def design_seq_primer(self):
        """
        Design sequencing primers to validate the mutations
        """
        seqprimerdesign = SEQprimer()
        seqprimers = []  # Store all sequencing primers in a list

        possible_primers = seqprimerdesign.all_fw_primers(gene_sequence=self.sequence_instance.sequence)  # Find all possible primers that fit the desired parameters
        
        count = 1  # Counter for the primer names
        for i in self.eblocks_design_instance.wt_eblocks:  
            length = i.end_index - i.start_index  # How many primers needed to sequence the whole gene block
            num_primers = math.ceil(length / seqprimerdesign.max_sequenced_region)
            
            size = int(length / num_primers)
            primer_start_indexes = [(i.start_index - seqprimerdesign.void_length) + (k * size) for k in range(num_primers)]  # Calculate the starting points of the sequencing primer
            primer_end_indexes = [i.start_index + (k * size) for k in range(1, num_primers + 1)]
            
            for si, ei in zip(primer_start_indexes, primer_end_indexes):
                closest_range, closest_primer, diff = seqprimerdesign.find_closest_primer(possible_primers, si)  # Find closest primer to starting position         
                if diff < seqprimerdesign.max_difference:
                    prim = SEQprimer(name=seqprimerdesign.get_name(count, i.block_number),
                                     sequence_5to3=closest_primer,
                                     is_forward=True,
                                     idx_start_seq=si + seqprimerdesign.void_length - diff,
                                     idx_end_seq=ei - diff,
                                     primer_idx_strt=closest_range[0],
                                     primer_idx_end=closest_range[1])
                    
                    start_idx, end_idx = SEQprimer.find_index_in_vector(str(self.sequence_instance.vector.seq).lower(), str(closest_primer).lower())
                    prim.idx_start = start_idx[0]
                    prim.idx_end = end_idx[0]

                    seqprimers.append(prim)
                    count += 1
                else:
                    print("No primer found within the desired range")

        mapped_primers = self.map_seqprimers_to_mutations(seqprimers)
        seqprimerdesign.mapped_primers = mapped_primers
        self.mapped_seqprimers_to_txt(seqprimerdesign.mapped_primers)  # Save mapped primers to mutations to file

        # Save primer information
        print(seqprimers)
        for i in seqprimers:
            print(i.sequence_5to3, type(i.sequence_5to3))
        #     print(i.name, i.idx_start, i.idx_end, i.5to3sequence)
        
        # Convert sequences to bytes
        df = seqprimerdesign.primers_to_dataframe(seqprimers)
        df.to_csv(os.path.join(self.output_dir, 'SEQprimers.csv'), index=False)
        
        return seqprimers, mapped_primers
        
    def map_seqprimers_to_mutations(self, primers: list):
        """Map the sequenced regions to the mutations that can be validated with the sequencing primers"""
        mapped_primers = {}
        for i in primers:
            mapped_primers[i.name] = []
            for mutation in self.mutation_instance.mutations:
                if i.idx_start_seq <= mutation.idx_dna[0] <= i.idx_end_seq:
                    mapped_primers[i.name].append(mutation.mutation)
        return mapped_primers
    
    def mapped_seqprimers_to_txt(self, primers: dict, filename='SEQprimers-mapped-mutations.txt'):
        with open(os.path.join(self.output_dir, filename), 'w') as f:
            for k, v in primers.items():
                f.write(f"{k}\n")
                for i in v:
                    if type(i) == list:
                        tmp = str(i)[1:-1]
                        tmp = tmp.replace("'", "")
                        f.write(f"\t{tmp}\n")
                    else:
                        i.replace("'", "")
                        f.write(f"\t{i}\n")


    
class Primer:
    def __init__(self,
                 name: str = None,
                 sequence_5to3: str = None,
                 complementarity_threshold: int = 4,
                 hairpin_threshold: int = 4,
                 is_forward: bool = False,
                 is_reverse: bool = False,
                 
                 idx_start: int = None,
                 idx_end: int = None):

        self.name = name
        self.sequence_5to3 = sequence_5to3
        self.complementarity_threshold = complementarity_threshold
        self.hairpin_threshold = hairpin_threshold
        self.is_forward = is_forward
        self.is_reverse = is_reverse
        self.idx_start = idx_start
        self.idx_end = idx_end

    def Tm(self, sequence):
        return round(primer3.calc_tm(str(sequence)), 2)
    
    def gc_content(self, primer):
        return round(100 * gc_fraction(primer, ambiguous="ignore"), 2)
    
    def calc_hairpin(self, sequence):
        """
        Calculate the hairpin formation in a given sequence
        """
        return primer3.calc_hairpin(str(sequence))
    
    @staticmethod
    def find_index_in_vector(vector, primer):
        circular_string = vector + vector
        start_index = 0
        start_indexes = []
        end_indexes = []

        while True:
            index = str(circular_string).find(primer, start_index)
            if index == -1:
                break
            if index >= len(vector):
                new_index = index - len(vector)
                if new_index not in start_indexes:
                    start_indexes.append(new_index)
                    end_indexes.append(new_index + len(primer))
            else:
                if index not in start_indexes:
                    start_indexes.append(index)
                    end_indexes.append(index + len(primer))
            start_index = index + len(primer)

        if len(start_indexes) > 1:
            # TODO What to do here? add them all?
            print(f"Multiple binding sites for primer {primer} in the vector sequence")
            sys.exit()

        return start_indexes, end_indexes
    
    @staticmethod
    def check_multiple_binding_sites(vector: str, sequence: str):
        # TODO Maybe use primerblast for this
        # TODO Write some checks for this function
        # TODO Check for exact cases where the sequence binds multiple times and almost exact cases where only few characters are different
        count = 0
        for i in range(len(vector) - len(sequence) + 1):
            sub = vector[i:i + len(sequence)]
            unique_chars = set(sub)
            if len(unique_chars) <= 2 or (len(unique_chars) == 3 and sub.count(sub[0]) == 2):
                count += 1
        if count > 1:
            print(f"Multiple binding sites for sequence {sequence} in the vector sequence")
        return count
    


class IVAprimer(Primer, DesignPrimers):
    def __init__(self,
                 
                 # Primer properties
                 name: str = None,
                 sequence_5to3: str = None,
                 complementarity_threshold: int = 4,
                 hairpin_threshold: int = 4,
                 is_forward: bool = False,
                 is_reverse: bool = False,
                 
                 template: str = None,
                 overhang: str = None,

                 idx_start: int = None,
                 idx_end: int = None,

                 # IVA primer properties
                 max_overhang_temp_IVA: int = 50,
                 max_template_temp_IVA: int = 60,
                 max_oh_length: int = 15,
                 init_size: int = 10,  # Initial size of the primer
                 vector_length: int = 0):  
        
        super().__init__(name=name, 
                         sequence_5to3=sequence_5to3, 
                         complementarity_threshold=complementarity_threshold, 
                         hairpin_threshold=hairpin_threshold, 
                         is_forward=is_forward,
                         is_reverse=is_reverse,
                         idx_start=idx_start,
                         idx_end=idx_end)
                         
        self.max_overhang_temp_IVA = max_overhang_temp_IVA
        self.max_template_temp_IVA = max_template_temp_IVA
        self.max_oh_length = max_oh_length
        self.init_size = init_size

        self.template  = template
        self.overhang = overhang

        self.vector_length = vector_length

        self.primers = []
        self.primers_df: pd.DataFrame = pd.DataFrame(columns=['eBlock', 
                                                              'Overhang', 
                                                              'Template', 
                                                              'direction'
                                                              'Tm Template',
                                                              'Tm Overhang']) 

    def Fw_name(self, n: int):
        return f"IVA_Fw_eBlock_{n}"
    
    def Rv_name(self, n: int):
        return f"IVA_Rv_eBlock_{n}"
    
    def Fw_overhang(self, block_end, fw_sequence, size=15):
        fw_oh = fw_sequence[Plasmid.circular_index(block_end-size, self.vector_length):Plasmid.circular_index(block_end, self.vector_length)]
        return fw_oh
        
    def Fw_template(self, block_end, fw_sequence, size=20):
        fw_template = fw_sequence[Plasmid.circular_index(block_end, self.vector_length):Plasmid.circular_index(block_end+size, self.vector_length)]
        return fw_template
    
    def Rv_overhang(self, block_begin, rv_sequence, size=15):
        begin = Plasmid.circular_index(block_begin, self.vector_length)
        end = Plasmid.circular_index(block_begin+size, self.vector_length)
        if begin < end:
            rv_oh = rv_sequence[begin:end]
        else:
            rv_oh = rv_sequence[begin:] + rv_sequence[:end]
        return rv_oh

    def Rv_template(self, block_begin, rv_sequence, size=20):
        rv_template = rv_sequence[Plasmid.circular_index(block_begin-size, self.vector_length):Plasmid.circular_index(block_begin, self.vector_length)]
        return rv_template
    
    def combine_Fw_primer(self, overhang, template):
        return overhang + template
    
    def combine_Rv_primer(self, template, overhang):
        return template + overhang
    
    def optimize_size(self, optimum, primer, pos, size, sequence, function):
        tm = self.Tm(primer)
        while (tm < (optimum)):
            size += 1
            primer = function(pos, sequence, size)
            tm = self.Tm(primer)
        while (tm > (optimum)):
            size -= 1
            primer = function(pos, sequence, size)
            tm = self.Tm(primer)
        return size
    
    def Tm_difference(self, primer1, primer2, threshold=3):
        """
        Check if the Tm of the primers are within the desired boundaries, otherwise optimize the primers
        """
        dTm_overhangs = abs(self.Tm(primer1.overhang) - self.Tm(primer2.overhang))
        dTm_templates = abs(self.Tm(primer1.template) - self.Tm(primer2.template))
        if dTm_overhangs > threshold:
            print(f"The overhang temperatures for {primer1.name} {primer2.name} exceed max Tm difference of {threshold} degrees")
        if dTm_templates > threshold:
            print(f"The template temperatures for {primer1.name} {primer2.name} exceed max Tm difference {threshold} degrees")
            
    def check_complementarity(self, primer1, primer2):
        overlap = self.get_overlap(primer1.sequence_5to3, primer2.sequence_5to3)
        if len(overlap) > self.complementarity_threshold:
            print(f"Complementarity between the primers {primer1.name} {primer2.name} exceeds threshold of {self.complementarity_threshold}")
        return overlap
    
    def get_overlap(self, s1, s2):
        s = difflib.SequenceMatcher(None, s1, s2)
        pos_a, _, size = s.find_longest_match(0, len(s1), 0, len(s2)) 
        return s1[pos_a:pos_a+size]
    
    def primers_to_dataframe(self, primers):
        self.primers_df.dropna(axis=1, inplace=True, how='all')
        for i in primers:
            new_row = pd.DataFrame({'eBlock': str(i.name),
                                    'Overhang': str(i.overhang),
                                    'Template': str(i.template),
                                    'direction': 'Forward' if i.is_forward else 'Reverse',
                                    'Tm Template': float(self.Tm(i.template)),
                                    'Tm Overhang': float(self.Tm(i.overhang))}, index=[0])
            self.primers_df = pd.concat([self.primers_df, new_row], ignore_index=True)
        return self.primers_df
    


class SEQprimer(Primer, DesignPrimers):
    def __init__(self,
                 
                 # Primer properties 
                 name: str = None,
                 sequence_5to3: str = None,
                 complementarity_threshold: int = 4,
                 hairpin_threshold: int = 4,
                 is_forward: bool = False,
                 is_reverse: bool = False,

                 idx_start: int = None,
                 idx_end: int = None,

                # SEQ primer properties
                 min_primer_length_seq: int = 18,
                 max_primer_length_seq: int = 24,
                 min_gc_content_primer_seq: int = 45,
                 max_gc_content_primer_seq: int = 55,
                 gc_clamp: bool = True,
                 void_length: int = 100,  # Approximate number of nucleotides that is skipped before sequencing starts
                 max_difference: int = 100,  
                 max_sequenced_region: int = 600,
                 idx_start_seq: int = -1,
                 idx_end_seq: int = -1,
                 primer_idx_strt: int = -1,
                 primer_idx_end: int = -1):  # Maximum number of nucleotides that can be sequenced using a single primer
        
        super().__init__(name=name,
                         sequence_5to3=sequence_5to3,
                         complementarity_threshold=complementarity_threshold,
                         hairpin_threshold=hairpin_threshold,
                         is_forward=is_forward,
                         is_reverse=is_reverse,
                         idx_start=idx_start,
                         idx_end=idx_end)
                         
        self.min_primer_length_seq = min_primer_length_seq
        self.max_primer_length_seq = max_primer_length_seq
        self.min_gc_content_primer_seq = min_gc_content_primer_seq
        self.max_gc_content_primer_seq = max_gc_content_primer_seq
        self.gc_clamp = gc_clamp
        self.void_length = void_length
        self.max_difference = max_difference
        self.max_sequenced_region = max_sequenced_region
        self.idx_start_seq = idx_start_seq
        self.idx_end_seq = idx_end_seq
        self.primer_idx_strt = primer_idx_strt
        self.primer_idx_end = primer_idx_end

        self.primers = []
        self.primers_df: pd.DataFrame = pd.DataFrame(columns=['eBlock', 
                                                              'sequence', 
                                                              'direction',
                                                              'Tm',
                                                              'GC content',
                                                              'begin position',
                                                              'end position'])
    
        self.mapped_primers = {}

    def get_name(self, n: int, number: int):
        return f"SEQ{n}_eBlock-{number}"

    def find_closest_primer(self, all_primers: dict, position: int):
        """
        Find the closest primer to a given position
        """
        closest_range = None
        closest_primer = None
        min_difference = float('inf')
        for prim_idx, prim in all_primers.items():
            difference = position - prim_idx[0]
            if abs(difference) < min_difference:
                min_difference = difference
                closest_range = prim_idx
                closest_primer = prim
                min_difference = difference
        return closest_range, closest_primer, min_difference

    def all_fw_primers(self, gene_sequence: str):
        result = {}
        for i in range(len(gene_sequence)):
            for j in range(i, len(gene_sequence)):
                option = gene_sequence[i:j]
                if self.min_primer_length_seq <= len(option) <= self.max_primer_length_seq:
                    if self.min_gc_content_primer_seq <= self.gc_content(option) <= self.max_gc_content_primer_seq:
                        if self.gc_clamp:
                            if option[-1].lower() == 'g' or option[-1].lower() == 'c':
                                result[(i, j)] = option
                        else:
                            result[(i, j)] = option
        return result
    
    def primers_to_dataframe(self, primers):
        self.primers_df.dropna(axis=1, inplace=True, how='all')
        for i in primers:
            new_row = pd.DataFrame({'eBlock': str(i.name),
                                    'sequence': str(i.sequence_5to3),
                                    'direction': 'Forward' if i.is_forward else 'Reverse',
                                    'Tm': float(self.Tm(str(i.sequence_5to3))),
                                    'GC content': float(self.gc_content(str(i.sequence_5to3))),
                                    'begin position': int(i.idx_start_seq),
                                    'end position': int(i.idx_end_seq)}, index=[0])
            self.primers_df = pd.concat([self.primers_df, new_row], ignore_index=True)
        return self.primers_df
        
    def all_rv_primers(self):
        # TODO in case FW does not work, implement this
        pass