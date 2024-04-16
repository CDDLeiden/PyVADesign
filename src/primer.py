import os
import sys
import math
import difflib
import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt, GC

from .mutation import Mutation
from .sequence import Plasmid
from .eblocks import EblockDesign
from .utils import Utils, SnapGene, OutputToFile


# TODO Think of a solution when the primers are designed at the very beginning of the gene
# TODO Add examples of complementary and hairpin structures to the tests directory


class DesignPrimers:
    """
    This class designs IVA primers to open-up the expression plasmid and Sanger sequencing primers for validation. 
    """

    def __init__(self,
                 eblocks_design_instance: EblockDesign,
                 mutation_instance: Mutation,
                 sequence_instance: Plasmid,
                 snapgene_instance = None,
                 output_dir: str = None,
                 verbose: bool = True):

        self.eblocks_design_instance = eblocks_design_instance
        self.mutation_instance = mutation_instance
        self.sequence_instance = sequence_instance
        self.output_dir = output_dir
        self.snapgene_instance = snapgene_instance
        self.verbose = verbose

    def run_design(self):
        """
        Run the design of the primers
        """

        primerinstance = Primer()
        primers = {}
        with OutputToFile(os.path.join(self.output_dir, 'primer-warnings.txt')):
            ivaprimers = self.design_iva_primer()
            for i in ivaprimers:
                primers[i.name] = i.sequence_5to3

            seqprimers = self.design_seq_primer()
            for i in seqprimers:
                primers[i.name] = i.sequence_5to3
            self.snapgene_instance.primers_to_fasta(primers=primers)

            # Check primers for hairpin formation and multiple binding sites
            for k, v in primers.items():
                max_hairpin, _, _ = primerinstance.check_hairpin(v)
                n_binding_sites = primerinstance.check_multiple_binding_sites(vector=self.sequence_instance.vector.seq, sequence=v)
                        
    def design_iva_primer(self):
        """
        Design IVA primers to open-up the expression plasmid
        """

        fw_sequence = self.sequence_instance.sequence
        rv_sequence = self.sequence_instance.reverse_complement(fw_sequence)
        
        ivaprimerdesign = IVAprimer()
        ivaprimers = []  # Store all IVA primers in a list

        for eblock in self.eblocks_design_instance.wt_eblocks:  # Loop over gene blocks and design IVA primers (starting with initial sequences that are optimized later on)

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
            rv_combined = ivaprimerdesign.combine_Rv_primer(final_rv_oh, final_rv_template)

            iva_fw_primer = IVAprimer(name=ivaprimerdesign.Fw_name(eblock.block_number),
                                      sequence_5to3=''.join(fw_combined),
                                      is_forward=True,
                                      template=''.join(final_fw_template),
                                      overhang=''.join(final_fw_oh))
            ivaprimers.append(iva_fw_primer)
            
            iva_rv_primer = IVAprimer(name=ivaprimerdesign.Rv_name(eblock.block_number),
                                      sequence_5to3=self.sequence_instance.invert_sequence(''.join(rv_combined)),
                                      is_reverse=True,
                                      template=''.join(final_rv_template),
                                      overhang=''.join(final_rv_oh))
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
                    seqprimers.append(prim)
                    count += 1
                else:
                    print("No primer found within the desired range")

        seqprimerdesign.primers = seqprimers
        seqprimerdesign.mapped_primers = self.map_seqprimers_to_mutations(seqprimers)
        self.mapped_seqprimers_to_txt(seqprimerdesign.mapped_primers)  # Save mapped primers to mutations to file

        # Save primer information
        df = seqprimerdesign.primers_to_dataframe(seqprimers)
        df.to_csv(os.path.join(self.output_dir, 'SEQprimers.csv'), index=False)
        
        return seqprimers
        
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
                 is_reverse: bool = False):

        self.name = name
        self.sequence_5to3 = sequence_5to3
        self.complementarity_threshold = complementarity_threshold
        self.hairpin_threshold = hairpin_threshold
        self.is_forward = is_forward
        self.is_reverse = is_reverse

    def Tm(self, sequence):
        return round(mt.Tm_NN(sequence), 2)
    
    def gc_content(self, primer):
        return round(GC(primer), 2)
    
    def check_hairpin(self, sequence: str):
        # TODO Check using this website http://biotools.nubic.northwestern.edu/OligoCalc.html
        max_hairpin = 0
        for i in range(0, len(sequence) +1):
            for j in range(1, len(sequence) +1):
                fragment = sequence[i:i+j]
                complementary = Plasmid.reverse_complement(fragment)
                complementary_inverted = complementary[::-1]
                if len(fragment) >= self.hairpin_threshold:
                    # Search for complementary regions in sequence
                    if complementary_inverted in sequence:
                        if len(fragment) > max_hairpin:
                            max_hairpin = len(fragment)
                            # print(f"Hairpin formation in sequence {sequence} exceeds threshold of {self.hairpin_threshold} ({len(fragment)}) with {fragment} and {complementary_inverted}")
        return max_hairpin, fragment, complementary_inverted
    
    @staticmethod
    def check_multiple_binding_sites(vector: str, sequence: str):
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

                 # IVA primer properties
                 max_overhang_temp_IVA: int = 50,
                 max_template_temp_IVA: int = 60,
                 max_oh_length: int = 15,
                 init_size: int = 10):  # Initial size of the primer
        
        super().__init__(name=name, 
                         sequence_5to3=sequence_5to3, 
                         complementarity_threshold=complementarity_threshold, 
                         hairpin_threshold=hairpin_threshold, 
                         is_forward=is_forward,
                         is_reverse=is_reverse)
                         
        self.max_overhang_temp_IVA = max_overhang_temp_IVA
        self.max_template_temp_IVA = max_template_temp_IVA
        self.max_oh_length = max_oh_length
        self.init_size = init_size

        self.template  = template
        self.overhang = overhang

        self.primers = []
        self.primers_df: pd.DataFrame = pd.DataFrame(columns=['eBlock', 
                                                              'Overhang', 
                                                              'Template', 
                                                              'direction'
                                                              'Tm Template',
                                                              'Tm Overhang']) 
                                                              # 'end position',
                                                              # 'begin position'])

    def Fw_name(self, n: int):
        return f"IVA_Fw_eBlock_{n}"
    
    def Rv_name(self, n: int):
        return f"IVA_Rv_eBlock_{n}"

    def Fw_overhang(self, block_end, fw_sequence, size=15):
        fw_oh = fw_sequence[block_end-size:block_end]
        return fw_oh
        
    def Fw_template(self, block_end, fw_sequence, size=20):
        fw_template = fw_sequence[block_end:block_end+size]
        return fw_template

    def Rv_template(self, block_begin, rv_sequence, size=20):
        rv_template = rv_sequence[block_begin-size:block_begin]
        return rv_template

    def Rv_overhang(self, block_begin, rv_sequence, size=15):
        rv_oh = rv_sequence[block_begin:block_begin+size]
        return rv_oh
    
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
            # TODO CHANGE
            print(f"The overhang temperatures for Fw and Rv primer of {primer1.name} exceed max Tm difference of {threshold} degrees")
        if dTm_templates > threshold:
            print(f"The template temperatures for Fw and Rv primer of {primer1.name} exceed max Tm difference {threshold} degrees")
            
    def check_complementarity(self, primer1, primer2):
        # Primer pairs should not have complementary regions
        overlap = self.get_overlap(primer1.sequence_5to3, primer2.sequence_5to3)
        if len(overlap) > self.complementarity_threshold:
            print(f"Complementarity between the IVA primers for {primer1.name} exceeds threshold of {self.complementarity_threshold}")
        return overlap
    
    def get_overlap(self, s1, s2):
        s = difflib.SequenceMatcher(None, s1, s2)
        pos_a, _, size = s.find_longest_match(0, len(s1), 0, len(s2)) 
        return s1[pos_a:pos_a+size]
    
    def primers_to_dataframe(self, primers):
        for i in primers:
            self.primers_df = self.primers_df.append({'eBlock': i.name,
                                                        'Overhang': i.overhang,
                                                        'Template': i.template,
                                                        'direction': 'Forward' if i.is_forward else 'Reverse',
                                                        'Tm Template': self.Tm(i.template),
                                                        'Tm Overhang': self.Tm(i.overhang)}, ignore_index=True)
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
                         is_reverse=is_reverse)
                         
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
                                                              'Sequence', 
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
                    if self.min_gc_content_primer_seq <= GC(option) <= self.max_gc_content_primer_seq:
                        if self.gc_clamp:
                            if option[-1].lower() == 'g' or option[-1].lower() == 'c':
                                result[(i, j)] = option
                        else:
                            result[(i, j)] = option
        return result
    
    def primers_to_dataframe(self, primers):
        for i in primers:
            self.primers_df = self.primers_df.append({'eBlock': i.name,
                                                      'Sequence': i.sequence_5to3,
                                                      'direction': 'Forward' if i.is_forward else 'Reverse',
                                                      'Tm': self.Tm(i.sequence_5to3),
                                                      'GC content': self.gc_content(i.sequence_5to3),
                                                      'begin position': i.idx_start_seq,
                                                      'end position': i.idx_end_seq}, ignore_index=True)
        return self.primers_df
    
    def all_rv_primers(self):
        # TODO in case FW does not work
        pass