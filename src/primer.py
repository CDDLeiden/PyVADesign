import os
import sys
import math
import difflib
import random
import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt, GC
from design_gene_blocks import DesignEblocks
from utils_old import load_pickle, extract_filename

from mutation import Mutation
from sequence import Plasmid
from eblocks import Eblocks, EblockDesign


class DesignPrimers:
    """
    This class designs IVA and sequencing primers
    """

    def __init__(self,
                 eblock_instance: Eblocks,
                 eblocks_design_instance: EblockDesign,
                 mutation_instance: Mutation,
                 sequence_instance: Plasmid,
                 output_dir: str = None):

        self.eblock_instance = eblock_instance
        self.eblocks_design_instance = eblocks_design_instance
        self.mutation_instance = mutation_instance
        self.sequence_instance = sequence_instance
        self.output_dir = output_dir

        # Primer design parameters
        self.max_overhang_temp_IVA: int = 50
        self.max_template_temp_IVA: int = 60
        self.max_oh_length: int = 15

        self.min_primer_length_seq: int = 18
        self.max_primer_length_seq: int = 24
        self.min_gc_content_primer_seq: int = 45
        self.max_gc_content_primer_seq: int = 55
        self.gc_clamp: bool = True # TODO Add GC clamp
        self.void_length: int = 100  # Approximate number of nucleotides that is skipped before sequencing starts
        self.max_sequenced_region: int = 600  # Maximum number of nucleotides that can be sequenced with a single primer

        self.complementarity_threshold: int = 4
        self.hairpin_threshold: int = 4

        # TODO Add examples of complementary and hairpin structures to the tests directory
        # TODO Check unique primer binding sequence in the plasmid
        # TODO For sequencing primer check CG lock
        # TODO Check multiple binding sites

        self.primers_IVA: dict = {}
        self.primers_IVA_df: pd.DataFrame = pd.DataFrame(columns=['Eblock', 'fw_sequence', 'rv_sequence', 'FW Overhang', 
                                                                  'FW Template', 'RV Overhang', 'RV Template', 'Tm FW Template',
                                                                  'Tm Rv Template', 'Tm FW Overhang', 'Tm RV Overhang', 'end position',
                                                                  'begin position'])

        # TODO Think about how to store the primers
        self.primers_SEQ: dict = {}

    def find_closest_primer(self, all_primers: dict, position: int):
        closest_key = None
        min_difference = float('inf')
        for key, value in all_primers.items():
            print(key, value)
            difference = abs(position - key[0])
            print(difference)
            if difference < min_difference:
                min_difference = difference
                closest_key = key
        return closest_key
    
    def run_SEQprimer(self):
        # TODO Calculate the how many mutations can be captured with a single primer (for now say 600 bp)
        # Calculate the length of the gene block and how many primers are needed to sequence the whole block
        # First cluster the mutations and calculate the distance between them (how many mutations can be sequenced with single primer?)
        # Then design the primers
        possible_primers = self.all_possible_fw_seqprimers()
        for key, value in possible_primers.items():
            print(key, value)
        for k, v in self.eblocks_design_instance.wt_eblocks.items():
            begin, end = DesignEblocks.gene_block_range(k)
            length = end - begin
            num_primers = math.ceil(length / self.max_sequenced_region)
            print(begin, end, num_primers)
            # Calculate the starting points of the sequencing primer
            size = int(length / num_primers)
            print(size)
            primer_start_index = begin - self.void_length
            print("primer_start_index", primer_start_index)
            end = begin + size
            for i in range(num_primers):
                # Find closest primer to starting position
                closest_primer = self.find_closest_primer(possible_primers, primer_start_index)
                print("closest distance", closest_primer)

    def all_possible_fw_seqprimers(self):
        result = {}
        for i in range(len(self.sequence_instance.sequence)):
            for j in range(i, len(self.sequence_instance.sequence)):
                option = self.sequence_instance.sequence[i:j]
                if self.min_primer_length_seq <= len(option) <= self.max_primer_length_seq:
                    gc_content = GC(option)
                    if self.min_gc_content_primer_seq <= gc_content <= self.max_gc_content_primer_seq:
                        if self.gc_clamp:
                            if option[-1] == 'G' or option[-1] == 'C':
                                result[(i, j)] = option
                        else:
                            result[(i, j)] = option
        return result
    

                
    # def design_fw_SEQprimer(self, eblock, start, end):
    #     # TODO Check if the primer is within the desired parameters
    #     # TODO Check if the primer is unique in the plasmid
    #     # For simplicity, make all primers forward primers
    #     # Start with a primer size of 10 and increase untill the Tm is within the desired range
    #     # If this cannot be found, primers can be shifted 50 bp to the right or the left 
    #     # If this does not work, make a rv primer 
    #     # Make sure to end with a GC clamp
    #     # If Tm too high randomly jump to the left or right and retry. Make sure to calculate the distance to the beginning mutation
        
    #     init_primer = self.sequence_instance.sequence[start-5:start+5]
    #     # Calculate GC content and Tm
    #     gc_content = self.calculate_gc_content(init_primer)
    #     tm = self.Tm(init_primer)
    #     count_left = 9
    #     count_right = 9
    #     results = {} # Store the results of the primer design
    #     tries = 0
    #     print(init_primer, gc_content)
    #     while (tries < 100):
    #         if (gc_content > self.min_gc_content_primer_seq and gc_content < self.max_gc_content_primer_seq):
    #             results['GC content'] = gc_content
    #             count_right += 1
    #             init_primer = self.sequence_instance.sequence[start-count_left:start+count_right]
    #             gc_content = self.calculate_gc_content(init_primer)
    #             tm = self.Tm(init_primer)
    #             print(f"Len {len(init_primer)} Tm {tm} GC {gc_content} Primer {init_primer}")
    #         else:
    #             # Shift starting position randomly to the left or right
    #             start = start + random.choice([-3, 3])
    #             init_primer = self.sequence_instance.sequence[start-count_left:start+count_right]
    #             gc_content = self.calculate_gc_content(init_primer)
    #             tm = self.Tm(init_primer)
    #         # print(tries)
    #         tries += 1
    #     sys.exit()
    #     # TODO Alternative approach > search for locations in the sequence where the GC content is within the desired range

    def calculate_gc_content(self, primer):
        return round(GC(primer), 2)

    def design_rv_SEQprimer(self, eblock, start, end):
        pass

    def run_IVAprimer(self):
        
        # Gene sequence
        # TODO Think of a solution when the primers are designed at the very beginning of the gene
        fw_sequence = self.sequence_instance.sequence
        rv_sequence = self.sequence_instance.reverse_complement(fw_sequence)

        # Loop over gene blocks and design primers
        for eblock, _ in self.eblocks_design_instance.wt_eblocks.items():
            begin_pos, end_pos = DesignEblocks.gene_block_range(eblock)

            # Create initial primers (later to be optimized)
            init_size = 10
            init_fw_oh = self.IVA_Fw_overhang(end_pos, fw_sequence, size=init_size)
            size = self.optimize_size(self.max_overhang_temp_IVA, init_fw_oh, end_pos, init_size, fw_sequence, self.IVA_Fw_overhang)
            final_fw_oh = self.IVA_Fw_overhang(end_pos, fw_sequence, size=size)

            init_fw_template = self.IVA_Fw_template(end_pos, fw_sequence, size=init_size)
            size = self.optimize_size(self.max_template_temp_IVA, init_fw_template, end_pos, init_size, fw_sequence, self.IVA_Fw_template)
            final_fw_template = self.IVA_Fw_template(end_pos, fw_sequence, size)

            init_rv_oh = self.IVA_Rv_overhang(begin_pos, rv_sequence, size=init_size)
            size = self.optimize_size(self.max_overhang_temp_IVA, init_rv_oh, begin_pos, init_size, rv_sequence, self.IVA_Rv_overhang)
            final_rv_oh = self.IVA_Rv_overhang(begin_pos, rv_sequence, size)

            init_rv_template = self.IVA_Rv_template(begin_pos, rv_sequence, size=init_size)
            size = self.optimize_size(self.max_template_temp_IVA, init_rv_template, begin_pos, init_size, rv_sequence, self.IVA_Rv_template)
            final_rv_template = self.IVA_Rv_template(begin_pos, rv_sequence, size)

            # Store primers and optimize 
            primerpair = self.parse_primerpair(eblock, fw_sequence, rv_sequence, final_fw_oh, final_fw_template, final_rv_oh, final_rv_template, end_pos, begin_pos)            

            # Combine primers
            fw_combined = self.combine_primers(final_fw_oh, final_fw_template)
            rv_combined = self.combine_primers(final_rv_oh, final_rv_template)

            primerpair["FW Primer (5>3)"] = ''.join(fw_combined)
            primerpair["RV Primer (5>3)"] = ''.join(rv_combined)

            # Check primers and make sure they are within the desired parameters
            self.Tm_difference(primerpair)
            overlap = self.check_complementarity(primerpair)
            max_hairpin_fw, _, _ = self.check_hairpin(primerpair["FW Primer (5>3)"])
            max_hairpin_rv, _, _ = self.check_hairpin(primerpair["RV Primer (5>3)"])
            n_binding_sites = self.check_multiple_binding_sites(primerpair["FW Primer (5>3)"])
            n_binding_siters = self.check_multiple_binding_sites(primerpair["RV Primer (5>3)"])

            primerpair['Max hairpin length'] = max(max_hairpin_fw, max_hairpin_rv)
            primerpair['Max complementary length'] = len(overlap)

            # Store primers in a dataframe
            df_row = pd.DataFrame([primerpair])
            self.primers_IVA_df = self.primers_IVA_df.append(df_row, ignore_index=True)

        # Save primers to file
        self.primers_IVA_df.to_csv(os.path.join(self.output_dir, 'primers_IVA.csv'), index=False)
        return self.primers_IVA_df

    def parse_primerpair(self, eblock, fw_sequence, rv_sequence, fw_oh, fw_template, rv_oh, rv_template, end, begin):
        primerpair = {}
        primerpair['Eblock'] = eblock
        primerpair['fw_sequence'] = ''.join(fw_sequence)
        primerpair['rv_sequence'] = ''.join(rv_sequence)
        primerpair['FW Overhang'] = ''.join(fw_oh)
        primerpair['FW Template'] = ''.join(fw_template)
        primerpair['RV Overhang'] = ''.join(rv_oh)
        primerpair['RV Template'] = ''.join(rv_template)
        primerpair['Tm FW Template'] = self.Tm(fw_template)
        primerpair['Tm Rv Template'] = self.Tm(rv_template)
        primerpair['Tm FW Overhang'] = self.Tm(fw_oh)
        primerpair['Tm RV Overhang'] = self.Tm(rv_oh)
        primerpair['end position'] = int(end)
        primerpair['begin position'] = int(begin)
        return primerpair

    def Tm(self, sequence):
        return round(mt.Tm_NN(sequence), 2)

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
    
    def IVA_Fw_overhang(self, block_end, fw_sequence, size=15):
        fw_oh = fw_sequence[block_end-size:block_end]
        return fw_oh
        
    def IVA_Fw_template(self, block_end, fw_sequence, size=20):
        fw_template = fw_sequence[block_end:block_end+size]
        return fw_template

    def IVA_Rv_template(self, block_begin, rv_sequence, size=20):
        rv_template = rv_sequence[block_begin-size:block_begin]
        return rv_template

    def IVA_Rv_overhang(self, block_begin, rv_sequence, size=15):
        rv_oh = rv_sequence[block_begin:block_begin+size]
        return rv_oh

    def combine_primers(self, overhang, template_binding):
        return overhang + template_binding
            
    def get_overlap(self, s1, s2):
        s = difflib.SequenceMatcher(None, s1, s2)
        pos_a, _, size = s.find_longest_match(0, len(s1), 0, len(s2)) 
        return s1[pos_a:pos_a+size]

    def Tm_difference(self, primerpair: dict, threshold=3):
        """
        Check if the Tm of the primers are within the desired boundaries, otherwise optimize the primers
        """
        dTm_overhangs = abs(primerpair['Tm FW Overhang'] - primerpair['Tm RV Overhang'])
        dTm_templates = abs(primerpair['Tm FW Template'] - primerpair['Tm Rv Template'])
        if dTm_overhangs > threshold:
            print(f"The overhang temperatures for Fw and Rv primer of {primerpair['Eblock']} exceed max Tm difference of {threshold} degrees")
        if dTm_templates > threshold:
            print(f"The template temperatures for Fw and Rv primer of {primerpair['Eblock']} exceed max Tm difference {threshold} degrees")

    def check_complementarity(self, primerpair: dict):
        # Primer pairs should not have complementary regions
        overlap = self.get_overlap(primerpair["FW Primer (5>3)"], primerpair["RV Primer (5>3)"])
        if len(overlap) > self.complementarity_threshold:
            print(f"Complementarity between the IVA primers for {primerpair['Eblock']} exceeds threshold of {self.complementarity_threshold}")
        return overlap

    def check_hairpin(self, primer: str):
        # TODO Check using this website http://biotools.nubic.northwestern.edu/OligoCalc.html
        max_hairpin = 0
        for i in range(0, len(primer) +1):
            for j in range(1, len(primer) +1):
                fragment = primer[i:i+j]
                complementary = Plasmid.reverse_complement(fragment)
                complementary_inverted = complementary[::-1]
                if len(fragment) >= self.hairpin_threshold:
                    # Search for complementary regions in primer
                    if complementary_inverted in primer:
                        if len(fragment) > max_hairpin:
                            max_hairpin = len(fragment)
                            # print(f"Hairpin formation in primer {primer} exceeds threshold of {self.hairpin_threshold} ({len(fragment)}) with {fragment} and {complementary_inverted}")
        return max_hairpin, fragment, complementary_inverted
    
    def check_multiple_binding_sites(self, primer: str):
        # TODO Write some checks for this function
        count = 0
        for i in range(len(self.sequence_instance.vector.seq) - len(primer) + 1):
            sub = self.sequence_instance.vector.seq[i:i + len(primer)]
            unique_chars = set(sub)
            if len(unique_chars) <= 2 or (len(unique_chars) == 3 and sub.count(sub[0]) == 2):
                count += 1
        if count > 1:
            print(f"Multiple binding sites for primer {primer} in the vector sequence")
        return count