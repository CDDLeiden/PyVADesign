import os
import sys
import difflib
import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt
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
                 sequence_instance: Plasmid):

        self.eblock_instance = eblock_instance
        self.eblocks_design_instance = eblocks_design_instance
        self.mutation_instance = mutation_instance
        self.sequence_instance = sequence_instance

        # Primer design parameters
        self.max_overhang_temp_IVA = 50
        self.max_template_temp_IVA = 60
        self.max_oh_length = 15

        self.min_primer_length_seq = 18
        self.max_primer_length_seq = 24
        self.min_gc_content_primer_seq = 0.45
        self.max_gc_content_primer_seq = 0.55
        self.gc_clamp = False # TODO Add GC clamp
        self.void_length = 100  # Approximate number of nucleotides that is skipped before sequencing

        self.complementarity_threshold = 4
        self.hairpin_threshold = 3
        # TODO Check self-complementarity
        # TODO Check hairpin formation
        # TODO Add examples of these to the tests directory

        self.primers_IVA: dict = {}
        self.primers_IVA_df: pd.DataFrame = pd.DataFrame(columns=['Eblock', 'fw_sequence', 'rv_sequence', 'FW Overhang', 
                                                                  'FW Template', 'RV Overhang', 'RV Template', 'Tm FW Template',
                                                                  'Tm Rv Template', 'Tm FW Overhang', 'Tm RV Overhang', 'end position',
                                                                  'begin position'])

        # TODO Think about how to store the primers
        self.primers_SEQ: dict = {}

    def run_SEQprimer(self):
        # TODO
        # First cluster the mutations and calculate the distance between them (how many mutations can be sequenced with single primer?)
        # Then design the primers
        pass

    def run_IVAprimer(self):
        
        # Gene sequence
        # TODO Think of a solution when the primers are designed at the very beginning of the gene
        fw_sequence = self.sequence_instance.sequence
        rv_sequence = self.sequence_instance.reverse_complement(fw_sequence)

        # Loop over gene blocks and design primers
        for eblock, eblock_seq in self.eblocks_design_instance.wt_eblocks.items():
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
            # primerpair = self.optimize_primerpair_overhang(primerpair)
            # primerpair = self.optimize_primerpair_template(primerpair)

            # Combine primers
            fw_combined = self.combine_primers(final_fw_oh, final_fw_template)
            rv_combined = self.combine_primers(final_rv_oh, final_rv_template)

            # TODO Add direction (5>3)
            primerpair['FW Primer'] = fw_combined
            primerpair['RV Primer'] = rv_combined

            # TODO Add to test
            # primerpair['FW Primer'] = 'GGGAAAATTCCAGGATCTAT'
            # primerpair['RV Primer'] = 'GGGAAAATTCCAGGATCTAT'

            # Check primers and make sure they are within the desired parameters
            self.Tm_difference(primerpair)
            self.check_complementarity(primerpair)
            self.check_hairpin(primerpair['FW Primer'])
            self.check_hairpin(primerpair['RV Primer'])

            # Store primers in a dataframe
            df_row = pd.DataFrame([primerpair])
            self.primers_IVA_df = self.primers_IVA_df.append(df_row, ignore_index=True)

            # TODO Save primers to file

        return self.primers_IVA_df

    def check_hairpin(self, primer: str):
        # TODO
        pass

    def parse_primerpair(self, eblock, fw_sequence, rv_sequence, fw_oh, fw_template, rv_oh, rv_template, end, begin):
        primerpair = {}
        primerpair['Eblock'] = eblock
        primerpair['fw_sequence'] = fw_sequence
        primerpair['rv_sequence'] = rv_sequence
        primerpair['FW Overhang'] = fw_oh
        primerpair['FW Template'] = fw_template
        primerpair['RV Overhang'] = rv_oh
        primerpair['RV Template'] = rv_template
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
        overlap = self.get_overlap(primerpair['FW Primer'], primerpair['RV Primer'])
        if len(overlap) > self.complementarity_threshold:
            print(f"Complementarity between the IVA primers for {primerpair['Eblock']} exceeds threshold of {self.complementarity_threshold}")

    # def create_primers_fname(self, fp, extension="_primers.txt"):
    #     fname = extract_filename(fp)
    #     fname_ext = fname + extension
    #     return fname_ext

    # def write_primers_to_file(self, df, fp, snapgene_file):
    #     """
    #     Write primers to TXT file that can be imported in snapgene

    #     Args:
    #         df (pd.DataFrame): _description_
    #         fp (str): _description_
    #     """    
    #     fname = self.create_primers_fname(snapgene_file)
    #     with open(os.path.join(fp, fname), 'w+') as out:
    #         for _, row in df.iterrows():
    #             name = row['Gene Block']
    #             fw_primer = row['Forward primer (5>3)']
    #             rv_primer = row['Reverse primer (5>3)']
    #             out.write(self.short_block_name(name) + '_Fw' + ';' + str(fw_primer) + '\n')
    #             out.write(self.short_block_name(name) + '_Rv' + ';' + str(rv_primer) + '\n')


    # def optimize_primerpair_template(self, primerpair: dict):
    #     """
    #     Optimize the primer pairs to make sure they are within the desired parameters
    #     """
    #     dt_tm = abs(primerpair['Tm FW Template'] - primerpair['Tm Rv Template'])
    #     lowest_dt_tm = dt_tm
    #     highest_tm = max(primerpair['Tm FW Template'], primerpair['Tm Rv Template'])
    #     while highest_tm < ((self.max_template_temp_IVA) + 0.5):
    #         if primerpair['Tm FW Template'] < primerpair['Tm Rv Template']:
    #             toprocess = primerpair['FW Template']
    #             direction = 'FW'
    #         else:
    #             toprocess = primerpair['RV Template']
    #             direction = 'RV'
    #         # Add 1 nucleotide to the primer with lowest Tm
    #         if direction == 'FW':
    #             primerpair['FW Template'] = self.IVA_Fw_template(primerpair['end position'], primerpair['fw_sequence'], size=len(toprocess)+1)
    #             primerpair['Tm FW Template'] = self.Tm(primerpair['FW Template'])
    #         elif direction == 'RV':
    #             primerpair['RV Template'] = self.IVA_Rv_template(primerpair['begin position'], primerpair['rv_sequence'], size=len(toprocess)+1)
    #             primerpair['Tm RV Template'] = self.Tm(primerpair['RV Template'])
    #         highest_tm = max(primerpair['Tm FW Template'], primerpair['Tm Rv Template'])
    #         dt_tm = abs(primerpair['Tm FW Template'] - primerpair['Tm Rv Template'])
    #         if dt_tm < lowest_dt_tm:
    #             lowest_dt_tm = dt_tm
    #         else:
    #             if direction == 'FW':
    #                 primerpair['FW Template'] = self.IVA_Fw_template(primerpair['end position'], primerpair['fw_sequence'], size=len(toprocess))
    #                 primerpair['Tm FW Template'] = self.Tm(primerpair['FW Template'])
    #             elif direction == 'RV':
    #                 primerpair['RV Template'] = self.IVA_Rv_template(primerpair['begin position'], primerpair['rv_sequence'], size=len(toprocess))
    #                 primerpair['Tm RV Template'] = self.Tm(primerpair['RV Template'])
    #             break
    #         highest_tm = max(primerpair['Tm FW Template'], primerpair['Tm Rv Template'])
    #         if highest_tm > ((self.max_template_temp_IVA) + 0.5): # revert back to previous primer
    #             if direction == 'FW':
    #                 primerpair['FW Template'] = self.IVA_Fw_template(primerpair['end position'], primerpair['fw_sequence'], size=len(toprocess))
    #                 primerpair['Tm FW Template'] = self.Tm(primerpair['FW Template'])
    #             elif direction == 'RV':
    #                 primerpair['RV Template'] = self.IVA_Rv_template(primerpair['begin position'], primerpair['rv_sequence'], size=len(toprocess))
    #                 primerpair['Tm RV Template'] = self.Tm(primerpair['RV Template'])
    #             break
    #     return primerpair

    # def optimize_primerpair_overhang(self, primerpair: dict):
    #     """
    #     Optimize the primer pairs to make sure they are within the desired parameters
    #     """
    #     dt_oh = abs(primerpair['Tm FW Overhang'] - primerpair['Tm RV Overhang'])
    #     lowest_dt_oh = dt_oh
    #     highest_tm = max(primerpair['Tm FW Overhang'], primerpair['Tm RV Overhang'])
    #     while highest_tm < ((self.max_overhang_temp_IVA) + 0.5):
    #         if primerpair['Tm FW Overhang'] < primerpair['Tm RV Overhang']:
    #             toprocess = primerpair['FW Overhang']
    #             direction = 'FW'
    #         else:
    #             toprocess = primerpair['RV Overhang']
    #             direction = 'RV'
    #         # Add 1 nucleotide to the primer with lowest Tm
    #         if direction == 'FW':
    #             primerpair['FW Overhang'] = self.IVA_Fw_overhang(primerpair['end position'], primerpair['fw_sequence'], size=len(toprocess)+1)
    #             primerpair['Tm FW Overhang'] = self.Tm(primerpair['FW Overhang'])
    #         elif direction == 'RV':
    #             primerpair['RV Overhang'] = self.IVA_Rv_overhang(primerpair['begin position'], primerpair['rv_sequence'], size=len(toprocess)+1)
    #             primerpair['Tm RV Overhang'] = self.Tm(primerpair['RV Overhang'])
    #         highest_tm = max(primerpair['Tm FW Overhang'], primerpair['Tm RV Overhang'])
    #         dt_oh = abs(primerpair['Tm FW Overhang'] - primerpair['Tm RV Overhang'])
    #         if dt_oh < lowest_dt_oh:
    #             lowest_dt_oh = dt_oh
    #         # else: # revert back to previous primer
    #         #     print("revert")
    #         #     if direction == 'FW':
    #         #         primerpair['FW Overhang'] = self.IVA_Fw_overhang(primerpair['end position'], primerpair['fw_sequence'], size=len(toprocess))
    #         #         primerpair['Tm FW Overhang'] = self.Tm(primerpair['FW Overhang'])
    #         #     elif direction == 'RV':
    #         #         primerpair['RV Overhang'] = self.IVA_Rv_overhang(primerpair['begin position'], primerpair['rv_sequence'], size=len(toprocess))
    #         #         primerpair['Tm RV Overhang'] = self.Tm(primerpair['RV Overhang'])
    #         #     break
    #     return primerpair
