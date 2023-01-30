import os
import sys
import difflib
import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt
from design_gene_blocks import DesignEblocks
from utils import load_pickle, extract_filename

class DesignPrimers:
    """
    """

    def __init__(self,
                 wt_gene_blocks_fp: str,
                 mut_gene_blocks_fp: str,
                 output_location: str,
                 input_gene_path: str,
                 snapgene_file: str):
    
        self.wt_gene_blocks_fp = wt_gene_blocks_fp
        self.mut_gene_blocks_fp = mut_gene_blocks_fp
        self.output_location = output_location
        self.input_gene_path = input_gene_path
        self.snapgene_file = snapgene_file

        self.overhang_temp = 50
        self.template_temp = 60

    def run(self):
        
        # Load input
        gene_blocks = load_pickle(self.wt_gene_blocks_fp)
        print(gene_blocks)
        fw_sequence = DesignEblocks.read_seq(self.input_gene_path)
        print(len(fw_sequence))
        sys.exit()
        rv_sequence = self.reverse_complement(fw_sequence)

        # Extract unique gene blocks
        # unique_gene_blocks = self.extract_unique_gene_blocks(gene_blocks)

        primers = {}
        for gb_name, gb in gene_blocks.items():
            
            begin_pos, end_pos = DesignEblocks.gene_block_range(gb_name)

            print(gb_name)
            
            # Design initial primers and optimize temperatures
            init_fw_oh = self.IVA_Fw_overhang(end_pos, fw_sequence)
            print(init_fw_oh)
            print("TEST")
            size = self.optimize_tm(self.overhang_temp, init_fw_oh, end_pos, 15, fw_sequence)
            final_fw_oh = self.IVA_Fw_overhang(end_pos, fw_sequence, size=size)

            init_fw_template = self.IVA_Fw_template(end_pos, fw_sequence)
            size = self.optimize_tm(self.template_temp, init_fw_template, end_pos, 15, fw_sequence)
            final_fw_template = self.IVA_Fw_template(end_pos, fw_sequence, size)

            init_rv_oh = self.IVA_Rv_overhang(begin_pos, rv_sequence)
            size = self.optimize_tm(self.overhang_temp, init_rv_oh, begin_pos, 15, rv_sequence)
            final_rv_oh = self.IVA_Rv_overhang(begin_pos, rv_sequence, size)
            final_rv_oh = self.invert_sequence(final_rv_oh)

            init_rv_template = self.IVA_Rv_template(begin_pos, rv_sequence)
            size = self.optimize_tm(self.template_temp, init_rv_template, begin_pos, 15, rv_sequence)
            final_rv_template = self.IVA_Rv_template(begin_pos, rv_sequence, size)
            final_rv_template = self.invert_sequence(final_rv_template)

            # Combine primers
            fw_combined = self.combine_primers(final_fw_oh, final_fw_template)
            rv_combined = self.combine_primers(final_rv_oh, final_rv_template)

            primers[gb] = [final_fw_oh, final_fw_template, final_rv_oh, final_rv_template, fw_combined, rv_combined]
            
        df = self.store_output_in_df(primers)

        # Check temperatures
        self.check_tms_within_bounds_of_each_other(df)

        # Check primer complementarity
        self.check_complementarity_primers(df)

        # Write to file
        self.df_to_csv(df, os.path.join(self.output_location))

        # Also write primers to file that snapgene can import
        if self.snapgene_file:
            self.write_primers_to_file(df, self.output_location, self.snapgene_file)
        
        print("Primers written to file")
        print("Make sure that primer binds nowhere else in sequence")

    def reverse_complement(self, sequence):
        pairs = {"a": "t", "c":"g", "t":"a", "g":"c"}
        reverse = ""
        for nucleotide in sequence:
            rev_nucl = pairs[nucleotide]
            reverse += rev_nucl
        return reverse

    def melting_temperature(self, sequence):
        return round(mt.Tm_NN(sequence), 2)

    def invert_sequence(self, sequence):
        return sequence[::-1]

    def df_to_csv(self, df, outpath, fname="IVA_primers.csv"):
        df.to_csv(os.path.join(outpath, fname))

    def optimize_tm(self, optimum, primer, pos, size, sequence, nsteps=20):

        best_tm = self.melting_temperature(primer)
        best_diff = round(abs(optimum - best_tm), 2)

        for i in range(nsteps):

            tm = self.melting_temperature(primer)
            diff = round(abs(optimum - tm), 2)

            if diff < best_diff:
                best_tm = tm
                best_diff = diff

            if tm > optimum:
                size -= 1
                primer = self.IVA_Fw_overhang(pos, sequence, size)
            elif tm < optimum:
                size += 1
                primer = self.IVA_Fw_overhang(pos, sequence, size)

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
        
    # def extract_unique_gene_blocks(self, gene_blocks):
    #     unique_gene_blocks = []
    #     for key, value in gene_blocks.items():
    #         print(key)
    #         if not key in unique_gene_blocks:
    #             unique_gene_blocks.append(value[0])
    #     return unique_gene_blocks

    def store_output_in_df(self, primers):
        
        gene_blocks = []
        final_fw_oh = []
        final_fw_oh_tm = []
        final_fw_template = []
        final_fw_template_tm = []
        final_rv_oh = []
        final_rv_oh_tm = []
        final_rv_template = []
        final_rv_template_tm = []
        fw_combined = []
        fw_combined_tm = []
        rv_combined = []
        rv_combined_tm = []

        for key, value in primers.items():
            gene_blocks.append(key)
            final_fw_oh.append(value[0])
            final_fw_template.append(value[1])
            final_rv_oh.append(value[2])
            final_rv_template.append(value[3])
            fw_combined.append(value[4])
            rv_combined.append(value[5])

        # Calculate melting temperatures
        final_fw_oh_tm = [self.melting_temperature(i) for i in final_fw_oh]
        final_fw_template_tm = [self.melting_temperature(i) for i in final_fw_template]
        final_rv_oh_tm = [self.melting_temperature(i) for i in final_rv_oh]
        final_rv_template_tm = [self.melting_temperature(i) for i in final_rv_template]
        fw_combined_tm = [self.melting_temperature(i) for i in fw_combined]
        rv_combined_tm = [self.melting_temperature(i) for i in rv_combined]

        # Calculate lengths
        final_fw_oh_len = [len(i) for i in final_fw_oh]
        final_fw_template_len = [len(i) for i in final_fw_template]
        final_rv_oh_len = [len(i) for i in final_rv_oh]
        final_rv_template_len = [len(i) for i in final_rv_template]
        fw_combined_len = [len(i) for i in fw_combined]
        rv_combined_len = [len(i) for i in rv_combined]

        # Store in data frame
        df = pd.DataFrame()
        df['Gene Block'] = gene_blocks
        df['Fw Overhang'] = final_fw_oh
        df['Fw Overhang Tm'] = final_fw_oh_tm
        df['Fw Overhang length'] = final_fw_oh_len
        df['Fw Template'] = final_fw_template
        df['Fw Template Tm'] = final_fw_template_tm
        df['Fw Template length'] = final_fw_template_len
        df['Rv Overhang'] = final_rv_oh
        df['Rv Overhang Tm'] = final_rv_oh_tm
        df['Rv Overhang length'] = final_rv_oh_len
        df['Rv Template'] = final_rv_template
        df['Rv Template Tm'] = final_rv_template_tm
        df['Rv Template length'] = final_rv_template_len
        df['Forward primer (5>3)'] = fw_combined
        df['Forward primer Tm'] = fw_combined_tm
        df['Forward primer length'] = fw_combined_len
        df['Reverse primer (5>3)'] = rv_combined
        df['Reverse primer Tm'] = rv_combined_tm
        df['Reverse primer length'] = rv_combined_len
        df['dTm Overhangs'] = abs(df['Fw Overhang Tm'] - df['Rv Overhang Tm']).apply(lambda x: round(x, 2))
        df['dTm Templates'] = abs(df['Fw Template Tm'] - df['Rv Template Tm']).apply(lambda x: round(x, 2))

        return df

    def get_overlap(self, s1, s2):
        s = difflib.SequenceMatcher(None, s1, s2)
        pos_a, _, size = s.find_longest_match(0, len(s1), 0, len(s2)) 
        return s1[pos_a:pos_a+size]

    def create_primers_fname(self, fp, extension="_primers.txt"):
        fname = extract_filename(fp)
        fname_ext = fname + extension
        return fname_ext

    def write_primers_to_file(self, df, fp, snapgene_file):
        """
        Write primers to TXT file that can be imported in snapgene

        Args:
            df (pd.DataFrame): _description_
            fp (str): _description_
        """    
        fname = self.create_primers_fname(snapgene_file)
        with open(os.path.join(fp, fname), 'w+') as out:
            for _, row in df.iterrows():
                name = row['Gene Block']
                fw_primer = row['Forward primer (5>3)']
                rv_primer = row['Reverse primer (5>3)']
                out.write(short_name(name) + '_Fw' + ';' + str(fw_primer) + '\n')
                out.write(short_name(name) + '_Rv' + ';' + str(rv_primer) + '\n')

    def check_complementarity_primers(self, df, threshold=4):
        # Primer pairs should not have complementary regions
        for _, row in df.iterrows():
            s1 = row['Forward primer (5>3)']
            s2 = row['Reverse primer (5>3)']
            overlap = self.get_overlap(s1, s2)
            if len(overlap) > threshold:
                print(f"Complementarity between the primers for {row['Gene Block']} exceeds threshold of {threshold}")

    def check_tms_within_bounds_of_each_other(self, df, threshold=4):
        # Primer pairs should have a Tm within 5°C of each other
        for _, row in df.iterrows():
            if row['dTm Overhangs'] > threshold:
                print(f"The overhang temperatures for Fw and Rv primer of {row['Gene Block']} exceed max Tm difference of {threshold} degrees")
            if row['dTm Templates'] > threshold:
                print(f"The template temperatures for Fw and Rv primer of {row['Gene Block']} exceed max Tm difference {threshold} degrees")