import os
import sys
import pickle
import difflib
import pandas as pd
from utils import extract_filename
from snapgene_output import short_name
from Bio.SeqUtils import MeltingTemp as mt
from design_gene_blocks import read_seq, gene_block_range

def reverse_complement(sequence):
    pairs = {"a": "t", "c":"g", "t":"a", "g":"c"}
    reverse = ""
    for nucleotide in sequence:
        rev_nucl = pairs[nucleotide]
        reverse += rev_nucl
    return reverse

def melting_temperature(sequence):
    return round(mt.Tm_NN(sequence), 2)

def invert_sequence(sequence):
    return sequence[::-1]

def df_to_csv(df, outpath, fname="IVA_primers.csv"):
    df.to_csv(os.path.join(outpath, fname))

def optimize_tm(optimum, primer, pos, size, sequence, nsteps=20):

    best_tm = melting_temperature(primer)
    best_diff = round(abs(optimum - best_tm), 2)

    for i in range(nsteps):

        tm = melting_temperature(primer)
        diff = round(abs(optimum - tm), 2)

        if diff < best_diff:
            best_tm = tm
            best_diff = diff

        if tm > optimum:
            size -= 1
            primer = IVA_Fw_overhang(pos, sequence, size)
        elif tm < optimum:
            size += 1
            primer = IVA_Fw_overhang(pos, sequence, size)

    return size

def IVA_Fw_overhang(block_end, fw_sequence, size=15):
    fw_oh = fw_sequence[block_end-size:block_end]
    return fw_oh
    
def IVA_Fw_template(block_end, fw_sequence, size=20):
    fw_template = fw_sequence[block_end:block_end+size]
    return fw_template

def IVA_Rv_template(block_begin, rv_sequence, size=20):
    rv_template = rv_sequence[block_begin-size:block_begin]
    return rv_template

def IVA_Rv_overhang(block_begin, rv_sequence, size=15):
    rv_oh = rv_sequence[block_begin:block_begin+size]
    return rv_oh

def combine_primers(overhang, template_binding):
    return overhang + template_binding
    
def extract_unique_gene_blocks(gene_blocks):
    unique_gene_blocks = []
    for _, value in gene_blocks.items():
        if not value[0] in unique_gene_blocks:
            unique_gene_blocks.append(value[0])
    return unique_gene_blocks

def store_output_in_df(primers):
    
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
    final_fw_oh_tm = [melting_temperature(i) for i in final_fw_oh]
    final_fw_template_tm = [melting_temperature(i) for i in final_fw_template]
    final_rv_oh_tm = [melting_temperature(i) for i in final_rv_oh]
    final_rv_template_tm = [melting_temperature(i) for i in final_rv_template]
    fw_combined_tm = [melting_temperature(i) for i in fw_combined]
    rv_combined_tm = [melting_temperature(i) for i in rv_combined]

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

def get_overlap(s1, s2):
    s = difflib.SequenceMatcher(None, s1, s2)
    pos_a, _, size = s.find_longest_match(0, len(s1), 0, len(s2)) 
    return s1[pos_a:pos_a+size]

def create_primers_fname(fp, extension="_primers.txt"):
    fname = extract_filename(fp)
    fname_ext = fname + extension
    return fname_ext

def write_primers_to_file(df, fp, snapgene_file):
    """
    Write primers to TXT file that can be imported in snapgene

    Args:
        df (pd.DataFrame): _description_
        fp (str): _description_
    """    
    fname = create_primers_fname(snapgene_file)
    with open(os.path.join(fp, fname), 'w+') as out:
        for _, row in df.iterrows():
            name = row['Gene Block']
            fw_primer = row['Forward primer (5>3)']
            rv_primer = row['Reverse primer (5>3)']
            out.write(short_name(name) + '_Fw' + ';' + str(fw_primer) + '\n')
            out.write(short_name(name) + '_Rv' + ';' + str(rv_primer) + '\n')

def check_complementarity_primers(df, threshold=4):
    # Primer pairs should not have complementary regions
    # TODO DOUBLE CHECK THIS FUNCTION
    for _, row in df.iterrows():
        s1 = row['Forward primer (5>3)']
        s2 = row['Reverse primer (5>3)']
        overlap = get_overlap(s1, s2)
        if len(overlap) > threshold:
            print(f"Complementarity between the primers for {row['Gene Block']} exceeds threshold of {threshold}")

def check_tms_within_bounds_of_each_other(df, threshold=4):
    # Primer pairs should have a Tm within 5°C of each other
    for _, row in df.iterrows():
        if row['dTm Overhangs'] > threshold:
            print(f"The overhang temperatures for Fw and Rv primer of {row['Gene Block']} exceed max Tm difference of {threshold} degrees")
        if row['dTm Templates'] > threshold:
            print(f"The template temperatures for Fw and Rv primer of {row['Gene Block']} exceed max Tm difference {threshold} degrees")

def load_pickle(fp):
    with open(fp, 'rb') as handle:
        obj = pickle.load(handle)
    return obj

def overhang_temp():
    return 50

def template_temp():
    return 60

def main(result_path, gene_path, output_path, snapgene_file):
    
    # Load input
    gene_blocks = load_pickle(result_path)
    fw_sequence = read_seq(gene_path)
    rv_sequence = reverse_complement(fw_sequence)

    unique_gene_blocks = extract_unique_gene_blocks(gene_blocks)

    primers = {}
    for gb in unique_gene_blocks:
        
        begin_pos, end_pos = gene_block_range(gb)
        
        # Design initial primers and optimize temperatures
        init_fw_oh = IVA_Fw_overhang(end_pos, fw_sequence)
        size = optimize_tm(overhang_temp(), init_fw_oh, end_pos, 15, fw_sequence)
        final_fw_oh = IVA_Fw_overhang(end_pos, fw_sequence, size=size)

        init_fw_template = IVA_Fw_template(end_pos, fw_sequence)
        size = optimize_tm(template_temp(), init_fw_template, end_pos, 15, fw_sequence)
        final_fw_template = IVA_Fw_template(end_pos, fw_sequence, size)

        init_rv_oh = IVA_Rv_overhang(begin_pos, rv_sequence)
        size = optimize_tm(overhang_temp(), init_rv_oh, begin_pos, 15, rv_sequence)
        final_rv_oh = IVA_Rv_overhang(begin_pos, rv_sequence, size)
        final_rv_oh = invert_sequence(final_rv_oh)

        init_rv_template = IVA_Rv_template(begin_pos, rv_sequence)
        size = optimize_tm(template_temp(), init_rv_template, begin_pos, 15, rv_sequence)
        final_rv_template = IVA_Rv_template(begin_pos, rv_sequence, size)
        final_rv_template = invert_sequence(final_rv_template)

        # Combine primers
        fw_combined = combine_primers(final_fw_oh, final_fw_template)
        rv_combined = combine_primers(final_rv_oh, final_rv_template)

        primers[gb] = [final_fw_oh, final_fw_template, final_rv_oh, final_rv_template, fw_combined, rv_combined]
        
    df = store_output_in_df(primers)

    # Check temperatures
    check_tms_within_bounds_of_each_other(df)

    # Check primer complementarity
    check_complementarity_primers(df)

    # Write to file
    df_to_csv(df, output_path)
    # Also write primers to file that snapgene can import
    if snapgene_file:
        write_primers_to_file(df, output_path, snapgene_file)
    
    print("Primers written to file")
    print("Make sure that primer binds nowhere else in sequence")