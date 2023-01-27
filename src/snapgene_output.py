import os
import sys
import pandas as pd
from Bio import SeqIO
from utils import extract_filename
from design_gene_blocks import DesignEblocks
from design_IVA_primers import load_pickle, extract_unique_gene_blocks, reverse_complement, invert_sequence

def read_snapgene_dna_file(fp):
    with open(fp, 'rb') as handle:
        for record in SeqIO.parse(handle, "snapgene"):
            return record

def create_gff3_fname(fp, extension="_features.gff3"):
    fname = extract_filename(fp)
    fname_ext = fname + extension
    return fname_ext

def len_sequence(sequence):
    return len(sequence)

def find_index(string, substring,start=None, stop=None):
    if start and stop:
        idx = string.find(substring, start, stop)
    else:
        idx = string.find(substring)
    return idx

def gff3_header(length_sequence):
    result = ["##gff-version 3.2.1", f"##sequence-region myseq 1 {str(length_sequence)}"]
    return result

def gff3_line(begin_pos, end_pos, name, hex_color):
    line = ['myseq', '.', 'gene', str(begin_pos), str(end_pos), '.', '.', '.', f"Name={name};color={hex_color}"]
    return line

def write_features_to_gff3(outpath, 
                           length_sequence,
                           snapgene_file, 
                           gene_blocks,
                           fw_overhangs,
                           rv_overhangs,
                           fw_templates,
                           rv_templates,
                           mutations):
    
    fname = create_gff3_fname(snapgene_file)
    header_lines = gff3_header(length_sequence)
    
    # Write header to file
    with open(os.path.join(outpath, fname), 'w+') as f:
        for i in header_lines:
            f.write(i + '\n')
        
        # Write all geneblocks to file (color = grey)
        for key, value in gene_blocks.items():
            result = gff3_line(value[0], value[1], DesignEblocks.short_name(key), '#939393')
            f.write('\t'.join(result) + '\n')
        
        # Write all overhangs to file (color = pink)
        for key, value in fw_overhangs.items():
            result = gff3_line(value[0], value[1], DesignEblocks.short_name(key) + '_Fw_OH', '#FF00D4')
            f.write('\t'.join(result) + '\n')
        for key, value in rv_overhangs.items():
            result = gff3_line(value[0], value[1], DesignEblocks.short_name(key) + "_Rv_OH", '#FF00D4')
            f.write('\t'.join(result) + '\n')
        
        # Write all template binding region to file (color = green)
        for key, value in fw_templates.items():
            result = gff3_line(value[0], value[1], DesignEblocks.short_name(key) + "_Fw_Template", '#06FF92')
            f.write('\t'.join(result) + '\n')
        for key, value in rv_templates.items():
            result = gff3_line(value[0], value[1], DesignEblocks.short_name(key) + "_Rv_Template", '#06FF92')
            f.write('\t'.join(result) + '\n')

        for key, value in mutations.items():
            result = gff3_line(value[0], value[1], key, value[2])
            f.write('\t'.join(result) + '\n')

def write_mutations_to_gff3(fp, gene_block_indexes, hex_color="#FF0000"):
    df = pd.read_csv(fp, sep='\t', header=0)
    results = {}  # d[mutation] = [begin, end]
    for idx, row in df.iterrows():
        for key, value in gene_block_indexes.items():
            if row['gene block name'] == key:
                begin = gene_block_indexes[key][0] + row['index mutation']
                end = gene_block_indexes[key][0] + row['index mutation'] + 2
                results[row['mutation']] = [begin, end, hex_color]
                continue
    return results

def gene_block_indexes_in_vector(vector_sequence, gene_blocks):
    gene_block_indexes = {}  # d[gene_block_name] =[begin position in vector, end position in vector]
    for key, value in gene_blocks.items():
        begin_idx = find_index(str(vector_sequence), value)  # Position in vector, NOT in target gene
        end_idx = begin_idx + len(value)
        gene_block_indexes[key] = [begin_idx, end_idx]
    return gene_block_indexes

def extract_primer_index(df, column, sequence, gene_block_indexes):
    primer_indexes = {}
    for idx, row in df.iterrows():
        name = row['Gene Block']
        primer = row[str(column)]
        if 'Rv' in column:
            primer = invert_sequence(primer)
            primer= reverse_complement(primer)
        [begin, end] = gene_block_indexes[name]
        idx = find_index(str(sequence), primer, begin-40, end+40) # Rv and Fw template binding region are outside of gene block range
        if not idx == -1:  # Index not found
            primer_indexes[name] = [idx, idx+len(primer)]
        else:
            print(f"{column}: {primer} not found in vector sequence")
            print(primer)
    return primer_indexes

def main(wt_gene_blocks_fp, mut_gene_blocks_fp, primers_fp, output_location, snapgene_file, gene_blocks_info_fp):
    
    # Read snapgene vector and extract length (=needed for gff3 file)
    record = read_snapgene_dna_file(snapgene_file)
    lentgh_seq = len_sequence(record.seq)

    # Read gene blocks and get the WT ones
    wt_gene_blocks = load_pickle(wt_gene_blocks_fp)
    mut_gene_blocks = load_pickle(mut_gene_blocks_fp)
    unique_gene_blocks = extract_unique_gene_blocks(mut_gene_blocks)
    gene_blocks = {k: wt_gene_blocks[k] for k in unique_gene_blocks}  # Only gene blocks that contain mutations
    
    # Start and end position of the gene block in the vector (=needed to create features later)
    gene_block_indexes = gene_block_indexes_in_vector(record.seq, gene_blocks)

    # Read primer data
    df = pd.read_csv(primers_fp)
    fw_overhangs = extract_primer_index(df, 'Fw Overhang', str(record.seq), gene_block_indexes)
    fw_templates = extract_primer_index(df, 'Fw Template', str(record.seq), gene_block_indexes)
    rv_overhangs = extract_primer_index(df, 'Rv Overhang', str(record.seq), gene_block_indexes)
    rv_templates = extract_primer_index(df, 'Rv Template', str(record.seq), gene_block_indexes)

    mutations = write_mutations_to_gff3(gene_blocks_info_fp, gene_block_indexes)

    # Write gene blocks and overhang/template regions to snapgene output files
    write_features_to_gff3(output_location, 
                           lentgh_seq,
                           snapgene_file,
                           gene_block_indexes, 
                           fw_overhangs, 
                           rv_overhangs, 
                           fw_templates, 
                           rv_templates,
                           mutations)