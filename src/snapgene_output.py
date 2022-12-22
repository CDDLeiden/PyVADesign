import os
import sys
import pandas as pd
from Bio import SeqIO
from design_IVA_primers import load_pickle, extract_unique_gene_blocks, reverse_complement, invert_sequence
from design_gene_blocks import gene_block_range

def read_snapgene_dna_file(fp):
    with open(fp, 'rb') as handle:
        for record in SeqIO.parse(handle, "snapgene"):
            return record

def len_sequence(sequence):
    return len(sequence)

def extract_filename(fp):
    base=os.path.basename(fp)
    fname = os.path.splitext(base)[0]
    return fname

def write_primers_to_file(df, fp):
    # SnapGene docs:
        # FOR PRIMER LIST:
        # Name <tab> sequence <tab> notes
        # Name;sequence;notes
        # Name;sequence;notes
    fname = create_primers_fname(fp)
    with open(os.path.join(fp, fname), 'w+') as out:
        for _, row in df.iterrows():
            name = row['Gene Block']
            fw_primer = row['Forward primer (5>3)']
            rv_primer = row['Reverse primer (5>3)']
            out.write(name+'_fw' + ';' + fw_primer + '\n')
            out.write(name+'_rv' + ';' + rv_primer + '\n')

def create_primers_fname(fp, extension="_primers.txt"):
    fname = extract_filename(fp)
    fname_ext = fname + extension
    return fname_ext

def create_gff3_fname(fp, extension="_features.gff3"):
    fname = extract_filename(fp)
    fname_ext = fname + extension
    return fname_ext

def find_index(string, substring,start=None, stop=None):
    if start and stop:
        idx = string.find(substring, start, stop)
    else:
        idx = string.find(substring)
    return idx

def write_gff3_header(length_sequence):
    result = ["##gff-version 3.2.1", f"##sequence-region myseq 1 {str(length_sequence)}"]
    return result

def write_features_to_gff3(outpath, 
                           length_sequence, 
                           gene_blocks,
                           fw_overhangs,
                           rv_overhangs,
                           fw_templates,
                           rv_templates):
    # TODO ADD MUTATIONS AS WELL?
    fname = create_gff3_fname(outpath)
    header_lines = write_gff3_header(length_sequence)
    # Write header to file
    with open(os.path.join(outpath, fname), 'w+') as f:
        for i in header_lines:
            f.write(i + '\n')
        # Write all geneblocks to file (grey)
        for key, value in gene_blocks.items():
            result = ['myseq', '.', 'gene', str(value[0]), str(value[1]), '.', '.', '.', f"Name={key};color=#939393"]
            f.write('\t'.join(result) + '\n')
        # Write all overhangs to file (pink)
        for key, value in fw_overhangs.items():
            result = ['myseq', '.', 'gene', str(value[0]), str(value[1]), '.', '.', '.', f"Name={key}_Fw_overhang;color=#FF00D4"]
            f.write('\t'.join(result) + '\n')
        for key, value in rv_overhangs.items():
            result = ['myseq', '.', 'gene', str(value[0]), str(value[1]), '.', '.', '.', f"Name={key}_Rv_overhang;color=#FF00D4"]
            f.write('\t'.join(result) + '\n')
        # Write all template binding region to file (green)
        for key, value in fw_templates.items():
            result = ['myseq', '.', 'gene', str(value[0]), str(value[1]), '.', '.', '.', f"Name={key}_Fw_template;color=#06FF92"]
            f.write('\t'.join(result) + '\n')
        for key, value in rv_templates.items():
            result = ['myseq', '.', 'gene', str(value[0]), str(value[1]), '.', '.', '.', f"Name={key}_Rv_template;color=#06FF92"]
            f.write('\t'.join(result) + '\n')

if __name__ == "__main__":

    # TODO NOT ONLY ADD A FEATURE LIST, ALSO ADD A PRIMER LIST FOR IMPORTING PRIMERS INTO SNAPGENE
    # TODO ALSO PUT MUTATIONS IN THE FEATURES FILE??
    # FIX OUTPUT FILENAME!!!! (NOT EXAMPLE )

    snapgene_file = r"C:\Users\Rosan\Documents\git\my_repositories\design_gene_blocks\example\example_data\snapgene_vector.dna"
    wt_blocks_file = r"C:\Users\Rosan\Documents\git\my_repositories\design_gene_blocks\example\example_output\wt_gene_blocks.npy"
    mut_gene_blocks_file = r"C:\Users\Rosan\Documents\git\my_repositories\design_gene_blocks\example\example_output\gene_blocks.npy"
    outpath = r"C:\Users\Rosan\Documents\git\my_repositories\design_gene_blocks\example\example_output"
    primers_path = r"C:\Users\Rosan\Documents\git\my_repositories\design_gene_blocks\example\example_output\IVA_primers.csv"

    record = read_snapgene_dna_file(snapgene_file)
    lentgh_seq = len_sequence(record.seq)

    # Read gene blocks and get the WT ones
    wt_gene_blocks = load_pickle(wt_blocks_file)
    mut_gene_blocks = load_pickle(mut_gene_blocks_file)
    unique_gene_blocks = extract_unique_gene_blocks(mut_gene_blocks)
    gene_blocks = {k: wt_gene_blocks[k] for k in unique_gene_blocks}
    
    # Write all overhangs to file (pink)
    df = pd.read_csv(primers_path)
    fw_ohs = df['Fw Overhang'].tolist()
    fw_overhangs = {}
    fw_templates = {}
    rv_overhangs = {}
    rv_templates = {}

    # Write gene blocks to file
    gene_block_indexes = {}
    for key, value in gene_blocks.items():
        
        begin_idx_gene_block_in_vector = find_index(str(record.seq), value)
        end_idx_gene_block_in_vector = begin_idx_gene_block_in_vector + len(value)
        gene_block_indexes[key] = [begin_idx_gene_block_in_vector, end_idx_gene_block_in_vector]

    # FW overhangs
    for idx, row in df.iterrows():
        name = row['Gene Block']
        overhang = row['Fw Overhang']
        [begin, end] = gene_block_indexes[name]
        idx = find_index(str(record.seq), overhang, begin, end)
        fw_overhangs[name] = [idx, idx+len(overhang)]

    # FW templates
    for idx, row in df.iterrows():
        name = row['Gene Block']
        template = row['Fw Template']
        [begin, end] = gene_block_indexes[name]
        idx = find_index(str(record.seq), template, begin, end+40)
        fw_templates[name] = [idx, idx+len(template)]

    # RV Overhang
    # TODO REVERSE TEMPLATE AND OVERHANG
    for idx, row in df.iterrows():
        name = row['Gene Block']
        template = row['Rv Overhang']
        template_invert = invert_sequence(template)
        rv_template = reverse_complement(template_invert)
        [begin, end] = gene_block_indexes[name]
        idx = find_index(str(record.seq), rv_template, begin, end)
        rv_overhangs[name] = [idx, idx+len(template)]

    # Rv template
    for idx, row in df.iterrows():
        name = row['Gene Block']
        template = row['Rv Template']
        template_invert = invert_sequence(template)
        rv_template = reverse_complement(template_invert)
        [begin, end] = gene_block_indexes[name]
        idx = find_index(str(record.seq), rv_template, begin-40, end)
        print(idx)
        rv_templates[name] = [idx, idx+len(template)]

    # Write primers to file
    write_primers_to_file(df, outpath)

    write_features_to_gff3(outpath, lentgh_seq, gene_block_indexes, fw_overhangs, rv_overhangs, fw_templates, rv_templates)


