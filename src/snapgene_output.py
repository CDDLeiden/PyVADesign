import os
import sys
import pandas as pd
from Bio import SeqIO
from design_IVA_primers import DesignPrimers
from utils_old import extract_filename, load_pickle


class SnapGeneOutput:
    """
    Class to write output to files that can be opened in SnapGene
    """

    def __init__(self,
                 wt_gene_blocks_fp,
                 mut_gene_blocks_fp,
                 primers_fp,
                 output_location,
                 snapgene_file,
                 gene_blocks_info_fp):
        
        self.wt_gene_blocks_fp = wt_gene_blocks_fp
        self.mut_gene_blocks_fp = mut_gene_blocks_fp
        self.primers_fp = primers_fp
        self.output_location = output_location
        self.snapgene_file = snapgene_file
        self.gene_blocks_info_fp = gene_blocks_info_fp

        self.type_mutation = "Mutation"
        self.type_insert = "Insert"
        self.type_deletion = "Deletion"
        self.type_combined = "Combined"

        self.color_mutation = '#FF0000'
        self.color_combination = '#D8FF00'
        self.color_insert = '#0017FF'
        self.color_deletion = '#FF5900'
        self.color_gene_block = '#939393'
        # self.color_primer = ''
        self.color_primer_OH = '#FF00D4'
        self.color_primer_template = '#06FF92'

    # TODO Move all +1/+2 to the designIVA class

    def run(self):
        
        # Read snapgene vector and extract length (=needed for gff3 file)
        record = self.read_snapgene_dna_file(self.snapgene_file)
        lentgh_seq = len(record.seq)

        # Read gene blocks and get the WT ones
        wt_gene_blocks = load_pickle(self.wt_gene_blocks_fp)
        mut_gene_blocks = load_pickle(self.mut_gene_blocks_fp)
        unique_gene_blocks = list(set([i[0] for i in mut_gene_blocks.values()]))
        gene_blocks = {k: wt_gene_blocks[k] for k in unique_gene_blocks}  # Only gene blocks that contain mutations
        
        # Start and end position of the gene block in the vector (=needed to create features later)
        gene_block_indexes = self.gene_block_indexes_in_vector(record.seq, gene_blocks)
        print(gene_block_indexes)

        # Read primer data
        df = pd.read_csv(self.primers_fp)
        fw_overhangs = self.extract_primer_index(df, 'Fw Overhang', str(record.seq), gene_block_indexes)
        fw_templates = self.extract_primer_index(df, 'Fw Template', str(record.seq), gene_block_indexes)
        rv_overhangs = self.extract_primer_index(df, 'Rv Overhang', str(record.seq), gene_block_indexes)
        rv_templates = self.extract_primer_index(df, 'Rv Template', str(record.seq), gene_block_indexes)

        mutations = self.write_mutations_to_gff3(self.gene_blocks_info_fp, gene_block_indexes)

        # Write gene blocks and overhang/template regions to snapgene output files
        self.write_features_to_gff3(self.output_location, 
                                    lentgh_seq,
                                    self.snapgene_file,
                                    gene_block_indexes, 
                                    fw_overhangs, 
                                    rv_overhangs, 
                                    fw_templates, 
                                    rv_templates,
                                    mutations)


    def create_gff3_fname(self, fp, extension="_features.gff3"):
        fname = extract_filename(fp)
        fname_ext = fname + extension
        return fname_ext

    def find_index(self, string, substring,start=None, stop=None):
        if start and stop:
            idx = string.find(substring, start, stop)
        else:
            idx = string.find(substring)
        return idx

    def gff3_header(self, length_sequence):
        result = ["##gff-version 3.2.1", f"##sequence-region myseq 1 {str(length_sequence)}"]
        return result

    def gff3_line(self, begin_pos, end_pos, name, hex_color):
        line = ['myseq', '.', 'gene', str(begin_pos), str(end_pos), '.', '.', '.', f"Name={name};color={hex_color}"]
        return line

    def write_features_to_gff3(self, 
                               outpath, 
                               length_sequence,
                               snapgene_file, 
                               gene_blocks,
                               fw_overhangs,
                               rv_overhangs,
                               fw_templates,
                               rv_templates,
                               mutations):
        
        fname = self.create_gff3_fname(snapgene_file)
        header_lines = self.gff3_header(length_sequence)
        
        # Write header to file
        with open(os.path.join(outpath, fname), 'w+') as f:
            for i in header_lines:
                f.write(i + '\n')
            
            # Write all geneblocks to file (color = grey)
            for key, value in gene_blocks.items():
                result = self.gff3_line(value[0], value[1], DesignPrimers.short_block_name(key), '#939393')
                f.write('\t'.join(result) + '\n')
            
            # Write all overhangs to file (color = pink)
            for key, value in fw_overhangs.items():
                result = self.gff3_line(value[0], value[1], DesignPrimers.short_block_name(key) + '_Fw_OH', '#FF00D4')
                f.write('\t'.join(result) + '\n')
            for key, value in rv_overhangs.items():
                result = self.gff3_line(value[0], value[1], DesignPrimers.short_block_name(key) + "_Rv_OH", '#FF00D4')
                f.write('\t'.join(result) + '\n')
            
            # Write all template binding region to file (color = green)
            for key, value in fw_templates.items():
                result = self.gff3_line(value[0], value[1], DesignPrimers.short_block_name(key) + "_Fw_Template", '#06FF92')
                f.write('\t'.join(result) + '\n')
            for key, value in rv_templates.items():
                result = self.gff3_line(value[0], value[1], DesignPrimers.short_block_name(key) + "_Rv_Template", '#06FF92')
                f.write('\t'.join(result) + '\n')

            for key, value in mutations.items():
                result = self.gff3_line(value[0], value[1], key, value[2])
                f.write('\t'.join(result) + '\n')

    def write_mutations_to_gff3(self, fp, gene_block_indexes):
        df = pd.read_csv(fp, sep='\t', header=0)
        results = {}  # d[mutation] = [begin, end, color]
        for idx, row in df.iterrows():
            for key, value in gene_block_indexes.items():
                if row['gene block name'] == key:
                    begin = gene_block_indexes[key][0] + row['index mutation'] -3
                    end = gene_block_indexes[key][0] + row['index mutation'] -1
                    print(begin, end)
                    type = row['type']
                    color = self.mutation_colors(type)
                    results[row['mutation']] = [begin, end, color]
                    continue
        return results

    def gene_block_indexes_in_vector(self, vector_sequence, gene_blocks):
        gene_block_indexes = {}  # d[gene_block_name] = [begin position in vector, end position in vector]
        for key, value in gene_blocks.items():
            begin_idx = self.find_index(str(vector_sequence), value)  # Position in vector, NOT in target gene
            end_idx = begin_idx + len(value)
            if (begin_idx == -1) or (end_idx == -1):
                print(f"Gene block {key} not found in vector sequence. Check whether your target gene is correct in your vector.")
                sys.exit()
            gene_block_indexes[key] = [begin_idx + 1, end_idx + 1]
        return gene_block_indexes

    def extract_primer_index(self, df, column, sequence, gene_block_indexes):
        primer_indexes = {}
        for idx, row in df.iterrows():
            name = row['Gene Block']
            primer = row[str(column)]
            if 'Rv' in column:
                primer = DesignPrimers.invert_sequence(primer)
                primer= DesignPrimers.reverse_complement(primer)
            [begin, end] = gene_block_indexes[name]
            idx = self.find_index(str(sequence), primer, begin-40, end+40) # Rv and Fw template binding region are outside of gene block range
            if not idx == -1:  # Index not found
                primer_indexes[name] = [idx, idx+len(primer)]
            else:
                print(f"{column}: {primer} not found in vector sequence")
                print(primer)
        return primer_indexes
    
    def mutation_colors(self, mutation_type):
        if mutation_type == self.type_mutation:
            return self.color_mutation
        elif mutation_type == self.type_insert:
            return self.color_insert
        elif mutation_type == self.type_deletion:
            return self.color_deletion
        elif mutation_type == self.type_combined:
            return self.color_combination
        
    @staticmethod
    def read_snapgene_dna_file(fp):
        with open(fp, 'rb') as handle:
            for record in SeqIO.parse(handle, "snapgene"):
                return record
        