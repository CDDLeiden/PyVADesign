import os
import sys
import pickle
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from design_gene_blocks import read_seq, gene_block_range

def reverse_complement(sequence):
    return sequence.reverse_complement()

def melting_temperature(sequence):
    return round(mt.Tm_NN(sequence), 2)

def invert_sequence(sequence):
    return sequence[::-1]

def load_pickle(fp):
    with open(fp, 'rb') as handle:
        obj = pickle.load(handle)
    return obj

def initial_Fw_IVA_primers(end_block, dna_seq_fw):
    # Initial fw primer
    primer_fw_oh = dna_seq_fw[end_block-(initial_oh_length()):end_block]
    primer_fw_template_binding = dna_seq_fw[end_block:end_block+(initial_template_binding_length())]
    return primer_fw_oh, primer_fw_template_binding

def initial_Rv_IVA_primers(begin_block, dna_seq_rv):
    # Initial Rv primer
    primer_rv_oh = dna_seq_rv[begin_block-(initial_oh_length()):begin_block]
    primer_rv_template_binding = dna_seq_rv[begin_block-(initial_template_binding_length() + initial_oh_length()):begin_block-(initial_oh_length())]
    return primer_rv_oh, primer_rv_template_binding

def combine_fw_primers(primer_oh, primer_template):
    result = str(primer_oh) + str(primer_template)  # 5>3
    return result

def combine_rv_primers(primer_oh, primer_template):
    result = str(primer_template) + str(primer_oh)  # 3>5
    result = invert_sequence(result)  # 5>3
    return result
        
def initial_oh_length():
    return 15

def initial_template_binding_length():
    return 20

def extract_unique_gene_blocks(gene_blocks):
    unique_gene_blocks = []
    for _, value in gene_blocks.items():
        if not value[0] in unique_gene_blocks:
            unique_gene_blocks.append(value[0])
    return unique_gene_blocks

if __name__ == "__main__":

    # TODO Design primers for insertion
    # TODO Design primers for deletion of gene fragment

    # Length of 18-24 bases
    # 40-60% G/C content
    # Start and end with 1-2 G/C pairs
    # Melting temperature (Tm) of 50-60°C
    # Primer pairs should have a Tm within 5°C of each other
    # Primer pairs should not have complementary regions

    # The rules set for finding primers are (a) the Tm value of sequences flanking the mutagenic codon is greater than 45; (b) The GC content of the primer is between 40 and 60%. 
    # https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html


    result_path = r"example_output\gene_blocks.npy"
    gene_path = r"example_data\mtb_DnaE1_seq.txt"
    output_path = r"example_output"
    
    gene_blocks = load_pickle(result_path)
    dna_seq_53 = read_seq(gene_path)
    dna_seq_35 = reverse_complement(dna_seq_53)

    unique_gene_blocks = extract_unique_gene_blocks(gene_blocks)

    # First make primers to delete the gene block
    with open(os.path.join(output_path, 'IVA_deletion_primers.txt'), 'w+') as out:
        header = ['gene_block', 'fw_oh', 'fw_oh_tm', 'fw_template_binding', 'fw_template_binding_tm', 'rv_oh', 'rv_oh_tm', 'rv_teplate_binding', 'rv_teplate_binding_tm', 'fw_primer', 'rv_primer']
        out.write('\t'.join(header) + '\n')
        for gb in unique_gene_blocks:
            
            print(gb)

            begin_pos, end_pos = gene_block_range(gb)
            
            # Design initial primers
            primer_fw_oh, primer_fw_template_binding = initial_Fw_IVA_primers(end_pos, dna_seq_53)
            primer_rv_oh, primer_rv_template_binding = initial_Rv_IVA_primers(begin_pos, dna_seq_35)  

            # Combine primers
            fw_deletion_combined = combine_fw_primers(primer_fw_oh, primer_fw_template_binding)
            rv_deletion_combined = combine_rv_primers(primer_rv_oh, primer_rv_template_binding)

            # Write to file
            # TODO ADD LENGTH OF PRIMERS TO OUTPUT FILE
            # TODO ADD EXPECTED OUTPUT SIZE TO FILE
            out.write(gb + '\t' + str(primer_fw_oh) + '\t' + str(melting_temperature(primer_fw_oh)) + '\t' + str(primer_fw_template_binding) + '\t' + str(melting_temperature(primer_fw_template_binding)) \
                          + '\t' + str(primer_rv_oh) + '\t' + str(melting_temperature(primer_rv_oh)) + '\t' + str(primer_rv_template_binding) + '\t' + str(melting_temperature(primer_rv_template_binding)) \
                          + '\t' + str(fw_deletion_combined) + '\t' + str(rv_deletion_combined) + '\n')

    # Next create the primers for subcloning consisting of a primer to amplify the insert and a primer to amplify the vector



            
        