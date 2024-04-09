import os
import sys
import pandas as pd

from .sequence import Plasmid


class Utils:

    def __init__(self):
            pass

    @staticmethod
    def natural_amino_acids():
        """
        This function returns a string of valid amino acids.
        """
        return "acdefghiklmnpqrstvwy"
    
    @staticmethod
    def DNA_codons():
        DNA_Codons = {
        "ATG": "start",
        "TAA": "stop", "TAG": "stop", "TGA": "stop",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "TGT": "C", "TGC": "C",
        "GAT": "D", "GAC": "D",
        "GAA": "E", "GAG": "E",
        "TTT": "F", "TTC": "F",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        "CAT": "H", "CAC": "H",
        "ATA": "I", "ATT": "I", "ATC": "I",
        "AAA": "K", "AAG": "K",
        "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "ATG": "M",
        "AAT": "N", "AAC": "N",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "TGG": "W",
        "TAT": "Y", "TAC": "Y"}
        return DNA_Codons
    
    @staticmethod
    def read_codon_usage(fp):
        codon_usage = {}
        codon_usage_no_U = {}
        df = pd.read_csv(fp, sep=';')
        for _, row in df.iterrows():
            codon_usage[row['Triplet']] = row['Number']
        for key, value in codon_usage.items():
            newkey = key.replace('U', 'T').lower()
            codon_usage_no_U[newkey] = value
        return codon_usage_no_U
    

class SnapGene:
    """
    Class for generating SnapGene files
    """

    def __init__(self,
                 sequence_instance: Plasmid,
                 output_dir: str = None):
            
            self.primer_filename = 'primers.fasta'
            self.eblocks_filename = 'eBlocks.gff3'

            self.output_dir = output_dir
            self.sequence_instance = sequence_instance
            
    def primers_to_fasta(self, primers: dict):
        """
        This function converts the primers to features that can be read by SnapGene.
        """
        self.make_file(self.primer_filename, header=False)
        with open(os.path.join(self.output_dir, self.primer_filename), 'a') as f:
            for k, v in primers.items():
                f.write(f">{k}\n")
                f.write(f"{v}\n")

    def eblocks_to_gff3(self, eblocks: dict):
        """
        This function converts the eBlocks to features that can be read by SnapGene.
        """
        self.make_file(self.eblocks_filename, header=True)
        with open(os.path.join(self.output_dir, self.eblocks_filename), 'a') as f:
            for k, v in eblocks.items():
                line = SnapGene.gff3_line(v[0], v[1], k, v[2])
                f.write('\t'.join(line) + '\n')

    def make_file(self, filename, header=False):
        try:
            with open(os.path.join(self.output_dir, filename), 'r') as f:
                pass
        except FileNotFoundError:
            with open(os.path.join(self.output_dir, filename), 'w') as f:
                if header:
                    f.write("\n".join(SnapGene.gff3_header(len(self.sequence_instance.vector.seq)+"\n")))
                else:
                    pass

    def remove_empty_lines(self):
        # TODO Remove last empty line from the file
        pass

    @staticmethod
    def gff3_header(length_sequence, version="3.2.1", sequence_name="myseq"):
        result = [f"##gff-version {version}", f"##sequence-region {sequence_name} 1 {len(length_sequence)}"]
        return result

    @staticmethod
    def gff3_line(begin_pos, end_pos, name, hex_color, type='gene'):
        # TODO Change feature, gene, CDS etc to the correct type
        # TODO Change Myseq to the correct sequence name?
        line = ['myseq', '.', f"{type}", str(begin_pos), str(end_pos), '.', '.', '.', f"Name={name};color={hex_color}"]
        return line
    
    @staticmethod
    def gff3_colors():
        colors = {
            'mutation': '#FF0000',
            'combination': '#D8FF00',
            'insert': '#0017FF',
            'deletion': '#FF5900',
            'eBlock': '#939393',
            'IVAprimer': '#FF00D4',
            'SEQprimer': '#06FF92'}
        return colors
    
    @staticmethod
    def count_substring_occurance(sequence: str, substring: str):
        return str(sequence).count(str(substring))
    

        # TODO Remove last empty line from the file
    

    # TODO Check if file exists, otherwise create it and append header gff3 file 
    
    # @staticmethod
    # def write_gff3_line(file, line):
    #     file.write('\t'.join(line) + '\n')

    # @staticmethod
    # TDOD remove last empty line in file