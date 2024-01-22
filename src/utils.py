import sys
import pandas as pd


class Utils:
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