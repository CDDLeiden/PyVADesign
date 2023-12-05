import sys
from Bio import SeqIO


class Sequence:
    @staticmethod
    def read_single_fasta(fp: str) -> str:
        """
        This function reads a single fasta file and returns the sequence.
        """
        for num, record in enumerate(SeqIO.parse(fp, "fasta")):
            sequence = record.seq
            if num > 0:
                print("Please provide a single sequence in FASTA format.")
                sys.exit()
        return sequence
    
    @staticmethod
    def check_sequence(sequence: str) -> str:
        """
        This function checks the input sequence.
        """
        if Sequence.is_dna(sequence) and Sequence.contains_start_stop_codon(sequence):
            return sequence

    @staticmethod
    def is_dna(sequence: str) -> bool:
        """
        This function checks if the sequence is DNA.
        """
        if set(sequence.upper()).issubset("ATGC"):
            return True
        else:
            print("Please provide a DNA sequence")

    @staticmethod
    def contains_start_stop_codon(sequence: str) -> bool:
        """
        This function checks if the sequence contains start and stop codons.
        """
        if sequence.startswith("ATG") and sequence.endswith(("TAA", "TAG", "TGA")):
            return True
        else:
            print("Please provide a sequence with start and stop codons")

    @staticmethod
    def translate_sequence(sequence: str) -> str:
        """
        Translate DNA sequence to protein sequence
        """
        return sequence.translate()