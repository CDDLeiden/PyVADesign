import sys
from Bio import SeqIO

# TODO Add vector here as well (so that mutations in N and C terminal can be made)

class Plasmid:
    """
    This class contains functions to parse and translate DNA sequences and vectors.
    """

    def __init__(self):
        self.sequence: str = None  # Gene sequence
        self.seqid = None
        self.organism = None
        self.vector = None
 
    def parse_sequence(self, fp: str) -> str:
        """
        This function parses the sequence from a text file and checks the input.
        """
        sequence, seqid = self.read_single_fasta(fp)
        self.sequence = sequence
        self.seqid = seqid

        result = self.check_sequence(sequence)
        return result
    
    def parse_vector(self, fp: str):
        """
        This function parses the vector from a DNA file.
        """
        vector = self.read_snapgene_dna_file(fp)
        self.vector = vector
        # TODO Add some checks here for the vector

    @staticmethod
    def read_snapgene_dna_file(fp: str):
        with open(fp, 'rb') as handle:
            for record in SeqIO.parse(handle, "snapgene"):
                return record
        
    @staticmethod
    def read_single_fasta(fp: str) -> str:
        """
        This function reads a single fasta file and returns the sequence.
        """
        for num, record in enumerate(SeqIO.parse(fp, "fasta")):
            sequence = record.seq
            seqid = record.id
            if num > 0:
                print("Please provide a single sequence in FASTA format.")
                sys.exit()
        return sequence, seqid
    
    @staticmethod
    def check_sequence(sequence: str) -> str:
        """
        This function checks the input sequence.
        """
        if Plasmid.is_dna(sequence) and Plasmid.contains_start_stop_codon(sequence):
            return 1
        else:
            print("Please provide a DNA sequence with start and stop codons.")
            sys.exit()

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
        if sequence.upper().startswith("ATG") and sequence.upper().endswith(("TAA", "TAG", "TGA")):
            return True
        else:
            print("Please provide a sequence with start and stop codons")

    @staticmethod
    def translate_sequence(sequence: str) -> str:
        """
        Translate DNA sequence to protein sequence
        """
        return sequence.translate()
    
    @staticmethod
    def invert_sequence(sequence):
        """
        Invert sequence
        """
        return sequence[::-1]
    
    @staticmethod
    def reverse_complement(sequence):
        """
        Reverse complement sequence
        """
        pairs = {"a": "t", "c":"g", "t":"a", "g":"c"}
        reverse = ""
        for nucleotide in sequence:
            rev_nucl = pairs[nucleotide]
            reverse += rev_nucl
        return reverse
    
