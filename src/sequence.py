import sys
from Bio import SeqIO


class Plasmid:
    """
    This class contains functions to parse and translate DNA sequences and vectors.
    """

    def __init__(self):
        self.sequence: str = None  # Gene sequence
        self.seqid: str = None
        self.organism: str = None
        self.vector: str = None  # Vector sequence with gene cloned into it
        self.color: str = "#d3d3d3"
 
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
    def find_index_in_vector(vector, sequence: str):
        idx_begin = str(vector).lower().find(str(sequence).lower())
        idx_end = idx_begin + len(sequence)
        if idx_begin == -1 or idx_end == -1:
            print(f"Gene block {sequence} not found in vector sequence. Check whether your target gene is correct in your vector.")
            sys.exit()
        return idx_begin, idx_end
        
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
        This function checks whether the sequence is DNA (only contains GCTA) and contains a start and stop codon.
        """
        if Plasmid.is_dna(sequence) and Plasmid.contains_start_stop_codon(sequence):
            return 1
        else:
            print("Please provide a DNA sequence with start and stop codons.")
            sys.exit()

    @staticmethod
    def check_vector(vector: str) -> str:
        """
        This function checks whether the vector is DNA (only contains GCTA) and contains the gene sequence.
        """
        # TODO Make this function
        pass

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