import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from biotite.sequence import NucleotideSequence


# .dna file can also be opened in benchling

#TODO Nucleotide sequence class vs protein sequence class


class Plasmid:
    """
    This class contains functions to parse and translate DNA sequences and vectors.
    """

    def __init__(self,
                 sequence: str = None,  # Gene sequence # TODO Rename to Gene sequence?
                 seqid: str = None,
                 organism: str = None,
                 vector: str = None,  # Vector sequence with gene cloned into it
                 vector_id: str = None,
                 gene_start_idx: int = -1,
                 gene_end_idx: int = -1,
                 color: str = "#d3d3d3"):  # TODO Obtain color from mutation colors
        
        self.sequence = sequence  # GoI
        self.seqid = seqid
        self.organism = organism
        self.vector = vector
        self.vector_id = vector_id
        self.gene_start_idx = gene_start_idx
        self.gene_end_idx = gene_end_idx
        self.color = color

    def save_vector(self, vector, output_dir, filename):
        """
        This function saves the vector to a file.
        """
        with open(os.path.join(output_dir, filename), "w") as output_handle:
            output_handle.write(f">{self.vector_id}\n")
            for i in range(0, len(vector), 60):
                output_handle.write(str(vector[i:i+60]) + "\n")

    def parse_vector(self, fp: str):
        """
        This function parses the vector from a DNA file.
        """
        if fp.endswith(".dna"):
            vector = self.read_snapgene_dna_file(fp)
        elif fp.endswith(".gb"):
            vector = self.read_genbank_file(fp)
        else:
            print("Please provide a SnapGene or GenBank file.")
            sys.exit()
        print(vector.id)
        self.vector = vector
        self.vector_id = vector.id
        # print(self.vector_id)
        # TODO Add some checks here for the vector
 
    def parse_sequence(self, fp: str) -> str:
        """
        This function parses the sequence from a text file and checks the input.
        """
        sequence, seqid = self.read_single_fasta(fp)
        self.sequence = sequence
        self.seqid = seqid
        self.gene_start_idx, self.gene_end_idx = self.find_index_in_vector(vector=self.vector.seq, sequence=self.sequence)
        self.check_dna_sequence(sequence)
        
    def mutate_vector(self, idx_start, idx_end, sequence, mutation_type):
        """
        This function mutates the input vector with a given sequence
        """
        # TODO WRITE A FUNCTION TO CHECK THE LENGTH OF THE MUTATED VECTOR
        # TODO FUND OUT WHAT IS WRONG WITH THIS FUNCTION!!
        if (idx_start < idx_end) and (idx_start >= 0):
            mutated_vector = self.vector.seq[:idx_start] + sequence + self.vector.seq[idx_end:]
        elif (idx_start > idx_end):
            restoend = len(self.vector.seq) - idx_start
            mutated_vector = sequence[restoend:] + self.vector.seq[idx_end:idx_start] + sequence[:restoend]
        else:
            print("Error in mutation")  # TODO 
            sys.exit()

        if mutation_type == "Mutation" or mutation_type == "Combined":  # Check whether vector length is correct
            if len(mutated_vector) != len(self.vector.seq):
                print("Error in mutation")
                sys.exit()
        return mutated_vector
        
    @staticmethod
    def circular_index(index, sequence_length):
        """
        This function returns the circular (=positive) index in a sequence.
        """
        return (index + sequence_length) % sequence_length
    
    @staticmethod
    def slice_circular_sequence(sequence, start_index, end_index):
        start_index = Plasmid.circular_index(start_index, len(sequence))
        if end_index < start_index:
            sliced_sequence = sequence[start_index:] + sequence[:end_index]
        else:
            sliced_sequence = sequence[start_index:end_index]
        return sliced_sequence
      
    @staticmethod
    def find_index_in_vector(vector, sequence: str):
        idx_begin = str(vector).lower().find(str(sequence).lower())
        idx_end = idx_begin + len(sequence)
        idx_begin = Plasmid.circular_index(idx_begin, len(vector))
        idx_end = Plasmid.circular_index(idx_end, len(vector))
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
    def read_genbank_file(fp: str):
        with open(fp, 'r') as handle:
            for record in SeqIO.parse(handle, "genbank"):
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
    def check_dna_sequence(sequence: str) -> bool:
        """
        This function checks whether the sequence is an actual nucleotide sequence and contains a start and stop codon.
        """
        if not (NucleotideSequence(sequence).is_valid() and Plasmid.contains_start_stop_codon(sequence)):
            raise ValueError("Invalid nucleotide sequence: Please provide a valid sequence containing start and stop codons.")
        return True

    @staticmethod
    def contains_start_stop_codon2(sequence: str) -> bool:
        """
        This function checks if the sequence contains start and stop codons.
        """
        stop_codons = ["TAA", "TAG", "TGA"]
        sequence = sequence.upper()
        if (sequence.startswith("ATG")) and (sequence[-3:] in stop_codons):
            return True
        else:
            return False

    @staticmethod
    def check_vector(vector: str) -> str:
        """
        This function checks whether the vector is DNA (only contains GCTA) and contains the gene sequence.
        """
        # TODO Make this function
        pass
    
    # @staticmethod
    # def invert_sequence(sequence):
    #     return sequence[::-1]