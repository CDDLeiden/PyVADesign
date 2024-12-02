#!/usr/bin/env python3

import os
import sys
import math
from Bio import SeqIO
import biotite.sequence as seq
from biotite.sequence import NucleotideSequence



class Gene:
    """
    This class represents a gene with a given sequence.

    Attributes:
    ----------
    sequence : str
        The nucleotide sequence of the gene.
    seqid : str
        The sequence identifier.
    stopcodon : bool
        Whether the sequence contains a stop codon.
    residues : dict
        The residues of the gene with the format: resnum: (residue, codon)
    """
    def __init__(self,
                 sequence: str = None,
                 seqid: str = None,
                 stopcodon: bool = True):
        
        self.sequence = sequence
        self.seqid = seqid
        self.stopcodon = stopcodon
        self.residues: dict = {} # Format: resnum: (residue, codon)

    def parse_sequence(self, fp: str) -> str:
        """This function parses the sequence from a .fasta file and checks the input."""
        self.sequence, self.seqid = self.read_single_fasta(fp)
        if self.stopcodon:
            if not (NucleotideSequence(self.sequence).is_valid() and self.contains_start_stop_codon(self.sequence)):
                raise ValueError("Invalid nucleotide sequence: Please provide a valid sequence containing start and stop codons.")
        else:
            if not NucleotideSequence(self.sequence).is_valid():
                raise ValueError("Invalid nucleotide sequence: Please provide a valid sequence.")
        self.get_residues()
            
    def get_residues(self):
        """This function returns the residues of the gene."""
        count = 1
        for i in range(0, math.ceil(len(self.sequence)), 3):
            codon = self.sequence[i:i+3]
            residue = seq.CodonTable.default_table()[str(codon).upper()]
            self.residues[count] = (residue, codon)
            count += 1
        return self.residues
        
    @staticmethod
    def read_single_fasta(fp: str):
        """This function reads a single fasta file and returns the sequence."""
        record = next(SeqIO.parse(fp, "fasta"))
        return record.seq, record.id
            
    @staticmethod
    def contains_start_stop_codon(sequence: str) -> bool:
        """This function checks if the sequence contains start and stop codons."""
        stop_codons = ["TAA", "TAG", "TGA"]
        sequence = sequence.upper()
        if (sequence.startswith("ATG")) and (sequence[-3:] in stop_codons):
            return True
        else:
            return False



class Vector:
    """
    This class represents a vector with a given sequence.

    Attributes:
    ----------
    gene : Gene
        The gene that is inserted into the vector.
    vector : SeqRecord
        The vector sequence.
    vector_id : str
        The vector identifier.
    organism : str
        The organism of the vector.
    gene_start_idx : int
        The start index of the gene in the vector.
    gene_end_idx : int
        The end index of the gene in the vector.
    length : int
        The length of the vector.
    color : str
        The color of the vector for visualization.
    """
    def __init__(self,
                 gene: Gene,
                 vector = None,
                 vector_id = None,
                 organism = None,
                 gene_start_idx = -1,
                 gene_end_idx = -1,
                 length = None,
                 color: str = "#d3d3d3"):

        self.vector = vector
        self.gene = gene
        self.vector_id = vector_id
        self.organism = organism
        self.gene_start_idx = gene_start_idx
        self.gene_end_idx = gene_end_idx
        self.color = color
        self.length = length

    def parse_vector(self, fp: str):
        """This function parses the vector from a DNA file in the SnapGene or GenBank format."""
        if fp.endswith(".dna"):
            vector = self.read_snapgene_dna_file(fp)
        elif fp.endswith(".gb"):
            vector = self.read_genbank_file(fp)
        else:
            raise ValueError("Please provide a SnapGene or GenBank file.")
        self.vector = vector
        self.vector_id = vector.id
        self.gene_start_idx, self.gene_end_idx = self.find_index_in_vector(self.gene.sequence)
        self.length = len(self.vector.seq)

    def find_index_in_vector(self, sequence: str):
        idx_begin = str(self.vector.seq).lower().find(str(sequence).lower())
        if idx_begin == -1:  # sequence overlapping with begin/end of vector
            new_vector = self.vector.seq + self.vector.seq
            idx = str(new_vector).lower().find(str(sequence).lower())
            idx_begin = idx - len(self.vector.seq)
            idx_begin = self.circular_index(idx_begin, len(self.vector))
        idx_end = idx_begin + len(sequence)
        idx_end = self.circular_index(idx_end, len(self.vector))
        if idx_end > len(self.vector):
            idx_end = idx_end - len(self.vector)
            idx_end = self.circular_index(idx_end, len(self.vector))
        if idx_begin == -1 or idx_end == -1:
            raise ValueError("Gene block not found in vector sequence.")
        return idx_begin, idx_end

    def save_vector(self, vector, output_dir, filename):
        """This function saves the vector to a file."""
        with open(os.path.join(output_dir, filename), "w") as output_handle:
            output_handle.write(f">{self.vector_id}\n")
            for i in range(0, len(vector), 60):
                output_handle.write(str(vector[i:i+60]) + "\n")

    def mutate_vector(self, idx_start, idx_end, sequence, mutation):
        """This function mutates the input vector with a given sequence."""
        if (idx_start < idx_end) and (idx_start >= 0):
            mutated_vector = self.vector.seq[:idx_start] + sequence + self.vector.seq[idx_end:]
        elif (idx_start > idx_end):
            restoend = len(self.vector.seq) - idx_start
            mutated_vector = sequence[restoend:] + self.vector.seq[idx_end:idx_start] + sequence[:restoend]
        else:
            raise ValueError("Invalid start and end indices.")
        self.check_length_mutated_vector(mutated_vector, mutation)
        return mutated_vector
    
    def check_length_mutated_vector(self, mutated_vector, mutation):
        """This function checks whether the mutated vector has the correct length."""
        if mutation.type == "Mutation" or mutation.type == "Combined":
            if len(mutated_vector) != len(self.vector.seq):
                raise ValueError("The mutated vector has a different length than the original vector.")
        elif mutation.type == "Insertion":
            if len(mutated_vector) != len(self.vector.seq) + mutation.length_insert:
                raise ValueError("The length of the mutated vector is incorrect.")
        elif mutation.type == "Deletion":
            if len(mutated_vector) != len(self.vector.seq) - mutation.length_deletion:
                raise ValueError("The length of the mutated vector is incorrect.")
                
    @staticmethod
    def circular_index(index, sequence_length):
        """This function returns the circular (=positive) index in a sequence."""
        return (index + sequence_length) % sequence_length
    
    @staticmethod
    def slice_circular_sequence(sequence, start_index, end_index):
        start_index = Vector.circular_index(start_index, len(sequence))
        if end_index < start_index:
            sliced_sequence = sequence[start_index:] + sequence[:end_index]
        else:
            sliced_sequence = sequence[start_index:end_index]
        return sliced_sequence
    
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