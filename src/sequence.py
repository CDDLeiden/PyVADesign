import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


# TODO Rename Plasmid


class Plasmid:
    """
    This class contains functions to parse and translate DNA sequences and vectors.
    """

    def __init__(self):
        self.sequence: str = None  # Gene sequence
        self.seqid: str = None
        self.organism: str = None
        
        self.vector: str = None  # Vector sequence with gene cloned into it
        self.vector_id: str = None

        self.gene_start_idx: int = -1
        self.gene_end_idx: int = -1
        self.color: str = "#d3d3d3"

    def to_gb(self, output_dir, filename):
        """
        This function saves the vector to a GenBank file.
        """
        # TODO Update this function 
        # Create a sequence record
        sequence = Seq("ATGCGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA")
        record = SeqRecord(sequence, id="test_sequence", description="Example DNA sequence")
        # Add molecule type to annotations
        record.annotations["molecule_type"] = "DNA"
        # Add organism information to annotations
        record.annotations["organism"] = "Example organism"
        # Add date to the locus line
        record.annotations["date"] = datetime.today().strftime('%d-%b-%Y').upper()
        # Add features to the sequence record
        features = [
            SeqFeature(FeatureLocation(1, 50), type="gene", qualifiers={"gene": "gene1"}),
            SeqFeature(FeatureLocation(51, 150), type="CDS", qualifiers={"gene": "gene1", "product": "protein1"}),
            # Add more features as needed
        ]
        record.features.extend(features)
        # Write the record to a GenBank file
        SeqIO.write(record, "test_sequence.gb", "genbank")

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
        vector = self.read_snapgene_dna_file(fp)
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
        result = self.check_sequence(sequence)
        return result
    
    def mutate_vector(self, idx_start, idx_end, sequence):
        """
        This function mutates the input vector with a given sequence
        """
        mutated_vector = self.vector.seq[:idx_start] + sequence + self.vector.seq[idx_end:]
        return mutated_vector
        
    @staticmethod
    def circular_index(index, sequence_length):
        """
        This function returns the circular index of a sequence.
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