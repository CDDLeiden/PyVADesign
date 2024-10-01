import os
import sys
import random
import tempfile
import itertools
import numpy as np
import pandas as pd
import biotite.sequence as seq
import biotite.database.entrez as entrez
import biotite.sequence.io.fasta as fasta
import biotite.sequence.io.genbank as gb
from biotite.sequence import AnnotatedSequence

from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

from .sequence import Vector, Gene



class SnapGene:
    """
    Class for handling SnapGene files and generating SnapGene features
    """
    def __init__(self,
                 vector_instance: Vector,
                 gene_instance: Gene,
                 output_dir: str = None):
        
            self.output_dir = output_dir
            self.vector_instance = vector_instance
            self.gene_instance = gene_instance

    @staticmethod
    def gff3_header(length_sequence, version="3.2.1", sequence_name="myseq"):
        result = [f"##gff-version {version}", f"##sequence-region {sequence_name} 1 {len(length_sequence)}"]
        return result

    @staticmethod
    def gff3_line(begin_pos, end_pos, name, hex_color, type):
        line = ['myseq', '.', f"{type}", str(begin_pos), str(end_pos), '.', '.', '.', f"Name={name};color={hex_color}"]
        return line
    
    @staticmethod
    def gff3_colors():
        colors = {
            'mutation': '#FF0000',
            'combination': '#D8FF00',
            'insert': '#0017FF',
            'deletion': '#FF5900',
            'IVAprimer': '#FF00D4',
            'SEQprimer': '#06FF92'}
        return colors
    

# TODO Describe in the tutorial how to obtain the genome IDs from NIH
# TODO Check licensing of this code that was adopted from ... (https://www.biotite-python.org/latest/examples/gallery/sequence/misc/codon_usage.html)
class CodonUsage:
    """
    Class for generating codon usage tables using biotite
    """    
    def __init__(self,
                genome_id: str = None,
                output_dir: str = None):

        self.genome_id = genome_id
        self.output_dir = output_dir

    def run(self):
        """
        Obtain the most common codon for each amino acid in a genome
        """
        genome = self.get_genome_features()
        codon_counter = CodonUsage.count_codons(genome)
        codon_counter = CodonUsage.get_codon_usage(genome, codon_counter)
        relative_frequencies = CodonUsage.get_relative_frequencies(genome, codon_counter)
        max_codons = CodonUsage.most_common_codon_per_aa(relative_frequencies)
        return max_codons

    def get_genome_features(self):
        """
        Get the CDS features of a genome
        """
        try:
            gb_file = gb.GenBankFile.read(
                entrez.fetch(self.genome_id, tempfile.gettempdir(), "gb", "nuccore", "gb"))
        except Exception as e:
            print(f"Error fetching genome with ID '{self.genome_id}': {e}")
            return None
        genome = gb.get_annotated_sequence(gb_file, include_only=["CDS"])
        if isinstance(genome, AnnotatedSequence):
            return genome
        else:
            raise ValueError("No CDS features found in genome")
    
    def relative_frequencies_to_csv(self, frequencies_dict):
        """
        Save relative frequencies in a CSV file
        """
        outpath = os.path.join(self.output_dir, f"{self.genome_id}_codon_usage.csv")
        with open(outpath, "w") as file:
            file.write("Amino Acid,Codon,Relative Frequency\n")
            for amino_acid, codons in frequencies_dict.items():
                for codon, freq in codons:
                    file.write(f"{amino_acid},{codon},{freq}\n")

    @staticmethod
    def count_codons(genome):
        """
        Count the occurrence of each codon in a genome
        The symbols [0 1 2 3] represent ['A' 'C' 'G' 'T'] respectively
        """
        codon_counter = {
            codon: 0
            for codon in itertools.product(*([range(len(genome.sequence.alphabet))] * 3))}
        return codon_counter
    
    @staticmethod
    def get_codon_usage(genome, codon_counter):
        """
        Get the codon usage of a genome
        """
        for feature in genome.annotation:
            cds = genome[feature]  # Get the coding sequence
            if len(cds) % 3 != 0:  # malformed CDS
                continue
            for i in range(0, len(cds), 3):  # Count the codons
                codon_code = tuple(cds.code[i:i+3])
                codon_counter[codon_code] += 1
        return codon_counter
    
    @staticmethod
    def get_relative_frequencies(genome, codon_counter):
        """
        Convert the total counts into relative frequencies
        """
        table = seq.CodonTable.default_table()
        relative_frequencies = {}
        for amino_acid_code in range(20):
            codon_codes_for_aa = table[amino_acid_code]
            total = 0  # Get the total amount of codon occurrences for the amino acid
            for codon_code in codon_codes_for_aa:
                total += codon_counter[codon_code]
            for codon_code in codon_codes_for_aa:
                codon_counter[codon_code] /= total
                amino_acid = seq.ProteinSequence.alphabet.decode(amino_acid_code)
                codon = genome.sequence.alphabet.decode_multiple(codon_code)
                codon = "".join(codon)
                freq = codon_counter[codon_code]
                if amino_acid not in relative_frequencies:  # Store relative frequencies in dictionary
                    relative_frequencies[amino_acid] = []
                relative_frequencies[amino_acid].append((codon, round(freq, 3)))
        return relative_frequencies
    
    @staticmethod
    def most_common_codon_per_aa(relative_frequencies):
        """
        Get the most common codon for each amino acid
        """
        max_freqs = {}
        for amino_acid, codons in relative_frequencies.items():
            max_codon = max(codons, key=lambda x: x[1])
            max_freqs[amino_acid] = max_codon
        return max_freqs