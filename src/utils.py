#!/usr/bin/env python3

import os
import tempfile
import itertools
import subprocess
import biotite.sequence as seq
import biotite.database.entrez as entrez
import biotite.sequence.io.genbank as gb
from biotite.sequence import AnnotatedSequence

def get_current_git_commit():
    """Get the current git commit hash."""
    try:
        commit_hash = subprocess.check_output(["git", "rev-parse", "HEAD"]).strip().decode('utf-8')
        return commit_hash
    except subprocess.CalledProcessError:
        return None

class CodonUsage:
    """
    Class for generating codon usage tables using biotite

    Attributes:
    ----------
    genome_id : str
        The genome ID
    output_dir : str
        The output directory
    """    
    def __init__(self,
                genome_id: str = None,
                output_dir: str = None):

        self.genome_id = genome_id
        self.output_dir = output_dir

    def run(self):
        """Obtain the most common codon for each amino acid in a genome."""
        genome = self.get_genome_features()
        codon_counter = CodonUsage.count_codons(genome)
        codon_counter = CodonUsage.get_codon_usage(genome, codon_counter)
        relative_frequencies = CodonUsage.get_relative_frequencies(genome, codon_counter)
        self.relative_frequencies_to_csv(relative_frequencies)
        max_codons = CodonUsage.most_common_codon_per_aa(relative_frequencies)
        return max_codons, relative_frequencies
    
    def codon_usage_exists(self):
        """Check if the codon usage CSV file exists."""
        path = os.path.join(self.output_dir, f"{self.genome_id}_codon-usage.csv")
        return os.path.exists(path)

    def get_genome_features(self):
        """Get the CDS features of a genome."""
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
        
    def load_relative_frequencies(self):
        """Load relative frequencies from a CSV file."""
        path = os.path.join(self.output_dir, f"{self.genome_id}_codon-usage.csv")
        frequencies_dict = {}
        with open(path, "r") as file:
            content = file.readlines()
            for line in content[1:]:
                amino_acid, codon, freq = line.strip().split(",")
                freq = float(freq)
                if amino_acid not in frequencies_dict:
                    frequencies_dict[amino_acid] = []
                frequencies_dict[amino_acid].append((codon, freq))
        max_codons = CodonUsage.most_common_codon_per_aa(frequencies_dict)
        self.relative_frequencies = frequencies_dict
        return max_codons, frequencies_dict

    def relative_frequencies_to_csv(self, frequencies_dict):
        """Save relative frequencies in a CSV file."""
        outpath = os.path.join(self.output_dir, f"{self.genome_id}_codon-usage.csv")
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
        """Get the codon usage of a genome."""
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
        """Convert the total counts into relative frequencies."""
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
        """Get the most common codon for each amino acid."""
        max_freqs = {}
        for amino_acid, codons in relative_frequencies.items():
            max_codon = max(codons, key=lambda x: x[1])
            max_freqs[amino_acid] = max_codon
        return max_freqs