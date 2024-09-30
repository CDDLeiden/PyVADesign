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


# TODO Rethink name of this class
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

    def make_dir(self, directory='clones'):
        newpath = os.path.join(self.output_dir, directory)
        if not os.path.exists(newpath):
            os.makedirs(newpath)
        self.output_dir = newpath
            
    # TODO Move this file to the primers file and primers class
    def primers_to_fasta(self, primers: dict, directory: str, filename='primers.fasta'):
        """
        This function converts the primers to features that can be read by SnapGene.
        """
        self.make_file(directory, filename, header=False)
        with open(os.path.join(self.output_dir, filename), 'a') as f:
            for k, v in primers.items():
                f.write(f">{k}\n")
                f.write(f"{v}\n")
    
    def make_file(self, directory, filename, header=False):
        try:
            with open(os.path.join(directory, filename), 'r') as f:
                pass
        except FileNotFoundError:
            with open(os.path.join(directory, filename), 'w') as f:
                if header:
                    f.write("\n".join(SnapGene.gff3_header(self.vector_instance.vector.seq)))
                    f.write("\n")
                else:
                    pass

    def eblocks_to_gff3(self, eblocks: dict, output_dir, type='gene', filename='eblocks.gff3', header=True):
        """
        This function converts the eBlocks to features that can be read by SnapGene.
        """
        self.make_file(output_dir, filename, header=header)
        with open(os.path.join(output_dir, filename), 'a') as f:
            for k, v in eblocks.items():
                line = self.gff3_line(v[0], v[1], k, v[2], type)
                f.write('\t'.join(line) + '\n')


    def eblocks_to_genbank(self, wtvector, mutvector, eblocks: dict, output_dir, type='gene', filename='eblocks.gb', header=True, max_filename_length=16):
        """
        This function saves a vector to a GenBank (gb) file
        """
        sequence = Seq(mutvector)
        record = SeqRecord(sequence, id=self.gene_instance.seqid, name=self.gene_instance.seqid, description="")
        record.annotations["molecule_type"] = "DNA"
        record.annotations["organism"] = self.vector_instance.organism
        record.annotations["date"] = datetime.today().strftime('%d-%b-%Y').upper()
        record.name = record.name + "_" + filename
        # Limit filename length characters
        if len(record.name) > max_filename_length:
            record.name = record.name[:max_filename_length]
        features = []  # Add eBlock and mutations as features
        for k, v in eblocks.items():
            if v[0] > v[1]: # Start index is larger than end index
                joint_location = CompoundLocation([FeatureLocation(v[0], len(mutvector)), FeatureLocation(0, v[1])])
                joint_feature = SeqFeature(joint_location, type="gene", qualifiers={"gene": k, "color": v[2]})
                features.append(joint_feature)
            else:
                feature = SeqFeature(FeatureLocation(v[0], v[1]), type=type, qualifiers={"gene": k, "color": v[2]})
                features.append(feature)
        record.features.extend(features)     
        outpath = os.path.join(output_dir, filename)
        SeqIO.write(record, outpath, "genbank")


    # TODO Move this function to the primers file and primers class
    def add_primers_to_genbank_file(self, genbank_file, primer, direction):
        """
        This function saves primer data to an existing GenBank file
        """
        seq_record = SeqIO.read(genbank_file, "genbank")
        if direction == "fw":
            site = {"start": int(primer.idx_start), "end": int(primer.idx_end), "sequence": str(primer.sequence_5to3)}
        elif direction == "rv":
            site = {"start": int(primer.idx_start), "end": int(primer.idx_end), "sequence": seq.NucleotideSequence(primer.sequence_5to3).reverse()}
        else:
            raise ValueError("Primer is neither forward (fw) nor reverse (rv).")
        # Add primer binding sites as features to the SeqRecord
        feature_location = FeatureLocation(start=site["start"], end=site["end"])
        feature = SeqFeature(location=feature_location, type="primer_bind")
        feature.qualifiers["note"] = primer.name
        feature.qualifiers["primer_sequence"] = site["sequence"]
        seq_record.features.append(feature)
        # for feature in seq_record.features:
        #     print(feature)
        SeqIO.write(seq_record, genbank_file, "genbank")

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
    

    
class OutputToFile:
    """
    Context manager for redirecting stdout to a file
    """
    def __init__(self, filepath):
        self.filename = filepath
        self.stdout = sys.stdout
        self.file = open(filepath, 'w')

    def __enter__(self):
        sys.stdout = self.file

    def __exit__(self, exc_type, exc_value, traceback):
        sys.stdout = self.stdout
        self.file.close()


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
        # self.relative_frequencies_to_csv(relative_frequencies)
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
    











