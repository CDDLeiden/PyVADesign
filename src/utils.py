import os
import sys
import random
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

# from .sequence import Plasmid
# TODO Make insert + deletion clones


class Utils:

    def __init__(self):
            pass
    
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
    
    @staticmethod
    def check_directory(directory, verbose):
        # Check if the directory exists
        if not os.path.exists(directory):
            raise FileNotFoundError(f"Directory {directory} does not exist.")
        # Check if the directory is empty
        if not os.listdir(directory):
            if verbose:
                print(f"Directory {directory} is empty.")
        else:
            if verbose:
                print(f"Directory {directory} is not empty. Files might get overwritten or appended to.")




class SnapGene:
    """
    Class for generating SnapGene files
    """

    def __init__(self,
                 sequence_instance,
                 output_dir: str = None):
        
            self.output_dir = output_dir
            self.sequence_instance = sequence_instance

    def make_dir(self, directory='clones'):
        newpath = os.path.join(self.output_dir, directory)
        if not os.path.exists(newpath):
            os.makedirs(newpath)
        self.output_dir = newpath
            
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
                    f.write("\n".join(SnapGene.gff3_header(self.sequence_instance.vector.seq)))
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


    def eblocks_to_genbank(self, wtvector, mutvector, eblocks: dict, output_dir, type='gene', filename='eblocks.gb', header=True):
        """
        This function saves a vector to a GenBank (gb) file
        """
        sequence = Seq(mutvector)
        record = SeqRecord(sequence, id=self.sequence_instance.seqid, description="")
        record.annotations["molecule_type"] = "DNA"
        record.annotations["organism"] = self.sequence_instance.organism
        record.annotations["date"] = datetime.today().strftime('%d-%b-%Y').upper()
        features = []  # Add eBlock and mutations as features
        for k, v in eblocks.items():
            print(k, v)
            if v[0] > v[1]: # Start index is larger than end index
                # location1 = FeatureLocation(v[0], len(vector))
                # location2 = FeatureLocation(0, v[1])
                # feature1 = SeqFeature(location1, type=type, qualifiers={"gene": k, "color": v[2]})
                # feature2 = SeqFeature(location2, type=type, qualifiers={"gene": k, "color": v[2]})
                print(len(mutvector))
                joint_location = CompoundLocation([FeatureLocation(v[0], len(mutvector)), FeatureLocation(0, v[1])])
                joint_feature = SeqFeature(joint_location, type="gene", qualifiers={"gene": k, "color": v[2]})
                features.append(joint_feature)
            else:
                feature = SeqFeature(FeatureLocation(v[0], v[1]), type=type, qualifiers={"gene": k, "color": v[2]})
                features.append(feature)
        record.features.extend(features)     
        outpath = os.path.join(output_dir, filename)
        SeqIO.write(record, outpath, "genbank")

    def add_primers_to_genbank_file(self, genbank_file, primer):
        """
        This function saves primer data to an existing GenBank file
        """
        seq_record = SeqIO.read(genbank_file, "genbank")
        print(seq_record.id)
        if primer.is_forward:
            site = {"start": int(primer.idx_start), "end": int(primer.idx_end), "sequence": str(primer.sequence_5to3)}
        elif primer.is_reverse:
            site = {"start": int(primer.idx_start), "end": int(primer.idx_end), "sequence": self.sequence_instance.invert_sequence(primer.sequence_5to3)}
        else:
            raise ValueError("Primer is neither forward nor reverse.")
        # Add primer binding sites as features to the SeqRecord
        feature_location = FeatureLocation(start=site["start"], end=site["end"])
        feature = SeqFeature(location=feature_location, type="primer_bind")
        feature.qualifiers["note"] = primer.name
        feature.qualifiers["primer_sequence"] = site["sequence"]
        seq_record.features.append(feature)
        # Print the features to identify any issues
        for feature in seq_record.features:
            print(feature)

        SeqIO.write(seq_record, genbank_file, "genbank")

    @staticmethod
    def gff3_header(length_sequence, version="3.2.1", sequence_name="myseq"):
        result = [f"##gff-version {version}", f"##sequence-region {sequence_name} 1 {len(length_sequence)}"]
        return result

    @staticmethod
    def gff3_line(begin_pos, end_pos, name, hex_color, type):
        # TODO Change feature, gene, CDS etc to the correct type
        # TODO Change Myseq to the correct sequence name?
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
    
    @staticmethod
    def count_substring_occurance(sequence: str, substring: str):
        return str(sequence).count(str(substring))
    
    @staticmethod
    def parse_codon_usage():
        pass



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