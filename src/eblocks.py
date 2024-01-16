import os
import sys
import math
import random
# import openpyxl
import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
# from dna_features_viewer import GraphicFeature, GraphicRecord
from utils import read_codon_usage, DNA_Codons, write_pickle, natural_amino_acids

script_dir = os.path.dirname(os.path.abspath(__file__))
template_path_96 = os.path.join(script_dir, 'data/eblocks-plate-upload-template-96.xlsx')
template_path_384 = os.path.join(script_dir, 'data/eblocks-plate-upload-template-384.xlsx')

# TODO Add @classmethod function to some of the functions
# TODO Make a separate class for plotting maybe?
# TODO Change filepath (str) to Path object

class Eblocks:
    def __init__(self):
        self.eblock_parameters = {"common_param": "shared_value"}
        self.block_sequences = []

    def set_parameters(self, parameters):
        self.eblock_parameters = parameters

    def set_block_sequences(self, block_sequences):
        self.block_sequences = block_sequences

    def display_parameters(self):
        print("Eblocks Parameters:", self.eblock_parameters)
        print("Block Sequences:", self.block_sequences)


class EblockDesign:
    def __init__(self, 
                 eblocks_instance: Eblocks,
                 # File paths for input files
                 mutations_file: str = None,
                 gene_file: str = None,
                 output_fp: str = None,
                # Parameters for design

                 ):
        self.eblocks_instance = eblocks_instance
        self.block_sequences = []

    def run_design_eblocks(self):
        """
        This function runs the design process for eblocks.
        """
        # mutations = [
        #     Mutation("SNP", 1, "ATG", "GTG", "eblock1"),
        #     Mutation("SNP", 2, "ATG", "GTG", "eblock2"),
        #     Mutation("SNP", 3, "ATG", "GTG", "eblock3"),
        # ]


        # Your design logic to generate block_sequences
        self.block_sequences = ["sequence1", "sequence2", "sequence3"]

        # Set the block_sequences in the Eblocks instance
        self.eblocks_instance.set_block_sequences(self.block_sequences)


class Utils:
    @staticmethod
    def natural_amino_acids():
        """
        This function returns a string of valid amino acids.
        """
        return "acdefghiklmnpqrstvwy"


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


class Mutation:
    # TODO Update attributes according to type and add list of mutations for combined
    # TODO Also when processing the mutations, update this class 
    # TODO - add mutation_type_colors
    # TODO After completion, write some tests for this class

    def __init__(self, 
                 type: str = None,
                 mutation: list = None,
                 insert: str = None,
                 deletion: str = None,
                 position: int = None, 
                 wt_residue: str = None, 
                 mut_residue: str = None):
        self.type = type
        self.mutation = mutation
        self.position = position
        self.wt_residue = wt_residue
        self.mut_residue = mut_residue
        self.mutations = []
            
    @classmethod
    def read_mutations(cls, fp: str):
        """
        This function reads the mutations from a text file and returns a list of Mutation objects.
        """
        mutations = []
        with open(fp, "r") as f:
            content = f.readlines()
            for line in content:
                line = line.strip()
                if len(line) == 1:
                    mutations.append(Mutation(type="Mutation", mutation=[line[0]]))
                elif (len(line) == 2) and (line[0] == "Combined"):
                    mutations.append(Mutation(type="Combined", mutation=line[1].split("-")[1]))
                elif (len(line) == 2) and (line[0] == "Deletion"):
                    mutations.append(Mutation(type="Deletion", mutation=line[1]))
                elif (len(line) == 2) and (line[0] == "Insert"):
                    mutations.append(Mutation(type="Insert", mutation=line[1]))
                else:
                    print(f"Please check format of mutation {line}")
                    sys.exit()
        return mutations

    @classmethod
    def perform_checks(cls, fp: str):
        """
        This function processes the mutations from a file and returns a list of Mutation objects.
        """
        mutations = cls.read_mutations(fp)
        Mutation.check_duplicates(mutations)
        Mutation.check_nonnatural_aas(mutations)
        Mutation.check_format(mutations)

    @classmethod
    def check_duplicates(cls, mutations: list):
        """
        This function checks whether there are any duplicate mutations.
        """
        mutations_list = [mutation.mutation for mutation in mutations]
        if len(mutations_list) != len(set(mutations_list)):
            print("Please check for duplicate mutations.")
            sys.exit()

    @classmethod
    def check_nonnatural_aas(cls, mutations: list):
        """
        This function checks whether there are any non-natural amino acids.
        """
        for mutation in mutations:
            if mutation.type == "Mutation":
                wt_residue = mutation.mutation[0][0].lower()
                mut_residue = mutation.mutation[0][-1].lower()
                if (wt_residue or mut_residue) not in Utils.natural_amino_acids():
                    print(f"Please check for non-natural amino acids in mutation {mutation.mutation}")
                    sys.exit()
            elif mutation.type == "Combined":
                for m in mutation.mutation:
                    wt_residue = m[0].lower()
                    mut_residue = m[-1].lower()
                    if (wt_residue or mut_residue) not in Utils.natural_amino_acids():
                        print(f"Please check for non-natural amino acids in mutation {mutation.mutation}")
                        sys.exit()
            elif mutation.type == "Deletion":
                start_res = mutation.mutation.split("-")[0][0].lower()
                end_res = mutation.mutation.split("-")[1][0].lower()
                if (start_res or end_res) not in Utils.natural_amino_acids():
                    print(f"Please check for non-natural amino acids in mutation {mutation.mutation}")
                    sys.exit()
            elif mutation.type == "Insert":
                start_res = mutation.mutation.split("-")[0][0].lower()
                insert = mutation.mutation.split("-")[1]
                if (start_res or insert) not in Utils.natural_amino_acids():
                    print(f"Please check for non-natural amino acids in mutation {mutation.mutation}")
                    sys.exit()

    @classmethod
    def check_format(cls, mutations: list):
        """
        This function checks the format of the mutations.
        """
        for mutation in mutations:
            if mutation.type == "Mutation":
                    if (mutation.mutation[0][0].lower() in Utils.natural_amino_acids()) \
                    and (mutation.mutation[0][-1].lower() in Utils.natural_amino_acids()) \
                    and (isinstance(int(mutation.mutation[1:-1]), int)):
                        continue
                    else:
                        print(f"Please check format of mutation {mutation.mutation}")
                        sys.exit()
            elif mutation.type == "Deletion":
                    if (mutation.mutation.split("-")[0][0].lower() in Utils.natural_amino_acids()) \
                    and (mutation.mutation.split("-")[1][0].lower() in Utils.natural_amino_acids()) \
                    and (isinstance(int(mutation.mutation.split("-")[0][1:]), int)) \
                    and (isinstance(int(mutation.mutation.split("-")[1][1:]), int)):
                        continue
                    else:
                        print(f"Please check format of mutation {mutation.mutation}")
                        sys.exit()
            elif mutation.type == "Insert":
                    if (mutation.mutation.split("-")[0][0].lower() in Utils.natural_amino_acids()) \
                    and (''.join(mutation.mutation.split("-")[1]).lower().issubset(Utils.natural_amino_acids())) \
                    and (isinstance(int(mutation.mutation.split("-")[0][1:]), int)):
                        continue
                    else:
                        print(f"Please check format of mutation {mutation.mutation}")
                        sys.exit()
            elif mutation.type == "Combined":
                for m in mutation.mutation:
                    if (m[0].lower() in Utils.natural_amino_acids()) \
                    and (m[-1].lower() in Utils.natural_amino_acids()) \
                    and (isinstance(int(m[1:-1]), int)):
                        continue
                    else:
                        print(f"Please check format of mutation {mutation.mutation}")
                        sys.exit()
            else:
                pass
