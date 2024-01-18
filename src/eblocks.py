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
# from utils import read_codon_usage, DNA_Codons, write_pickle, natural_amino_acids

from mutation import Mutation
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
                 mutation_instance: Mutation,
                 # File paths for input files
                 gene_file: str = None,
                 output_fp: str = None,
                # Parameters for design

                 ):
        self.eblocks_instance = eblocks_instance
        self.mutation_instance = mutation_instance
        self.block_sequences = []

    def run_design_eblocks(self):
        """
        This function runs the design process for eblocks.
        """
        # Your design logic to generate block_sequences
        self.block_sequences = ["sequence1", "sequence2", "sequence3"]

        # Set the block_sequences in the Eblocks instance
        self.eblocks_instance.set_block_sequences(self.block_sequences)
