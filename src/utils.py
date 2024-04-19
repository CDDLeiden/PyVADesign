import os
import sys
import random
import pandas as pd
from datetime import datetime

from .sequence import Plasmid


class Utils:

    def __init__(self):
            pass
    
    def make_random_mutations(self,
                              nucleotide_sequence: str,
                              n_single_mutations: int, 
                              n_multiple_mutations: int, 
                              n_deletions: int, 
                              n_insertions: int,
                              output_dir: str):
        """
        Make random mutations in a protein sequence.
        """
        protein_sequence = Plasmid.translate_sequence(nucleotide_sequence)
        residues = [i + str(j) for i, j in zip(protein_sequence, range(1, len(protein_sequence) + 1))]

        random_mutations = {}

        # List containing all natural amino acids and options for our mutations
        amino_acids = Utils.natural_amino_acids()

        # Randomly select single mutations
        selected_single_mutants = self.random_single_mutation(residues, amino_acids, n_single_mutations)
        random_mutations['single_mutations'] = selected_single_mutants
        print(f"Generated {n_single_mutations} single mutations: ", selected_single_mutants)

        # Randomly select double mutations
        selected_double_mutants = self.random_multiple_mutation(residues, amino_acids, n_multiple_mutations)
        random_mutations['double_mutations'] = selected_double_mutants
        print(f"Generated {n_multiple_mutations} paired mutations: ", selected_double_mutants)

        # Randomly select insertions
        selected_insertions = self.random_insert(residues, amino_acids, n_insertions)
        random_mutations['insertions'] = selected_insertions
        print(f"Generated {n_insertions} insertions: ", selected_insertions)

        # Randomly select deletions
        selected_deletions = self.random_deletion(residues, n_deletions)
        random_mutations['deletions'] = selected_deletions
        print(f"Generated {n_deletions} deletions: ", selected_deletions)

        n_mutations = n_single_mutations + n_multiple_mutations + n_insertions + n_deletions
        print(f"Total number of mutations: {n_mutations}")

        self.random_mutations_to_file(output_dir, n_mutations, random_mutations)
    
    def random_single_mutation(self, residues, aas, n):
        """
        Randomly select single mutations in a protein sequence.

        Parameters
        ----------
        residues : list
            List of residues in the protein sequence
        choices : list
            List of all natural amino acids
        n : int
            Number of mutations to sample
        """
        res = random.sample(residues, n)
        mut = random.choices(aas, k=n)
        mutants = [i + j for i, j in zip(res, mut)]
        return mutants

    def random_multiple_mutation(self, residues, aas, n, max_distance_between_mutants=50, max_number_mutations=5):
        """
        Randomly select multiple mutations in a protein sequence that will be combined in one mutant.

        Parameters
        ----------
        residues : list
            List of residues in the protein sequence
        aas : list
            List of all natural amino acids
        n : int
            Number of mutations to sample
        max_distance_between_mutants : int
            Maximum distance (in residues) between the two mutations, default is set to 10
        """
        mutants = []  # List to store the mutants
        region_to_be_sampled = residues[(max_distance_between_mutants + 1):-max_distance_between_mutants]
        res1 = random.sample(region_to_be_sampled, n)
        var1 = random.choices(aas, k=n) 
        mut1 = [i + j for i, j in zip(res1, var1)]
        for i in mut1:
            temp_muts = []
            num_mutations = random.sample(range(1, max_number_mutations), 1)[0]
            pos_x = random.sample(range(1, max_distance_between_mutants), num_mutations)
            vars_x = random.choices(aas, k=num_mutations)
            for j, k in zip(pos_x, vars_x):
                index = residues.index(i[0:-1])
                mut_x = residues[index + j] + k
                temp_muts.append(mut_x)
            mutants.append(i + '-' + '-'.join(temp_muts))
        return mutants
    
    def random_deletion(self, residues, n, max_length_deletion=10):
        """
        Randomly generate deletions in a protein sequence.

        Parameters
        ----------
        residues : list
            List of residues in the protein sequence
        n : int
            Number of deletions to sample
        max_length_deletion : int
            Maximum length of the deletion, default is set to 10
        """
        deletions = []
        for i in range(n):
            len_deletion = random.sample(range(1, max_length_deletion), 1)[0]
            res_b = random.sample(residues, 1)[0]
            res_e = residues[residues.index(res_b) + len_deletion]
            deletion = res_b + '-' + res_e
            deletions.append(deletion)
        return deletions
    
    def random_insert(self, residues, aas, n, max_length_insertion=10):
        """
        Randomly generate insertions in a protein sequence.

        Parameters
        ----------
        residues : list
            List of residues in the protein sequence
        choices : list
            List of all natural amino acids
        n : int
            Number of insertions to sample
        max_length_insertion : int 
            Maximum length of the insertion, default is set to 10
        """
        inserts = []  # List to store the inserts
        for i in range(n):
            len_insertion = random.sample(range(1, max_length_insertion), 1)[0]
            insertion = random.choices(aas, k=len_insertion)
            insertion = ''.join(insertion)
            residue = random.sample(residues, 1)[0]
            insert = residue + '-' + insertion
            inserts.append(insert)
        return inserts
    
    def random_mutations_to_file(output_dir, n, mutations):
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d")
        mutationsfile = os.path.join(output_dir, f'mutations_random_N={n}_{dt_string}.txt')

        with open(mutationsfile, 'w') as f:
            for i in mutations['single_mutations'].values():
                f.write(i.upper() + '\n')
            for i in mutations['double_mutations'].values():
                f.write('Combined ' + i.upper() + '\n')
            for i in mutations['insertions'].values():
                f.write('Insert ' + i.upper() + '\n')
            for i in mutations['deletions'].values():
                f.write('Deletion ' + i.upper() + '\n')

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
    


class SnapGene:
    """
    Class for generating SnapGene files
    """

    def __init__(self,
                 sequence_instance: Plasmid,
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
    
    def remove_empty_lines(self):
        # TODO Remove last empty line from the file
        pass



class OutputToFile:
    """
    Context manager for redirecting stdout to a file.
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