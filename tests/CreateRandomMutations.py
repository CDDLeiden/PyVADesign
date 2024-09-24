import sys
import os
import random
from datetime import datetime

sys.path.append(os.path.join(os.path.dirname(__file__), '../src'))

from sequence import Plasmid
from utils import Utils


# python src/CLI.py -g tutorial/files/A0QX55.fasta -v tutorial/files/vector.dna -s data/codon_usage/Mycobacterium_smegmatis.csv -m tests/randominput/1_mutations_random_N=120_2024-04-22.txt -o tests/randomoutput/1



class RandomMutations:

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
        protein_sequence = RandomMutations.translate_sequence(nucleotide_sequence)
        residues = [i + str(j) for i, j in zip(protein_sequence, range(1, len(protein_sequence) + 1))]

        random_mutations = {}

        # List containing all natural amino acids and options for our mutations
        amino_acids = Utils.aas()

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

    def random_mutations_to_file(self, output_dir, n, mutations):
        now = datetime.now()
        dt_string = now.strftime("%Y-%m-%d")
        mutationsfile = os.path.join(output_dir, f'mutations_random_N={n}_{dt_string}.txt')

        with open(mutationsfile, 'w') as f:
            for i in mutations['single_mutations']:
                f.write(i.upper() + '\n')
            for i in mutations['double_mutations']:
                f.write('Combined ' + i.upper() + '\n')
            for i in mutations['insertions']:
                f.write('Insert ' + i.upper() + '\n')
            for i in mutations['deletions']:
                f.write('Deletion ' + i.upper() + '\n')

    @staticmethod
    # TODO Change this function if possible
    def translate_sequence(sequence: str) -> str:
        return sequence.translate()


if __name__ == "__main__":

    # TODO Check for different gene lengths and different mutations what are feasible input numbers

    # sequence
    sequence_file = r"tutorial\files\A0QX55.fasta"
    vector_file = r"tutorial/files/vector.dna"
    sequence_instance = Plasmid()
    sequence_instance.parse_vector(vector_file)
    sequence_instance.parse_sequence(sequence_file)

    # output directory
    output_dir = r"tests/randominput"
    n_single_mutations = 100
    n_multiple_mutations = 50
    n_deletions = 30
    n_insertions = 30

    random_mutations = RandomMutations()
    random_mutations.make_random_mutations(nucleotide_sequence=sequence_instance.sequence,
                                          n_single_mutations=n_single_mutations, 
                                          n_multiple_mutations=n_multiple_mutations,
                                          n_deletions=n_deletions,
                                          n_insertions=n_insertions, 
                                          output_dir=output_dir)