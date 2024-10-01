import sys
import os
import random
from datetime import datetime

sys.path.append(os.path.join(os.path.dirname(__file__), '../src'))

from sequence import Vector, Gene


class RandomMutations:

    def __init__(self):
        pass

    def make_random_mutations(self,
                                nucleotide_sequence: str,
                                n_single_mutations: int, 
                                n_multiple_mutations: int, 
                                n_deletions: int, 
                                n_insertions: int,
                                output_dir: str,
                                filename):
        """
        Make random mutations in a protein sequence.
        """
        protein_sequence = RandomMutations.translate_sequence(nucleotide_sequence)
        residues = [i + str(j) for i, j in zip(protein_sequence, range(1, len(protein_sequence) + 1))][:-1]

        random_mutations = {}

        # List containing all natural amino acids and options for our mutations
        amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

        # Randomly select single mutations
        selected_single_mutants = self.random_single_mutation(residues, amino_acids, n_single_mutations)
        random_mutations['single_mutations'] = selected_single_mutants
        # print(f"Generated {n_single_mutations} single mutations: ", selected_single_mutants)

        # Randomly select double mutations
        selected_double_mutants = self.random_multiple_mutation(residues, amino_acids, n_multiple_mutations)
        random_mutations['double_mutations'] = selected_double_mutants
        # print(f"Generated {n_multiple_mutations} paired mutations: ", selected_double_mutants)

        # Randomly select insertions
        selected_insertions = self.random_insert(residues, amino_acids, n_insertions)
        random_mutations['insertions'] = selected_insertions
        # print(f"Generated {n_insertions} insertions: ", selected_insertions)

        # Randomly select deletions
        selected_deletions = self.random_deletion(residues, n_deletions)
        random_mutations['deletions'] = selected_deletions
        # print(f"Generated {n_deletions} deletions: ", selected_deletions)

        n_mutations = n_single_mutations + n_multiple_mutations + n_insertions + n_deletions
        # print(f"Total number of mutations: {n_mutations}")

        self.random_mutations_to_file(output_dir, n_mutations, random_mutations, filename)


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
            # len_deletion = random.sample(range(1, max_length_deletion), 1)[0]
            # res_b = random.sample(residues, 1)[0]
            # res_e = residues[residues.index(res_b) + len_deletion]
            # deletion = res_b + '-' + res_e
            # deletions.append(deletion)
                    # Ensure len_deletion is at least 1 and at most max_length_deletion
            len_deletion = random.randint(1, min(max_length_deletion, len(residues)))

            # Select a residue from the list
            res_b = random.choice(residues)

            # Get the starting index of the selected residue
            start_index = residues.index(res_b)

            # Check if we can safely calculate the end index
            if start_index + len_deletion < len(residues):
                res_e = residues[start_index + len_deletion]
                deletion = res_b + '-' + res_e
                deletions.append(deletion)
            else:
                # If out of range, handle it, e.g., skip or create a deletion with available residues
                # Here, we'll just create a deletion that ends at the last residue
                deletion = res_b + '-' + residues[-1]
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

    def random_mutations_to_file(self, output_dir, n, mutations, filename=None):
        if not filename:
            now = datetime.now()
            dt_string = now.strftime("%Y-%m-%d")
            mutationsfile = os.path.join(output_dir, f'mutations_random_N={n}_{dt_string}.txt')
        else:
            mutationsfile = os.path.join(output_dir, filename)

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
    def translate_sequence(sequence: str) -> str:
        return sequence.translate()


if __name__ == "__main__":

    current_dir = os.path.dirname(__file__)
    base_dir = os.path.abspath(os.path.join(current_dir, '..'))

    # sequence
    gene = 'A0QX55'
    sequence_file = os.path.join(base_dir, "example_data", "Msmegmatis_DnaE1", "A0QX55.fasta")
    vector_file = os.path.join(base_dir, "example_data", "Msmegmatis_DnaE1", "vector.dna")
        
    # Output directory
    output_dir = os.path.join(base_dir, "tests", "randominput")

    gene_instance = Gene()
    gene_instance.parse_sequence(sequence_file)

    vector_instance = Vector(gene=gene_instance)
    vector_instance.parse_vector(vector_file)

    # output directory
    for i in range(1, 200):
        # randomly choose number of mutations
        n_single_mutations = random.randint(1, 100)
        n_multiple_mutations = random.randint(1, 50)
        n_deletions = random.randint(1, 50)
        n_insertions = random.randint(1, 50)

        random_mutations = RandomMutations()
        random_mutations.make_random_mutations(nucleotide_sequence=gene_instance.sequence,
                                            n_single_mutations=n_single_mutations, 
                                            n_multiple_mutations=n_multiple_mutations,
                                            n_deletions=n_deletions,
                                            n_insertions=n_insertions, 
                                            output_dir=output_dir,
                                            filename=f"{i}_{gene}_mutations.txt")