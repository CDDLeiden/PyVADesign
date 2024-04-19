import sys

sys.path.append('../src')

from sequence import Plasmid
from utils import Utils

# TODO Check for different gene lengths and different mutations what are feasible input numbers

# output directory
output_dir = "files/randominput"
n_single_mutations = 100
n_multiple_mutations = 5
n_deletions = 10
n_insertions = 5


# sequence
sequence_file = '../tutorial/files/A0QX55.fasta'
vector_file = '../tutorial/files/vector.dna'
sequence_instance = Plasmid()
sequence_instance.parse_vector(vector_file)
sequence_instance.parse_sequence(sequence_file)


# Create random mutations
utils = Utils()
utils.create_random_mutations(nucleotide_sequence=sequence_instance.sequence,
                              n_single_mutations=n_single_mutations, 
                              n_multiple_mutations=n_multiple_mutations,
                              n_deletions=n_deletions,
                              n_insertions=n_insertions, 
                              output_dir="output")