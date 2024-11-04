import sys
import itertools
import numpy as np


class Mutation:
    """
    This class represents a mutation.
    """
    def __init__(self, 
                 type: str = None,
                 input: str = None,
                 name: str = None,
                 mutation: list = None,
                 n_mutants: int = None,
                 insert: str = None,
                 idx_dna: list = None,
                 idx_dna_deletion_begin: int = None,
                 idx_dna_deletion_end: int = None,
                 length_deletion: int = None,
                 length_insert: int = None,
                 is_singlemutation = False,
                 is_insert = False,
                 is_deletion = False,
                 is_multiplemutation = False):
        
        self.type = type
        self.input = input
        self.name = name
        self.mutation = mutation
        self.n_mutants = n_mutants
        self.insert = insert
        self.idx_dna = idx_dna
        self.idx_dna_deletion_begin = idx_dna_deletion_begin
        self.idx_dna_deletion_end = idx_dna_deletion_end
        self.length_deletion = length_deletion
        self.length_insert = length_insert
        self.is_singlemutation = is_singlemutation
        self.is_insert = is_insert
        self.is_deletion = is_deletion
        self.is_multiplemutation = is_multiplemutation
        
        self.mutations = []
        self.silent_mutations = []
        self.unprocessed_mutations = []
        self.unvalid_mutations = []
        self.colors = {'Mutation': '#000000', 'Insert': '#FF0000', 'Deletion': '#0000FF', 'Combined': '#00ff00', 'Silent': '#FFA500'}

    def parse_mutations(self, fp: str):
        """
        This function parses the mutations from a text file and checks the input
        """
        mutations = self.read_mutations(fp)
        Mutation.check_nonnatural_aas(mutations)
        Mutation.check_format(mutations)
                
    def read_mutations(self, fp: str):
        """
        This function reads the mutations from a text file and returns a list of Mutation objects.
        """
        mutations = []
        with open(fp, "r") as f:
            content = f.readlines()
            
            for line in content:
                stripped_line = line.strip()
                if not stripped_line or stripped_line.startswith("#"):
                    continue

                try:
                    str_spl_line = line.strip().split()

                    # Single mutation
                    if len(str_spl_line) == 1:
                        idx = int(str_spl_line[0][1:-1]) * 3
                        mut = Mutation(
                            name=str_spl_line[0],
                            input=line,
                            mutation=[str_spl_line[0]],
                            idx_dna=[idx],
                            is_singlemutation=True,
                            type = "Mutation")
                        mutations.append(mut)

                    # Combined mutation
                    elif (len(str_spl_line) == 2) and (str_spl_line[0] == "Combined"):
                        muts = str_spl_line[1].split("-")
                        idxs = list()
                        for m in muts:
                            idx = int(m[1:-1]) * 3
                            idxs.append(idx)

                        mut = Mutation(
                            name='-'.join(muts),
                            input=line,
                            mutation=muts,
                            idx_dna=idxs,
                            is_multiplemutation=True,
                            type = "Combined")
                        mutations.append(mut)

                    # Deletion
                    elif (len(str_spl_line) == 2) and (str_spl_line[0] == "Deletion"):
                        mut_begin = int(str_spl_line[1].split("-")[0][1:])
                        mut_end = int(str_spl_line[1].split("-")[1][1:])
                        mut_length = int(mut_end - mut_begin)
                        idx_end = (mut_begin * 3) + (mut_length * 3)
                        idxs = list(range((mut_begin * 3), idx_end, 3))

                        mutation = Mutation(
                            name=str_spl_line[1],
                            input=line,
                            mutation=str_spl_line[1],
                            idx_dna_deletion_begin=int(mut_begin) *3,
                            idx_dna_deletion_end=int(mut_end) *3,
                            idx_dna=idxs,
                            length_deletion=mut_length * 3,
                            is_deletion=True,
                            type = "Deletion")
                        mutations.append(mutation)
                        
                    # Insertion
                    elif (len(str_spl_line) == 2) and (str_spl_line[0] == "Insert"):
                        idx = int(str_spl_line[1].split("-")[0][1:]) * 3

                        mutation = Mutation(
                            name=str_spl_line[1],
                            input=line,
                            mutation=str_spl_line[1],
                            idx_dna=[idx],
                            insert=str_spl_line[1].split("-")[1],
                            length_insert=len(str_spl_line[1].split("-")[1]) * 3,
                            is_insert=True,
                            type = "Insert")
                        mutations.append(mutation)

                    else:
                        raise ValueError(f"Please check format of mutation {line}")
                except:
                    raise ValueError(f"Please check format of mutation {line}")
                
        valid_mutations = self.validate_mutations(mutations)
        if len(valid_mutations) != len(mutations):
            invalids = self.find_unique_items(mutations, valid_mutations)
            invalid_names = [i.name for i in invalids]
            raise ValueError(f"Please check the format of the mutation {invalid_names}.")
        mutations = self.sort_mutations(mutations)  # sort mutation by index
        self.mutations = mutations
        self.n_mutants = len(mutations)
        return mutations
        
    def find_unique_items(self, list1, list2):
        set1 = set(list1)
        set2 = set(list2)
        unique_items = set1.symmetric_difference(set2)
        return list(unique_items)
    
    def validate_mutations(self, mutations):
        valid = [m for m in mutations if isinstance(m.idx_dna, list) and m.idx_dna]
        return valid 
    
    def sort_mutations(self, mutations: list):
        """
        This function sorts the mutations by the position of the mutation.
        """
        sorted_mutations = sorted(mutations, key=lambda x: x.idx_dna[0])
        return sorted_mutations
    
    def remove_index(self, idx: int):
        self.unprocessed_mutations.append(self.mutations[idx])
        del self.mutations[idx]
    
    def print_mutations(self, padding: int = 10):
        """
        This function prints the mutations.
        """
        print("The selected mutations are:")
        for mut in self.mutations:
            if type(mut.mutation) == list:  
                tmp = str(mut.mutation)[1:-1]  # Remove list brackets
                tmp = tmp.replace("'", "")
                print(f"\t{mut.type.ljust(padding)}\t{tmp.ljust(padding)}")
            else:
               mutation = mut.mutation.replace("'", "")
               print(f"\t{mut.type.ljust(padding)}\t{mutation.ljust(padding)}")

    def extract_indices(self):
        """
        This function extracts the indices of the mutations.
        """
        indices = []  # Store all indices of the mutations
        constraints = []  # Store constraints 
        for mutation in self.mutations:
            self._process_mutation_extraction(mutation, indices, constraints)
        indices = np.asarray(indices).reshape(-1, 1).flatten()
        constraints_indices = self._generate_constraint_indices(indices, constraints)
        return indices, constraints_indices

    def _process_mutation_extraction(self, mutation, indices, constraints):
        """Process a single mutation and update indices and constraints."""
        if mutation.is_singlemutation:
            indices.append(mutation.idx_dna[0])
        elif mutation.is_insert:
            self._handle_insert_extraction(mutation, indices, constraints)
        elif mutation.is_deletion:
            self._handle_deletion_extraction(mutation, indices, constraints)
        elif mutation.is_multiplemutation:
            self._handle_multiple_mutation_extraction(mutation, indices, constraints)

    def _handle_insert_extraction(self, mutation, indices, constraints):
        """Handle insert mutations."""
        start_insert = mutation.idx_dna[0]
        end_insert = start_insert + mutation.length_insert
        constraints.append((start_insert, end_insert))
        indices.extend([start_insert, end_insert])

    def _handle_deletion_extraction(self, mutation, indices, constraints):
        """Handle deletion mutations."""
        indices.append(mutation.idx_dna_deletion_begin)
        indices.append(mutation.idx_dna_deletion_end)
        constraints.append((mutation.idx_dna_deletion_begin, mutation.idx_dna_deletion_end))
        # print(f"Deletion indices: {mutation.idx_dna_deletion_begin}, {mutation.idx_dna_deletion_end}")

    def _handle_multiple_mutation_extraction(self, mutation, indices, constraints):
        """Handle multiple mutations."""
        idxs = list(mutation.idx_dna)
        indices.extend(idxs)
        constraints.append(tuple(idxs))
        # print(f"Multiple mutation indices: {idxs}")
        # print(f"Multiple mutation constraints: {constraints}")

    def _generate_constraint_indices(self, indices, constraints):
        """Generate constraint indices from the given constraints."""
        constraint_indices = []
        for con in constraints:
            for pair in itertools.combinations(con, 2):
                idx_a = np.where(indices == pair[0])[0]
                idx_b = np.where(indices == pair[1])[0]
                if idx_a.size > 0 and idx_b.size > 0:
                    constraint_indices.append((int(idx_a[0]), int(idx_b[0])))
        return constraint_indices

    @classmethod
    def check_nonnatural_aas(cls, mutations: list):
        """
        This function checks whether there are any non-natural amino acids in the mutations.
        """
        for mutation in mutations:
            if mutation.is_singlemutation:
                wt_residue = mutation.mutation[0][0].lower()
                mut_residue = mutation.mutation[0][-1].lower()
                if (wt_residue or mut_residue) not in Mutation.aas():
                    raise ValueError(f"Please check for non-natural amino acids in mutation {mutation.input}")
            elif mutation.is_multiplemutation:
                for m in mutation.mutation:
                    wt_residue = m[0].lower()
                    mut_residue = m[-1].lower()
                    if (wt_residue or mut_residue) not in Mutation.aas():
                        raise ValueError(f"Please check for non-natural amino acids in mutation {mutation.input}")
            elif mutation.is_deletion:
                start_res = mutation.mutation.split("-")[0][0].lower()
                end_res = mutation.mutation.split("-")[1][0].lower()
                if (start_res or end_res) not in Mutation.aas():
                    raise ValueError(f"Please check for non-natural amino acids in mutation {mutation.input}")
            elif mutation.is_insert:
                start_res = mutation.mutation.split("-")[0][0].lower()
                insert = mutation.mutation.split("-")[1]
                if (start_res or insert) not in Mutation.aas():
                    raise ValueError(f"Please check for non-natural amino acids in mutation {mutation.input}")

    @classmethod
    def check_format(cls, mutations: list):
        """
        This function checks the format of the mutations.
        """
        for mutation in mutations:
            if mutation.is_singlemutation:
                    if (mutation.mutation[0][0].lower() in Mutation.aas()) \
                    and (mutation.mutation[0][-1].lower() in Mutation.aas()) \
                    and (isinstance(int(mutation.mutation[0][1:-1]), int)):
                        continue
                    else:
                        raise ValueError(f"Please check format of mutation: {mutation.input}")
            elif mutation.is_deletion:
                    if (mutation.mutation.split("-")[0][0].lower() in Mutation.aas()) \
                    and (mutation.mutation.split("-")[1][0].lower() in Mutation.aas()) \
                    and (isinstance(int(mutation.mutation.split("-")[0][1:]), int)) \
                    and (isinstance(int(mutation.mutation.split("-")[1][1:]), int)):
                        continue
                    else:
                        raise ValueError(f"Please check format of mutation: {mutation.input}")
            elif mutation.is_insert:
                    if (mutation.mutation.split("-")[0][0].lower() in Mutation.aas()) \
                    and (set(''.join(mutation.mutation.split("-")[1]).lower())).issubset(Mutation.aas()) \
                    and (isinstance(int(mutation.mutation.split("-")[0][1:]), int)):
                        continue
                    else:
                        raise ValueError(f"Please check format of mutation: {mutation.input}")
            elif mutation.is_multiplemutation:
                for m in mutation.mutation:
                    if (m[0].lower() in Mutation.aas()) \
                    and (m[-1].lower() in Mutation.aas()) \
                    and (isinstance(int(m[1:-1]), int)):
                        continue
                    else:
                        raise ValueError(f"Please check format of mutation: {mutation.input}")
            else:
                pass
    
    @staticmethod
    def aas():
        return "acdefghiklmnpqrstvwy"