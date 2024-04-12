import sys
import numpy as np
from .utils import Utils


class Mutation:
    """
    This class represents the information for the selected mutations.
    """

    def __init__(self, 
                 type: str = None,
                 input: str = None,
                 mutation: list = None,
                 n_mutants: int = None,
                 insert: str = None,
                 deletion: str = None,
                 idx_dna: list = None,
                #  idx_start: int = None,
                #  idx_end: int = None,
                 idx_dna_deletion_begin: int = None,
                 idx_dna_deletion_end: int = None,
                 length_deletion: int = None,
                 length_insert: int = None,
                 position: int = None, 
                 wt_residue: str = None, 
                 mut_residue: str = None,
                 
                is_singlemutation = False,
                is_insert = False,
                is_deletion = False,
                is_multiplemutation = False):
        
        self.type = type
        self.input = input
        self.mutation = mutation
        self.position = position
        self.idx_dna = idx_dna
        self.paired = []
        self.insert = insert
        self.deletion = deletion
        self.length_insert = length_insert
        self.idx_dna_deletion_begin = idx_dna_deletion_begin
        self.idx_dna_deletion_end = idx_dna_deletion_end
        # self.idx_start = idx_start
        # self.idx_end = idx_end
        self.length_deletion = length_deletion
        self.wt_residue = wt_residue
        self.wt_residue = wt_residue
        self.mut_residue = mut_residue
        self.n_mutantsn = n_mutants
        self.mutations = []
        self.colors = {'Mutation': 'black', 'Insert': 'red', 'Deletion': 'blue', 'Combined': 'green'}

        self.is_singlemutation = is_singlemutation
        self.is_insert = is_insert
        self.is_deletion = is_deletion
        self.is_multiplemutation = is_multiplemutation

    def parse_mutations(self, fp: str):
        """
        This function parses the mutations from a text file and checks the input.
        """
        mutations = self.read_mutations(fp)
        result = self.perform_checks(mutations)
        return result
                
    def read_mutations(self, fp: str):
        """
        This function reads the mutations from a text file and returns a list of Mutation objects.
        """
        mutations = []
        with open(fp, "r") as f:
            content = f.readlines()
            count = 0
            for line in content:
                try:
                    str_spl_line = line.strip().split()
                    if len(str_spl_line) == 1:
                        count += 1
                        idx = int(str_spl_line[0][1:-1]) * 3
                        mutations.append(Mutation(type="Mutation", 
                                                input=line, 
                                                mutation=[str_spl_line[0]],
                                                idx_dna=[idx],
                                                is_singlemutation =True))
                    elif (len(str_spl_line) == 2) and (str_spl_line[0] == "Combined"):
                        count += 1
                        muts = str_spl_line[1].split("-")
                        idxs = list()
                        for m in muts:
                            idx = int(m[1:-1]) * 3
                            idxs.append(idx)
                        mutations.append(Mutation(type="Combined", 
                                                input=line, 
                                                mutation=muts,
                                                idx_dna=idxs,
                                                is_multiplemutation=True))
                    elif (len(str_spl_line) == 2) and (str_spl_line[0] == "Deletion"):
                        count += 1
                        mut_begin = int(str_spl_line[1].split("-")[0][1:])
                        mut_end = int(str_spl_line[1].split("-")[1][1:])
                        mut_length = int(mut_end - mut_begin)
                        idx_end = (mut_begin * 3) + (mut_length * 3)
                        idxs = list(range((mut_begin * 3), idx_end, 3))
                        mutations.append(Mutation(type="Deletion", 
                                                input=line, 
                                                mutation=str_spl_line[1],
                                                idx_dna_deletion_begin=int(mut_begin) *3,
                                                idx_dna_deletion_end=int(mut_end) *3,
                                                idx_dna=idxs,
                                                length_deletion=mut_length * 3,
                                                is_deletion=True))
                    elif (len(str_spl_line) == 2) and (str_spl_line[0] == "Insert"):
                        count += 1
                        idx = int(str_spl_line[1].split("-")[0][1:]) * 3
                        mutations.append(Mutation(type="Insert", 
                                                  input=line, 
                                                  mutation=str_spl_line[1],
                                                  idx_dna=[idx],
                                                  insert=str_spl_line[1].split("-")[1],
                                                  length_insert=len(str_spl_line[1].split("-")[1]) * 3,
                                                  is_insert=True))
                    else:
                        print(f"Please check format of mutation {line}")
                        sys.exit()
                except:
                    print(f"Please check format of mutation {line}")
                    sys.exit()

        self.mutations = mutations
        self.n_mutants = count
        return mutations
    
    def print_mutations(self):
        """
        This function prints the mutations.
        """
        print("The selected mutations are:")
        for mut in self.mutations:
            if type(mut.mutation) == list:  
                tmp = str(mut.mutation)[1:-1]  # Remove list brackets
                tmp = tmp.replace("'", "")
                print(f"\t{mut.type}\t{tmp}")
            else:
               mutation = mut.mutation.replace("'", "")
               print(f"\t{mut.type}\t{mutation}")

    def mean_idxs_dna(self) -> list:
        """
        This function returns the first index of each mutation.
        """
        mean_idxs = []
        for mutation in self.mutations:
            if len(mutation.idx_dna) > 1:
                mean_idxs.append(float(np.mean(mutation.idx_dna)))
            elif len(mutation.idx_dna) == 1:
                mean_idxs.append(mutation.idx_dna[0])
        return mean_idxs

    @classmethod
    def perform_checks(cls, mutations: list):
        """
        This function processes the mutations from a file and returns a list of Mutation objects.
        """
        Mutation.check_nonnatural_aas(mutations)
        Mutation.check_format(mutations)
        return 1

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
                    print(f"Please check for non-natural amino acids in mutation {mutation.input}")
                    sys.exit()
            elif mutation.type == "Combined":
                for m in mutation.mutation:
                    wt_residue = m[0].lower()
                    mut_residue = m[-1].lower()
                    if (wt_residue or mut_residue) not in Utils.natural_amino_acids():
                        print(f"Please check for non-natural amino acids in mutation {mutation.input}")
                        sys.exit()
            elif mutation.type == "Deletion":
                start_res = mutation.mutation.split("-")[0][0].lower()
                end_res = mutation.mutation.split("-")[1][0].lower()
                if (start_res or end_res) not in Utils.natural_amino_acids():
                    print(f"Please check for non-natural amino acids in mutation {mutation.input}")
                    sys.exit()
            elif mutation.type == "Insert":
                start_res = mutation.mutation.split("-")[0][0].lower()
                insert = mutation.mutation.split("-")[1]
                if (start_res or insert) not in Utils.natural_amino_acids():
                    print(f"Please check for non-natural amino acids in mutation {mutation.input}")
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
                    and (isinstance(int(mutation.mutation[0][1:-1]), int)):
                        continue
                    else:
                        print(f"Please check format of mutation: {mutation.input}")
                        sys.exit()
            elif mutation.type == "Deletion":
                    if (mutation.mutation.split("-")[0][0].lower() in Utils.natural_amino_acids()) \
                    and (mutation.mutation.split("-")[1][0].lower() in Utils.natural_amino_acids()) \
                    and (isinstance(int(mutation.mutation.split("-")[0][1:]), int)) \
                    and (isinstance(int(mutation.mutation.split("-")[1][1:]), int)):
                        continue
                    else:
                        print(f"Please check format of mutation: {mutation.input}")
                        sys.exit()
            elif mutation.type == "Insert":
                    if (mutation.mutation.split("-")[0][0].lower() in Utils.natural_amino_acids()) \
                    and (set(''.join(mutation.mutation.split("-")[1]).lower())).issubset(Utils.natural_amino_acids()) \
                    and (isinstance(int(mutation.mutation.split("-")[0][1:]), int)):
                        continue
                    else:
                        print(f"Please check format of mutation: {mutation.input}")
                        sys.exit()
            elif mutation.type == "Combined":
                for m in mutation.mutation:
                    if (m[0].lower() in Utils.natural_amino_acids()) \
                    and (m[-1].lower() in Utils.natural_amino_acids()) \
                    and (isinstance(int(m[1:-1]), int)):
                        continue
                    else:
                        print(f"Please check format of mutation: {mutation.input}")
                        sys.exit()
            else:
                pass

        