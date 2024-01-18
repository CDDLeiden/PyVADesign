import sys
from src.eblocks import Utils


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
                 length_deletion: int = None,
                 position: int = None, 
                 wt_residue: str = None, 
                 mut_residue: str = None):
        
        self.type = type
        self.input = input
        self.mutation = mutation
        self.position = position
        self.idx_dna = idx_dna
        # self.idx_start = idx_start
        # self.idx_end = idx_end
        self.length_deletion = length_deletion
        self.wt_residue = wt_residue
        self.mut_residue = mut_residue
        self.n_mutantsn = n_mutants
        self.mutations = []

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
                                                idx_dna=[idx]))
                    elif (len(str_spl_line) == 2) and (str_spl_line[0] == "Combined"):
                        count += 1
                        mutations = str_spl_line[1].split("-")
                        idxs = []
                        for m in mutations:
                            idx = int(m[1:-1]) * 3
                            idxs.append(idx)
                        mutations.append(Mutation(type="Combined", 
                                                input=line, 
                                                mutation=mutations,
                                                idx_dna=idxs))
                    elif (len(str_spl_line) == 2) and (str_spl_line[0] == "Deletion"):
                        count += 1
                        mut_begin = int(str_spl_line[1].split("-")[0][1:])
                        mut_end = int(str_spl_line[1].split("-")[1][1:])
                        mut_length = int(mut_end - mut_begin)
                        idx_end = (mut_begin * 3) + (mut_length * 3)
                        idxs = list(range((mut_begin * 3), idx_end, 3))
                        mutations.append(Mutation(type="Deletion", 
                                                input=line, 
                                                mutation=str_spl_line[1]),
                                                idx_dna=idxs)
                    elif (len(str_spl_line) == 2) and (str_spl_line[0] == "Insert"):
                        count += 1
                        idx = int(str_spl_line[1].split("-")[0][1:-1]) * 3
                        mutations.append(Mutation(type="Insert", 
                                                input=line, 
                                                mutation=str_spl_line[1],
                                                idx_dna=[idx]))
                    else:
                        print(f"Please check format of mutation {line}")
                        sys.exit()
                except:
                    print(f"Please check format of mutation {line}")
                    sys.exit()

        self.mutations = mutations
        self.n_mutants = count
        return mutations
    
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

    @staticmethod
    def mutation_colors():
        """
        This function returns a dictionary of colors for each mutation type.
        """
        return {'Mutation': 'black', 'Insert': 'red', 'Deletion': 'blue', 'Combined': 'green'}