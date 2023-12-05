import sys
from src.eblocks import Utils


class Mutation:
    """
    This class represents a mutation.
    """

    def __init__(self, 
                 type: str = None,
                 input: str = None,
                 mutation: list = None,
                 insert: str = None,
                 deletion: str = None,
                 position: int = None, 
                 wt_residue: str = None, 
                 mut_residue: str = None):
        self.type = type
        self.input = input
        self.mutation = mutation
        self.position = position
        self.wt_residue = wt_residue
        self.mut_residue = mut_residue
        self.mutations = []

    def parse_mutations(self, fp: str):
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
            for line in content:
                str_spl_line = line.strip().split()
                if len(str_spl_line) == 1:
                    mutations.append(Mutation(type="Mutation", input=line, mutation=[str_spl_line[0]]))
                elif (len(str_spl_line) == 2) and (str_spl_line[0] == "Combined"):
                    mutations.append(Mutation(type="Combined", input=line, mutation=str_spl_line[1].split("-")))
                elif (len(str_spl_line) == 2) and (str_spl_line[0] == "Deletion"):
                    mutations.append(Mutation(type="Deletion", input=line, mutation=str_spl_line[1]))
                elif (len(str_spl_line) == 2) and (str_spl_line[0] == "Insert"):
                    mutations.append(Mutation(type="Insert", input=line, mutation=str_spl_line[1]))
                else:
                    print(f"Please check format of mutation {line}")
                    sys.exit()

        self.mutations = mutations
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
