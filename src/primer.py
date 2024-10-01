import os
import sys
import math
import primer3
import difflib
import pandas as pd
import biotite.sequence as seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import MeltingTemp as mt

from .mutation import Mutation
from .sequence import Vector, Gene
from .eblocks import EblockDesign
from .utils import SnapGene

from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation


# TODO Add examples of complementary and hairpin structures to the tests directory

class Primer:
    def __init__(self,
                 name: str = None,
                 sequence_5to3: str = None,
                 is_forward: bool = False,
                 is_reverse: bool = False,
                 idx_start: int = None,
                 idx_end: int = None,
                 tm: float = None,
                 gc_content: float = None,):

        self.name = name
        self.sequence_5to3 = sequence_5to3
        self.is_forward = is_forward
        self.is_reverse = is_reverse
        self.idx_start = idx_start
        self.idx_end = idx_end
        self.tm = tm
        self.gc_content = gc_content

    # @staticmethod
    # def check_multiple_binding_sites(vector: str, sequence: str):
    #     # TODO Maybe use primerblast for this
    #     # TODO Write some checks for this function
    #     # TODO Check for exact cases where the sequence binds multiple times and almost exact cases where only few characters are different
    #     count = 0
    #     for i in range(len(vector) - len(sequence) + 1):
    #         sub = vector[i:i + len(sequence)]
    #         unique_chars = set(sub)
    #         if len(unique_chars) <= 2 or (len(unique_chars) == 3 and sub.count(sub[0]) == 2):
    #             count += 1
    #     if count > 1:
    #         print(f"Multiple binding sites for sequence {sequence} in the vector sequence")
    #     return count


class Primer3Config:
    """
    This class contains the configuration for the primer3 design
    """
    # TODO Add some other parameters as well from https://primer3.org/manual.html
    @staticmethod
    def global_args_primer_pair():
        return {
            # TODO Add task
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_PICK_INTERNAL_OLIGO': 1,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_MIN_SIZE': 12,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_MAX_POLY_X': 100,
            'PRIMER_INTERNAL_MAX_POLY_X': 100,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8,
            'PRIMER_PAIR_MAX_DIFF_TM': 5.0,
        }

    @staticmethod
    def global_args_sequencing_primer():
        # TODO Add GC clamp
        return {
            'PRIMER_TASK': 'pick_sequencing_primers',
            'PRIMER_OPT_SIZE': 16,
            'PRIMER_MIN_SIZE': 16,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 58.0,
            'PRIMER_MIN_TM': 55.0,
            'PRIMER_MAX_TM': 60.0,
            'PRIMER_MIN_GC': 40.0,
            'PRIMER_MAX_GC': 60.0,
            'PRIMER_MAX_POLY_X': 5,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 8,
            'PRIMER_MAX_SELF_END': 3,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0
        }

    @staticmethod
    def seq_args_primer_pair(sequence, product_start=None, product_end=None, force_left=None, force_right=None):
        args = {
            'SEQUENCE_ID': 'SEQ0001',
            'SEQUENCE_TEMPLATE': sequence,
            'SEQUENCE_INCLUDED_REGION': [0, len(sequence)],
        }
        if product_start is not None and product_end is not None:
            args['PRIMER_PRODUCT_SIZE_RANGE'] = [[product_start, product_end]]
        if force_left is not None:
            args['SEQUENCE_FORCE_LEFT_START'] = force_left
        if force_right is not None:
            args['SEQUENCE_FORCE_RIGHT_START'] = force_right
        return args
    
    @staticmethod
    def seq_args_sequencing_primer(sequence, start, end):
        return {
            'SEQUENCE_ID': 'SEQ0001',
            'SEQUENCE_TEMPLATE': sequence,
            'SEQUENCE_TARGET': [start, end],
        }



class DesignPrimers:
    """
    This class designs primers to open-up the expression plasmid and Sequencing primers for validation. 
    """

    def __init__(self,
                 eblocks_design_instance: EblockDesign,
                 mutation_instance: Mutation,
                 vector_instance: Vector,
                 output_dir: str = None,
                 minimum_overlap: int = 15,
                 void_length: int = 100,
                 verbose: bool = True):

        self.eblocks_design_instance = eblocks_design_instance
        self.mutation_instance = mutation_instance
        self.vector_instance = vector_instance
        self.output_dir = output_dir
        self.minimum_overlap = minimum_overlap
        self.void_length = void_length
        self.verbose = verbose

    def run_design(self):
        """
        Run the design of the primers
        """
        min_oh = self.check_overlap()  # Check if the overlap is large enough
        
        print("Designing primer pairs ...")
        primerpairs = self.run_design_primerpair(min_oh)  # store primers: eblock.name : [fw_primer, rv_primer]

        # Add primers to Genbank file of clones
        for mut, eblock in self.eblocks_design_instance.eblocks.items():
            for k, v in primerpairs.items():
                if k == eblock.name:
                    self.add_primers_to_genbank_file(genbank_file=os.path.join(self.output_dir, 'clones', f"{mut.name}", f"{mut.name}.gb"), primer=v[0])
                    self.add_primers_to_genbank_file(genbank_file=os.path.join(self.output_dir, 'clones', f"{mut.name}", f"{mut.name}.gb"), primer=v[1])

        for k, v in primerpairs.items():  # Save primers to fasta file
            self.primers_to_fasta(name=v[0].name, seq=v[0].sequence_5to3, directory=self.output_dir, filename='primers.fasta')
            self.primers_to_fasta(name=v[1].name, seq=v[1].sequence_5to3, directory=self.output_dir, filename='primers.fasta')

        print("Designing sequencing primers ...")
        seqprimers = self.run_design_seqprimer()

        # Save primers to fasta file (both pairs and sequenging primers)
        for k, v in seqprimers.items():
            self.primers_to_fasta(name=v[0].name, seq=v[0].sequence_5to3, directory=self.output_dir, filename='primers.fasta')
            self.primers_to_fasta(name=v[1].name, seq=v[1].sequence_5to3, directory=self.output_dir, filename='primers.fasta')

        for mut, eblock in self.eblocks_design_instance.eblocks.items():
            for k, v in seqprimers.items():
                if k == eblock.name:
                    self.add_primers_to_genbank_file(genbank_file=os.path.join(self.output_dir, 'clones', f"{mut.name}", f"{mut.name}.gb"), primer=v[0])
                    self.add_primers_to_genbank_file(genbank_file=os.path.join(self.output_dir, 'clones', f"{mut.name}", f"{mut.name}.gb"), primer=v[1])
        
        print("Finished designing primers.")

        
    def run_design_seqprimer(self):
        """
        For each eBlock design a Forward and Reverse sequencing primer using primer3
        """
        primers = {}  # store primers: d[eblock.name] = [fw, rv]
        for eblock in self.eblocks_design_instance.wt_eblocks:
            start = self.vector_instance.circular_index(eblock.start_index - self.void_length, self.vector_instance.length)
            end = self.vector_instance.circular_index(start + 2*self.void_length, self.vector_instance.length)
            if start > end:  # change the 0-point of the sequence
                sequence = self.vector_instance.vector.seq[-1000:] + self.vector_instance.vector.seq[0:-1000]
                idx_eblock = sequence.find(eblock.sequence)
                start = self.vector_instance.circular_index(idx_eblock - self.void_length, self.vector_instance.length)
                end = self.vector_instance.circular_index(start + 2*self.void_length, self.vector_instance.length)
            else:
                sequence = self.vector_instance.vector.seq

            # Obtain FW primer
            result  = self.design_sequencing_primer(sequence=sequence, start=start, end=end)
            fw_result = self.parse_primer3_result(result, eblock, type='seq', direction='forward')
            primers[eblock.name] = [fw_result]
            # Obtain RV primer
            start = self.vector_instance.circular_index(eblock.end_index - self.void_length , self.vector_instance.length)
            end = self.vector_instance.circular_index(start + 2*self.void_length, self.vector_instance.length)
            result  = self.design_sequencing_primer(sequence=self.vector_instance.vector.seq, start=start, end=end)

            # Find the closest start index that is higher than the end index
            possible_starts = []
            for i in range(len(result) - 1):
                try:
                    possible_starts.append(int(result['PRIMER_RIGHT'][i]['COORDS'][0]))
                except:
                    break
                
            index = self.find_closest_higher_index(possible_starts, end)
            rv_results = self.parse_primer3_result(result, eblock, type='seq', direction='reverse', index=index)
            primers[eblock.name].append(rv_results)
        return primers
        
    def run_design_primerpair(self, min_oh):
        primers = {}
        for eblock in self.eblocks_design_instance.wt_eblocks:

            mid_index = (eblock.start_index + eblock.end_index) // 2
            length_product = len(self.vector_instance.vector.seq) - len(eblock.sequence)
            len_start = length_product - 75 # TODO Approximation
            len_end = length_product + 75

            sequence_template = self.vector_instance.vector.seq[mid_index:] + self.vector_instance.vector.seq[:mid_index]  # Linearize plasmid to allow for primer design
            half_eblock = len(eblock.sequence) // 2

            fw_end_range = self.vector_instance.circular_index(half_eblock, len(self.vector_instance.vector.seq))
            fw_start_range = self.vector_instance.circular_index(half_eblock - min_oh, len(self.vector_instance.vector.seq))
        
            rv_start_range = self.vector_instance.circular_index(len(self.vector_instance.vector.seq) - half_eblock, len(self.vector_instance.vector.seq))
            rv_end_range = self.vector_instance.circular_index(len(self.vector_instance.vector.seq) - half_eblock + min_oh, len(self.vector_instance.vector.seq))

            result = self.find_primerpair(sequence_template, fw_start_range, fw_end_range, rv_start_range, rv_end_range, len_start, len_end)
            fw, rv = self.parse_primer3_result(result, eblock, type='pair')
            primers[eblock.name] = [fw, rv]
        return primers

    def design_sequencing_primer(self, sequence, start, end):
        try:
            result = primer3.bindings.design_primers(seq_args=Primer3Config.seq_args_sequencing_primer(sequence, start, end),
                                                     global_args=Primer3Config.global_args_sequencing_primer())
            # Check if a sequencing primer was found
            if 'PRIMER_LEFT_0_SEQUENCE' or 'PRIMER_RIGHT_0_SEQUENCE' in result:
                    seqprimer = result['PRIMER_LEFT_0_SEQUENCE']
            return result
                
        except KeyError:  # If no suitable primer is found
            raise ValueError("No sequencing primer found.")
        
    def find_closest_higher_index(self, numbers, target):
        """Find closest number in list that is higher than target"""
        higher_numbers = [num for num in numbers if num > target]
        if not higher_numbers:
            return None
        closest_higher = min(higher_numbers)
        return numbers.index(closest_higher)
    
    def make_fasta_file(self, directory, filename, header=False):
        try:
            with open(os.path.join(directory, filename), 'r') as f:
                pass
        except FileNotFoundError:
            with open(os.path.join(directory, filename), 'w') as f:
                if header:
                    f.write("\n".join(SnapGene.gff3_header(self.vector_instance.vector.seq)))
                    f.write("\n")
                else:
                    pass
    
    def primers_to_fasta(self, name, seq, directory, filename='primers.fasta'):
        """
        This function converts the primers to features that can be read by SnapGene.
        """
        self.make_fasta_file(directory, filename, header=False)
        with open(os.path.join(self.output_dir, filename), 'a') as f:
            f.write(f">{name}\n")
            f.write(f"{seq}\n")

    def parse_primer3_result(self, primer3output, eblock, type, index=0, direction=None,):
        if type == "pair":
            # Forward (left) primer
            start_idx, end_idx = self.vector_instance.find_index_in_vector(str(primer3output[f"PRIMER_LEFT_{index}_SEQUENCE"]))
            fw_primer = Primer()
            fw_primer.name = eblock.name + "_fw"
            fw_primer.sequence_5to3 = primer3output[f"PRIMER_LEFT_{index}_SEQUENCE"]
            fw_primer.is_forward = True
            fw_primer.idx_start = start_idx
            fw_primer.idx_end = end_idx
            fw_primer.tm = round(primer3output[f'PRIMER_LEFT_{index}_TM'], 2)
            fw_primer.gc_content = round(primer3output[f'PRIMER_LEFT_{index}_GC_PERCENT'], 2)
            # Reverse (right) primer
            tmpseq = seq.NucleotideSequence(str(primer3output[f'PRIMER_RIGHT_{index}_SEQUENCE'])).complement().reverse()
            start_idx, end_idx = self.vector_instance.find_index_in_vector(str(tmpseq))
            rv_primer = Primer()
            rv_primer.name = eblock.name + "_rv"
            rv_primer.sequence_5to3 = primer3output[f'PRIMER_RIGHT_{index}_SEQUENCE']
            rv_primer.is_reverse = True
            rv_primer.idx_start = start_idx
            rv_primer.idx_end = end_idx
            rv_primer.tm = round(primer3output[f'PRIMER_RIGHT_{index}_TM'], 2)
            rv_primer.gc_content = round(primer3output[f'PRIMER_RIGHT_{index}_GC_PERCENT'], 2)
            return fw_primer, rv_primer
        elif type == 'seq':
            # TODO Does direction matter here?
            if direction == 'forward':
                start_idx, end_idx = self.vector_instance.find_index_in_vector(str(primer3output[f'PRIMER_LEFT_{index}_SEQUENCE']))
                seq_primer = Primer()
                seq_primer.name = eblock.name + f"_{direction}_seq"
                seq_primer.sequence_5to3 = primer3output[f'PRIMER_LEFT_{index}_SEQUENCE']
                seq_primer.is_forward = True
                seq_primer.idx_start = start_idx
                seq_primer.idx_end = end_idx
                seq_primer.tm = round(primer3output[f'PRIMER_LEFT_0_TM'], 2)
                seq_primer.gc_content = round(primer3output[f'PRIMER_LEFT_{index}_GC_PERCENT'], 2)
                return seq_primer
            elif direction == 'reverse':
                start_idx, end_idx = self.vector_instance.find_index_in_vector(str(primer3output[f'PRIMER_RIGHT_{index}_SEQUENCE']))
                seq_primer = Primer()
                seq_primer.name = eblock.name + f"_{direction}_seq"
                seq_primer.sequence_5to3 = primer3output[f'PRIMER_RIGHT_{index}_SEQUENCE']
                seq_primer.is_reverse = True
                seq_primer.idx_start = start_idx
                seq_primer.idx_end = end_idx
                seq_primer.tm = round(primer3output[f'PRIMER_RIGHT_{index}_TM'], 2)
                seq_primer.gc_content = round(primer3output[f'PRIMER_RIGHT_{index}_GC_PERCENT'], 2)
                return seq_primer
        else:
            raise ValueError("Invalid type. Please specify 'pair' or 'seq'.")
    
    def find_primerpair(self, sequence: str, fw_range_start: int, fw_range_end: int, rv_range_start: int, rv_range_end: int, i_start: int, i_end: int):
        
        left = fw_range_start
        right = rv_range_start 

        while True:
            try:
                result = primer3.bindings.design_primers(seq_args=Primer3Config.seq_args_primer_pair(sequence,
                                                                                                     product_start=i_start,
                                                                                                     product_end=i_end,
                                                                                                     force_left=left,
                                                                                                     force_right=right),
                                                        global_args=Primer3Config.global_args_primer_pair())
                check = result['PRIMER_LEFT_0_SEQUENCE']
                return result

            except KeyError: # Exhaust all possibilities in the desired range
                if left < fw_range_end:
                    left += 1
                elif right < rv_range_end:
                    left = fw_range_start  # Reset left to the start range
                    right += 1
                else:
                    raise ValueError("Primer pair not found.")

    def add_primers_to_genbank_file(self, genbank_file, primer):
        """
        This function saves primer data to an existing GenBank file
        """
        # TODO : Check direction of the primer in benchling and see whether this is correct
        seq_record = SeqIO.read(genbank_file, "genbank")
        feature = None
        site = {}

        if primer.is_forward == True:
            direction = "forward"
            if primer.idx_start < primer.idx_end:
                # Case where primer doesn't overlap with the end of the plasmid
                site = {"start": int(primer.idx_start), "end": int(primer.idx_end), "sequence": str(primer.sequence_5to3)}
                feature_location = FeatureLocation(start=site["start"], end=site["end"])
            else:
                # Case where primer overlaps with the end of the plasmid
                size = self.vector_instance.length - primer.idx_start
                site1 = {"start": int(primer.idx_start), "end": int(self.vector_instance.length), "sequence": seq.NucleotideSequence(primer.sequence_5to3[:size])}
                site2 = {"start": 0, "end": int(primer.idx_end), "sequence": seq.NucleotideSequence(primer.sequence_5to3[size:])}
                feature_location = CompoundLocation([FeatureLocation(start=site1["start"], end=site1["end"]), FeatureLocation(start=site2["start"], end=site2["end"])])
                site["sequence"] = site1["sequence"] + site2["sequence"]

        elif primer.is_reverse == True:
            direction = "reverse"
            if primer.idx_start < primer.idx_end:
                site = {"start": int(primer.idx_start), "end": int(primer.idx_end), "sequence": seq.NucleotideSequence(primer.sequence_5to3).reverse()}
                feature_location = FeatureLocation(start=site["start"], end=site["end"])
            else:
                size = self.vector_instance.length - primer.idx_start
                site1 = {"start": int(primer.idx_start), "end": int(self.vector_instance.length), "sequence": seq.NucleotideSequence(primer.sequence_5to3[:size]).reverse()}
                site2 = {"start": 0, "end": int(primer.idx_end), "sequence": seq.NucleotideSequence(primer.sequence_5to3[size:]).reverse()}
                feature_location = CompoundLocation([FeatureLocation(start=site1["start"], end=site1["end"]), FeatureLocation(start=site2["start"], end=site2["end"])])
                site["sequence"] = site1["sequence"] + site2["sequence"]
        else:
            raise ValueError("Primer is neither forward (fw) nor reverse (rv).")
        # Add primer binding sites as features to the SeqRecord
        feature = SeqFeature(location=feature_location, type="primer_bind")
            # feature_location1 = FeatureLocation(start=site1["start"], end=site1["end"])
            # feature_location2 = FeatureLocation(start=site2["start"], end=site2["end"])
            # feature = SeqFeature(location=feature_location1, type="primer_bind")
            # feature = SeqFeature(location=feature_location2, type="primer_bind")
        feature.qualifiers["note"] = "Primer binding site"
        feature.qualifiers["label"] = primer.name + "_" + direction
        feature.qualifiers["primer_sequence"] = site["sequence"]
        feature.qualifiers["direction"] = direction
        seq_record.features.append(feature)
        # for feature in seq_record.features:
        #     print(feature)
        SeqIO.write(seq_record, genbank_file, "genbank")
        
    def check_overlap(self):
        min_oh = self.eblocks_design_instance.min_overlap - self.minimum_overlap
        if min_oh < 0:
            raise ValueError("Overhang is too small.")
        else:
            return min_oh