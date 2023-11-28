import os
import pickle
import pandas as pd

# For plasmid viewing
import biotite.sequence as seq
import biotite.sequence.io.genbank as gb
import biotite.sequence.graphics as graphics
import biotite.database.entrez as entrez
import matplotlib.pyplot as plt

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
    "TAT": "Y", "TAC": "Y",
}

def natural_amino_acids():
    return ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', 'R', 'H', 'K', 'D', 'E']


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

def extract_filename(fp):
    base=  os.path.basename(fp)
    fname = os.path.splitext(base)[0]
    return fname

def write_pickle(obj, out_fp, fname="mut_gene_blocks.npy"):
    with open(os.path.join(out_fp, fname), 'wb') as handle:
        pickle.dump(obj, handle)

def load_pickle(fp):
    with open(fp, 'rb') as handle:
        obj = pickle.load(handle)
    return obj



def read_plasmid(dna_fp, organism, gene_name):
    # TODO Read file and parse into (plasmid mapping)
    # TODO Read .DNA file using biopython ans parse into (plasmid mapping)
    # TODO How to Extract the features from a .dna file?
    # TODO Find a way to parse the .dna file (See snapgene output for this)
    annotation = seq.Annotation([
        seq.Feature(
            "source",
            [seq.Location(0, 1500)],
            {"organism": f"{organism}"}
        ),

        # # Ori
        # seq.Feature(
        #     "rep_origin",
        #     [seq.Location(600, 700, seq.Location.Strand.REVERSE)],
        #     {"regulatory_class": "promoter", "note": "MyProm"}
        # ),

        # # Promoter
        # seq.Feature(
        #     "regulatory",
        #     [seq.Location(1000, 1060)],
        #     {"regulatory_class": "promoter", "note": "MyProm"}
        # ),
        # seq.Feature(
        #     "protein_bind",
        #     [seq.Location(1025, 1045)],
        #     {"note": "repr"}
        # ),

        # Gene A
        seq.Feature(
            "regulatory",
            [seq.Location(1070, 1080)],
            {"regulatory_class": "ribosome_binding_site"}
        ),
        seq.Feature(
            "CDS",
            [seq.Location(1091, 1150)],
            {"product": "geneA"}
        ),


        # Terminator
        # seq.Feature(
        #     "regulatory",
        #     [seq.Location(310, 350)],
        #     {"regulatory_class": "terminator", "note": "MyTerm"}
        # ),

        # Primers
        # The labels will be too long to be displayed on the map
        # If you want to display them nevertheless, set the
        # 'omit_oversized_labels' to False
        # seq.Feature(
        #     "primer_bind",
        #     [seq.Location(1385, 1405)],
        #     {"note": "geneC"}
        # ),
        # seq.Feature(
        #     "primer_bind",
        #     [seq.Location(345, 365, seq.Location.Strand.REVERSE)],
        #     {"note": "geneC_R"}
        # ),

        # # Terminator
        # seq.Feature(
        #     "regulatory",
        #     [seq.Location(310, 350)],
        #     {"regulatory_class": "terminator", "note": "MyTerm"}
        # ),
    ])


    fig = plt.figure(figsize=(8.0, 8.0))
    ax = fig.add_subplot(111, projection="polar")
    graphics.plot_plasmid_map(
        ax, 
        annotation, 
        plasmid_size=1500, 
        label="My plasmid",
        label_properties={"fontsize": 8})

    ticks = ax.get_xticks()
    labels = ax.get_xticklabels()

    fig.tight_layout()
    # plt.show()
    return fig, ax



# def log_to_file_and_console(logfile, *messages):
#     with open(logfile, 'a') as log_file:
#         for message in messages:
#             # Write each message to the log file
#             log_file.write(message)
#     # Print all messages to the console
#     print(*messages, end='')

# def create_or_clear_file(file_path):
#     try:
#         # Check if the file exists
#         with open(file_path, 'r') as file:
#             # Close the file
#             pass
#     except FileNotFoundError:
#         # If the file doesn't exist, create it
#         with open(file_path, 'w'):
#             pass