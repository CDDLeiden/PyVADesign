import os
import pickle
import pandas as pd

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

def log_to_file_and_console(logfile, *messages):
    with open(logfile, 'a') as log_file:
        for message in messages:
            # Write each message to the log file
            log_file.write(message)
    # Print all messages to the console
    print(*messages, end='')

def create_or_clear_file(file_path):
    try:
        # Check if the file exists
        with open(file_path, 'r') as file:
            # Close the file
            pass
    except FileNotFoundError:
        # If the file doesn't exist, create it
        with open(file_path, 'w'):
            pass