import os
import pickle
import pandas as pd
import matplotlib.pyplot as plt

# For plasmid viewing
import biotite.sequence as seq
import biotite.sequence.io.genbank as gb
import biotite.sequence.graphics as graphics
import biotite.database.entrez as entrez
from biotite.sequence import Feature, Location, Annotation


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


def extract_snapgene_features(vector, organism=None):
    """
    Extract features from snapgene file and return them as a list of biotite sequence features.

    Parameters
    ----------
    vector : biotite sequence object
        The vector sequence
    organism : str
        The source organism name (optional)
    """
    vector_features = []
    for i in vector.features:
        feature_type = i.type
        quals = i.qualifiers
        label_qualifier = i.qualifiers.get('label', '')
        locations = [Location(int(i.location.start), int(i.location.end))]

        if label_qualifier:
            qual = {"label": ' '.join(label_qualifier)}
        else:
            qual = {}

        if (i.type == "rep_origin") or (i.type == "protein_bind") or (i.type == "terminator"):
            if quals['note']:
                qual = {"note": ' '.join(quals['label'])}

        elif (i.type == "CDS") or (i.type == "gene"):
            if quals['label']:
                qual = {"product": ' '.join(quals['label'])}
        elif (i.type == "regulatory"):
            if quals['label']:
                qual = {"label": ' '.join(quals['label'])}
        elif (i.type == "misc_feature"):
            if quals['label']:
                qual = {"note": ' '.join(quals['label'])}

        vector_features.append(Feature(feature_type, locations, qual))
    # Add organism source
    if organism:
            vector_features.append(Feature("source", [Location(0, len(vector.seq))], {"organism": f"{organism}"}))
    return vector_features


def plot_vector(vector, features, figsize=(8,8), fontsize=10):
    """
    Plot a plasmid map of the vector sequence.
    """
    annotation = Annotation(features)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="polar")
    graphics.plot_plasmid_map(
        ax, 
        annotation, 
        plasmid_size=len(vector.seq), 
        label=f"{vector.name}",
        label_properties={"fontsize": fontsize})
    ticks = ax.get_xticks()
    # Only show a few ticks
    ax.set_xticks(ticks[::2])
    fig.tight_layout()
    return ax, fig