import pandas as pd

DNA_Codons = {
    # 'M' - START, '_' - STOP
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
    # "TAA": "_", "TAG": "_", "TGA": "_"
}

def read_codon_usage(fp="codon_usage_smegmatis.csv"):
    # Obtained from http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=246196
    codon_usage = {}
    codon_usage_no_U = {}
    df = pd.read_csv(fp, sep=';')
    for idx, row in df.iterrows():
        codon_usage[row['Triplet']] = row['frequency_number']
    for key, value in codon_usage.items():
        newkey = key.replace('U', 'T').lower()
        codon_usage_no_U[newkey] = value
    return codon_usage_no_U