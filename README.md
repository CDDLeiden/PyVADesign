# PyVADesign

## Description
PyVADesign does the design of dsDNA fragments and primers for the introduction of mutations in a gene of interest

## Installation Instructions

### Cloning the Repository

Clone the repository using the following command:

```bash
   git clone https://github.com/CDDLeiden/PyVADesign
```

### Creating a Conda Environment

To create a Conda environment using the provided ```environment.yaml``` file, follow these steps:

Create the environment from the ```environment.yaml``` file:

```bash
   conda env create -f environment.yaml
```

Activate the Environment: After the environment is created, you can activate it using:

```bash
   conda activate PyVADesign
```

## Usage

After repository cloning, you should be able to run:

```python
   python CLI.py -g tutorial-data/A0A0H2ZYQ2.fasta -v /tutorial-data/pACE_mmpL3-Mav.dna -m tutorial-data/mutations.txt -o tutorial-output 
```

The ```-g``` specifies the gene sequence, ```-v``` specifies the vector sequence with the gene already cloned into it. ```-o``` is the output directory and ```-m``` describes the mutations to be made.

The output directory ```tutorial-output``` will contain a ```clones``` directory that contains a clone in GenBank format of each mutant. A file named ```primers.fasta``` contains primers to open-up the destination plasmid (fw/rv-eBlock-X) as well as sequencing primers (seq-fw/rv-eBlock-X).

Also see a more detailed [tutorial]() for the example of how to use PyVADesign


