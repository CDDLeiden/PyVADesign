# Design_Gene_Blocks

## Description

eBlocks can be used to clone a set of mutations. 
This program will optimize the design of eBlocks so that the cost will be as low as possible.

For running this program you will need:
- A txt file with mutations, with one mutation per line (e.g L342K, Y425A)
- The gene sequence

The program will give you the sequence of the gene blocks as well PCR primers to open up the destination plasmid.
Additionaly, the program can create files with the eblocks, mutations and primers in a file that can directly be imported into snapgene.

![gene blocks]("https://gitlab.com/rcmkuin/design_gene_blocks/-/blob/main/doc/eblocks.png")

## Installation

For installation of the required packages conda is needed.

First create a conda environment

`conda create -n eblocks python=3.8`

Next, activate the conda environment using

`conda activate eblocks`

Finally, install the packages using

`conda install --file requirements.txt`

## Usage

A folder with example data has been included in the 'example' directory

You can run this example by typing 

`python3 .\src\main.py -i example\example_data\mtb_DnaE1_seq.txt -m example\example_data\mutations.txt -o example\example_output -s example\example_data\snapgene_vector.dna`

This example will search for optimal eblocks for the provided mutations and provides snapgene-readable files.

This script will use the codon usage of _M. smegmatis_, if you wish to use the codon usage of another organism, you should add this to the data directory 

## Tutorial

Describe tutorial


## Mutations input format (.txt) file

This section outlines the format for specifying selected mutations using a text file (.txt). The format ensures clarity and precision in describing single mutations, inserts, deletions, and double mutations.

### Single mutations

For single mutations, each mutation should be specified on a separate line. The format for a single mutation is exemplified below:

G432E
R436Q
I451A
A484S

### Inserts

To specify an insert, start with "Insert," followed by the residue after which the new residues should be inserted, and then list the residues to be inserted. For instance, if you wish to insert 3 residues (PLR) after residue 'A770,' use the following format:

Insert A770-PLR

### Deletions

To specify a deletion, begin with "Deletion," followed by the start and end residues to be deleted. For example, if you want to delete residues from I537 to K562, use the following format:

Deletion I537-K562

### Double mutations

For combined mutants, use "Combined" followed by the mutations you want to combine. If you, for instance, want to combine I447D and A484S in a single eBlock, use the format below:

Combined I477D-A484S

The order of the input does not matter, so an example file could be

G432E
R436Q
Combined I477D-A484S
A484S
Deletion I537-K562
A540L
Insert A770-PLR

Mixing inserts, deletions, single and double mutations


## Codon usage table

DECRIBE HOW TO OBTAIN


