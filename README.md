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
