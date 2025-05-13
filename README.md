# PyVADesign

## Description

The generation and analysis of diverse mutants of a protein is a powerful tool to understand the functioning of a protein. However, the generation of such mutants can be time-consuming while the commercial option of buying a series of mutant plasmids can be expensive. In contrast, the insertion of a synthesized double stranded DNA (dsDNA) fragment into a plasmid is a fast and low-cost method to generate a large library of mutants with one or more point mutations, insertions, or deletions. To aid in the design of these DNA fragments present PyVADesign: a Python package that makes the design and ordering of dsDNA fragments straight-forward and cost-effective. In PyVADesign, the mutations of interest are clustered in different cloning groups for efficient exchange into the target plasmid. Additionally, primers that prepare the target plasmid for insertion of the dsDNA fragment, as well as primers for sequencing, are automatically designed within the same program.

## workflow

User’s input includes a list of desired mutations and a plasmid sequence with gene of interest (blue). PyVADesign clusters the mutations on the gene of interest without explicit user input. Next, the dsDNA fragments and primers are designed to prepare the plasmid for dsDNA fragment insertion (shown in green, magenta and red colours). To linearize the input plasmid for insertion of the dsDNA fragment, PCR reactions are performed using the PyVADesign primers. Finally, the dsDNA fragment can be inserted in the linearized plasmid in a single step.

![overview_v3](https://github.com/user-attachments/assets/77969f38-d03a-4f5e-886c-56776c49b8c9)

## Installation Instructions

### Cloning the Repository

Clone the repository using the following command:

```bash
git clone https://github.com/CDDLeiden/PyVADesign.git
```

### Creating a Conda Environment

To create a Conda environment using the provided ```environment.yaml``` file, follow these steps:

Navigate to the project directory

```bash
cd PyVADesign
```

Create the environment from the ```environment.yaml``` file:

```bash
conda env create -f environment.yaml
```

After the environment is created, you can activate it using:

```bash
conda activate PyVADesign
```

## Usage

To quickly get started with PyVADesign and test its functionality, you can run the provided test files. This will allow you to verify that everything is working as expected without needing to configure your own data.

```python
python pyvadesign.py --gene tutorial-data/A0A0H2ZYQ2.fasta --vector tutorial-data/pACE_mmpL3-Mav.dna --mutations tutorial-data/mutations.txt --output tutorial-output/amount-optimization
```

The ```--gene``` specifies the gene sequence in .fasta format, ```--vector``` specifies the vector sequence with the gene already cloned into it in .dna format. ```--output``` is the output directory and ```--mutations``` describes the mutations to be made.

The output directory ```tutorial-output/amount-optimization``` will contain a ```clones``` directory that contains a clone in GenBank format of each variant. A file named ```primers.fasta``` contains primers to open-up the vector (fw/rv-DNABlock-X) as well as sequencing primers (seq-fw/rv-DNABlock-X).

Also see a more detailed [tutorial]() for the example of how to use PyVADesign

## License

This project is licensed under the [GNU General Public License v3 (GPL-3.0)](LICENSE), and includes code from the following projects:

- [DnaFeaturesViewer](https://edinburgh-genome-foundry.github.io/DnaFeaturesViewer/) (MIT License)
- [Primer3](https://github.com/primer3-org/primer3) (GPL License)

### GNU General Public License (GPL)
The full text of the GPL-3.0 License is available in the `LICENSE` file or can be viewed here: [https://www.gnu.org/licenses/gpl-3.0.html](https://www.gnu.org/licenses/gpl-3.0.html).

### MIT License
The full text of the MIT License is also included below:

MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.




