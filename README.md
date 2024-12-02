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




