# VCF processing

The code available in this directory will download the raw VCF file from [Zenodo](https://doi.org/10.5281/zenodo.2850876) and performed the different steps described in the manuscript to generate a phased VCF used in the subsequent analyzes. This will produce:
* a filtered and phased VCF,
* a TSV file with the information related to the scored variants of *SmSULT-OR*. This file was later hand edited to produce the first supplementary table of the manuscript.

## Prerequisites

To run the script properly, you will need the following tools:
* [BEDtools](https://bedtools.readthedocs.io/en/latest/)
* [VCFtools](https://vcftools.github.io/index.html)
* the [EMBOSS package](https://github.com/kimrutherford/EMBOSS)
* [vt](https://github.com/atks/vt)
* java
* [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html)

You also probably need to update the `vcf-processing.sh` script with your custom path to:
* the v5 genome of *Schistosoma mansoni*,
* the Beagle jar file.

## Usage

Run the `./vcf-processing.sh` from this folder. No other script needs to be run.
