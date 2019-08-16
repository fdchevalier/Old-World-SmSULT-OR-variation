# Alignment

The code available in this directory performed the alignement steps from the fastq files to the final BAM files used to perform variant calling.

## Prerequisites

To run the scripts properly, the following tools are needed:
* [bwa](http://bio-bwa.sourceforge.net/)
* [SAMtools](http://www.htslib.org/)
* [GATK](https://software.broadinstitute.org/gatk/download/)
* [picard](https://broadinstitute.github.io/picard/)
* [BEDtools](https://bedtools.readthedocs.io/en/latest/)

The subfolders from `exome_pipeline` need to be in the system `PATH`. To include them, use this command: `export PATH=$PATH:$PWD/exome_pipeline:$PWD/exome_pipeline/exome_scripts:$PWD/exome_pipeline/lib_stat_scripts`.

Paired fastq files of a given library must be in a dedicated library folder. All library folders must be in a data folder.

## Usage

Run the `exome_pipeline/all_exome_pipeline.sh -h` from this folder to get directions. This script is designed to run on a [Sun Grid Engine](https://en.wikipedia.org/wiki/Oracle_Grid_Engine) cluster. It launches the `exome_pipeline/exome_pipeline.sh` for each library on a node and generate library statistics at the end.

After alignment, to obtain information regarding read depth and coverage of exons of *SmSULT-OR*, run the `SmSULT-OR locus stats.sh` script from the `SmSULT-OR locus` folder. This script (along with the R script) must be at the same level as the data folder that contains subdirectories of each library. If the data folder is not named `data`, both scripts in this folder need to be edited.
