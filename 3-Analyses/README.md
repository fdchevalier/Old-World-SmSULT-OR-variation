## 1-SmSULT-OR haplotypes

This folder contains the data and scripts to generate the haplotype sequences for the genetic analysis. The fasta file with haplotype sequences that will be generated can be loaded in [DnaSP](http://www.ub.edu/dnasp/). Alternatively the nexus file (with the aligned sequences and population information) can be used with DnaSP.

### Prerequisites

To run the scripts properly, the following tools are needed:
* [Clustal Omega](http://www.clustal.org/omega/)
* [SAMtools](http://www.htslib.org/)
* the [EMBOSS package](https://github.com/kimrutherford/EMBOSS) 

### Usage

If needed, run the following commands to generate the haplotype sequences (steps are detailed in the script):
```bash
cd "1-SmSULT-OR haplotypes"
./hap-processing.sh
cd ..
```


## 2-Prediction of impact of mutations

This folder contains the data files, scripts and results of the prediction from the Rosetta package.

### Prerequisites

To run the prediction script properly, the following tools are needed:
* [Openbabel](http://openbabel.org/wiki/Main_Page)
* the [Rosetta package](https://www.rosettacommons.org/)

### Usage

To run the prediction (be aware that this step takes a while), run the following commands (steps are detailed in the script):
```bash
cd "2-Prediction of impact of mutations/ddg_analysis"
./mutation-prediction.sh
cd ../..
```

To generate the graph:
```bash
cd "2-Prediction of impact of mutations"
./ddg_plot.R
cd ..
```


## 3-Resistance allele analysis

This folder contains the scripts to perform the haplotype analysis of *SmSULT-OR* and flanking regions as well as the analysis of OXA-R allele frequency.


### Prerequisites

To run the `OXA-R_map.R` script properly, the following data are needed:
* the [small scale cultural data](https://www.naturalearthdata.com/downloads/110m-cultural-vectors/) from the Natural Earth website
* the [small scale physical data](https://www.naturalearthdata.com/downloads/110m-physical-vectors/) from the Natural Earth website

To avoid modifying the script, the data needs to be reached from this path: `$HOME/data/natural_earth_data/`.

### Usage

To run the haplotype analysis, run the following commands:
```bash
cd "3-Resistance allele analysis"
./Haplotype_seq.R
cd ..
```

To run the analysis on frequency of OXA-R alleles and worms, run the following commands:
```bash
cd "3-Resistance allele analysis"
./OXA-R_map.R
cd ..
```


## 4-Selection simulation

This folder contains the data files and script to perform the selection simulation and plot the results.

To run the simulation, run the following commands (steps are detailed in the script):
```bash
cd "4-Selection simulation"
#rm -R data/    # Uncomment this step to rerun the simulation instead of only plotting the graph
./Selection.R
cd ..
```

## R prerequisites

To run the R scripts, specific packages are needed. To install them, run the following command:
```bash
R <<EOF
install.packages(c("ape", "doParallel", "gridExtra", "igraph", "magrittr", "mapplots", "parallel", "phangorn", "plotrix", "plyr", "poppr", "RColorBrewer", "reshape", "rgdal", "seqinr", "strataG", "vcfR", "viridis")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("rtracklayer"))
EOF
```
