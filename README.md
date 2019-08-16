# Oxamniquine resistance alleles are widespread in Old World *Schistosoma mansoni* and predate drug deployment

## Abstract
Do mutations required for adaptation occur *de novo*, or are they segregating within populations as standing genetic variation? This question is key to understanding adaptive change in nature, and has important practical consequences for the evolution of drug resistance. We provide evidence that alleles conferring resistance to oxamniquine (OXA), an antischistosomal drug, are widespread in natural parasite populations under minimal drug pressure and predate OXA deployment. OXA has been used since
the 1970s to treat *Schistosoma mansoni* infections in the New World where *S. mansoni* established during the slave trade. Recessive loss-of-function mutations within a parasite sulfotransferase (SmSULT-OR) underlie resistance, and several verified resistance mutations, including a deletion (p.E142del), have been identified in the New World. Here we investigate sequence variation in *SmSULT-OR* in *S. mansoni* from the Old World, where OXA has seen minimal usage. We sequenced exomes of 204 *S.
mansoni* parasites from West Africa, East Africa and the Middle East, and scored variants in *SmSULT-OR* and flanking regions. We identified 39 non-synonymous SNPs, 4 deletions, 1 duplication and 1 premature stop codon in the *SmSULT-OR* coding sequence, including one confirmed resistance deletion (p.E142del). We expressed recombinant proteins and used an in vitro OXA activation assay to functionally validate the OXA-resistance phenotype for four predicted OXA-resistance mutations. Three
aspects of the data are of particular interest: (i) segregating OXA-resistance alleles are widespread in Old World populations (4.29 – 14.91% frequency), despite minimal OXA usage, (ii) two OXA-resistance mutations (p.W120R, p.N171IfsX28) are particularly common (>5%) in East African and Middle-Eastern populations, (iii) the p.E142del allele has identical flanking SNPs in both West Africa and Puerto Rico, suggesting that parasites bearing this allele colonized the New World during the slave
trade and therefore predate OXA deployment. We conclude that standing variation for OXA resistance is widespread in *S. mansoni*.

Full article is available [here](https://doi.org/10.1101/657056)

## Code

The code used to generate and analyze the data presented in the mansucript is available in this repository.

The code is divided in different sections:
* [1-Alignment](1-Alignment/README.md): Scripts in this section were used to align the fastq data against the reference genome.
* [2-VCF processing](2-VCF%20processing/README.md): Scripts in this section will download and process (filtering and phasing) the raw VCF file containing variants of the 3 first Mb of the chromosome 6 of *S. mansoni*.
* [3-Analyses](3-Analyses/README.md): Scripts in this section will analyze the information from the phased variants, generate haplotype sequences for genetic analysis and run simulation.


## Prerequisites

To run the scripts correctly, the *S. mansoni* genome and the corresponding annotation file are needed:
* [v5 *S. mansoni* reference genome](ftp://ftp.sanger.ac.uk/pub/pathogens/Schistosoma/mansoni/Archive/S.mansoni/genome/Assembly-v5/sma_v5.0.chr.fa.gz) (default location set in scripts: `$HOME/data/sm_genome/sma_v5.0.chr.fa`)
* [v5 *S. mansoni* genome annotation](ftp://ftp.sanger.ac.uk/pub/pathogens/Schistosoma/mansoni/Archive/S.mansoni/genome/Gene_models/ARCHIVE/v5.19.05.11.chado.raw.gff) (default location set in scripts: `$HOME/data/sm_Gene_table/v5.07.08.12.chado.raw.gff`)

Other dependencies are specified in the readme file of each section.
