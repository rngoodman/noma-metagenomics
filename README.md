# Noma Metagenomics

This repository includes the code used to run the analysis in a recent study currently posted as a preprint on bioRxiv:

[**Shotgun metagenomic analysis of the oral microbiomes of children with noma reveals a novel disease-associated organism**](https://doi.org/10.1101/2025.06.24.661267)

Michael Olaleye, Angus M O'Ferrall, Richard N. Goodman, Deogracia W Kabila, Miriam Peters, Gregoire Falq, Joseph Samuel, Donal Doyle, Diana Gomez, Gbemisola Oloruntuyi, Shafiu Isah, Adeniyi S Adetunji, Elise N. Farley, Nicholas J Evans, Mark Sherlock, Adam P. Roberts, Mohana Amirtharajah, Stuart Ainsworth

*bioRxiv* 2025.06.24.661267; doi: [https://doi.org/10.1101/2025.06.24.661267](https://doi.org/10.1101/2025.06.24.661267)

The code used in the analysis is linked below for the **taxonomic based metagenomic analysis** of noma samples and healthy dataset (Table 1, Figures 1-3 and Supplementary Figures S1-S4), the **recovery and analysis of *Treponema* MAGs from noma samples** (Figure 4) and  **AMR profiling of noma metagenomes** (Figure 5).

## Part 1 - Taxonomic based metagenomic analysis

### Part 1.1 - [R markdown - Within noma dataset analysis (Table 1B, Figure 1A, Figure S1)](https://rngoodman.github.io/noma-metagenomics/code/Noma_swab_vs_saliva.html)
* Relative abundance
* non-parametric statistical tests

### Part 1.2 - [R markdown - Within healthy dataset analysis](https://rngoodman.github.io/noma-metagenomics/code/Healthy_vs_healthy.html)
* non-parametric statistical tests 

### Part 1.3 - [R markdown - Noma vs healthy dataset analysis (Table 1A, Figures 1B, 2 and 3, Figures S2-S4)](https://rngoodman.github.io/noma-metagenomics/code/Noma_vs_healthy.html)
* non-parametric statistical tests 
* Relative abundance
* Differential analysis
* Machine learning and multivariate statistical analyses

## Part 2 - Recovery and analysis of *Treponema* MAGs from noma samples

### Part 2.1 - [R script - Descriptive plot of *Treponema* MAGs (Figure 4A)](https://github.com/rngoodman/noma-metagenomics/blob/main/code/Treponema_MAGs.R)

### Part 2.2 - [R script - Tree of high-quality *Treponema* MAGs in context of previously characterised *Treponema* RefSeq genomes (Figure 4B)](https://github.com/rngoodman/noma-metagenomics/blob/main/code/tree.R)

### Part 2.3 - [Python script - ANI matrices for novel *Treponema* MAGs (Figure 4C)](https://github.com/rngoodman/noma-metagenomics/blob/main/code/ANI_visualisation_multi.py)

## Part 3 - AMR profiling of noma metagenomes

### Part 3.1 - [R script - AMR plots (Figure 5)](https://github.com/rngoodman/noma-metagenomics/blob/main/code/ARGs.R)
