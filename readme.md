CUBACR
======
Scripts for the codon usage bias analysis of Conserved and Variable regions from core genome proteins.  
CUBACR scripts are programmed to process multifasta files containing ortholog genes of the desired core-gene set (or corresponding to proteins with a defined expression level), make an amino acid guided codon alignment for each gene (using TranslatorX software), and to output the conserved and variable regions of all genes on the core-genome set.  
This CR (conserved) and VR (variable) sequences can then be used for:  

1. Calculation of modal codon usage frequencies (using G. Olsen software) for CR and VR.
2. Generate CR and VR sequence concatenates.
3. Correspondence analysis of RSCU for genes and modal frequencies, including CR and VR.
4. Calculate bootstrapped sequences using G. Olsen software.
4. Generate bootstrapped sequences for concatenated CR and VR, and calculate G. Olsen distance (all vs. all)
5. Calculate a neighbor joining tree using G. Olsen `freqs_2_nj_tree_linux` script (using Type 2 distance), and generate heatmap of modal codon usage frequences guided by the previously generated tree.
6. Heatmaps of the difference of modal codon usage between Highly/lowly expressed and Conserved/Variable sequences.
7. Calculate GC3, s-tAI, and Nc plots for all the genes and modal sequences.

# Installation
In order for the scripts to run the following programs must be installed on the system, and included in the linux $PATH.

## Requirements
- [TranslatorX](http://www.translatorx.co.uk/) (included with the scripts). Alignment software must be installed in order for TranslatorX to work. Muscle, clustal, mafft, t-coffee can be installed from Ubuntu (Linux) repository (sudo apt install ...).
- The software needs G. Olsen software which can be found in [link](http://www.life.illinois.edu/gary/programs/codon_usage.html). The installation instructions are clearly provided by the author. The software must be in the linux $PATH to work correctly. Additionally, in order to run in ubuntu some lines of the code must be modified (See CUBES page. Modified files are provided with CUBES)  
- a local installation of [codonw](http://codonw.sourceforge.net/).
- [phylip](http://evolution.genetics.washington.edu/phylip.html) package for the tRNA tree generation. Can be installed from the linux repository.
- [R-project] (https://www.r-project.org/) software with the following libraries: stringr, plyr, tidyr, ggplot2, ggrepel, ggthemes, ggpmisc, ggsignif, ggpubr, data.table, naturalsort, tidyr, agricolae, pheatmap, ape, tAI, doParallel, doRNG, ade4, phytools. 
  
The main script of this package is **cCvsE.sh**

# cCvsE.sh
To display usage instructions type: cCvsE.sh -h  
```
Usage: cCvsE.sh [-o] [-c] [-h] [-i] [-b] [-m] [-n] [-a] [-t] [-v]

        -o Skips G.Olsen Bootstrapped distances calculation
        -c Skips AA guided codon alingment calculation
        -h Show this Help
        -i n for bootstrap of the optimization algorithm (50 default)
        -b Skips Bootstrapped codon distances calculation
        -m Skips modal sequences calculations
        -n Skips CR and VR modal sequences calculations
        -a Skips CA analysis
        -t Skips Conservation table generation
        -v Defines variability cutoff, must be 0<v<1

```
Generates an amino acid guided codon alignment and calcualtes codon usage for conserved protein sequences, both on higly and lowly expressed proteins.

Then calculates the distance between modal and bootstraped sequences.


## Input
The input for cCvsE.sh requires the following folder structure:
From the root, which has all the scripts
```
/─┬─SP1_folder/┬─ CC/
  │            ├─ HEP/┬─ locus_tag1.fa
  │            │      ├─ locus_tag2.fa
  │            │      ¦
  │            │      ├─ locus_tagi.fa
  │            │      ¦
  │            │      ¦
  │            │      └─ locus_tagN.fa
  │            ├─ LEP/┬─ locus_tag1'.fa
  │            │      ├─ locus_tag2'.fa
  │            │      ¦
  │            │      ├─ locus_tagi'.fa
  │            │      ¦
  │            │      ¦
  │            │      └─ locus_tagN'.fa
  │            ├─ PHE/
  │            └─ trna.txt, s_opts_DCBS_GNM.txt
  │
  │
  ├─SP2_folder/
  ¦
  ¦
  ├─SPi_folder/
  ¦
  ¦
  ├─SP4_folder/
  ¦
  ¦
```

`trna.txt` -> relative frequence of tRNA gene copies, as described in CUBES readme.md.  
`s_opts_DCBS_GNM.txt` -> Sij obtained with CUBES using all the genome.

Content of the folders:

* HEP and LEP folders:
  `locus_tagi.fa`: these files are multifasta files which will be used by TranslatorX for the amino acid based codon alignment. The first ortholog sequence in the file, must be the corresponding to SPi, the representative specie of the order/family. The requirement for the file to be named with the locus tag can be fulfilled by running the `bulk_rename.sh` script, on the `SPi_folder`. This is required because after the alignment the sequence order in the multifasta can be changed.
* CC: this folder must contain C1.fa -> Cn.fa (multifasta files, each with all the CDS sequences of the core-genome (C1->n) of representative specie -SPi- of the set), GNM.fa (genome sequence of SPi) and single.fa (singletons of the set). Folder names are important, since from them the species names will be taken for some of the plots.
* PHE: this folder must contain a PHE.fa file with the codin sequences of the putative/predicted highly expressed proteins (usually ribosomal proteins, tRNA ligases, etc. It can also be user defined)

## Instructions
Just run the script!
./cCvsE.sh -o -b -v 0.5

The distance calculation takes time, so if you don't need it, just skip it with -o -b.
-v 0.5 command tells the script to define variable positions as the positions in which the mean amino acid frequency (not counting gaps and over the total amino acids types present) is less than 0.5
- At least one fo HEP or LEP folders must exist. HEP and LEP are for Highly or Lowly expressed proteins, however, if we only want to analyze conservation (there is no expression data), the multifasta files must be located on the HEP folder, regardless of the expression level.

## Relevant Outputs
* Sequences representing the modal frequency of highly/lowly expressed conserved and variable regions.
  `SPi_HEP_CR.modal_freq.fas, SPi_HEP_VR.modal_freq.fas, SPi_LEP_CR.modal_freq.fas, SPi_LEP_VR.modal_freq.fas`.
* Conservation table and plot
  `SPi_CE.table`, a table with the proportion of conserved and variable regions for each gene.
  `coverage_plot.svg`, a density plot of the length of conserved and variable regions for HEP and LEP.
* Correspondence analysis plots and tables
  `CA_genes.svg, CA_codons.svg`.

# GC3.sh
```
Usage: GC3.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help
```
This script scans for singletons, HEP (VR - CR), LEP (VR - CR) and PHE fasta files, and computes the GC percent of the third (GC3) codon base. Generates a GC3 vs gene set plot.

## Input
Requires the outputs of cCvsE.sh 

## Outputs

`GC3_CR_VR.svg` Bar plot of GC3 for the modal sequences of singletons, HEP (VR - CR), LEP (VR - CR) and PHE.

# tAi_Modal.sh
```
Usage: tAi_Modal.sh [-y] [-h]

        -y Skips command confirmation
        -h Show this Help

```
This script scans all the child directories for C1->n, singleton,  HEP (VR - CR), LEP (VR - CR), PHE, GNM .modal_freq files and computes de adaptation index s-tAi. Generates a tAi vs Distance plot.

## Input
Requires the outputs of cCvsE.sh 
Requires the following files present on the SPi_folder:
`trna.txt` -> relative frequency of tRNA gene copies, as used in CUBES.
`s_opts_DCBS_GNM.txt` -> Sij obtained with CUBES using all the genome.

## Output

`tAi_CR_VR.svg` Bar plot of s-tAI for the modal sequences of singletons, HEP (VR - CR), LEP (VR - CR) and PHE.

# calculate_dist_tree_heatmap.sh
```
Usage: calculate_dist_tree_heatmap.sh [-h]]

        -h Show this Help
```
Generates a NJ distances tree and heatmap of CUF for Ci, Singletons, and CR, VR, HEP and LEP sets.


## Input
Requires the outputs of cCvsE.sh 

## Output
`tree.svg, tree_heetmap.svg` a neighbor joining tree calculated using G. Olsen `freqs_2_nj_tree_linux` script, and a heatmap of modal codon usage frequences that uses the NJ tree to order the rows.


# Nc_plots.sh 
```
Usage: Nc_plots.sh [-h]

        -h Show this Help
```
This script calculates Nc using codonw, and generates different plots and tables.
 
## Input
Requires the outputs of cCvsE.sh 

## Output
`modal_seqs_NC_GC3s.tab, NC_GC3s.tab` -> Tables with Nc information for all the genes in the genome, and for modal sequences of HEP/LEP CR/VR, singletons, GNM and C1->n sets.
`Nc_genes.svg, Nc_genes_vs_GC3s.svg, Nc_modal.svg, Nc_modal_vs_GC3s.svg` Bar Plots of Nc for the different sequences, and the corresponding Nc plots.


# CRdelta_Heatmap.sh
```
Usage: CRdelta_heatmap.sh [-h]

        -h Show this Help
```
Calculates RSCU difference between HEP-CR, HEP and PHE, and plots different heatmaps.


## Input
Requires the outputs of cCvsE.sh 
Requires the outputs of Nc_plots.sh 

## Output
The output of this script is a set of heatmaps with the difference of modal codon usage frequency between HEP-CR/VR, LEP-CR/VR, HEP and LEP. The interpretation of the heatmap can be useful to discover the weight of efficiency and accuracy on the choice of codons for each set of genes.
Files: `delta-heatmap.svg, delta-heatmap2.svg, delta-heatmap3.svg`.
To integrate the results for all species analyzed, the `make_heatmap_fig.sh` can be run.
The output of this script are three Heatmap(n).svg files.
