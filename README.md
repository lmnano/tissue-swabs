# Comparison of RAD Loci generated from minimally invasive swab samples and invasive tissue samples

## Intro

This is a project about comparing tissue and skin swab DNA samples from Proteus anguinus.

## Scripts

### Work in progress

First script () is for making a blast databases from catalogs and genome to use as reference when blasting.

Following this are the scripts for blasting catalogs against different databases (tissue, genome). There were three different options used for this BLAST (), local Bowtie2 () and global Bowtie2 ().

The results that matched were run through an R script.

The results that didn't match were blasted against the BLAST nt database (). This was done only for the best result of the three options (BLAST, local and global Bowtie2). These result were transformed into the right format (BLAST outfmt 0), and used as MEGAN6 input.



### Figures

The output from MEGAN6 was used as the input for making figures.

In the R script scripts/tissue_swab_result_analysisFINAL.R is code used for producing certain tables and figures in the corresponding paper.

The first part of the script loads data and produces a table (table S2) with proportions of reads mapped with different alogirhms and in different combinations, and two figures with proportions of matches.

The second part produces pie charts used for making figure 5, with proportions of different taxa in exogenous DNA, and data for a table (table S5) with top ten most commonly found species in the exogenous DNA with Blast.



