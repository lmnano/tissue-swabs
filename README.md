# Comparison of RAD Loci generated from minimally invasive swab samples and invasive tissue samples

## Intro

This is a project about comparing tissue and skin swab DNA samples from Proteus anguinus.

## Scripts

### Work in progress, scripts not yet uploaded

First scripts (01_makeBlastDB.sh, 01_bowtie2build.sh) are for making a blast databases from catalogs and genome to use as reference when blasting and building a bowtie2 index, respectively.

Following this are the scripts for blasting catalogs against different databases (tissue, genome) and bowte2 indexes. There were three different options used for this BLAST (02_blastAll.sh), local Bowtie2 (02_bowtie2Local.sh) and global Bowtie2 (02_bowtie2Global.sh). The Bowtie2 scripts are followed by the 03_sam2bam.sh script which transforms the sam output from bowtie2 to (sorted) bam and adds bai index file. The 04_run_tissue_swab_top_hits.R.sh filters the blast results.

These results are analysed with the 05_tissue_swab_result_analysisFINAL2.Rmd R markdown script. This script is written in a way that it takes all combinatios of results from Blast and Bowtie2: swabs against tissue and vice versa (AllSvT and AllTvS), split swab sequences and split tissue sequences, and vice versa (SplitSvT and SplitTvS), split swab against tissue (SplitAllSvT), split tissue against swab (SplitAllTvS). The abbriveations are used in scripts to describe different combinations. The tables produced in the script produce more data than used in the paper. 

The results that didn't match were blasted against the BLAST nt database (06_blastBt2_unalignedNT.sh). This was done only for the local Bowtie2 results, which were considered the best for the swabs vs. tissue combination (AllSvT). These result were transformed into the right format (BLAST outfmt 0), and used as MEGAN6 input.

Unmapped results for individual samples were also blasted against NT database (07_individual_blastBt2_unalignedNT.sh)

Output from MEGAN6 was used to determine ratios of different taxa. For whole catalogs this was done using the 09_tissue_swabs_taxa_analysisFINAL2.Rmd markdown and for the individual samples the tables used for producing the figures were made with 10_tissue_swabs_indidvidual_taxa_analysis_MeganFINAL.Rmd, the latter markdown uses the megan_individual_taxa_df.R script to work. Script 08_tissue_swabs_taxa_analysis_megan_resultsFINAL.Rmd was used to produce the top 10 most commons species table directly from the blast results, not MEGAN6 output. 
