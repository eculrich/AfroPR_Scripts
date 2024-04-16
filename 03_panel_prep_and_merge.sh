#!/bin/bash

# Script written by E. Ulrich on 2.27.24 adapted from scripts previously written by E. Ulrich and M. Nieves-Colón
    # E. Ulrich: ulriche227@gmail.com
    # M. Nieves-Colón: mnievesc@umn.edu

# Script used to QC filter and prepare genomic plink files and then merge with AfroPR sample plink data
# Input plink files are assumed to be in the format PANEL.fam, PANEL.bim, PANEL.bed. 
    # The names "PANEL" and "SAMPLE" are used throughout this script as placeholders for a provided panel/sample population file name. 
    # Files are assumed to be in the same directory as this script

# Software used: PLINK 1.9, Python 3.6.8, snpflip 0.0.6

################
### PANEL QC ###
################

# for use on an HPC/supercomputer
module load plink

# 1. Make a list of populations to keep from the panel. This step should also filter out all related individuals if necessary
# Specific steps will vary depending on sample info format and the criteria used to include/exclude individuals
# Include individuals to keep in "keeplist" file, one individual per row, with columns for family ID and within-family ID

# 2. Filter the input plink files to keep only individuals from desired populations
plink --bfile PANEL --allow-extra-chr --keep keeplist --out PANEL.sex --make-bed

# 3. For analyses using only autosomal data, limit panel to only chromosomes 1-22
plink --bfile PANEL.sex --allow-extra-chr --chr 1-22 --out PANEL.aut --make-bed

# 4. Filter out non-SNP variants (indels, etc.)
plink --bfile PANEL.aut --snps-only --out PANEL.snp --make-bed

# 5. Remove duplicate variants
awk '{print $2}' PANEL.snp.bim | sort | uniq -d > dupslist
plink --bfile PANEL.snp --exclude dupslist --out PANEL.snp.rmd --make-bed

# 6. Filter out SNPs (>5% missing) and individuals (>10% missing) with missing data
plink --bfile PANEL.snp.rmd --mind 0.1 --geno 0.05 --out PANEL.snp.rmd.mg --make-bed

# 7. Apply HWE filter
plink --bfile PANEL.snp.rmd.mg --hwe 1e-5 --out PANEL.snp.rmd.mg.hwe --make-bed

# 8. Rename SNPs in chr:pos format
paste -d" " <(awk '{print $2}' PANEL.snp.rmd.mg.hwe.bim) \
<(awk -v OFS=':' '{print $1,$4}' PANEL.snp.rmd.mg.hwe.bim) > update_snpids
plink --bfile PANEL.snp.rmd.mg.hwe --update-name update_snpids --out PANEL.snprnm --make-bed

# 9. Run snpflip to check strand direction. Remove ambiguous SNPs and reverse strands as needed.
# Assumes a local installation of snpflip and access to the human genome GRCh37 fasta file human_g1k_v37.fasta
./snpflip -f human_g1k_v37.fasta -b PANEL.snprnm -o PANEL.flip
mkdir snpflip_out
mv PANEL.flip* snpflip_out

# 9.1. Remove ambiguous SNPs 
plink --bfile PANEL.snprnm --exclude snpflip_out/*ambiguous --out PANEL.flip.amb --make-bed

# 9.2. Flip reversed strands
plink --bfile PANEL.flip.amb --flip snpflip_out/*reverse --out PANEL.final --make-bed

###############
### MERGING ###
###############

# Number of merges varied based on analyses, as plink only permits 2-way merges
# Only one merge here included as an example; assumes QC'd sample data in SAMPLE.final plink format

# 1. Extract common SNPs in both panels
plink --bfile SAMPLE.final --extract PANEL.final.bim --out SAMPLE_PANEL.int --make-bed
plink --bfile PANEL.final --extract SAMPLE.final.bim --out PANEL_SAMPLE.int --make-bed

# 2. Attempt to merge
plink --bfile SAMPLE_PANEL.int --bmerge PANEL_SAMPLE.int.bed PANEL_SAMPLE.int.bim PANEL_SAMPLE.int.fam --out SAMPLE_PANEL.merge --make-bed

# 3. If step 2 fails, exlude multi-allelic SNPs and re-merge (skip if step 2 does not fail)
plink --bfile SAMPLE_PANEL.int --exclude SAMPLE_PANEL.merge-merge.missnp --make-bed --out SAMPLE_PANEL.int.rm
plink --bfile PANEL_SAMPLE.int --exclude SAMPLE_PANEL.merge-merge.missnp --make-bed --out PANEL_SAMPLE.int.rm
plink --bfile SAMPLE_PANEL.int --bmerge PANEL_SAMPLE.int.rm.bed PANEL_SAMPLE.int.rm.bim PANEL_SAMPLE.int.rm.fam --out SAMPLE_PANEL.merge --make-bed

# 4. Perform a final filter for missing data on the merged panel
plink --bfile SAMPLE_PANEL.merge --mind 0.1 --geno 0.05 --out SAMPLE_PANEL.final --make-bed