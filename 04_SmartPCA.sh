#!/bin/bash

# Script for performing PCA with SmartPCA on merged genomewide SNP data
# Written by M. Nieves Colon on 11.5.16 and modified by E. Ulrich on 4.20.22
    # E. Ulrich: ulriche227@gmail.com
    # M. Nieves-Colón: mnievesc@umn.edu

# Software used: smartpca version 13050, PLINK 1.9
# Usage: 03_SmartPCA.sh $1=input $2=output
    # $1=../smartpca_data/AfroPR-1KGP-Puno.merge.mg
    # $2=AfroPR-1KGP

smartpca=/home/mnievesc/shared/programs/bin/smartpca 
module load plink

# 1. LD trim merged dataset
plink --bfile ${1} --indep-pairwise 50 10 0.1 --out ${2}
plink --bfile ${1} --extract ${2}.prune.in --geno 0.1 --make-bed --out ${2}.LDpruned

# 2. Save as .ped .map format for input to smartpca. Output missing phenotype as 1.
plink --bfile ${2}.LDpruned --recode --out ${2}.LDpruned --output-missing-phenotype 1

# 3. Make .par file for smpartpca
echo -e genotypename: ${2}.LDpruned.ped"\n"snpname: ${2}.LDpruned.map"\n"indivname: ${2}.LDpruned.ped"\n"evecoutname: ${2}.LDpruned.evec"\n"evaloutname: ${2}.LDpruned.eval"\n"numoutlieriter: 0 > ${2}.LDpruned.par

# 4. Run smartpca given directory structure above
$smartpca -p ${2}.LDpruned.par > ${2}.smartpca.log
