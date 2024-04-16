#!/bin/bash

# Script for using ADMIXTURE to determine ancestry from genomewide SNP data
# Written by M. Nieves Colon on 5.13.15 (based on ADMIXTURE manual)
    # E. Ulrich: ulriche227@gmail.com
    # M. Nieves-Colón: mnievesc@umn.edu

# Modified 4.20.19 to run with bootstrap

# Input files should be in binary plink (bed/bim/fam) or .ped/.map plink and merged for the sample and reference panels
# Usage: ./04_ADMIXTURE.sh  $1=merged_dataset_basename

# Unsupervised analysis with bootstrap and using 4 threads.
# Run ADMIXTURE for K 2-10 with the crossvalidation parameter activated such
# that we can see the CV errors. The ideal K is the one with the lowest cv errors.
# The tee command is for taking the standard output and writing it to a file.

# Working directory
cd /home/mnievesc/shared/projects/AfroPR/analyses/4_globanc_K3/ADMIXTURE+WilcoxTest
admixture=/home/mnievesc/shared/programs/bin/admixture

# Run ADMIXTURE in parallel for all .bed files and for K 2 to 10.
for K in {2..10}; do $admixture --cv -B -j5 ../Smartpca/$1.LDpruned.bed $K | tee $1.log${K}.out; done

# 4. Use grep to see the cv parameter (-h is for supressing file name)
grep -h CV $1.log*.out
grep -h CV $1.log*.out > $1.Kcomparison.txt