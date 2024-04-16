#!/bin/bash

# Pipeline for processing raw plink files from array genotyping on Kexol.
# Written by M. Nieves-Colon with input from C. Quinto - November 2017   
    # mnievesc@umn.edu

# Modified with commands used for analysis of AfroPR project dataset.
# Input is raw plink binary files for dataset: AfroPR.sex.bed .bim .fam
# Usage: 2_dataqc_kexol.sh | tee dataqc+snpflip.log

# Working Directory:
WD=/data/users/mnievesc/AfroPRProject/analyses/2.DataQC


##################### Per-marker QC #####################

# 1. Update physical coordinates, chromosome assigmnents and SNP names.
# With this step we make sure that the SNPids (rsids), SNP positions and chromosome assignments
# for each variant are up to date. For the later we want to make sure that any variant listed
# at chromosome 0 is updated or removed if we cannot ascertain where it maps to.
# In this case I suspect all the variants with rsids that include "SNP" map to chr M.
# All others will be excluded

# 1A. Update chromosome code for MT variants listed as being on chr 0.
# From plink documentation: --update-chr, --update-cm, --update-map, and --update-name update variant chromosomes,
# centimorgan positions, base-pair positions, and IDs, respectively. By default, the new value is read from column 2 and
# the (old) variant ID from column 1, but you can adjust these positions with #the second and third parameters.
awk '{ if ($1=="0") print $2" MT"}' ../1_TXTtoPlink/AfroPR.sex.bim | grep -E "SNP" > update_chr.txt
plink --bfile ../1_TXTtoPlink/AfroPR.sex --update-chr update_chr.txt --make-bed --out AfroPR.sex.t1
echo " "

# Sanity Check - Check for changes in chromosome assignments
grep "105SNP6248RT_C" AfroPR.sex.t1.bim
grep "105SNP6248RT_C" ../1_TXTtoPlink/AfroPR.sex.bim
echo " "


# 1B. Remove the sites that remain mapped to chr 0 using plink command --not-chr which excludes
# variants on specified chromosome
plink --bfile AfroPR.sex.t1 --not-chr 0 --make-bed --out AfroPR.sex.t2
echo " "


# 1C. Update physical coordinates for chrMT sites. All sites are listed as being on position 0 at the
# moment. This is a problem because they will appear as duplicates later on (sites with different ID but same position).
# To make list of sites to update coordinates for we must provide the variant and then the new position.
# Use awk to print out all SNPs mapping to chr MT/26. The position of the SNP is contained in notation like this
# "105SNP6248RT_C". So use awk to strip out all text before and including SNP. Then use sed to exclude all
# letters and the underscore character.
paste <(awk '{ if ($1=="26") print $2}' AfroPR.sex.t2.bim) \
<(awk '{ if ($1=="26") print $2}' AfroPR.sex.t2.bim | awk 'BEGIN {FS="SNP";}{print $2}' |  sed 's/[a-zA-Z]*//g' |  sed 's/_//g') > update_map.txt
plink --bfile AfroPR.sex.t2 --update-map update_map.txt --make-bed --out AfroPR.sex.t3

# Sanity check
grep "105SNP6248RT_C" AfroPR.sex.t2.bim
grep "105SNP6248RT_C" AfroPR.sex.t3.bim
echo " "


echo "########################################### "
raw=`wc -l ../1_TXTtoPlink/AfroPR.sex.bim | cut -f1 -d" "`
echo "Number of variants before processing: $raw"
t3=`wc -l AfroPR.sex.t3.bim | cut -f1 -d" "`
nomap=`echo "$raw-$t3" | bc`
echo "Number of variants that did not map to any chromosomes: $nomap"
echo "Total number of variants after updating names, names and positions: $t3"
echo "########################################### "
echo " "



# 2. Remove duplicate SNPs (requires R)

# 2A. Transform binary file to transposed format for input into R
plink --bfile AfroPR.sex.t3 --recode --transpose --out AfroPR.sex.t3
echo " "

# 2B. Run Rscript to remove duplicates
Rscript ~/AfroPRProject/scripts/remove_duplicateSNP_v2.R AfroPR.sex.t3.tped AfroPR.sex.t4.tped

# 2C. Make new tfam file (back in UNIX terminal)
cp AfroPR.sex.t3.tfam AfroPR.sex.t4.tfam

# 2D. Do a second check of duplicates using plink. Make lisf of duplicates.
# Option ids-only and suppres-first only keep first id of the list of two duplicates.
plink --tfile AfroPR.sex.t4 --list-duplicate-vars ids-only suppress-first
wc -l plink.dupvar #
echo " "

# 4E. Remove duplicates again using plink. This time transform output back to binary format.
plink --tfile AfroPR.sex.t4 --exclude plink.dupvar --make-bed --out AfroPR.sex.t5
echo " "

# Sanity check
echo "########################################### "
dup=`wc -l AfroPR.sex.t3.bim | cut -f1 -d" "`
dedup=`wc -l AfroPR.sex.t4.tped | cut -f1 -d" "`
dedup2=`wc -l AfroPR.sex.t5.bim | cut -f1 -d" "`
echo "Number of variants before duplicate removal with R: $dup"
echo "Number of variants after duplicate removal with R: $dedup"
echo "Number of variants after duplicate removal with plink: $dedup2"
echo "########################################### "
echo " "



# 3. Exclude structural variants

# 3A. Make list of indel variants by selecting those that have I or D as alleles
awk '$5 ~/I/ || $5 ~/D/ || $6 ~/I/ || $6 ~/D/ ' AfroPR.sex.t5.bim | cut -f2 > indels.txt

# 3B. Remove indels with plink
plink --bfile AfroPR.sex.t5 --exclude indels.txt --make-bed --out  AfroPR.sex.t6
echo " "

# Sanity check
indels=`wc -l indels.txt | cut -f1 -d" "`
noindels=`wc -l AfroPR.sex.t6.bim | cut -f1 -d" "`
echo "########################################### "
echo "Number of indels identified: $indels"
echo "Number of SNPs after removing indels: $noindels"
echo "########################################### "
echo " "



# 4. Flip SNPs to forward strand (requires snpflip and perl)

# 4A. Rename chromosome codes for bim file. Make temporary bim file with new codes just to get list
# of SNPs for flipping with actual file in plink.
perl -p -e 's/^23/X/' AfroPR.sex.t6.bim | perl -p -e 's/^24/Y/' | perl -p -e 's/^25/X/' | perl -p -e 's/^26/M/' > temp.t6.bim

# 4B. Run snpflip and store output in directory
snpflip -b temp.t6.bim -f /data/programs/LiftOver/human_g1k_v37.fasta -o snpflip_output.t6
mkdir snpflip_output.t6
mv snpflip_output.t6* snpflip_output.t6/

# 4C. Organize list of reverse and ambiguous SNPs
sort snpflip_output.t6/snpflip_output.t6.reverse | uniq > snpflip_output.t6/snpflip_output.t6.reverse.uniq
sort snpflip_output.t6/snpflip_output.t6.ambiguous | uniq > snpflip_output.t6/snpflip_output.t6.ambiguous.uniq

# Sanity check
wc -l snpflip_output.t6/snpflip_output.t6.*
echo " "

# 4D. Exclude all SNPs listed as ambiguous
plink --bfile AfroPR.sex.t6 --exclude snpflip_output.t6/snpflip_output.t6.ambiguous.uniq --make-bed --out AfroPR.sex.t7

# 4E. Flip all SNPs listed by snpflip as reversed.
plink --bfile AfroPR.sex.t7 --flip snpflip_output.t6/snpflip_output.t6.reverse.uniq --make-bed --out AfroPR.flip

# Sanity check
echo "########################################### "
rev=`wc -l snpflip_output.t6/snpflip_output.t6.reverse.uniq| cut -f1 -d" "`
amb=`wc -l snpflip_output.t6/snpflip_output.t6.ambiguous.uniq| cut -f1 -d" "`
postflip=`wc -l AfroPR.flip.bim | cut -f1 -d" "`
echo "Number of SNPs to exclude (ambiguous): $amb"
echo "Number of SNPs to flip (reverse): $rev"
echo "Number of SNPs after excluding ambiguous and flipping: $postflip"
echo "########################################### "
echo " "


##################### Per-individual QC #####################

# 5. Check for discordant sex information
# Use plink --check-sex command to calculate the mean homozygosity rate across
# X-chromosome markers for each individual.

# 5A. Check sex
plink --bfile AfroPR.flip --check-sex --out AfroPR.flip
grep "PROBLEM" AfroPR.flip.sexcheck
echo " "



# 7. Check rates of missing genotypes per sample.
# 7A.Get list of missing genotypes (imiss)
plink --bfile AfroPR.flip --missing --out AfroPR.flip
echo " "

# 7B.  Sort ismiss file to look at individuals with most missing data.
# Tail everything except for first line to get rid of header.
tail -n +2 AfroPR.flip.imiss | sort --key 6 --numeric-sort
echo " "

# 7C. Get list of samples with >10% missing call rate
tail -n +2 AfroPR.flip.imiss | awk '$6 > 0.10' | wc -l
tail -n +2 AfroPR.flip.imiss | awk '$6 > 0.10'
echo " "

# 7D. Sort lmiss file to look at individuals with most missing data.
# Tail everything except for first line to get rid of header.
tail -n +2 AfroPR.flip.lmiss | sort --key 6 --numeric-sort | awk '$6 > 0.10'
echo " "

# 7E. Also estimate heterozygosity
plink --bfile AfroPR.flip --het --out AfroPR.flip
echo " "


# 8. Plot figures to illustrate missinges (requires R)
mkdir QCplots
Rscript ~/AfroPRProject/scripts/missingness_qcplots_v3.R AfroPR.flip pre AfroPR-project
echo " "


# 9. Check IBS to find related individuals.

# 9A. Generate list of pairwise IBD
plink --bfile AfroPR.flip --genome --out AfroPR.flip
echo " "

# 9B. Use script from Anderson et al 2010 paper to identify individuals that fail
# threshold of 0.5 (first degree relatives. Output is fail-IBD-QC.txt list of individuals
# with high IBD. Requires perl.
perl ~/AfroPRProject/scripts/run-IBD-QC_mod82316.pl AfroPR.flip


# 9C. Remove problematic individuals. Choose individuals after looking at figures and
# plink analyses. Create two datasets one with all chromosomes and another restricted
# to just autosomes. Make list externally. DID NOT RUN THIS YET need to evaluate manually
#plink --bfile AfroPR.flip --remove inds-remove.txt --make-bed --out AfroPR
#plink --bfile AfroPR.flip --remove inds-remove.txt --autosome --make-bed --out AfroPR.aut


# 10. Re-do missingess estimates and generate new figures (requires R)
# plink --bfile AfroPR.aut --missing --out AfroPR.aut
# plink --bfile AfroPR.aut --het --out AfroPR.aut
# Rscript missingess_qcplots.R AfroPR.aut post




##################### Additional filtering #####################

# 11. Additional filtering.

# 11A. Make two datasets - one restricted to autosomes and one with all chromosomes
plink --bfile AfroPR.flip --autosome --make-bed --out AfroPR.flip.aut

# 11B. Missing call rates for SNPs and samples.
# --geno filters out all variants with missing call rates exceeding the provided value
# (default 0.1) to be removed, while --mind does the same for samples.
plink --bfile AfroPR.flip.aut --mind 0.1 --geno 0.05 --make-bed --out AfroPR.flip.aut.mg
echo " "
plink --bfile AfroPR.flip --mind 0.1 --geno 0.05 --make-bed --out AfroPR.flip
echo " "


# 11C. HWE threshold filtering
plink --bfile AfroPR.flip.aut --hwe 1e-5 --make-bed --out AfroPR.aut.final
echo " "
plink --bfile AfroPR.flip --hwe 1e-5 --make-bed --out AfroPR.final
echo " "




##################### Rename SNPs chr:pos  #####################

# 12. Rename SNP IDS from rsid to chr:pos format. This will allow us to increase overlap With
# reference panels that may have same variants annotated with different positions.

# 12A. Make list of new SNP names
paste -d" " <(awk '{print $2}' AfroPR.aut.final.bim) <(awk -v OFS=':' '{print $1,$4}' AfroPR.aut.final.bim) > AfroPR.aut.final.update_snpids.txt
paste -d" " <(awk '{print $2}' AfroPR.final.bim) <(awk -v OFS=':' '{print $1,$4}' AfroPR.final.bim) > AfroPR.final.update_snpids.txt

# 12B. Use plink to update SNP IDs
plink --bfile AfroPR.aut.final --update-name AfroPR.aut.final.update_snpids.txt --make-bed --out AfroPR.aut.final.snprnm
plink --bfile AfroPR.final --update-name AfroPR.final.update_snpids.txt --make-bed --out AfroPR.final.snprnm

echo " "
echo "Done! Final files have suffix X.aut.final.bim .bed .fam"



# 13. Clean up
mkdir intermediate_files
mv temp* intermediate_files/
mv *txt intermediate_files/
mv plink* intermediate_files/
mv AfroPR.sex* intermediate_files/
mv AfroPR.flip* intermediate_files/
mv snpflip* intermediate_files/
tar -czvf intermediate_files.tar.gz intermediate_files/
rm -r intermediate_files/


# References
# 1. https://www.gnu.org/software/sed/manual/html_node/Regular-Expressions.html
# 2. https://stackoverflow.com/questions/25447324/how-to-use-cut-with-multiple-character-delimiter-unix
