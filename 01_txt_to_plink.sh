#!/bin/bash

# Script to convert all text files to one map and ped file.
# Written for AfroPR project data by M. Nieves Colon - 10/03/2018
  # mnievesc@umn.edu

# 1. Create map file, all samples have same variants so any sample will do
# mapfile: chr rsid pos bp (space delimited)
awk 'NR>11 {print $3, $2, 0, $4}' AC36C6EG5V-GRC171426543.txt  > AfroPR.map


# 2. Create pedfile including data for all individuals (space delimited)
# pedfile: famID indID father mother sex phenotype variant1:A1 A2 variant2:A1 A2, etc...
# Each row is an individual. Use for loop to process all 30 samples.
for SAMPLE in `ls *txt | grep -v Batch1547_NG_samplesfile.txt`;
  do
    #echo $SAMPLE
    ID=$(basename $SAMPLE)
    NAME=$(echo $ID | cut -d "-" -f2 | cut -d"." -f1)
    echo "Processing sample $NAME"


    # For all textfiles, row 12 has the first variant. Name population then export the first info field
    # which has sample ID and set father, mother, sex and phenotype to 0. Next xport fields 5 and 6 which correspond
    # to the first diploid genotype (these will be columns 7 and 8 in ped).
    awk 'NR==12  {print "AfroPR", $1, 0, 0, 0, 0, $5, $6}' $ID  > $NAME.temp.ped

    # Transpose fields 5 and 6 of the text file for all subsequent rows in textfile (row 12+).
    # These will be columns 9+ with information on diploid genotypes for this individual.
    # Must use gsub to make this space delimited.
    awk 'NR>12 {print $0}' $ID | cut -f5,6 | paste -s | awk 'gsub("\t"," ")' > $NAME.temp.row

    # Paste the two temp files together to make the final pedfile for this individual.
    # Must use -d" " to make space delimited.
    paste -d" " $NAME.temp.ped $NAME.temp.row > $NAME.ped

    # Remove temp files
    rm *temp*

    # Sanity check
    # Should only have one row and (711059 SNPs *2 diploid)=1422118 variants + 6 info fields = 1422124 columns
    echo "Rows: `wc -l $NAME.ped`"  # count rows
    echo "Genotype Columns: `awk -F " " '{print NF; exit}' $NAME.ped`" # count columns
    echo " "
  done


# 3. Cat the pedfiles together to make a single file including all 30 individuals.
# Recode this file in plink to ensure proper input.
cat *ped > AfroPR.ped

# Sanity check
cut -d" " -f1-14 AfroPR.ped

plink --file AfroPR --allow-extra-chr 0 --missing-genotype '-' --recode --out AfroPR.raw
# 711059 variants and 30 people pass filters and QC.


# 4. Update sex information.
# --update-sex expects a file with FIDs and IIDs in the first two columns, and sex information (1 or M = male, 2 or ## F = female, 0 = missing) in the (n+2)th column.
awk 'NR>1 {print $0}' Batch1547_NG_samplesfile.txt  | awk '{print "AfroPR", $2, $6}' > updatesex.txt
plink --file AfroPR.raw --allow-extra-chr 0 --missing-genotype '-' \
 --update-sex updatesex.txt --make-bed --out AfroPR.sex
# Total genotyping rate is 0.996202.
#711059 variants and 30 people pass filters and QC.


# 5. Move data files to working directory and remove individual .ped files
mv AfroPR* ../../analysis/1_TXTtoPlink/
rm *ped


# References
# 1. https://www.thelinuxrain.com/articles/transposing-rows-and-columns-3-methods
# 2. https://www.unix.com/unix-for-dummies-questions-and-answers/101733-converting-space-delimited-file-tab-delimited-file.html
# 3. http://www.theunixschool.com/2012/09/examples-how-to-change-delimiter-of-file-Linux.html
