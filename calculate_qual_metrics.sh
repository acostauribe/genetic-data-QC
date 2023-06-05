#! /user/bin/bash

## I. Data description 
bcftools stats $1.vcf.gz > $1.vchk.txt
  
## II. 1 Individual missingness
vcftools --gzvcf $1.vcf.gz --missing-indv --out $1

# Generates a file reporting the missingness on a per-individual basis. 
# The output file has the suffix ".imiss".
# Individuals whose missingness is >10% should be added to 'flagged_samples' file
  
## II. 2. Individual depth 
vcftools --gzvcf $1.vcf.gz --depth --out $1

# Generates a file containing the mean depth per individual. 
# This file has the suffix ".idepth".
# Individuals whose mean depth is <20 should be added to 'flagged_samples' file
  
## II. 3. Individual heterozygosity 
vcftools --gzvcf $1.vcf.gz --het --out $1

# Inbreeding coefficient, F, is estimated for each individual using a method of moments. 
# The resulting file has the suffix ".het"
# Individuals whose heterozygosity deviated more than 3 SD from the main should be added to 'flagged_samples' file
  
##III. 1. Site missingness 
vcftools --gzvcf $1.vcf.gz --missing-site --out $1

# Generates a file reporting the missingness on a per-site basis. 
# The file has the suffix ".lmiss".
  
## III. 2.Site depth
vcftools --gzvcf $1.vcf.gz --site-mean-depth --out $1

# Generates a file containing the mean depth per site across all individuals. 
# This output file has the suffix ".ldepth.mean"
  
## III. 3 VQRS filtering
# General distribution of depth, missingness, heterozygosity
vcftools --gzvcf $1.vcf.gz --FILTER-summary --out $1
# Generates a summary of the number of SNPs and Ts/Tv ratio for each FILTER category. 
# The output file has the suffix ".FILTER.summary"
