#!/bin/bash
# This is a short script for the data preprocessing to harmonize the different types of data from the same cohort. 
# It is also good to do it before the quality control process
# Juliana Acosta-Uribe 2023

# Specify file name and sequencing type
RAW_dataset='ADFTD.vqsr.snp.indel' # Name of the original file
seq='EXOME' # Type of sequencing (GENOME, EXOME or ARRAY)

# Provide additional files
redlat_samples='redlat_samples.txt' # List of the samples to be retained from the original file. Take a look at the sample IDs in the raw file, as these may be different from the RedLat IDs. For array data in plink files you must give two column ID. See plink documentation 
new_ids='redlat_ids.txt' # List of the new names to be given to the samples. You need to create a file (e.g. ids.txt) where every line has two columns: 'old sample ID' 'new sample ID'. For array data in plink files you must give a four column file. old sample FID' old sample IID' 'new sample FID''new sample IID' See plink documentation.  
targets='xgen-exome-hyb-panel-v2-targets-hg38.bed' # If your data is an EXOME, you need to provide the assay targets. make sure the vcf dile and the bed have the chromosomes using the same style

# Provide a path to your reference genome
fasta_file='/home/acostauribe/public_html/Utilities/hg38.fa.gz'


if [ $seq == 'ARRAY' ]; then
  # I. Extract ReDLat samples from the original file
  plink --bfile $RAW_dataset --keep $redlat_samples --make-bed --out $RAW_dataset.redlat
  # Note: I initially obtained basic qc metrics in the entire data and identified duplicated samples,and those are not included in $redlat_samples
  
  # II. Rename samples according to ReDLat sequence IDs
  plink --bfile $RAW_dataset.redlat --update-ids $new_ids --make-bed --out $RAW_dataset.redlat.id
    
  # III. Recode variants according to the reference genome using PLINK2
  plink2 --bfile $RAW_dataset.redlat.id --ref-from-fa force --fa $fasta_file --real-ref-alleles --make-bed --out $RAW_dataset.redlat.id.ref
  
  # IV. Change the . in the bim files to a numerical sequence
  plink --bfile $RAW_dataset.redlat.id.ref --set-missing-var-ids @:#[hg38]\$1,\$2 --keep-allele-order --make-bed --out $RAW_dataset.redlat.id.ref.var
  awk 'BEGIN {count=1}; {if ($2 ~ /\./) {sub(/\./,"INDEL"(count++));print} else {print} }' $RAW_dataset.redlat.id.ref.var.bim > $RAW_dataset.redlat.id.ref.var.bim2
  mv $RAW_dataset.redlat.id.ref.var.bim2 $RAW_dataset.redlat.id.ref.var.bim 

elif [[ $seq == 'EXOME' || $seq == 'GENOME' ]]; then
  
  # I. Extract targets (if it's an Exome) and ReDLat samples from the original file
  if [ $seq == 'EXOME' ]; then
    vcftools --gzvcf $RAW_dataset.vcf.gz --bed $targets --keep $redlat_samples --recode --recode-INFO-all --out $RAW_dataset.redlat
    # WARNING: Chromosomes in the VCF and in the Target file must have the same encoding (e.g. chr#) If you need to fix this, look at Step IV
  else 
    vcftools --gzvcf $RAW_dataset.vcf.gz --keep $redlat_samples --recode --recode-INFO-all --out $RAW_dataset.redlat
  fi
  
  mv $RAW_dataset.redlat.recode.vcf $RAW_dataset.redlat.vcf
  bgzip --threads 4 $RAW_dataset.redlat.vcf
  # you can add --threads # to make it faster
  
  # II. Rename samples according to ReDLat sequence IDs
  bcftools reheader --samples $new_ids $RAW_dataset.redlat.vcf.gz > $RAW_dataset.redlat.id.vcf.gz
  # https://samtools.github.io/bcftools/bcftools.html#reheader
  
  # III. bgzip and index vcf
  tabix -p vcf $RAW_dataset.redlat.id.vcf.gz #File should be indexed with Tabix
  
  # IV. Recode according to the reference genome alleles and normalize INDELs
  bcftools norm --check-ref ws --fasta-ref $fasta_file --output-type z $RAW_dataset.redlat.id.vcf.gz > $RAW_dataset.redlat.id.ref.vcf.gz
  # --check-ref warn (w), exclude (x), or set/fix (s)
  # --output-type compressed VCF (z) will bgzip the output
  # fasta file was downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/
  # The index file fasta.fai was created using http://www.htslib.org/doc/samtools-faidx.html
  # WARNING: Chromosomes must be encoded as chr#. 
  # If you need to change chromosome designation from # to chr# you need to create a file (e.g. chromosomes.txt) where every line has two columns: 'old name' 'new name' (e.g. 1 chr1)
  # bcftools annotate --rename-chrs chromosomes.txt --output-type z file.vcf.gz > file.chr.vcf.gz
fi


  
