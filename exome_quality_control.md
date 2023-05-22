---
title: "Exome Quality Control Pipeline"
author: "Juliana Acosta-Uribe"
date: '2022-12-23'
output: 
  html_document: 
    toc: yes
    fig_caption: yes
    number_sections: yes
keep_md: true
---

# Quality control tutorial

This script/tutorial is an overview of the quality control (QC) processes that were done for the ReD-Lat Exomes.\
Developed by Juliana Acosta-Uribe December 2022

You will need: - List of samples to exclude: 'excluded_samples_exome.txt' - fasta_file of the reference genome: '\~/public/utilities/hg38.fa.gz' - a list with the redlat ids and the libary id new_ids'exome_to_redlatIDs.txt' - targets of the exome panel 'xgen-exome-hyb-panel-v2-targets-hg38.bed' - a dataframe with disclosed sex of samples. sex should be encoded 1=male, 2=female 'exome_sample_sex.txt'

## 1. Set up environment and data for QC

### 1.a Set up R

It may be necessary to start R before, as you need to have knitr to run this markdown properly

```{r markdown-setup}
## Install and load required packages
#install.packages("psych")
library(psych)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("knitr")
library(knitr)
#install.packages("dplyr")
library(dplyr)
#install.packages("kinship2")
library(kinship2)

## Set your working directory
knitr::opts_chunk$set(root.dir = "/home/acostauribe/Genetic-Sequencing_raw-data/HudsonAlpha_J.Nicholas.Cochran/ReDLat/Exome_Seq/Joint_call_fixed",
                      tidy=TRUE,
                      engine.path = list(plink = '/home/acostauribe/bin/plink',
                                         king = '/home/acostauribe/bin/king'))
PREFIX="redlat_exomes"
```

### 1.b Data preprocessing

We will use **vcftools** and **bcftools** to clean our "starting dataset".

Since these are exomes, we first need to extract the 'exome sequencing targets' from the raw dataset. Targets are a '.bed' file downloaded from IDT website for xGen™ Exome Hybridization Panel v2 hg38 [the sequencing platform that was used]. This .bed file has three columns: *chr*, *start bp*, *end bp*. Additionally, we will remove individuals that are not part of the ReDLat project Exomes, but are in the same file.

**Note**: VCFtools output can be compressed directly into a .gz file by adding `–stdout | gzip -c > newname.gz`, but it wont generate a .log file.

Then, we will use bcftools to perform Left- alignment and normalization of indels, as well as a check to verify that the REF/ALT alleles are correct. However, it will NOT fix strand issues in your VCF (make sure VCF is properly aligned).

**Keep in mind**: Unlike R, all bash code chunks in a markdown are independent in memory, and the variables created in previous chunks will not be available in latter chunks. you have to define your variables/alias at each chunk

```{bash extract-targets-align-reference-normalize, eval=FALSE, cache=FALSE, include=FALSE}
RAW_dataset='ADFTD.vqsr.snp.indel'
seq='exome'
targets='xgen-exome-hyb-panel-v2-targets-hg38.bed'
RAW_dataset_exclude_samples='excluded_samples_exome.txt'
fasta_file='~/public/utilities/hg38.fa.gz'
new_ids='exome_to_redlatIDs.txt'

# I. Extract targets and ReDLat samples from the original file
vcftools --gzvcf $RAW_dataset.vcf.gz 
--bed $targets 
--remove $RAW_dataset_exclude_samples 
--recode 
--recode-INFO-all 
--out $RAW_dataset.redlat

mv $RAW_dataset.redlat.recode.vcf $RAW_dataset.redlat.vcf
bgzip $RAW_dataset.redlat.vcf

# II. Rename samples according to ReDLat sequence IDs
# https://samtools.github.io/bcftools/bcftools.html#reheader

bcftools reheader --samples $new_ids $RAW_dataset.redlat.vcf.gz > redlat_$seq.vcf

# III. bgzip and index vcf
bgzip redlat_$seq.vcf
tabix -p vcf redlat_$seq.vcf.gz #File should be indexed with Tabix

# IV. Check reference allele and normalize INDELs
bcftools norm --check-ref ws --fasta-ref $fasta_file --output-type z redlat_$seq.vcf.gz > redlat_$seq.temp.vcf
# --check-ref warn (w), exclude (x), or set/fix (s)
# --output-type compressed VCF (z) 
# fasta file was downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/
# The index file fasta.fai was created using http://www.htslib.org/doc/samtools-faidx.html

mv redlat_$seq.temp.vcf.gz redlat_$seq.vcf.gz
tabix -p vcf redlat_$seq.vcf.gz #File should be indexed with Tabix

#Tidy up

mkdir 1.data_preprocessing
mv $RAW_dataset* 1.data_preprocessing
mv $RAW_dataset_exclude_samples 1.data_preprocessing
mv $new_ids 1.data_preprocessing
mv $targets 1.data_preprocessing
```

## 2: General report of data quality and statistics previous to QC

For the following steps we are going to split the data into autosomes and the sex chromosomes. The reason we do this is because quality metrics like *Depth* or *Missingness* are very different between autosomes (where all individuals are diploid) and X and Y (where individuals can be diploid or haploid)

### 2.a Autosomes quality report

#### i. Extract autosomes

```{bash retain-autosomes, eval=FALSE, cache=FALSE, include=FALSE}
PREFIX='redlat_exomes'

## Split into autosomes and sex chromosomes
## It is important to know how the chromosomes numbers are encoded in your vcf, e.g. they can be '1' or 'chr1'. IN this case is chr1
vcftools --gzvcf $PREFIX.vcf.gz --chr chr1 --chr chr2 --chr chr3 --chr chr4 --chr chr5 --chr chr6 --chr chr7 --chr chr8 --chr chr9 --chr chr10 --chr chr11 --chr chr12 --chr chr13 --chr chr14 --chr chr15 --chr chr16 --chr chr17 --chr chr18 --chr chr19 --chr chr20 --chr chr21 --chr chr22 --recode --recode-INFO-all --out $PREFIX.autosomes

mv $PREFIX.autosomes.recode.vcf $PREFIX.autosomes.vcf
bgzip $PREFIX.autosomes.vcf
```

#### ii. Calculate quality metrics

```{bash autosome-preQC-get-stats, eval=FALSE, cache=FALSE, include=FALSE}
PREFIX='redlat_exomes'

## I. Data description 
bcftools stats $PREFIX.autosomes.vcf.gz > $PREFIX.autosomes.vchk.txt
# General distribution of depth, missingness, heterozygosity
vcftools --gzvcf $PREFIX.autosomes.vcf.gz --FILTER-summary --out $PREFIX.autosomes
# Generates a summary of the number of SNPs and Ts/Tv ratio for each FILTER category. 
# The output file has the suffix ".FILTER.summary"

## II. Individual missingness
vcftools --gzvcf $PREFIX.autosomes.vcf.gz --missing-indv --out $PREFIX.autosomes
# Generates a file reporting the missingness on a per-individual basis. 
# The output file has the suffix ".imiss".
# Individuals whose missingness is >10% should be added to 'flagged_samples' file

## III. Individual depth 
vcftools --gzvcf $PREFIX.autosomes.vcf.gz --depth --out $PREFIX.autosomes
# Generates a file containing the mean depth per individual. 
# This file has the suffix ".idepth".
# Individuals whose mean depth is <20 should be added to 'flagged_samples' file

## VI. Individual heterozygosity 
vcftools --gzvcf $PREFIX.autosomes.vcf.gz --het --out $PREFIX.autosomes
# Inbreeding coefficient, F, is estimated for each individual using a method of moments. 
# The resulting file has the suffix ".het"
# Individuals whose heterozygosity deviated more than 3 SD from the main should be added to 'flagged_samples' file

## V. Site missingness 
vcftools --gzvcf $PREFIX.autosomes.vcf.gz --missing-site --out $PREFIX.autosomes
# Generates a file reporting the missingness on a per-site basis. 
# The file has the suffix ".lmiss".

## VI. Site depth
vcftools --gzvcf $PREFIX.autosomes.vcf.gz --site-mean-depth --out $PREFIX.autosomes
# Generates a file containing the mean depth per site across all individuals. 
# This output file has the suffix ".ldepth.mean"
```

**Checkpoint:** At the end of the "General report" you should have at least the 7 following files: I. Sample based: 1. missingness per sample (.imiss) 2. depth per sample (.idepth) 3. heterozugosity per sample (.het) II. Site based 1. Missingness per site (.lmiss) 2. Mean depth per site (.ldepth.mean) 3. VQSR quality (vqsr) III. General stats (.vchk)

#### iii. Plot quality metrics

##### I. Sample based Metrics:

###### - Missingness per sample

```{r autosome-preQC-stats-sample-missingess}
imiss = read.delim(paste0(PREFIX,".autosomes.imiss"), header = T, sep = "")
# Get some stats of missingness rate per sample (given as a fraction)
imiss_F_MISS = describe(imiss$F_MISS) #describe function from the "Psych" package gives basic stats
rownames(imiss_F_MISS) = c("sample_missingness_F")
# Individuals whose missingness is >10% should be identified
high_missingness = filter(imiss, F_MISS>0.1)
# Save it. It may come useful for debugging
write.table(high_missingness, 
            (paste0(PREFIX,".autosomes.FLAGGED_sample_high_missingness.preQC.txt")),
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
# Missingness counts per sample (given as an integer where N = variants)
imiss_N_MISS = describe(imiss$N_MISS)
rownames(imiss_N_MISS) = c("sample_missingness_N")

missingness_sample = bind_rows(imiss_F_MISS,
                               imiss_N_MISS)
                               
# Plot Missingness rate per sample
png(file=(paste0(PREFIX,".autosomes.sample_missingness.preQC.png")), 
    width=600, height=350)
hist(imiss$F_MISS,
     xlab="Missingness rate",
     ylab="Samples", 
     main="Missingness rate per sample in - Autosomes preQC", 
     col="paleturquoise3")
dev.off()
```

###### - Depth per sample

```{r autosome-preQC-stats-sample-depth}
idepth = read.delim((paste0(PREFIX,".autosomes.idepth")), header = T, sep = "")

# Get some site depth stats:
depth_sample = describe(idepth$MEAN_DEPTH)
rownames(depth_sample) = c("sample_depth")

# Individuals whose mean depth is <20 should be identified
low_mean_depth = filter(idepth, MEAN_DEPTH<20)
write.table(low_mean_depth, 
            (paste0(PREFIX,".autosomes.FLAGGED_sample_low_mean_depth.preQC.txt")),
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')

# Plot depth per sample
png(file=(paste0(PREFIX,".autosomes.sample_depth.preQC.png")), 
    width=600, height=350)
hist(idepth$MEAN_DEPTH,
     xlab="Mean Depth ",
     ylab="Samples", 
     main="Mean Depth per sample - Autosomes preQC", 
     col="paleturquoise3",
     breaks=50)
dev.off()
```

###### - Individual heterozygosity

```{r autosome-preQC-stats-sample-heterozygosity}
het = read.delim((paste0(PREFIX,".autosomes.het")), header = T, sep = "")

# Get some stats:
heterozygosity_sample = describe(het$F)
rownames(heterozygosity_sample) = c("sample_heterozygosity_F")
# Identify the values for +3 and -3 standard deviations
heterozygosity_low_limit = mean(het$F)-(3*(sd(het$F)))
heterozygosity_high_limit = mean(het$F)+(3*(sd(het$F)))

# Plot heterozygosity per sample
png(file=(paste0(PREFIX,".autosomes.sample_heterozygosity.preQC.png")), 
    width=600, height=350)
hist(het$F,  
     freq=TRUE, 
     xlab="Heterozygosity F coefficient",  
     ylab="Samples", 
     main="Heterozygosity rate per sample - Autosomes preQC",
     col="paleturquoise3",
     breaks=50)
abline(v = (heterozygosity_low_limit), col="red")
abline(v = (heterozygosity_high_limit), col="red")
abline(v = (mean(het$F)), col="blue") 
legend("topleft",
       c("+/-3 SD","mean"),
       col=c("red","blue"),
       pch=16)
dev.off()

# Individuals whose heterozygosity deviated more than 3 SD from the mean should be identified
het_outlier_low = filter(het, F<heterozygosity_low_limit)
het_outlier_high = filter(het, F>heterozygosity_high_limit)
het_outlier_both = bind_rows(het_outlier_low,
                             het_outlier_high)
write.table(het_outlier_both, 
            (paste0(PREFIX,".autosomes.FLAGGED_sample_heterozygosity_outliers.preQC.txt")),
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
```

###### - Data distribution of sample based metrics

```{r autosome-preQC-stats-sample-summary}
# Generate a file with  descriptive statistics and distribution of the data
stats_sample = bind_rows(missingness_sample,
                         depth_sample,
                         heterozygosity_sample)
write.table(stats_sample,
           (paste0(PREFIX,".autosomes.sample_descriptive_statistics_preQC.txt")),
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
print(stats_sample)

# Create a dataframe  with the samples that failed the autosome quality thresholds
high_missingness_id = select(high_missingness, INDV)
low_mean_depth_id = select(low_mean_depth, INDV)
het_outlier_low_id = select(het_outlier_low, INDV)
het_outlier_high_id = select(het_outlier_high, INDV)
flagged_samples_autosomes = bind_rows(high_missingness_id,
                                      low_mean_depth_id,
                                      het_outlier_low_id, 
                                      het_outlier_high_id)

# Save it. It may come useful for debugging
write.table(flagged_samples_autosomes, 
            (paste0(PREFIX,".autosomes.FLAGGED_sample_allstats.preQC.txt")),
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
```

###### - Quality statistics per sample

This step is absolutely optional, but I like to generate a dataframe where I annotate all samples with their different QC metrics. It is very useful when you come back and check the reason why a sample was dropped.

```{r autosome-preQC-sample-quality}
# Take the information of individuals from the .idepth file. 
# We are changing the names of the columns to specify these values come from raw data
sample_metrics = rename(idepth, n_sites_raw= N_SITES, mean_depth_raw = MEAN_DEPTH)

# Using the "match" function, we will create a new column in 'sample_metrics' with the missingness per sample.
# Technically, imiss and idepth should have the same samples in the same order, but using the "match" function will be useful when we start dropping samples
sample_metrics$missingness_raw = imiss$F_MISS[match(sample_metrics$INDV, imiss$INDV)]

## Using the "match" function, we will create a new column in 'sample_metrics' with the heterozygosity per sample 
sample_metrics$heterozygosity_F_raw = het$F[match(sample_metrics$INDV, het$INDV)]

## Save as a file (optional)
write.table(sample_metrics,
           (paste0(PREFIX,".sample_based_all_statistics.preQC.txt")),
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE)
```

We will add additional columns in the X and Y quality report `sample_metrics` as a file

##### II Site based Metrics:

###### - Missingness per site

```{r autosome-preQC-stats-plots-site-missingness}
lmiss = read.delim((paste0(PREFIX,".autosomes.lmiss")), header = T, sep = "")

# Get missingness rate per site stats (given as a fraction)
lmiss_F_MISS = describe(lmiss$F_MISS)
rownames(lmiss_F_MISS) = c("site_missingness_F")
# Missingness counts (given as an integer)
lmiss_N_MISS = describe(lmiss$N_MISS)
rownames(lmiss_N_MISS) = c("site_missingness_N")
missingness_site = bind_rows(lmiss_F_MISS,
                             lmiss_N_MISS)
# Plot Missingness per site
png(file=(paste0(PREFIX,".autosomes.site_missingness.preQC.histogram.png")), 
    width=600, height=350)
hist(lmiss$F_MISS,
        xlab="Missingness rate",
        ylab="Number of sites", 
        main="Missingness rate per site - Autosomes preQC", 
        col="paleturquoise3",
        breaks=50)
dev.off()

png(file=(paste0(PREFIX,".autosomes.site_missingness_preQC.box-plot.png")), 
    width=600, height=350)
boxplot(lmiss$F_MISS,
        ylab="Missingness rate",
        xlab="Raw dataset", 
        main="Missingness rate per site - Autosomes preQC", 
        col="paleturquoise3")
dev.off()
```

###### - Depth per site

```{r autosome-preQC-stats-plots-site-depth}
ldepth.mean = read.delim((paste0(PREFIX,".autosomes.ldepth.mean")), header = T, sep = "")

# Get basic stats
depth_site = describe(ldepth.mean$MEAN_DEPTH)
rownames(depth_site) = c("site_mean_depth")

# Plot
png(file=(paste0(PREFIX,".autosomes.site_depth.preQC.png")), 
    width=600, height=350)
hist(ldepth.mean$MEAN_DEPTH,
     xlab="Mean depth",
     ylab="Sites", 
     main="Mean depth per site - Autosomes preQC", 
     col="paleturquoise3")
dev.off()

# Take a zoom at the lower end
png(file=(paste0(PREFIX,".autosomes.site.depth_preQC.low_end.png")), 
    width=600, height=350)
boxplot(ldepth.mean$MEAN_DEPTH,
        ylab="Mean depth per variant",
        xlab="Raw dataset", 
        main="Mean depth per site in Autosomes - preQC - 80X and lower",
        col = c("paleturquoise3"),
        ylim = c(0, 80))
dev.off()
```

###### - Data distribution of site based metrics

```{r autosome-preQC-stats-plots-site-summary}
stats_sites = bind_rows(missingness_site,
                         depth_site)
write.table(stats_sites,
           (paste0(PREFIX,".autosomes.site_descriptive_statistics.preQC.txt")),
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
print(stats_sites)
```

###### - VQSR (Variant Quality Score Recalibration)

```{r autosome-preQC-stats-plots-site-VQSR}
filter_raw = read.delim((paste0(PREFIX,".autosomes.FILTER.summary")), header = T, sep = "")
print(filter_raw)

# Plot it
filter_raw %>%
  filter(!is.na(N_VARIANTS)) %>%
  arrange(N_VARIANTS) %>%
  mutate(FILTER=factor(FILTER, FILTER)) %>%
  ggplot( aes(x=FILTER, y=N_VARIANTS, label = N_VARIANTS) ) +
    geom_segment( aes(x=FILTER ,xend=FILTER, y=0, yend=N_VARIANTS), color="grey") +
    geom_point(size=3, color="#69b3a2") +
    geom_text(vjust=-1, size = 3) +
    coord_flip() +
    theme_minimal() +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position="none") +
    scale_y_continuous(name = "Number of variants") +
    labs(title = "VQSR in Autosomes - preQC")
ggsave(paste0(PREFIX,".autosomes.site_VQSR.pre-QC.png"),
      width = 7,
      height = 4)  
```

Taking a close look to the descriptive statistics files can help decide which are the best thresholds for QC.

and tidy up!

```{bash clean-up-2autosomes}
PREFIX='redlat_exomes'
mkdir 2.autosomes_quality_report
mv $PREFIX.autosomes.* 2.autosomes_quality_report
```

### 2.b Sex chromosomes quality report

Since this data set is an Exome, we wanted to explore the distribution of markers prior to the quality control and determination of chromosomal sex. You will need a dataframe/table with the reported phenotypic sex of the included samples `exome_sample_sex.txt` sex should be encoded: 1=male, 2=female

#### i. Extract sex chromosomes

Before extracting the sex chromosomes check the file to see how are these encoded, a good idea is to do `grep X` or `grep x`. For example, in the ReDLat file, these are encoded as 'chrX' and 'chrY'

```{bash retain-sex-chromosomes, eval=FALSE, cache=FALSE, include=FALSE}
PREFIX='redlat_exomes'
for chr in X Y
do
vcftools --gzvcf $PREFIX.vcf.gz --chr chr$chr --recode --recode-INFO-all --out $PREFIX.$chr
mv $PREFIX.$chr.recode.vcf $PREFIX.$chr.vcf
bgzip $PREFIX.$chr.vcf
done
```

#### ii. Get quality metrics for X and Y

It is useful to know the number of expected males and females

```{r count-males-females}
## Import a dataframe with disclosed sex of samples
## sex should be encoded 1=male, 2=female
sample_sex = read.delim(file = "exome_sample_sex.txt", header=T, sep="")

## Count number of females
females = sum(sample_sex$sex == 2) 
print(paste0("Number of reported females in dataset: ", females))

## Count number of males
males = sum(sample_sex$sex == 1) 
print(paste0("Number of reported males in dataset: ", males))

## Annotate the sample_metrics dataframe 
sample_metrics$reported_sex = sample_sex$sex[match(sample_metrics$INDV, sample_sex$sample)]
```

Calculate quality statistics

```{bash XY-preQC-check, eval=FALSE, cache=FALSE, include=FALSE}
PREFIX='redlat_exomes'
# Make a list of only males 
awk '{ if($2 == 1) {print $1}}' exome_sample_sex.txt > male.samples.txt
# Make a list of only females 
awk '{ if($2 == 2) {print $1}}' exome_sample_sex.txt > female.samples.txt

# These metrics are discussed in section 2.a
for chr in X
do
bcftools stats $PREFIX.$chr.vcf.gz > $PREFIX.$chr.vchk
vcftools --gzvcf $PREFIX.$chr.vcf.gz --FILTER-summary --out $PREFIX.$chr
vcftools --gzvcf $PREFIX.$chr.vcf.gz --missing-indv --out $PREFIX.$chr
vcftools --gzvcf $PREFIX.$chr.vcf.gz --depth --out $PREFIX.$chr
vcftools --gzvcf $PREFIX.$chr.vcf.gz --missing-site --out $PREFIX.$chr
  # For site based metrics we will split the analyses into male and female
  for sex in female male
  do
  vcftools --gzvcf $PREFIX.$chr.vcf.gz --keep $sex.samples.txt --missing-site --out $PREFIX.$chr.$sex
  vcftools --gzvcf $PREFIX.$chr.vcf.gz --keep $sex.samples.txt --site-mean-depth --out $PREFIX.$chr.$sex
  done
done
```

You should end with 14 files. For each chromosome ( X and Y) you get: I. Sample based: 1. missingness per sample (.imiss) 2. depth per sample (.idepth) II. Site based 1. Missingness per site a. Females (female.lmiss) b. Males (male.lmiss) 2. Mean depth per site a. Females (female.ldepth.mean) b. Males (female.ldepth.mean) 3. VQSR quality III. General stats (.vchk)

#### iii. Plot quality metrics for X

Plot quality metrics of all samples, females and males. This script will also generate a .txt file with the descriptive statistics

##### I. Sample based Metrics:

###### - Missingness per sample

```{r X-preQC-check-plots-X-missigness}
imiss_X = read.delim((paste0(PREFIX,".X.imiss")), header = T, sep = "")
imiss_X$sex = sample_sex$sex[match(imiss_X$INDV,sample_sex$sample)]
imiss_X_female = filter(imiss_X, sex == 2)
imiss_X_male = filter(imiss_X, sex == 1)

# Missingness rate (given as a fraction)
imiss_X_F_MISS_all = describe(imiss_X$F_MISS)
rownames(imiss_X_F_MISS_all) = c("sample_F_missingness_X_all")
imiss_X_F_MISS_female = describe(imiss_X_female$F_MISS)
rownames(imiss_X_F_MISS_female) = c("sample_F_missingness_X_female")
imiss_X_F_MISS_male = describe(imiss_X_male$F_MISS)
rownames(imiss_X_F_MISS_male) = c("sample_F_missingness_X_male")

# Missingness counts (given as an integer where N = samples x 2)
imiss_X_N_MISS_all = describe(imiss_X$N_MISS)
rownames(imiss_X_N_MISS_all) = c("sample_N_missingness_X_all")
imiss_X_N_MISS_female = describe(imiss_X_female$N_MISS)
rownames(imiss_X_N_MISS_female) = c("sample_N_missingness_X_female")
imiss_X_N_MISS_male = describe(imiss_X_male$N_MISS)
rownames(imiss_X_N_MISS_male) = c("sample_N_missingness_X_male")

missingness_sample_X_chromosome = bind_rows(imiss_X_F_MISS_all,
                                          imiss_X_F_MISS_female,
                                          imiss_X_F_MISS_male,
                                          imiss_X_N_MISS_all,
                                          imiss_X_N_MISS_female,
                                          imiss_X_N_MISS_male)

# Plot Missingness per sample
# Missingness counts helps to see distribution
png(file=(paste0(PREFIX,".X.sample_missingness.preQC.png")), 
    width=600, height=350)
boxplot(imiss_X$F_MISS, imiss_X_female$F_MISS, imiss_X_male$F_MISS,
        ylab="Missingness",
        xlab="Raw dataset", 
        main="Missingness rate per sample in X chromosome - preQC",
        sub = "Distribution according to disclosed sex",
        col = c("lavender","paleturquoise3", "lightgoldenrod1"),
        names = c("all", "Female", "Male"))
dev.off()
```

###### - Depth per sample

```{r X-preQC-check-plots-X-depth}
idepth_X = read.delim((paste0(PREFIX,".X.idepth")), header = T, sep = "")
idepth_X$sex = sample_sex$sex[match(idepth_X$INDV,sample_sex$sample)]
idepth_X_female = filter(idepth_X, sex == 2)
idepth_X_male = filter(idepth_X, sex == 1)

# Get some site depth stats:
# Mean depth per sample
idepth_X_stats_all = describe(idepth_X$MEAN_DEPTH)
rownames(idepth_X_stats_all) = c("sample_depth_X_all")
idepth_X_stats_female = describe(idepth_X_female$MEAN_DEPTH)
rownames(idepth_X_stats_female) = c("sample_depth_X_female")
idepth_X_stats_male = describe(idepth_X_male$MEAN_DEPTH)
rownames(idepth_X_stats_male) = c("sample_depth_X_male")

depth_sample_X_chromosome = bind_rows(idepth_X_stats_all,
                                      idepth_X_stats_female,
                                      idepth_X_stats_male)
# Plot depth per sample
png(file=(paste0(PREFIX,".X.sample_depth.preQC.png")), 
    width=600, height=350)
boxplot(idepth_X$MEAN_DEPTH, idepth_X_female$MEAN_DEPTH, idepth_X_male$MEAN_DEPTH,
        ylab="Mean sample depth",
        xlab="Raw dataset", 
        main="Mean sample depth in X chromosome - preQC",
        sub = "Distribution according to disclosed sex",
        col = c("lavender","paleturquoise3", "lightgoldenrod1"),
        names = c("all", "Female", "Male"))
dev.off()
```

###### - Data distribution of sample based metrics

```{r X-preQC-stats-sample-metrics}
stats_sample_X_chromosome = bind_rows(missingness_sample_X_chromosome,
                                      depth_sample_X_chromosome)
write.table(stats_sample_X_chromosome,
           (paste0(PREFIX,".X.sample_descriptive_statistics.preQC.txt")),            
           col.names = TRUE, 
           row.names = TRUE, 
           quote = FALSE,
           sep = '\t')

print(stats_sample_X_chromosome)
```

##### II Site based Metrics:

###### - Missingness per site

```{r X-preQC-stats-plots-site-missingness}
lmiss_X = read.delim((paste0(PREFIX,".X.lmiss")), header = T, sep = "")
# a. Females (female.lmiss)
lmiss_X_female = read.csv((paste0(PREFIX,".X.female.lmiss")), header = T, sep = "")
# b. Males (male.lmiss)
lmiss_X_male = read.csv((paste0(PREFIX,".X.male.lmiss")), header = T, sep = "")

# Get some site missingness stats:
# Missingness rate (given as a fraction)
lmiss_X_F_MISS = describe(lmiss_X$F_MISS)
rownames(lmiss_X_F_MISS) = c("site_F_missingness_X_all")
lmiss_X_female_F_MISS = describe(lmiss_X_female$F_MISS)
rownames(lmiss_X_female_F_MISS) = c("site_F_missingness_X_female")
lmiss_X_male_F_MISS = describe(lmiss_X_male$F_MISS)
rownames(lmiss_X_male_F_MISS) = c("site_F_missingness_X_male")

# Missingness counts (given as an integer
lmiss_X_N_MISS = describe(lmiss_X$N_MISS)
rownames(lmiss_X_N_MISS) = c("site_N_missingness_X_all")
lmiss_X_female_N_MISS = describe(lmiss_X_female$N_MISS)
rownames(lmiss_X_female_N_MISS) = c("site_N_missingness_X_female")
lmiss_X_male_N_MISS = describe(lmiss_X_male$N_MISS)
rownames(lmiss_X_male_N_MISS) = c("site_N_missingness_X_male")

missingness_site_X_chromosome = bind_rows(lmiss_X_F_MISS,
                                          lmiss_X_female_F_MISS,
                                          lmiss_X_male_F_MISS,
                                          lmiss_X_N_MISS,
                                          lmiss_X_female_N_MISS,
                                          lmiss_X_male_N_MISS)

# Plot Missingness per site
# Missingness counts helps to see distribution
png(file=(paste0(PREFIX,".X.site_missingness.preQC.boxplot.png")), 
    width=600, height=350)
boxplot(lmiss_X_female$N_MISS, lmiss_X_male$N_MISS,
        ylab="Missingness counts",
        xlab="Raw dataset", 
        main="Missingness counts per site in X chromosome - preQC",
        sub = "Distribution according to disclosed sex",
        col = c("paleturquoise3", "lightgoldenrod1"),
        names = c("Female", "Male"))
dev.off()

# Take a closer look 
png(file=(paste0(PREFIX,".X.site_missingness.preQC.boxplot.low_end.png")), 
    width=600, height=350)
boxplot(lmiss_X_female$F_MISS, lmiss_X_male$F_MISS,
        ylab="Missingness rate",
        xlab="Raw dataset", 
        main="Missingness rate per site in X chromosome - preQC",
        sub = "Distribution according to disclosed sex",
        col = c("paleturquoise3", "lightgoldenrod1"),
        names = c("Female", "Male"),
        ylim = c(0, 0.05))
dev.off()
```

###### - Depth per site

```{r x-preQC-stats-plots-site-depth}
# a. Females (female.ldepth.mean)
ldepth_X_female = read.csv((paste0(PREFIX,".X.female.ldepth.mean")), header = T, sep = "")
# b. Males (female.ldepth.mean)
ldepth_X_male = read.csv((paste0(PREFIX,".X.male.ldepth.mean")), header = T, sep = "")

# Get some stats
ldepth_X_female_depth = describe(ldepth_X_female$MEAN_DEPTH)
rownames(ldepth_X_female_depth ) = c("site_mean-depth_X_female")
ldepth_X_male_depth = describe(ldepth_X_male$MEAN_DEPTH)
rownames(ldepth_X_male_depth ) = c("site_mean-depth_X_male")

mean_site_depth_X_chromosome = bind_rows(ldepth_X_female_depth,
                                         ldepth_X_male_depth)

# Plot as boxplots
png(file=(paste0(PREFIX,".X.site_depth.preQC.png")), 
    width=600, height=350)
boxplot(ldepth_X_female$MEAN_DEPTH, ldepth_X_male$MEAN_DEPTH,
        ylab="Mean depth per site",
        xlab="Raw dataset", 
        main="Mean depth per site in X chromosome - preQC",
        sub = "Distribution according to disclosed sex",
        col = c("paleturquoise3", "lightgoldenrod1"),
        names = c("Female", "Male"))
dev.off()

# Take a zoom at the lower end
png(file=(paste0(PREFIX,".X.site_depth.preQC.low_end.png")), 
    width=600, height=350)
boxplot(ldepth_X_female$MEAN_DEPTH, ldepth_X_male$MEAN_DEPTH,
        ylab="Mean depth per variant",
        xlab="Raw dataset", 
        main="Mean depth per site in X chromosome - preQC - 80X and lower",
        sub = "Distribution according to disclosed sex",
        col = c("paleturquoise3", "lightgoldenrod1"),
        names = c("Female", "Male"),
        ylim = c(0, 80))
dev.off()
```

###### - Data distribution of site based metrics

```{r X-preQC-stats-plots-site-summary}
stats_site_X_chromosome = bind_rows(missingness_site_X_chromosome,
                                    mean_site_depth_X_chromosome)
write.table(stats_site_X_chromosome,
            (paste0(PREFIX,".X.site_descriptive_statistics.preQC.txt")),
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
print(stats_site_X_chromosome)
```

###### - VQSR (Variant Quality Score Recalibration)

```{r x-preQC-stats-plots-site-VQSR}
filter_X_raw = read.delim((paste0(PREFIX,".X.FILTER.summary")), header = T, sep = "")
print(filter_X_raw)

# Plot it
filter_X_raw %>%
  filter(!is.na(N_VARIANTS)) %>%
  arrange(N_VARIANTS) %>%
  mutate(FILTER=factor(FILTER, FILTER)) %>%
  ggplot( aes(x=FILTER, y=N_VARIANTS, label = N_VARIANTS) ) +
    geom_segment( aes(x=FILTER ,xend=FILTER, y=0, yend=N_VARIANTS), color="grey") +
    geom_point(size=3, color="#69b3a2") +
    geom_text(vjust=-1, size = 3) +
    coord_flip() +
    theme_minimal() +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position="none") +
    scale_y_continuous(name = "Number of variants") +
    labs(title = "VQSR in X chromosome - preQC")
ggsave(paste0(PREFIX,".X.site_VQSR.pre-QC.png"),
      width = 7,
      height = 4)  
```

In the ReDLat dataset you can observe that all samples (regardless chromosomal sex) have a mean sample depth above 20X and the missingness rate per sample for X is below 1%.

#### iv. Plot quality metrics for Y

##### I. Sample based Metrics:

###### - Missingness per sample

```{r y-preQC-stats-sample-missingess}
imiss_Y = read.delim((paste0(PREFIX,".Y.imiss")), header = T, sep = "")
imiss_Y$sex = sample_sex$sex[match(imiss_Y$INDV,sample_sex$sample)]
imiss_Y_female = filter(imiss_Y, sex == 2)
imiss_Y_male = filter(imiss_Y, sex == 1)

# Get some site missingness stats:
# Missingness rate (given as a fraction)
imiss_Y_F_MISS_female = describe(imiss_Y_female$F_MISS)
rownames(imiss_Y_F_MISS_female) = c("sample_F_missingness_Y_female")
imiss_Y_F_MISS_male = describe(imiss_Y_male$F_MISS)
rownames(imiss_Y_F_MISS_male) = c("sample_F_missingness_Y_male")

# Missingness counts (given as an integer)
imiss_Y_N_MISS_female = describe(imiss_Y_female$N_MISS)
rownames(imiss_Y_N_MISS_female) = c("sample_N_missingness_Y_female")
imiss_Y_N_MISS_male = describe(imiss_Y_male$N_MISS)
rownames(imiss_Y_N_MISS_male) = c("sample_N_missingness_Y_male")

missingness_sample_Y_chromosome = bind_rows(imiss_Y_F_MISS_female,
                                          imiss_Y_F_MISS_male,
                                          imiss_Y_N_MISS_female,
                                          imiss_Y_N_MISS_male)
# Plot Missingness per sample
# Missingness counts helps to see distribution
png(file=(paste0(PREFIX,".Y.sample_missingness.preQC.png")), 
    width=600, height=350)
boxplot(imiss_Y_female$F_MISS, imiss_Y_male$F_MISS,
        ylab="Missingness",
        xlab="Raw dataset", 
        main="Missingness rate per sample in Y chromosome - preQC",
        sub = "Distribution according to disclosed sex",
        col = c("paleturquoise3", "lightgoldenrod1"),
        names = c("Female", "Male"))
dev.off()
```

###### - Depth per sample

```{r y-preQC-stats-sample-depth}
idepth_Y = read.delim((paste0(PREFIX,".Y.idepth")), header = T, sep = "")
idepth_Y$sex = sample_sex$sex[match(idepth_Y$INDV,sample_sex$sample)]
idepth_Y_female = filter(idepth_Y, sex == 2)
idepth_Y_male = filter(idepth_Y, sex == 1)

# Get some site depth stats:
# Mean depth per sample
idepth_Y_stats_female = describe(idepth_Y_female$MEAN_DEPTH)
rownames(idepth_Y_stats_female) = c("sample_depth_Y_female")
idepth_Y_stats_male = describe(idepth_Y_male$MEAN_DEPTH)
rownames(idepth_Y_stats_male) = c("sample_depth_Y_male")

depth_sample_Y_chromosome = bind_rows(idepth_Y_stats_female,
                                      idepth_Y_stats_male)
# Plot depth per sample
png(file=(paste0(PREFIX,".Y.sample_depth.preQC.png")), 
    width=600, height=350)
boxplot(idepth_Y_female$MEAN_DEPTH, idepth_Y_male$MEAN_DEPTH,
        ylab="Mean sample depth",
        xlab="Raw dataset", 
        main="Mean sample depth in Y chromosome - preQC",
        sub = "Distribution according to disclosed sex",
        col = c("paleturquoise3", "lightgoldenrod1"),
        names = c("Female", "Male"))
dev.off()
```

###### - Quality statistics per sample

```{r y-preQC-sample-quality}
stats_sample_Y_chromosome = bind_rows(missingness_sample_Y_chromosome,
                                      depth_sample_Y_chromosome)
write.table(stats_sample_X_chromosome,
           (paste0(PREFIX,".Y.sample_descriptive_statistics.preQC.txt")),            
           col.names = TRUE, 
           row.names = TRUE, 
           quote = FALSE,
           sep = '\t')
print(stats_sample_Y_chromosome)
```

##### II Site based Metrics:

###### - Missingness per site

```{r y-preQC-stats-plots-site-missingness}
# a. Females (female.lmiss)
lmiss_Y_female = read.csv((paste0(PREFIX,".Y.female.lmiss")), header = T, sep = "")
# b. Males (male.lmiss)
lmiss_Y_male = read.csv((paste0(PREFIX,".Y.male.lmiss")), header = T, sep = "")

# Get some site missingness stats:
# Missingness rate (given as a fraction of samples missing the variant)
lmiss_Y_female_F_MISS = describe(lmiss_Y_female$F_MISS)
rownames(lmiss_Y_female_F_MISS) = c("site_F_missingness_Y_female")
lmiss_Y_male_F_MISS = describe(lmiss_Y_male$F_MISS)
rownames(lmiss_Y_male_F_MISS) = c("site_F_missingness_Y_male")

# Missingness counts (given as an integer where N = samples x 2)
lmiss_Y_female_N_MISS = describe(lmiss_Y_female$N_MISS)
rownames(lmiss_Y_female_N_MISS) = c("site_N_missingness_Y_female")
lmiss_Y_male_N_MISS = describe(lmiss_Y_male$N_MISS)
rownames(lmiss_Y_male_N_MISS) = c("site_N_missingness_Y_male")

missingness_site_Y_chromosome = bind_rows(lmiss_Y_female_F_MISS,
                                     lmiss_Y_male_F_MISS,
                                     lmiss_Y_female_N_MISS,
                                     lmiss_Y_male_N_MISS)

# Plot Missingness per site
# Missingness counts helps to see distribution
png(file=(paste0(PREFIX,".Y.site_missingness.preQC.boxplot.png")), 
    width=600, height=350)
boxplot(lmiss_Y_female$N_MISS, lmiss_Y_male$N_MISS,
        ylab="Missingness counts",
        xlab="Raw dataset", 
        main="Missingness counts per site in Y chromosome - preQC",
        sub = "Distribution according to disclosed sex",
        col = c("paleturquoise3", "lightgoldenrod1"),
        names = c("Female", "Male"))
dev.off()
# Take a closer look 
png(file=(paste0(PREFIX,".Y.site_missingness.preQC.boxplot.low_end.png")), 
    width=600, height=350)
boxplot(lmiss_Y_female$F_MISS, lmiss_Y_male$F_MISS,
        ylab="Missingness rate",
        xlab="Raw dataset", 
        main="Missingness rate per site in Y chromosome - preQC",
        sub = "Distribution according to disclosed sex",
        col = c("paleturquoise3", "lightgoldenrod1"),
        names = c("Female", "Male"))
dev.off()
```

###### - Depth per site

```{r y-preQC-stats-plots-site-depth}
# a. Females (female.ldepth.mean)
ldepth_Y_female = read.csv((paste0(PREFIX,".Y.female.ldepth.mean")), header = T, sep = "")
# b. Males (female.ldepth.mean)
ldepth_Y_male = read.csv((paste0(PREFIX,".Y.male.ldepth.mean")), header = T, sep = "")

# Get some stats
ldepth_Y_female_depth = describe(ldepth_Y_female$MEAN_DEPTH)
rownames(ldepth_Y_female_depth ) = c("site_mean-depth_Y_female")
ldepth_Y_male_depth = describe(ldepth_Y_male$MEAN_DEPTH)
rownames(ldepth_Y_male_depth ) = c("site_mean-depth_Y_male")

mean_site_depth_Y_chromosome = bind_rows(ldepth_Y_female_depth,
                                         ldepth_Y_male_depth)
# Plot as boxplots
png(file=(paste0(PREFIX,".Y.site_depth.preQC.png")), 
    width=600, height=350)
boxplot(ldepth_Y_female$MEAN_DEPTH, ldepth_Y_male$MEAN_DEPTH,
        ylab="Mean depth per site",
        xlab="Raw dataset", 
        main="Mean depth per site in Y chromosome - preQC",
        sub = "Distribution according to disclosed sex",
        col = c("paleturquoise3", "lightgoldenrod1"),
        names = c("Female", "Male"))
dev.off()
```

###### - Data distribution of site based metrics

```{r y-preQC-stats-plots-site-summary}
stats_site_Y_chromosome = bind_rows(missingness_site_Y_chromosome,
                                    mean_site_depth_Y_chromosome)
write.table(stats_site_X_chromosome,
            (paste0(PREFIX,".Y.site_descriptive_statistic.preQC.txt")),
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
print(stats_site_Y_chromosome)
```

###### - VQSR (Variant Quality Score Recalibration)

```{r y-preQC-stats-plots-site-VQSR}
filter_Y_raw = read.delim((paste0(PREFIX,".Y.FILTER.summary")), header = T, sep = "")
print(filter_Y_raw)

# Plot it
filter_Y_raw %>%
  filter(!is.na(N_VARIANTS)) %>%
  arrange(N_VARIANTS) %>%
  mutate(FILTER=factor(FILTER, FILTER)) %>%
  ggplot( aes(x=FILTER, y=N_VARIANTS, label = N_VARIANTS) ) +
    geom_segment( aes(x=FILTER ,xend=FILTER, y=0, yend=N_VARIANTS), color="grey") +
    geom_point(size=3, color="#69b3a2") +
    geom_text(vjust=-1, size = 3) +
    coord_flip() +
    theme_minimal() +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position="none") +
    scale_y_continuous(name = "Number of variants") +
    labs(title = "VQSR in Y chromosome - preQC")
ggsave(paste0(PREFIX,".Y.site_VQSR.pre-QC.png"))
```

Y chromosome sample statistics can also give us a hint of chromosomal vs disclosed/phenotypic sex discordances. Notice how there are only 83 markers in the Y chromosome exome!

#### v. Check chromosomal sex in plink

Since this is exome data, and the exomic portion of the Y chromosome is so small, we will can only use plink's algorithm to check for sex using the X chromosome.

```{bash check-sex-plink, eval=FALSE, cache=FALSE, include=FALSE}
PREFIX='redlat_exomes'
sample_sex=exome_sample_sex.txt
plink=~/bin/plink

# I. Import files into plink and add sex information more info in https://www.cog-genomics.org/plink/1.9/input#vcf
# half calls will be set as haploid
$plink --vcf $PREFIX.X.vcf.gz --keep-allele-order --vcf-half-call h --double-id --make-bed --out $PREFIX.X.plink

awk '{print $1, $1, $2}' $sample_sex > sex_plink.txt

$plink --bfile $PREFIX.X.plink --update-sex sex_plink.txt --make-bed --out $PREFIX.X.plink.sex

# II.Remove X pseudoautosomal region (if your data is hg19, use 'hg19')
$plink --bfile $PREFIX.X.plink.sex --split-x hg38 --make-bed --out $PREFIX.X.plink.sex.split-x

# III. Check if variants have an id in the second column of the bim file.
# Assign IDs to variants if they dont have it (maximum lenght of variant IDs is 20, which means longer insertions will be left without a variant ID)
$plink --bfile $PREFIX.X.plink.sex.split-x --set-missing-var-ids '@:#' --new-id-max-allele-len 20 --make-bed --out $PREFIX.X.plink.sex.split-x.id

# IV. Prune for Linkage Disequilibrium 
#(make sure that the variants have an ID in the bim file)
$plink --bfile $PREFIX.X.plink.sex.split-x.id --indep-pairphase 20000 2000 0.5
# this produces plink.prune.in and plink.prune.out

# V. Retain independent markers (in linkage equilibrium)
$plink --bfile $PREFIX.X.plink.sex.split-x.id --extract plink.prune.in --make-bed --out $PREFIX.X.plink.sex.split-x.id.LD

# VI. Check sex 
$plink --bfile $PREFIX.X.plink.sex.split-x.id.LD --check-sex 0.3 0.7 --out $PREFIX.X.plink.sex.split-x.id.LD.Xsex
# more info on https://www.cog-genomics.org/plink/1.9/basic_stats#check_sex
```

Take a look at the results

```{r plot-check-sex-plink, cache = TRUE}
#install.packages("ggbeeswarm")
library(ggbeeswarm)
#install.packages("scales")
#library(scales)

# Load X chromosome F coefficients calculated by Plink
X_file = read.delim(file = paste0(PREFIX,".X.plink.sex.split-x.id.LD.Xsex.sexcheck"), 
                    header = TRUE, 
                    sep = '')
# Plot these coefficients comparing males vs. females
# Roughly you expect Female to have an F coefficient < 0.2-0.3 and males have an F coefficient > 0.7-0.8
ggplot(X_file, 
       aes(x = factor(PEDSEX,
                      labels = c("Male", "Female")),
           y = F,
           color = PEDSEX)) +
    geom_quasirandom(alpha = 0.7,
                   size = 1.5) + 
    labs(title = "Chromosomal sex assignement in samples",
         x = "Disclosed sex",
         y = "F coefficient X chromosome") +
  theme_minimal() +
  theme(legend.position = "none")
ggsave((paste0(PREFIX,".xy.chromosomal_sex_assignement.png")), width = 8, height = 5)
```

We can add the missingness of Y chromosome as another column to select the samples that will be taken out. Individuals that fail sex check should be added to 'flagged_samples' file

```{r determine-sex-check-fails, cache = TRUE}
X_file = read.delim(file = paste0(PREFIX,".X.plink.sex.split-x.id.LD.Xsex.sexcheck"), 
                    header = TRUE, 
                    sep = '')

# Create a new column in 'sample_metrics' with the X chromosome F coefficients
sample_metrics$X_Fcoeff = X_file$F[match(sample_metrics$INDV, X_file$IID)]
sample_metrics$X_STATUS = X_file$STATUS[match(sample_metrics$INDV, X_file$IID)]
sample_metrics$Y_F_MISS = imiss_Y$F_MISS[match(sample_metrics$INDV, imiss_Y$INDV)]

# Generate a file up to this point
write.table(sample_metrics,
           (paste0(PREFIX,".sample_based_all_statistics.preQC.txt")),
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE)

# I am taking a combination of X F coefficient and missingness in Y chromosome to filter samples that DEFINITELY fail check sex
sex_fail = filter(sample_metrics, 
                        (reported_sex == "1" & X_Fcoeff < 0.3 & Y_F_MISS > 0.1 |
                        reported_sex == "2" & X_Fcoeff > 0.7 & Y_F_MISS < 0.7)) 
 
# Samples that were flagged by plink but dont meet the Y chromosome condition are labeled as Warning
sex_warning = filter(sample_metrics, 
                        (reported_sex == "1" & X_Fcoeff < 0.7 & Y_F_MISS < 0.1 |
                        reported_sex == "2" & X_Fcoeff > 0.3 & Y_F_MISS > 0.7)) 

# Add the sex-fail samples to the flagged sample list
flagged_samples_sex = select(sex_fail, INDV)
```

All the samples of `sex_fail` and `sex_warning` should be double checked with the submitting center to determine if there was a sex labeling error.

Organize your files

```{bash clean-up-2xy}
PREFIX='redlat_exomes'

mkdir 2.xy_quality_report
mkdir 2.xy_quality_report/x
mkdir 2.xy_quality_report/y
mv $PREFIX.X.* 2.xy_quality_report/x
mv $PREFIX.Y.* 2.xy_quality_report/y
mv $PREFIX.xy* 2.xy_quality_report/
```

## 3: Genotype quality control

### 3.a Autosomes

Retain genotypes with Depth \>= 10 & Quality \>= 20 and check how this affected the missingness rates. Flagged samples in Autosome quality report and those with discrepancies in the expected vs. observed chromosomal sex should be removed.

```{bash autosomes-genotype-QC, eval=FALSE, cache=FALSE, include=FALSE}
PREFIX='redlat_exomes'
DP=10
GQ=20
FLAGGED="samples_to_remove.txt"

# I. Retain genotypes with Depth >= 20 & Quality >= 20
# Other genotypes will be set to missing
vcftools --gzvcf $PREFIX.autosomes.vcf.gz --minDP $DP --minGQ $GQ --keep $FLAGGED --recode --recode-INFO-all --out $PREFIX.autosomes.DP$DP.GQ$GQ
#--minDP <float> Includes only genotypes greater than or equal to the "--minDP"value. DP <10 are set as missing
#This option requires that the "DP" FORMAT tag is specified for all sites.
#--minGQ <float> Excludes all genotypes with a quality below the threshold specified. 
#This option requires that the "GQ" FORMAT tag is specified for all sites.

mv $PREFIX.autosomes.DP$DP.GQ$GQ.recode.vcf $PREFIX.autosomes.DP$DP.GQ$GQ.vcf
bgzip $PREFIX.autosomes.DP$DP.GQ$GQ.vcf
```

Get the statistics for the new dataset

```{bash autosome-GT-stats, eval=FALSE, cache=FALSE, include=FALSE}
PREFIX='redlat_exomes'
DP=10
GQ=20

## I. Data description 
bcftools stats $PREFIX.autosomes.DP$DP.GQ$GQ.vcf.gz > $PREFIX.autosomes.DP$DP.GQ$GQ.vchk.txt
# General distribution of depth, missingness, heterozygosity
vcftools --gzvcf $PREFIX.autosomes.DP$DP.GQ$GQ.vcf.gz --FILTER-summary --out $PREFIX.autosomes.DP$DP.GQ$GQ
# Generates a summary of the number of SNPs and Ts/Tv ratio for each FILTER category. 
# The output file has the suffix ".FILTER.summary"

## II. Individual missingness
vcftools --gzvcf $PREFIX.autosomes.DP$DP.GQ$GQ.vcf.gz --missing-indv --out $PREFIX.autosomes.DP$DP.GQ$GQ
# Generates a file reporting the missingness on a per-individual basis. 
# The output file has the suffix ".imiss".
# Individuals whose missingness is >10% should be added to 'flagged_samples' file

## III. Individual depth 
vcftools --gzvcf $PREFIX.autosomes.DP$DP.GQ$GQ.vcf.gz --depth --out $PREFIX.autosomes.DP$DP.GQ$GQ
# Generates a file containing the mean depth per individual. 
# This file has the suffix ".idepth".
# Individuals whose mean depth is <20 should be added to 'flagged_samples' file

## VI. Individual heterozygosity 
vcftools --gzvcf $PREFIX.autosomes.DP$DP.GQ$GQ.vcf.gz --het --out $PREFIX.autosomes.DP$DP.GQ$GQ
# Inbreeding coefficient, F, is estimated for each individual using a method of moments. 
# The resulting file has the suffix ".het"
# Individuals whose heterozygosity deviated more than 3 SD from the main should be added to 'flagged_samples' file

## V. Site missingness 
vcftools --gzvcf $PREFIX.autosomes.DP$DP.GQ$GQ.vcf.gz --missing-site --out $PREFIX.autosomes.DP$DP.GQ$GQ
# Generates a file reporting the missingness on a per-site basis. 
# The file has the suffix ".lmiss".

## VI. Site depth
vcftools --gzvcf $PREFIX.autosomes.DP$DP.GQ$GQ.vcf.gz --site-mean-depth --out $PREFIX.autosomes.DP$DP.GQ$GQ
# Generates a file containing the mean depth per site across all individuals. 
# This output file has the suffix ".ldepth.mean"
```

Plot missingness rates per site and per sample comparing raw dataset and after genotype quality control. Declare the same GQ and DP values

##### I. Sample based Metrics:

###### - Missingness per sample

```{r autosome-GT-stats-sample-missingess}
DP_autosomes = 10
GQ_autosomes = 20

imiss_GT = read.delim((paste0(PREFIX,".autosomes.DP",DP_autosomes,".GQ",GQ_autosomes,".imiss")), header = T, sep = "")
# Get stats of missingness rate per sample (given as a fraction)
imiss_GT_F_MISS = describe(imiss_GT$F_MISS) 
rownames(imiss_GT_F_MISS) = c("sample_missingness_F_GT")
# Individuals whose missingness is >10% should be identified
high_missingness = filter(imiss_GT, F_MISS>0.1)
# Save it. It may come useful for debugging
write.table(high_missingness, 
            (paste0(PREFIX,".autosomes.FLAGGED_sample_high_missingness.preQC.GT.txt")),
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')

# Missingness counts per sample (given as an integer where N = variants)
imiss_GT_N_MISS = describe(imiss_GT$N_MISS) 
rownames(imiss_GT_N_MISS) = c("sample_missingness_N_GT")

missingness_sample = bind_rows(imiss_GT_F_MISS,
                               imiss_GT_N_MISS)

# Plot Missingness rate per sample
png(file=(paste0(PREFIX,".autosomes.sample_missingness_preQC.GT.png")), 
    width=600, height=350)
boxplot(imiss$F_MISS, imiss_GT$F_MISS,
        ylab = "Missingness rate",
        names = c("Raw dataset", "After genotype QC"),
        main = "Missingness rate per sample - Autosomes genotype QC", 
        col = c("paleturquoise3", "lavender"))
dev.off()
```

###### - Depth per sample

```{r autosome-GT-stats-sample-depth}
idepth_GT = read.delim(paste0(PREFIX,".autosomes.DP",DP_autosomes,".GQ",GQ_autosomes,".idepth"), header = T, sep = "")

# Get some site depth stats:
depth_sample_GT = describe(idepth_GT$MEAN_DEPTH)
rownames(depth_sample) = c("sample_depth")
# Individuals whose mean depth is <20 should be identified
low_mean_depth = filter(idepth_GT, MEAN_DEPTH<20)
write.table(low_mean_depth, 
            (paste0(PREFIX,".autosomes.FLAGGED_sample_low_mean_depth.preQC.GT.txt")),
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')

# Plot depth per sample
png(file=(paste0(PREFIX,".autosomes.sample_depth.preQC.GT.png")), 
    width=600, height=350)
boxplot(idepth$MEAN_DEPTH, idepth_GT$MEAN_DEPTH,
     xlab="Mean Depth",
     ylab="Samples", 
     names = c("Raw dataset", "After genotype QC"),
     main="Mean Depth per sample - Autosomes genotype QC", 
     col = c("paleturquoise3", "lavender"),
     breaks=50)
dev.off()
```

###### - Individual heterozygosity

```{r autosome-GT-stats-sample-heterozygosity}
het_GT = read.delim((paste0(PREFIX,".autosomes.DP",DP_autosomes,".GQ",GQ_autosomes,".het")), header = T, sep = "")

### Get some stats:
heterozygosity_sample = describe(het_GT$F)
rownames(heterozygosity_sample) = c("sample_heterozygosity_F")
### Identify the values for +3 and -3 standard deviations
heterozygosity_low_limit = mean(het_GT$F)-(3*(sd(het_GT$F)))
heterozygosity_high_limit = mean(het_GT$F)+(3*(sd(het_GT$F)))

### Plot heterozygosity per sample [we wont be comparing pre vs podt GT here]
png(file=(paste0(PREFIX,".autosomes.sample_heterozygosity.preQC.GT.png")), 
    width=600, height=350)
hist(het_GT$F,  
     freq=TRUE, 
     xlab="Heterozygosity F coefficient",  
     ylab="Samples", 
     main="Heterozygosity rate per sample - Autosomes genotype QC",
     col="lavender",
     breaks=50)
abline(v = (heterozygosity_low_limit), col="red")
abline(v = (heterozygosity_high_limit), col="red")
abline(v = (mean(het$F)), col="blue") 
legend("topleft",
       c("+/-3 SD","mean"),
       col=c("red","blue"),
       pch=16)
dev.off()

### Individuals whose heterozygosity deviated more than 3 SD from the mean should be identified
het_outlier_low = filter(het_GT, F<heterozygosity_low_limit)
het_outlier_high = filter(het_GT, F>heterozygosity_high_limit)
het_outlier_both = bind_rows(het_outlier_low,
                             het_outlier_high)
write.table(het_outlier_both, 
            (paste0(PREFIX,".autosomes.FLAGGED_sample_heterozygosity_outliers.preQC.GT.txt")),
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
```

###### - Data distribution of sample based metrics

```{r autosome-preQC-stats-sample-metrics}
stats_sample = bind_rows(missingness_sample,
                         depth_sample,
                         heterozygosity_sample)
write.table(stats_sample,
            (paste0(PREFIX,".autosomes.sample_descriptive_statistics_preQC.GTtxt")),
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
print(stats_sample)

## Create a dataframe  with the samples that failed the autosome quality thresholds
high_missingness_id = select(high_missingness, INDV)
low_mean_depth_id = select(low_mean_depth, INDV)
het_outlier_low_id = select(het_outlier_low, INDV)
het_outlier_high_id = select(het_outlier_high, INDV)
flagged_samples_autosomes_GT = bind_rows(high_missingness_id,
                                      low_mean_depth_id,
                                      het_outlier_low_id, 
                                      het_outlier_high_id)

### Save it. It may come useful for debugging
write.table(flagged_samples_autosomes_GT, 
            (paste0(PREFIX,".autosomes.FLAGGED_sample_allstats.preQC.GT.txt")),
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
```

###### - Quality statistics per sample

```{r autosome-GT-sample-quality}
# Annotate the sample_metrics file with the new Missingness values
sample_metrics$missingness_GT = imiss_GT$F_MISS[match(sample_metrics$INDV, imiss_GT$INDV)]
# Generate a file
write.table(sample_metrics,
           (paste0(PREFIX,".sample_based_all_statistics.preQC.GT.txt")),
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE)

# Save Autosomes - genotype QC stats
stats_missingness_autosomeGT = bind_rows(imiss_GT_F_MISS,
                                         lmiss_GT_F_MISS)
write.table(stats_missingness_autosomeGT,
            "stats_missingness_autosomeGT.txt", 
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
print(stats_missingness_autosomeGT)
write.table(stats_sites,
           (paste0(PREFIX,".autosomes.site_descriptive_statistics.preQC.txt")),
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')

#Flag samples that are missing more than 10%
high_missingness_GT = filter(imiss_GT, F_MISS>0.1)
## Add the sex-fail samples to the flagged sample list
flagged_samples_autosomes_GT= select(high_missingness_GT, INDV)
```

##### II Site based Metrics:

###### - Missingness per site

```{r autosome-GT-stats-plots-site-missingness}
lmiss_GT = read.delim((paste0(PREFIX,".autosomes.DP",DP_autosomes,".GQ",GQ_autosomes,".lmiss")), header = T, sep = "")

# Get missingness rate per site stats (given as a fraction
lmiss_GT_F_MISS = describe(lmiss$F_MISS)
rownames(lmiss_GT_F_MISS) = c("site_missingness_F_GT")

# Plot Missingness per site
png(file=(paste0(PREFIX,".autosomes.site_missingness_preQC.GT.png")), 
    width=600, height=350)
boxplot(lmiss$F_MISS, lmiss_GT$F_MISS,
        ylab ="Missingness rate",
        names = c("Raw dataset", "After genotype QC"),
        main = "Missingness rate per site - Autosomes genotype QC", 
        col = c("paleturquoise3", "lavender"))
dev.off()

png(file=(paste0(PREFIX,".autosomes.site_missingness_preQC.GT.low_end.png")), 
    width=600, height=350)
boxplot(lmiss$F_MISS, lmiss_GT$F_MISS,
        ylab ="Missingness rate",
        names = c("Raw dataset", "After genotype QC"),
        main = "Missingness rate per site < 0.01 - Autosomes genotype QC", 
        col = c("paleturquoise3", "lavender"),
        ylim = c(0,0.01))
```

###### - Depth per site

```{r autosome-GT-stats-plots-site-depth}
ldepth.mean = read.delim((paste0(PREFIX,".autosomes.DP",DP_autosomes,".GQ",GQ_autosomes,".lmiss")), header = T, sep = "")

### Get basic stats
depth_site = describe(ldepth.mean$MEAN_DEPTH)
rownames(depth_site) = c("site_mean_depth")

### Plot
png(file=(paste0(PREFIX,".autosomes.site_depth.preQC.GT.png")), 
    width=600, height=350)
hist(ldepth.mean$MEAN_DEPTH,
     xlab="Mean depth",
     ylab="Sites", 
     main="Mean depth per site - Autosomes preQC", 
     col="paleturquoise3")
dev.off()

### Take a zoom at the lower end
png(file=(paste0(PREFIX,".autosomes.site.depth_preQC.low_end.png")), 
    width=600, height=350)
boxplot(ldepth.mean$MEAN_DEPTH,
        ylab="Mean depth per variant",
        xlab="Raw dataset", 
        main="Mean depth per site in Autosomes - preQC - 80X and lower",
        col = c("paleturquoise3"),
        ylim = c(0, 80))
dev.off()
```

```{bash clean-up-autosomes}
PREFIX='redlat_exomes'
mkdir 2.autosomes_quality_report
mv $PREFIX.autosomes.* 2.autosomes_quality_report
```

### 3.b X Chromosome

When we checked the missingness rates and depth rates in x chromosomes, the mean values were similar to those of autosomes (both for male and female samples), therefore it is reasonable to perform genotype quality control using GQ = 20 and testing for DP=10 and DP=20

```{bash X-genotype-QC, eval=FALSE, cache=FALSE, include=FALSE}
PREFIX='redlat_exomes'
DP=10
GQ=20

for DP in 10 20
do
## I. Retain genotypes with Quality >= 20
vcftools --gzvcf $PREFIX.X.vcf.gz --minDP $DP --minGQ $GQ --recode --recode-INFO-all --out $PREFIX.X.DP$DP.GQ$GQ
mv $PREFIX.X.DP$DP.GQ$GQ.recode.vcf $PREFIX.X.DP$DP.GQ$GQ.vcf
bgzip $PREFIX.X.DP$DP.GQ$GQ.vcf
## II. Check how this affected the missingness rates
## Per sample
vcftools --gzvcf $PREFIX.X.DP$DP.GQ$GQ.vcf.gz --missing-indv --out $PREFIX.X.DP$DP.GQ$GQ
# Per site
vcftools --gzvcf $PREFIX.X.DP$DP.GQ$GQ.vcf.gz --missing-site --out $PREFIX.X.DP$DP.GQ$GQ
#Stratified per sex
  for sex in female male
  do
  vcftools --gzvcf $PREFIX.X.DP$DP.GQ$GQ.vcf.gz --keep $sex.samples.txt --missing-site --out $PREFIX.X.$sex.DP$DP.GQ$GQ
  #vcftools --gzvcf $PREFIX.X.DP$DP.GQ$GQ.vcf.gz --keep $sex.samples.txt --site-mean-depth --out $PREFIX.X.$sex.DP$DP.GQ$GQ
  done
done
```

Make some comparative plots

```{r compare-missingness-rates-X-GT, chache = TRUE}
#Declare the same GQ and DP values as in the filtering
DP_1 = 10
DP_2 = 20
GQ_X = 20

# We had previously generated an X.lmiss and X.imiss files of raw dataset

## I. Missingness per sample. Stratified by male-female and the two depth filters
#DP_1
imiss_GT_X1 = read.delim((paste0(PREFIX,".X.DP",DP_1,".GQ",GQ_X,".imiss")), header = T, sep = "")
imiss_GT_X1$sex = sample_sex$sex[match(imiss_GT_X1$INDV,sample_sex$sample)]
imiss_GT_X1_female = filter(imiss_GT_X1, sex == 2)
imiss_GT_X1_male = filter(imiss_GT_X1, sex == 1)
#DP_2
imiss_GT_X2 = read.delim((paste0(PREFIX,".X.DP",DP_2,".GQ",GQ_X,".imiss")), header = T, sep = "")
imiss_GT_X2$sex = sample_sex$sex[match(imiss_GT_X2$INDV,sample_sex$sample)]
imiss_GT_X2_female = filter(imiss_GT_X2, sex == 2)
imiss_GT_X2_male = filter(imiss_GT_X2, sex == 1)

# Get stats
#DP_1
imiss_GT_X1_F_MISS_all = describe(imiss_GT_X1$F_MISS) 
rownames(imiss_GT_X1_F_MISS_all) = c(paste0("sample_chrX_missingness_F_GT.DP",DP_1,"_all"))
imiss_GT_X1_F_MISS_female = describe(imiss_GT_X1_female$F_MISS) 
rownames(imiss_GT_X1_F_MISS_female) = c(paste0("sample_chrX_missingness_F_GT.DP",DP_1,"_female"))
imiss_GT_X1_F_MISS_male = describe(imiss_GT_X1_male$F_MISS) 
rownames(imiss_GT_X1_F_MISS_male) = c(paste0("sample_chrX_missingness_F_GT.DP",DP_1,"_male"))
#DP_2
imiss_GT_X2_F_MISS_all = describe(imiss_GT_X2$F_MISS) 
rownames(imiss_GT_X2_F_MISS_all) = c(paste0("sample_chrX_missingness_F_GT.DP",DP_1,"_all"))
imiss_GT_X2_F_MISS_female = describe(imiss_GT_X2_female$F_MISS) 
rownames(imiss_GT_X2_F_MISS_female) = c(paste0("sample_chrX_missingness_F_GT.DP",DP_1,"_female"))
imiss_GT_X2_F_MISS_male = describe(imiss_GT_X2_male$F_MISS) 
rownames(imiss_GT_X2_F_MISS_male) = c(paste0("sample_chrX_missingness_F_GT.DP",DP_1,"_male"))

# Save Chromosome X - genotype QC stats
stats_missingness_sample_X_chr_GT = bind_rows(imiss_X_F_MISS_all,
                                          imiss_X_F_MISS_female,
                                          imiss_X_F_MISS_male,
                                          imiss_GT_X1_F_MISS_all,
                                          imiss_GT_X1_F_MISS_female,
                                          imiss_GT_X1_F_MISS_male,
                                          imiss_GT_X2_F_MISS_all,
                                          imiss_GT_X2_F_MISS_female,
                                          imiss_GT_X2_F_MISS_male)
write.table(stats_missingness_sample_X_chr_GT,
            "stats_missingness_sample_X_chr_GT.txt", 
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
print(stats_missingness_sample_X_chr_GT)

### Plot Missingness rate per sample
boxplot(imiss_X$F_MISS, imiss_GT_X1$F_MISS,imiss_GT_X2$F_MISS,
        imiss_X_female$F_MISS, imiss_GT_X1_female$F_MISS,imiss_GT_X2_female$F_MISS,
        imiss_X_male$F_MISS, imiss_GT_X1_male$F_MISS,imiss_GT_X2_male$F_MISS,
        ylab = "Missingness rate",
        names = c("Raw\nall", paste0("DP",DP_1,"\nall"), paste0("DP",DP_2,"\nall"),
                  "Raw\nfemale", paste0("DP",DP_1,"\nfemale"), paste0("DP",DP_2,"\nfemale"),
                  "Raw\nmale", paste0("DP",DP_1,"\nmale"), paste0("DP=",DP_2,"\nmale")),
        main = "Missingness rate per sample - Chromosome X genotype QC", 
        col = c("paleturquoise3", "lavender", "lightgoldenrod1"),
        cex.axis=0.8)


## II. Missingness per site. Stratified by male-female and the two depth filters
#DP_1
lmiss_GT_X1 = read.delim((paste0(PREFIX,".X.DP",DP_1,".GQ",GQ_X,".lmiss")), header = T, sep = "")
lmiss_GT_X1_female = read.delim((paste0(PREFIX,".X.female.DP",DP_1,".GQ",GQ_X,".lmiss")), header = T, sep = "")
lmiss_GT_X1_male = read.delim((paste0(PREFIX,".X.male.DP",DP_1,".GQ",GQ_X,".lmiss")), header = T, sep = "")

#DP_2
lmiss_GT_X2 = read.delim((paste0(PREFIX,".X.DP",DP_2,".GQ",GQ_X,".lmiss")), header = T, sep = "")
lmiss_GT_X2_female = read.delim((paste0(PREFIX,".X.female.DP",DP_2,".GQ",GQ_X,".lmiss")), header = T, sep = "")
lmiss_GT_X2_male = read.delim((paste0(PREFIX,".X.male.DP",DP_2,".GQ",GQ_X,".lmiss")), header = T, sep = "")

# Get stats
#DP_1
lmiss_GT_X1_F_MISS_all = describe(lmiss_GT_X1$F_MISS) 
rownames(lmiss_GT_X1_F_MISS_all) = c(paste0("site_chrX_missingness_F_GT.DP",DP_1,"_all"))
lmiss_GT_X1_F_MISS_female = describe(lmiss_GT_X1_female$F_MISS) 
rownames(lmiss_GT_X1_F_MISS_female) = c(paste0("site_chrX_missingness_F_GT.DP",DP_1,"_female"))
lmiss_GT_X1_F_MISS_male = describe(lmiss_GT_X1_male$F_MISS) 
rownames(lmiss_GT_X1_F_MISS_male) = c(paste0("site_chrX_missingness_F_GT.DP",DP_1,"_male"))
#DP_2
lmiss_GT_X2_F_MISS_all = describe(lmiss_GT_X2$F_MISS) 
rownames(lmiss_GT_X2_F_MISS_all) = c(paste0("site_chrX_missingness_F_GT.DP",DP_1,"_all"))
lmiss_GT_X2_F_MISS_female = describe(lmiss_GT_X2_female$F_MISS) 
rownames(lmiss_GT_X2_F_MISS_female) = c(paste0("site_chrX_missingness_F_GT.DP",DP_1,"_female"))
lmiss_GT_X2_F_MISS_male = describe(lmiss_GT_X2_male$F_MISS) 
rownames(lmiss_GT_X2_F_MISS_male) = c(paste0("site_chrX_missingness_F_GT.DP",DP_1,"_male"))

# Save Chromosome X - genotype QC stats
stats_missingness_site_X_chr_GT = bind_rows(lmiss_X_F_MISS,
                                          lmiss_X_female_F_MISS,
                                          lmiss_X_male_F_MISS,
                                          lmiss_GT_X1_F_MISS_all,
                                          lmiss_GT_X1_F_MISS_female,
                                          lmiss_GT_X1_F_MISS_male,
                                          lmiss_GT_X2_F_MISS_all,
                                          lmiss_GT_X2_F_MISS_female,
                                          lmiss_GT_X2_F_MISS_male)
write.table(stats_missingness_site_X_chr_GT,
            "stats_missingness_site_X_chr_GT.txt", 
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
print(stats_missingness_site_X_chr_GT)

# Plot
boxplot(lmiss_X$F_MISS, lmiss_GT_X1$F_MISS,lmiss_GT_X2$F_MISS,
        lmiss_X_female$F_MISS, lmiss_GT_X1_female$F_MISS,lmiss_GT_X2_female$F_MISS,
        lmiss_X_male$F_MISS, lmiss_GT_X1_male$F_MISS,lmiss_GT_X2_male$F_MISS,
        ylab = "Missingness rate",
        names = c("Raw\nall", paste0("DP",DP_1,"\nall"), paste0("DP",DP_2,"\nall"),
                  "Raw\nfemale", paste0("DP",DP_1,"\nfemale"), paste0("DP",DP_2,"\nfemale"),
                  "Raw\nmale", paste0("DP",DP_1,"\nmale"), paste0("DP=",DP_2,"\nmale")),
        main = "Missingness rate per site - Chromosome X genotype QC", 
        col = c("paleturquoise3", "lavender", "lightgoldenrod1"),
        cex.axis=0.8)

boxplot(lmiss_X$F_MISS, lmiss_GT_X1$F_MISS,lmiss_GT_X2$F_MISS,
        lmiss_X_female$F_MISS, lmiss_GT_X1_female$F_MISS,lmiss_GT_X2_female$F_MISS,
        lmiss_X_male$F_MISS, lmiss_GT_X1_male$F_MISS,lmiss_GT_X2_male$F_MISS,
        ylab = "Missingness rate",
        names = c("Raw\nall", paste0("DP",DP_1,"\nall"), paste0("DP",DP_2,"\nall"),
                  "Raw\nfemale", paste0("DP",DP_1,"\nfemale"), paste0("DP",DP_2,"\nfemale"),
                  "Raw\nmale", paste0("DP",DP_1,"\nmale"), paste0("DP=",DP_2,"\nmale")),
        main = "Missingness rate per site - Chromosome X genotype QC <0.1", 
        col = c("paleturquoise3", "lavender", "lightgoldenrod1"),
        cex.axis=0.8,
        ylim = c(0,0.1))


# Annotate the sample_metrics file with the new Missingness values
sample_metrics$missingness_Xchrom_GT_DP1 = imiss_GT_X1$F_MISS[match(sample_metrics$INDV, imiss_GT_X1$INDV)]
names(sample_metrics)[names(sample_metrics) == 'missingness_Xchrom_GT_DP1'] = paste0('missingness_Xchrom_GT_DP', DP_1)
sample_metrics$missingness_Xchrom_GT_DP2 = imiss_GT_X2$F_MISS[match(sample_metrics$INDV, imiss_GT_X2$INDV)]
names(sample_metrics)[names(sample_metrics) == 'missingness_Xchrom_GT_DP2'] = paste0('missingness_Xchrom_GT_DP', DP_2)

## Generate a file
write.table(sample_metrics,
           (paste0(PREFIX,".autosomal.xy.sample_based_all_statistics.preQC.GT.X.txt")),
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE,
            sep = '\t')

```

### 3.c Y Chromosome

We are testing two depth values, 10 and 15. The filtering is done in the vcf including both males and females, but missingness stats are based on males only

```{bash Y-genotype-QC, eval=FALSE, cache=FALSE, include=FALSE}
PREFIX='redlat_exomes'
#DP=10
GQ=20

for DP in 10 15
do
## I. Retain genotypes with Depth >= x & Quality >= 20
vcftools --gzvcf $PREFIX.Y.vcf.gz --minDP $DP --minGQ $GQ --recode --recode-INFO-all --out $PREFIX.Y.DP$DP.GQ$GQ
mv $PREFIX.Y.DP$DP.GQ$GQ.recode.vcf $PREFIX.Y.DP$DP.GQ$GQ.vcf
bgzip $PREFIX.Y.DP$DP.GQ$GQ.vcf
## II. Check how this affected the missingness rates
# Per site
vcftools --gzvcf $PREFIX.Y.DP$DP.GQ$GQ.vcf.gz --keep male.samples.txt --missing-site --out $PREFIX.Y.DP$DP.GQ$GQ
## Per sample
vcftools --gzvcf $PREFIX.Y.DP$DP.GQ$GQ.vcf.gz --keep male.samples.txt --missing-indv --out $PREFIX.Y.DP$DP.GQ$GQ
done
```

Make some comparative plots

```{r compare-missingness-rates-Y-GT}
#Declare the same GQ and DP values as in the filtering
DP_1 = 10
DP_2 = 15
GQ_Y = 20

# We had previously generated an Y.lmiss and Y.imiss files of raw dataset
## I. Missingness per sample
imiss_GT_Y1 = read.delim((paste0(PREFIX,".Y.DP",DP_1,".GQ",GQ_Y,".imiss")), header = T, sep = "")
imiss_GT_Y2 = read.delim((paste0(PREFIX,".Y.DP",DP_2,".GQ",GQ_Y,".imiss")), header = T, sep = "")

# Get stats
imiss_GT_Y1_F_MISS = describe(imiss_GT_Y1$F_MISS) 
rownames(imiss_GT_Y1_F_MISS) = c(paste0("sample_chrY_missingness_F_GT.DP",DP_1))
imiss_GT_Y2_F_MISS = describe(imiss_GT_Y2$F_MISS) 
rownames(imiss_GT_Y2_F_MISS) = c(paste0("sample_chrY_missingness_F_GT.DP",DP_2))

# Plot
boxplot(imiss_Y_male$F_MISS, imiss_GT_Y1$F_MISS,imiss_GT_Y2$F_MISS,
        ylab = "Missingness rate",
        names = c("Raw dataset", paste0("Genotype QC DP=",DP_1), paste0("Genotype QC DP=",DP_2)),
        main = "Missingness rate per sample - Chromosome Y genotype QC", 
        col = c("paleturquoise3", "lavender", "lightgoldenrod1"))

## II. Missingness per site
lmiss_GT_Y1 = read.delim((paste0(PREFIX,".Y.DP",DP_1,".GQ",GQ_Y,".lmiss")), header = T, sep = "")
lmiss_GT_Y2 = read.delim((paste0(PREFIX,".Y.DP",DP_2,".GQ",GQ_Y,".lmiss")), header = T, sep = "")
# Get stats
lmiss_GT_Y1_F_MISS = describe(lmiss_GT_Y1$F_MISS) 
rownames(lmiss_GT_Y1_F_MISS) = c(paste0("site_chrY_missingness_F_GT.DP",DP_1))
lmiss_GT_Y2_F_MISS = describe(lmiss_GT_Y2$F_MISS) 
rownames(lmiss_GT_Y2_F_MISS) = c(paste0("site_chrY_missingness_F_GT.DP",DP_2))

# plot
boxplot(lmiss_Y_male$F_MISS, lmiss_GT_Y1$F_MISS, lmiss_GT_Y2$F_MISS,
        ylab ="Missingness rate",
        names = c("Raw dataset", paste0("Genotype QC DP=",DP_1), paste0("Genotype QC DP=",DP_2)),
        main = "Missingness rate per site - Chromosome Y genotype QC", 
        col = c("paleturquoise3", "lavender", "lightgoldenrod1"))

# Annotate the sample_metrics file with the new Missingness values
sample_metrics$missingness_Ychrom_GT_DP1 = imiss_GT_Y1$F_MISS[match(sample_metrics$INDV, imiss_GT_Y1$INDV)]
names(sample_metrics)[names(sample_metrics) == 'missingness_Ychrom_GT_DP1'] = paste0('missingness_Ychrom_GT_DP', DP_1)
sample_metrics$missingness_Ychrom_GT_DP2 = imiss_GT_Y2$F_MISS[match(sample_metrics$INDV, imiss_GT_Y2$INDV)]
names(sample_metrics)[names(sample_metrics) == 'missingness_Ychrom_GT_DP2'] = paste0('missingness_Ychrom_GT_DP', DP_2)
## Generate a file
write.table(sample_metrics,
           (paste0(PREFIX,".autosomal.xy.sample_based_all_statistics.preQC.GT.X.Y.txt")),
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE,
            sep = '\t')

# Save Chromosome Y - genotype QC stats
stats_missingness_Y_GT = bind_rows(imiss_GT_Y1_F_MISS,
                                   imiss_GT_Y2_F_MISS,
                                   lmiss_GT_Y1_F_MISS,
                                   lmiss_GT_Y2_F_MISS)
write.table(stats_missingness_Y_GT,
            "stats_missingness_Y_GT.txt", 
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
print(stats_missingness_Y_GT)
```

## 4: Individual Quality Control

In the previous steps we already have identified individuals that need to be removed\
In Autosome QC we identified those with missingness \>10%, mean depth \<20 and heterozygosity outliers `flagged_samples_autosomes`\
In Sex-chromosome qc we identified those where the disclosed sex didnt match the genotypic sex `flagged_samples_sex`\
After genotype QC we identified the samples that now had missingness \>10% `flagged_samples_autosomes_GT`

We can use R to merge all of these samples as a single file: `flagged_samples.raw-data.GT.txt`

```{r individual-qc-flagged-samples}
flagged_samples = bind_rows(flagged_samples_autosomes,
                            flagged_samples_sex,
                            flagged_samples_autosomes_GT)
# in case some IDs are repeated, retain only unique values
flagged_samples_unique = unique(flagged_samples) 
# Save is as a .txt file
write.table(flagged_samples_unique,
            "flagged_samples.raw-data.GT.txt", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)
```

### 4.a Check Relatedness and duplicates

#### i. Import file into plink

To check for relatedness and duplicate samples we will import the vcf file with genotype qc into plink. We will exclude the samples on the `flagged_samples.raw-data.GT.txt` file.

*Note*: the IDs of the samples to exclude need to be formated according to plink's `fam_id ind_id` requirement, which in this case is the vcf ID twice.

```{bash import-autosomes-plink, eval=FALSE, cache=FALSE, include=FALSE}
PREFIX='redlat_exomes'
DP=10
GQ=20

## I. Format samples to exclude according to plink's fam_id, ind_id requirement
awk '{print $1, $1}' flagged_samples.raw-data.GT.txt > flagged_samples.raw-data.GT.plink.txt

## II. Import data into plink
~/bin/plink --vcf $PREFIX.autosomes.DP$DP.GQ$GQ.vcf.gz --keep-allele-order --vcf-half-call m --double-id --remove flagged_samples.raw-data.GT.plink.txt  --geno 0.1 --make-bed --out $PREFIX.autosomes.GT
# vcf half calls will be set as missing
# the new id will be the vcf sample id for Family_id and Individual_id
```

#### ii. Update the sample ID, family ID and ID of the parents.

The `--double-id` flag takes the vcf sample id and uses it as Family_id and Individual_id in the .fam file

`–update-ids` expects input with the following four columns per row:\
Old family ID (which in this case is the sample id of the vcf)\
Old within-family ID(which in this case is the sample id of the vcf)\
New family ID\
New within-family ID

Its important that the individuals that are **known** to be related have the same family ID.

Ideally, all the samples in the vcf should be represented in the `new_ids.txt`. If this file has individuals that are not in the vcf, plink will give a warning, but the command will still run.

`–update-parents` expects the following four columns per row:\
Family ID\
Within-family ID\
New paternal within-family ID\
New maternal within-family ID

```{bash update-sample-id, eval=FALSE, cache=FALSE, include=FALSE}
PREFIX='redlat_exomes'
ids='exome_new_ids.txt'
parents='exome_new_parents.txt'


# Update individual id and family id
~/bin/plink --bfile $PREFIX.autosomes.GT --update-ids $ids --make-bed --out $PREFIX.autosomes.GT.id

# Update code for parents. 
~/bin/plink --bfile $PREFIX.autosomes.GT.id --update-parents $parents --make-bed --out $PREFIX.autosomes.GT.id.parents
```

#### iii. Detect duplicate samples

We will use [king](https://www.kingrelatedness.com/) for the following steps. King can take both vcf and plink files and input, but using plink files makes it substantially faster.

```{bash check-duplicates, eval=FALSE, cache=FALSE, include=FALSE}
PREFIX='redlat_exomes'

# Unlike in plink, you have to put the .bed extension in the king prompt
~/bin/king -b $PREFIX.autosomes.GT.id.parents.bed --duplicate --rplot
```

This will create the file `king.con` which has the following columns:\
*FID1,* Family 1 ID\
*ID1,* Individual 1 ID\
*FID1,* Family 2 ID\
*ID2,* Individual 2 ID\
*N*,Number of variants\
*N_IBS0,* Number of variants not Identical by descent (IBD) *N_IBS1* Number of variants with 1 IBD allele\
*N_IBS2,* Number of variants with 2 IBD alleles\
*Concord,* Concordance of genotypes\
*HomConc,* Concordance of homozygous genotypes\
*HetConc,* Concordance of heterozygous genotypes\

It will also generate an R script `king_duplicateplot.R` that depicts the duplicate clusters (in case than several samples are identical)

```{r plot-duplicates}
source("king_duplicateplot.R", local = knitr::knit_global())
# This script will save a pdf file with the image bellow. 
# If you also want to see the image in the markdown you can run the following command
plot(g, vertex.label.dist=0.3,vertex.label.degree=pi/2,vertex.size=1,vertex.label.cex=0.5, edge.arrow.size=0.5, edge.arrow.width=1,
  layout=layout_with_fr, asp=0)
```

For the ReDLat dataset we have 1 pair of samples. Next step is to determine if this duplicate is the same sample sequenced twice (which is often done if the first sequencing was suspected to have low quality), or if it is supposed to represent two different samples.  If **both samples are known to correspond to the same individual** look into the other quality metrics (missingness or depth) and pick the one with the best quality data.\
If **both samples are supposed to correspond to the different individuals** you will have to discard both genomes/exomes as you cannot be sure which is the individual that this genetic data really represents.\
(Note: there are some ways in which this puzzle be solved, like comparing to known relatives of one of the possible individuals, or other methods of 'DNA profiling' approaches which wont be covered here)

For the ReDLat dataset we have reviewed the `king.con` file and decided to remove all the samples flagged as duplicates.

```{r duplicates-dataframe}
duplicate_df = read.delim("king.con", header = T)
ID_1 = select(duplicate_df, FID1, ID1)
ID_1_rename = rename(ID_1, FID = FID1, ID = ID1)
ID_2 = select(duplicate_df, FID2, ID2)
ID_2_rename = rename(ID_2, FID = FID2, ID = ID2)
duplicate_df2 = bind_rows(ID_1_rename , ID_2_rename )
flagged_samples_duplicates_plink = unique(duplicate_df2) 
write.table(flagged_samples_duplicates_plink,
            "flagged_samples_duplicates_plink.txt", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)
```

Remove duplicates from the plink file

```{bash remove-dupicates, cache = TRUE}
PREFIX='redlat_exomes'

~/bin/plink --bfile $PREFIX.autosomes.GT.id.parents --remove flagged_samples_duplicates_plink.txt --make-bed --out $PREFIX.autosomes.GT.id.parents.no-dups
```

#### iv. Check Relatedness

Given that this are exomes we will only check fisrt and second degree

```{bash king-relatedness, cache = TRUE}
PREFIX='redlat_exomes'

~/bin/king -b $PREFIX.autosomes.GT.id.parents.no-dups.bed --related --degree 2
```

King will generate two files:\
`king.kin:` Within-family kinship data, which is relatedness between individuals that share the same Family ID `king.kin0:` Between-family relatives, which is relatedness between individuals that *don't* share the same Family ID

These files should be checked manually to detect errors in the known relationships and samples with cryptic relatedness.

**NOTE:** Exomes of founder populations may have significant similarity as a result of background relatedness.\
I suggest to verify *InfType* up to third degree in `king.kin0` for cryptic relatedness.

For the ReDLat cohort two samples with cryptic relatedness were identified and manually saved in the file `flagged_samples_cryptic-rel_plink.txt`

### 4.b Organize all the flagged samples in the same file

We have already generated these r objects where each line corresponds to a vcf sample id\
*flagged_samples_autosomes*\
*flagged_samples_sex*\
*flagged_samples_autosomes_GT*\

We also generated an r object from the king duplicates. IDs in this file are in plink format\
*flagged_samples_duplicates_plink*

We will import the samples flagged for cryptic relatedness into R and use the `new_ids.txt` file to get the vcf id equivalency. Then we will generate a single file with all the IDs that need to be removed from the vcf.

```{r individual-qc-flagged-samples-all}

# King output for cryptic relationships
if (file.exists("flagged_samples_cryptic-rel_plink.txt")) {
  flagged_samples_cryptic_rel_plink = read.delim("flagged_samples_cryptic-rel_plink.txt", header = FALSE)
  colnames(flagged_samples_cryptic_rel_plink) = c('FID', 'ID') 
  # Merge flagged samples for cryptic relatedness and duplicates                                         
  flagged_samples_king_plinkformat = bind_rows(flagged_samples_cryptic_rel_plink,
                                              flagged_samples_duplicates_plink)
} else {
  flagged_samples_king_plinkformat = flagged_samples_duplicates_plink
}


# Read plink vcf id equivalencies
plink_ids = read.delim("exome_new_ids.txt", header = FALSE)
colnames(plink_ids) = c('vcf_id', 'vcf_id2', 'plink_fam', 'plink_iid') 

# Add the vcf id
flagged_samples_king_plinkformat$INDV = plink_ids$vcf_id[match(flagged_samples_king_plinkformat$ID, plink_ids$plink_iid)]  

# Retain the vcf ids of the samples flagged during the King analyses
flagged_samples_king = select(flagged_samples_king_plinkformat, INDV)

flagged_samples_all = bind_rows(flagged_samples_autosomes,
                                flagged_samples_sex,
                                flagged_samples_autosomes_GT,
                                flagged_samples_king )

# Save it
write.table(flagged_samples_all,
            "flagged_samples_all.txt", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)
```

### 4.c Remove flagged samples from the genotype QC'd vcf

```{bash autosomes-individual-QC, cache=FALSE}
PREFIX='redlat_exomes'
DP=10
GQ=20

vcftools --gzvcf $PREFIX.autosomes.DP$DP.GQ$GQ.vcf.gz --remove flagged_samples_all.txt --recode --recode-INFO-all --out $PREFIX.autosomes.DP$DP.GQ$GQ.ind

mv $PREFIX.autosomes.DP$DP.GQ$GQ.ind.recode.vcf $PREFIX.autosomes.DP$DP.GQ$GQ.ind.vcf
bgzip $PREFIX.autosomes.DP$DP.GQ$GQ.ind.vcf

```

## 5: Site Based Quality Control

Remove sites with Missingness \>= 10% and those that are below the desired SNP VQSLOD tranche

**NOTE**: 99.9% is the recommended default VQSLOD cutoff for SNPs in human genomic analysis

```{bash site-qc, cache = TRUE}
PREFIX='redlat_exomes'
DP=20
GQ=20

vcftools --gzvcf $PREFIX.autosomes.DP$DP.GQ$GQ.ind.vcf.gz --max-missing 0.9 --remove-filtered VQSRTrancheSNP99.00to99.90 --recode --recode-INFO-all --out $PREFIX.autosomes.DP$DP.GQ$GQ.ind.site
# --max-missing float is a 0-1range. 
# 0 allows sites that are completely missing and 1 indicates no missing data allowed

mv $PREFIX.autosomes.DP$DP.GQ$GQ.ind.site.recode.vcf $PREFIX.autosomes.DP$DP.GQ$GQ.ind.site.vcf
bgzip $PREFIX.autosomes.DP$DP.GQ$GQ.ind.site.vcf
```

We have already performed Genotype filtering so we will only retain sites where our sample has \>95% of individuals with sequencign information for that given variant. We will remove variants with missingness higher than 5%.

Notice how "--max-missing" in vcf tools excludes sites on the basis of the proportion of missing data (defined to be between 0 and 1, where 0 allows sites that are completely missing and 1 indicates no missing data allowed). This is the oposite of how plink "--geno" does it!!

```{bash autosomes-site-QC, cache=FALSE}
PREFIX='redlat_exomes'
DP=10
GQ=20
vcftools --gzvcf $PREFIX.autosomes.DP$DP.GQ$GQ.vcf.gz --max-missing 0.95 --recode --recode-INFO-all --out $PREFIX.autosomes.DP$DP.GQ$GQ.site

```
