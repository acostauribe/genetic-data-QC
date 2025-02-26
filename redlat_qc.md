---
output:
  html_document:
    keep_md: yes
---


# Quality control pipeline for the ReD-Lat genomic data

#### Developed by Juliana Acosta-Uribe for the ReD-Lat Consortium 2023

This pipeline is designed to be run as an [R markdown](https://rmarkdown.rstudio.com/lesson-1.html) file in R Studio. This way you can run it in a step-by-step mode. However, you could also run it directly from the r command line if you already have the `sample_data.txt` and the `problematic_relatedness.txt` files in your workspace.

```         
library(rmarkdown)
render("path/to/your/file.Rmd")
```

For this pipeline you will need the following:

### Data:

A plink formatted .*bed*, .*fam*, .*bim* set, or a bgzipped *.vcf.gz* file.

**Always Take a look at your files before beginning:** 

If you start with a plink dataset: \ 
- Do the Individual IDs (column 2 in *.fam*) match what you were expecting? \ 
- Have the families been given a Family ID (column 1 in *.fam*)? \ 
- Do the Individuals have their sex assigned (column 5 in *.fam*? \ 
- Do the variants in the *.bim* have an identifier (column 2 in *.bim*)? Some analyses will require this information, and we may have to incorporate it to the *.fam*/*.bim* if its not already there. \
Make sure your files are properly aligned and the alleles are being called from the correct strand. INDELs should be [left aligned and normalized](https://samtools.github.io/bcftools/bcftools.html#norm).

**Sample Data File** 

If you are starting with a *.vcf*, or your *.fam* does not have the sex /family of the samples already specified please provide an additional file with a header as follows:

**IID** ID of the sample as it is in the plink IID or in the VCF\
**SEX** should be specified (1=male, 2=female, 0=no data)\
**FID** Family ID of the sample\
**PID** Paternal ID of the sample. If this individual is also in the dataset, the PID should be identical to the father's IID\
**MID** Maternal ID of the sample. If this individual is also in the dataset, the MID should be identical to the mother's IID\
**PHENO** Optional, if you are interested in any downstream analyses involving the phenotype. (1=control, 2=case, -9=missing)

⚠️ The **FID** given to the parents needs to match the same **FID** given to that individual for their genome.

A trio would look like this (order of the columns does not matter)

| FID  | IID | PID | MID | SEX | PHENO |
|------|-----|-----|-----|-----|-------|
| FID1 | SON | DAD | MOM | 1   | 2     |
| FID1 | MOM | 0   | 0   | 2   | 1     |
| FID1 | DAD | 0   | 0   | 1   | 1     |

R is expecting a tab delimited file. If your file is delimited by spaces you can fix it with the following bash command `sed -i 's/ /\t/g'  sample_data.txt`


### Software:

-[R](https://www.r-project.org/)\
-[RStudio](https://posit.co/download/rstudio-desktop/)\
-[plink](https://www.cog-genomics.org/plink2/)\
-[king](https://www.kingrelatedness.com/)

If you are working with Exome or Genome data, you will also need:

-[vcftools](https://vcftools.github.io/man_latest.html)\
-[bcftools](https://samtools.github.io/bcftools/bcftools.html)


### Contents Table:

[1. Set up your environment](#section-1)\
[2. Customize the quality control process](#section-2)\
[3. Start Quality Control process](#section-3)\
[4. Genotype Quality control](#section-4)\
[5. Filter for missingness](#section-5)\
[6. Calculate heterozygosity](#section-6)\
[7. Check sex](#section-7)\
[8. Identify duplicates](#section-8)\
[9. Calculate relatedness](#section-9)\
[10. Remove variants & Individuals and with missingness ≥5%](#section-10)


## 1. Set up your environment {#section-1}

In this first step you will install the packages that R needs to run the quality control process and you will define the paths to your working directory and software.

```{r environment-setup}
# Install required R packages:
if (!require("knitr", quietly = TRUE))
install.packages("knitr")
library(knitr)

# Set your working directory:
knitr::opts_chunk$set(root.dir = "~/Genetic-Sequencing_raw-data/HudsonAlpha_J.Nicholas.Cochran/redlat/Exome_Seq/Joint_call_fixed",
                      dev = "png",
                      dpi = 300,
                      echo = FALSE,
                      cache = TRUE)
setwd("~/Genetic-Sequencing_raw-data/HudsonAlpha_J.Nicholas.Cochran/redlat/Exome_Seq/Joint_call_fixed")


# Set up path to software:
Sys.setenv(plink='/home/acostauribe/bin/plink')
Sys.setenv(king='/home/acostauribe/bin/king')
Sys.setenv(vcftools='/usr/bin/vcftools')
Sys.setenv(bcftools='/usr/bin/bcftools')

# Give the name of your starting file without the .bed or .vcf.gz extension
prefix="redlat_exome"
Sys.setenv(prefix=prefix)

# Define the type of data you will be working with. 
# Answers can be "GENOME", "EXOME", "ARRAY"
data_type="EXOME"
Sys.setenv(data_type=data_type)

# Specify your reference genome
# Answers can be 'hg19', 'hg38'
genome_alignment="hg38"
Sys.setenv(genome_alignment=genome_alignment)

# Import the text file that contains the information for your samples. 
# Look at Sample Data File specifications in the introduction. File needs to be tab delimited, with a header
sample_data_file="exome_sample_data.txt"

if (file.exists(sample_data_file)) {
    sample_data = read.delim(sample_data_file, header = TRUE, sep = '\t')
    print(paste(".fam file will be annotated with:", sample_data_file))
  } else { 
    print("No sample data was provided, assuming .fam file is already annotated") 
  }
Sys.setenv(sample_data_file=sample_data_file)
```

## 2. Customize the quality control process {#section-2}

Choose what you want to do in the analysis:

```{r customize-pipeline}
# Do you want to do filter variant calls for genotype depth (DP) and Genotype Quality (GQ)? 
## This can only be done if you start with a VCF file and requires that the "DP and 'GQ' FORMAT tag is included for each genotype.
## TRUE Includes only sites with DP/GQ (over all included individuals) greater than or equal to the given value 
genotype_filter=TRUE
# Define the values you want to use as threshold
DP=10
GQ=20

# Do you want to keep only the variants with PASS in the VQSR filtering?
## TRUE Removes all sites with a FILTER flag other than PASS.
vqsr_PASS=FALSE

# Do you want to check for known and cryptic relatedness among samples?
## For this analyses you will need to provide information of known families as described before.
## TRUE uses KING to perform relatedness check.
check_relatedness=FALSE

# Do you want to create directories and to organize your data as you go?
## TRUE generates directories and organizes your files as you run the pipeline
tidy_up=TRUE

# Make these values into system environment variables
Sys.setenv(genotype_filter=genotype_filter)
Sys.setenv(DP=DP)
Sys.setenv(GQ=GQ)
Sys.setenv(vqsr_PASS=vqsr_PASS)
Sys.setenv(check_relatedness=check_relatedness)
Sys.setenv(tidy_up=tidy_up)
```

## 3. Start Quality Control process {#section-3}

Before starting, take some time to look at your data and calculate some basic quality metrics. If you are starting with a vcf file from Exome or Genome data, after running the following chunk you should have at least the 6 following files:

**I. General statistics** [.vchk](https://samtools.github.io/bcftools/bcftools.html#stats)\
**II. Sample based metrics** [metrics](https://vcftools.sourceforge.net/man_latest.html#OUTPUT%20OPTIONS)\
1. missingness per sample (.imiss)\
2. depth per sample (.idepth)\
**III. Site based metrics**\
1. Missingness per site (.lmiss)\
2. Mean depth per site (.ldepth.mean)\
3. VQSR quality [vqsr](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-)


```{bash get-preQC-stats-2}
if [[ "$data_type" == "EXOME" || "$data_type" == "GENOME" ]]; then
    # Data description
    bcftools stats ${prefix}.vcf.gz > ${prefix}.preqc.vchk.txt
    
    # Individual missingness
    vcftools --gzvcf ${prefix}.vcf.gz --missing-indv --out ${prefix}.preqc
    # Generates a file reporting the missingness on a per-individual basis.
    # The output file has the suffix ".imiss".
    # Individuals whose missingness is >10% should be added to 'flagged_samples' file
    
    # Individual depth
    vcftools --gzvcf ${prefix}.vcf.gz --depth --out ${prefix}.preqc
    # Generates a file containing the mean depth per individual.
    # This file has the suffix ".idepth".
    #  Individuals whose mean depth is <20 should be added to 'flagged_samples' file
    
    # Site missingness
    vcftools --gzvcf ${prefix}.vcf.gz --missing-site --out ${prefix}.preqc
    # Generates a file reporting the missingness on a per-site basis.
    # The file has the suffix ".lmiss".
    
    # Remove the 'chr' prefix if chromosomes are labeled 'chr1'
    sed -i 's/chr//g' ${prefix}.preqc.lmiss
    
    # Site depth
    vcftools --gzvcf ${prefix}.vcf.gz --site-mean-depth --out ${prefix}.preqc
    # Generates a file containing the mean depth per site across all individuals.
    # This output file has the suffix ".ldepth.mean"
    
    # Remove the 'chr' prefix if chromosomes are labeled 'chr1'
    sed -i 's/chr//g' ${prefix}.preqc.ldepth.mean
    
    # VQRS filtering
    vcftools --gzvcf ${prefix}.vcf.gz --FILTER-summary --out ${prefix}.preqc
    # Gives the information on the number of variants that passed VQSR filtering, or are in specific tranches
    # The output file has the suffix ".FILTER.summary"
fi
```

If your starting dataset is a plink file from a snp array, you will first have to edit the .fam with the samples sex and the only statistics you get are missingness per sample and per site. \
If you don't add sex to the fam before calculating missingness, plink wont calculate missingness values for Y chromosome. 

```{r edit-fam-1}
if(data_type=="ARRAY"){
  fam = read.delim(paste0(prefix,".fam"),
                 sep=' ',
                 header=FALSE)

  # If your .fam doesn't have the sex already specified, it will add the SEX from sample_data
  if(sum(fam$V5) > length(fam$V5)*0.25) {
    fam$V5 = sample_data$SEX[match(fam$V2, sample_data$IID, nomatch = 0)]
    write.table(fam,
                paste0(prefix,".fam"), 
                col.names = FALSE, 
                row.names = FALSE, 
                quote = FALSE)
  }
}
```

Then, use plink for Site and individual missingness in SNP Array. \
It will produce .lmiss and .imiss files      
```{bash get-preQC-stats-1}
if [[ "$data_type" == "ARRAY" ]]; then
  $plink --bfile ${prefix} --missing --out ${prefix}.preqc
  rm *.hh
fi
```


You can then proceed to get some plots and statistics from your initial files.

**Sample based metrics**

```{r plot-preQC-stats-ind}
if (!require("psych", quietly = TRUE))
install.packages("psych")
library(psych)
if (!require("dplyr", quietly = TRUE))
install.packages("dplyr")
library(dplyr)

# Individual Missingness 
imiss = read.delim(file = paste0(prefix,".preqc.imiss"), 
                    header = TRUE, 
                    sep = '')
# Generate basic summary statistics
imiss_F = describe(imiss$F_MISS)
rownames(imiss_F) = c("sample_missingness")
print(imiss_F)
# Plot Missingness rate per sample
imiss_hist = hist(imiss$F_MISS,
             xlab="Missingness rate",
             ylab="Samples", 
             main="Missingness rate per sample preQC",
             col="cyan3",
             breaks=20)

# Depth per Individual (Only for Exome or Genome data)
if(data_type=="EXOME" || data_type=="GENOME"){
  # Read depth calculations
  idepth = read.delim((paste0(prefix,".preqc.idepth")), header = T, sep = "")
  
  # Generate basic summary statistics
  idepth_mean = describe(idepth$MEAN_DEPTH)
  rownames(idepth_mean) = c("sample_depth")
  print(idepth_mean)
  
  # Individuals whose mean depth is <20 should be identified
  low_mean_depth = filter(idepth, MEAN_DEPTH<20)
  if (nrow(low_mean_depth) > 0) {
      write.table(low_mean_depth, 
              "samples_low_mean_depth_raw.txt",
              col.names = FALSE, 
              row.names = FALSE, 
              quote = FALSE,
              sep = 't')
    } else {
      print("All samples have a mean depth higher than 20. 'samples_low_mean_depth_raw.txt' wont be generated")
  }
  
  # Plot depth per sample
  idepth_hist = hist(idepth$MEAN_DEPTH,
                     xlab="Mean Depth ",
                     ylab="Samples", 
                     main="Mean Depth per sample preQC", 
                     col="cyan3",
                     breaks=50)
}
```
![Mean Depths Per Sample](https://github.com/acostauribe/genetic-data-QC/blob/apr17/redlat_result/3_MeanDepthPerSample.png?raw=true)
![Missingness Per Sample](https://github.com/acostauribe/genetic-data-QC/blob/apr17/redlat_result/3_MissPerSample.png?raw=true)

Summarize your sample statistics

```{r generate-sample-statistics}
# Take the information of individuals from the .idepth file. 
# We are changing the names of the columns to specify these values come from raw data

if (data_type == "ARRAY"){
  sample_statistics = subset(imiss, select = c(IID, F_MISS, N_GENO))
  sample_statistics = rename(sample_statistics, Initial_missingness = F_MISS, Initial_n_sites = N_GENO)
} else if (exists("idepth")){ 
  sample_statistics = rename(idepth, IID = INDV, Initial_n_sites = N_SITES, Initial_mean_depth = MEAN_DEPTH)
  # Using the "match" function, we will create a new column in 'sample_statistics' with the missingness per sample.
  # imiss and idepth should have the same samples in the same order, but using the "match" function will be useful when   we start dropping samples
  sample_statistics$Initial_missingness = imiss$F_MISS[match(sample_statistics$IID, imiss$INDV)]
} else {
  sample_statistics = subset(imiss, select = c(INDV, F_MISS, N_DATA))
  sample_statistics = rename(imiss, IID = INDV, Initial_missingness = F_MISS, Initial_n_sites = N_DATA)
}
# Save as a file
write.table(sample_statistics,
            "sample_statistics.txt", 
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE)
```

**Site based metrics**

```{r plot-preQC-stats-site}
if (!require("ggplot2", quietly = TRUE))
install.packages("ggplot2")
library(ggplot2)

# Site Missingness
lmiss = read.delim(file = paste0(prefix,".preqc.lmiss"), 
                    header = TRUE, 
                    sep = '')

# Generate basic summary statistics
lmiss_F = describe(lmiss$F_MISS)
rownames(lmiss_F) = c("site_missingness")
print(lmiss_F)

# Plot Missingness rate site per chromosome

sorted_chr = unique(lmiss$CHR) #organize the chromosomes in ascending order
lmiss$CHR = factor(lmiss$CHR , levels=sorted_chr)
  
lmiss_box = boxplot(F_MISS ~ CHR,
                    data=lmiss,
                    xlab="Chromosome",
                    ylab="Missingness rate",
                    cex.axis = 0.5 ,
                    main="Missingness rate per site preQC",
                    col="cyan3")

# Depth per site - chromosome (Only for Exome or Genome data)
if(data_type=="EXOME" || data_type=="GENOME"){
  # Read depth calculations
  ldepth = read.delim((paste0(prefix,".preqc.ldepth.mean")), header = T, sep = "")

  # Generate basic summary statistics
  ldepth_mean = describe(ldepth$MEAN_DEPTH)
  rownames(ldepth_mean) = c("site_depth")
  print(ldepth_mean)
  ### Plot depth per site
  
  sorted_chr = unique(ldepth$CHR) #organize the chromosomes in ascending order
  ldepth$CHR = factor(ldepth$CHR , levels=sorted_chr)
  
  ldepth_box = boxplot(MEAN_DEPTH ~ CHR,
                       data=ldepth,
                       xlab="Chromosome",
                       ylab="Mean Depth per site",
                       cex.axis = 0.5 ,
                       main="Mean Depth per site preQC",
                       col="cyan3")
}

# Variant Quality Score Recalibration ranks (Only for Exome or Genome data)
if(data_type=="EXOME" || data_type=="GENOME"){
  vqsr = read.delim((paste0(prefix,".preqc.FILTER.summary")), header = T, sep = "")
  print(vqsr)
  #Plot it
  vqsr %>%
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
      labs(title = "Variant Quality Score Recalibration ranks preQC")
  #ggsave("VQSR-preQC.png")  
}
```
![Mean Depths Per Site](https://github.com/acostauribe/genetic-data-QC/blob/apr17/redlat_result/3_MeanDepthPerSite.png?raw=true)
![Missingness Per Site](https://github.com/acostauribe/genetic-data-QC/blob/apr17/redlat_result/3_MissPerSite.png?raw=true)
![VQSR](https://github.com/acostauribe/genetic-data-QC/blob/apr17/redlat_result/3_VQSR.png?raw=true)

Summarize statistical measures

```{r generate-statistical-measures}
if(data_type=="EXOME" || data_type=="GENOME"){
dataset_statistics = bind_rows(imiss_F,
                                 idepth_mean,
                                 lmiss_F,
                                 ldepth_mean)
} else {
  dataset_statistics = bind_rows(imiss_F,
                                  lmiss_F)
}

write.table(dataset_statistics,
            "dataset_statistics.txt", 
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE)

print(dataset_statistics)
```

Organize all the generated files

```{bash organize-files-1}
if [ "$tidy_up" == "TRUE" ]; then
  mkdir Pre-QC_stats
  mv ${prefix}.preqc.* ./Pre-QC_stats
fi
```

## 4. Genotype Quality control {#section-4}

This step is only possible to do when data is in a VCF file that has been properly annotated with genotype depth (DP) and genotype quality (GQ). Notice that no sites will be eliminated, as this only affects individual calls. For example if an individuals genotype for a given variant is below desired thresholds, it's turned into a missing value in that individual.

```{bash genotype-qc}
if [[ "$data_type" == "EXOME"  || "$data_type" == "GENOME" ]] && [[ "$genotype_filter" == "TRUE" ]]; then
  if [[ "$vqsr_PASS" == "TRUE" ]]; then
    $vcftools --gzvcf ${prefix}.vcf.gz --minDP $DP --minGQ $GQ --remove-filtered-all --recode --recode-INFO-all --out ${prefix}.DP$DP.GQ$GQ
  else
    $vcftools --gzvcf ${prefix}.vcf.gz --minDP $DP --minGQ $GQ --recode --recode-INFO-all --out ${prefix}.DP$DP.GQ$GQ
  fi
  mv ${prefix}.DP$DP.GQ$GQ.recode.vcf  ${prefix}.DP$DP.GQ$GQ.vcf
  bgzip ${prefix}.DP$DP.GQ$GQ.vcf
elif [[ "$genotype_filter" == 'FALSE' ]]; then
  echo "genotype_filter was set to FALSE, this step is skipped"
else
  echo "Plink data cannot be checked for depth or genotype quality"
fi
```

**Note**: VCFtools output can be compressed directly into a .gz file by adding `–stdout | gzip -c > newname.gz`, but it won't generate a .log file, which are useful for debugging.

Organize files

```{bash organize-files-2}
if [ "$tidy_up" == "TRUE" ]; then
  mkdir File_evolution
  if [[ "$data_type" == "EXOME" || "$data_type" == "GENOME" ]] && [[ "$genotype_filter" == "TRUE" ]]; then
    mv ${prefix}.vcf.gz ./File_evolution
  fi
fi
```

Update prefix variable

```{r update-prefix-1}
if((data_type=="EXOME" || data_type=="GENOME") && (genotype_filter=="TRUE")){
  prefix=(paste0(prefix,".DP",DP,".GQ",GQ))
  Sys.setenv(prefix=prefix)
}
```

## 5. Filter for missingness {#section-5}

Identify samples and variants missing more than 10% of data and remove them. The 10% or 0.1 threshold can be modified as desired.

```{bash filter-missingness}
if [[ "$data_type" == "EXOME" || "$data_type" == "GENOME" ]]; then
  # Identify individuals with high missingness
  $vcftools --gzvcf ${prefix}.vcf.gz --missing-indv --out ${prefix}.mind
  awk '{if ($5 > 0.1) {print $1} }' ${prefix}.mind.imiss > ${prefix}.mind.irem
  awk '{if ($5 <= 0.1) {print $1} }' ${prefix}.mind.imiss > ${prefix}.mind.ikeep
  
  # Remove variants and sites missing more than 10% of data
  $vcftools --gzvcf ${prefix}.vcf.gz --keep ${prefix}.mind.ikeep --max-missing 0.9 --recode --recode-INFO-all --out ${prefix}.miss
  mv ${prefix}.miss.recode.vcf ${prefix}.miss.vcf
  bgzip ${prefix}.miss.vcf
  
else 
  # Remove samples missing more than 10% of data
  $plink --bfile ${prefix} --mind 0.1 --make-bed --out ${prefix}.miss
  rm *.hh
fi
```

> Individuals missing more than 10% of the data get removed in this step. Removed individuals will be written to *prefix.mind.irem*

Update your sample_statistics file

```{r update-sample-metrics-1}
if(data_type=="EXOME" || data_type=="GENOME"){
  imiss_qc = read.delim(paste0(prefix,".mind.imiss"))
  sample_statistics$GQDP_missingness = imiss_qc$F_MISS[match(sample_statistics$IID, imiss_qc$INDV)]
}
```

Tidy up

```{bash organize-files-3}
if [ "$tidy_up" == "TRUE" ]; then
  if [[ "$data_type" == "GENOME" || "$data_type" == "EXOME" ]]; then
    if [[ "$genotype_filter" == "TRUE" ]]; then
      mkdir genotype_qc
      mv ${prefix}.mind.* ./genotype_qc
    else
      mv ${prefix}.mind.* ./Pre-QC_stats
    fi
    mv ${prefix}.vcf.gz ./File_evolution
  else
    mv ${prefix}.bed ./File_evolution
    mv ${prefix}.bim ./File_evolution
    mv ${prefix}.fam ./File_evolution
    mv ${prefix}.log ./File_evolution
    mv ${prefix}.miss.irem ./File_evolution
  fi
fi
```

Update prefix variable

```{r update-prefix-2}
prefix=(paste0(prefix,".miss"))
Sys.setenv(prefix=prefix)
```

## 6. Calculate heterozygosity {#section-6}

This is a recommended step when you are analyzing cohorts with either SNP array or Whole genome sequencing.

We prefer plink to calculate sample heterozygosity. Plink uses a *sliding window* approach to identify variants in linkage disequilibrium. There are many options to modify the behavior or this approach in [plink's documentation](https://www.cog-genomics.org/plink/1.9/ld#indep). The LD pruning requires that the *.bim* file has variant IDs in the second column. If no variants have been assigned, you could do a preliminary step using [--set-missing-var-ids](https://www.cog-genomics.org/plink/1.9/data#set_missing_var_ids).

```{bash import-vcf, echo=FALSE}
if [ "$data_type" == "GENOME" ] || [ "$data_type" == "EXOME" ]; then
# Import VCF and assign IDs to SNV
  $plink --vcf ${prefix}.vcf.gz --keep-allele-order --double-id --vcf-half-call m --set-missing-var-ids @:#\$1,\$2 --new-id-max-allele-len 1 --make-bed --out ${prefix}
# Assign IDs to indels
mv ${prefix}.bim ${prefix}.bim-original
awk 'BEGIN {count=1}; {if ($2 ~ /\./) {sub(/\./,"INDEL"(count++));print} else {print} }' ${prefix}.bim-original > ${prefix}.bim
fi
```
```{bash calculate-heterozygosity, echo=FALSE}
# I have this as an if statement because im still deciding if this is a good measure for exomes. 
if [ "$data_type" == "GENOME" ] || [ "$data_type" == "ARRAY" ] || [ "$data_type" == "EXOME" ]; then
  # Retain variants with MAF > 10% and individuals with low missing %
  $plink --bfile ${prefix} --maf 0.1 --make-bed --out ${prefix}.maf
  
  # Calculate LD
  #--indep-pairwise <window size>['kb'] <step size (variant ct)> <r^2 threshold>
  $plink --bfile ${prefix}.maf --indep-pairwise 50 10 0.2 
  
  # Retain variants not in LD (independent markers)
  $plink --bfile ${prefix}.maf --extract plink.prune.in --make-bed --out ${prefix}.maf.ld
  
  # Check heterozygosity
  $plink --bfile ${prefix}.maf.ld --het --out ${prefix}.het
  
  # Calculate missingness rates
  $plink --bfile ${prefix} --missing --out ${prefix}
  
  rm ${prefix}.maf.*
  rm plink.*
fi
```

**Identify heterozygosity outliers**

```{r plot-heterozygosity}
# Load data into R
het = read.delim(file = paste0(prefix,".het.het"), 
                    header = TRUE, 
                    sep = '')

# Generate basic summary statistics
heterozygosity_sample = describe(het$F)
rownames(heterozygosity_sample) = c("Sample_heterozygosity")
print(heterozygosity_sample)

# Calculate limits for excluding samples [3 standard deviations]
## Low threshold
heterozygosity_low_limit = mean(het$F)-(3*(sd(het$F)))
print(paste("Mean heterozygosity F value -3 Standard Deviations =", heterozygosity_low_limit))

## High threshold
heterozygosity_high_limit = mean(het$F)+(3*(sd(het$F)))
print(paste("Mean heterozygosity F value +3 Standard Deviations =", heterozygosity_high_limit))

# Individuals whose heterozygosity deviated more than 3 SD from the mean should be identified
het_outlier_low = filter(het, F<heterozygosity_low_limit)
het_outlier_high = filter(het, F>heterozygosity_high_limit)
het_outlier_both = bind_rows(het_outlier_low,
                             het_outlier_high)
het_outlier_id = select(het_outlier_both, FID, IID)

# Save your list. Useful for debugging
write.table(het_outlier_id,
            "heterozygosity_outliers.txt", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)
print(paste(nrow(het_outlier_id), "samples are heterozygosity outliers. IDs are written in heterozygosity_outliers.txt"))

# Plot heterozygosity per sample
heterozygosity_hist = hist(het$F,  
                           freq=TRUE, 
                           xlab="Heterozygosity F coefficient",  
                           ylab="Samples", 
                           main="Heterozygosity rate per sample",
                           col="cyan3",
                           breaks=100)
                    abline(v = (heterozygosity_low_limit), col="red")
                    abline(v = (heterozygosity_high_limit), col="red")
                    abline(v = (mean(het$F)), col="blue") 
                    legend("topleft",
                           c("+/-3 SD","mean"),
                           col=c("red","blue"),
                           pch=16)
```
![Heterozygosity Rate (Histogram)](https://github.com/acostauribe/genetic-data-QC/blob/apr17/redlat_result/6_Het_Hist.png?raw=true)

> Heterozygosity outliers will be removed along with those that fail sex-check on next step.

It is useful to visualize heterozygosity F statistic vs. missingness per sample

```{r heterozygosity-vs-missing}
# Load your data:
miss_imiss = read.delim(file = paste0(prefix,".imiss"), 
                    header = TRUE, 
                    sep = '')

plot(miss_imiss$F_MISS, het$F,
     xlab="Missingness",  
     ylab="Heterozygosity", 
     main="Heterozygosity rate per sample - Data before QC")
abline(h = (heterozygosity_low_limit), col="red")
abline(h = (heterozygosity_high_limit), col="red")
abline(h = (mean(het$F)), col="blue") 
legend("bottomright",
       legend = c("+/-3 SD", "mean"),
       col = c("red", "blue"),
       pch = 1)
```
![Heterozygosity Rate (Scatter)](https://github.com/acostauribe/genetic-data-QC/blob/apr17/redlat_result/6_Het_Scatter.png?raw=true)

Update your Sample Metrics dataframe

```{r update-sample-metrics-2}
sample_statistics$heterozygosity = het$F[match(sample_statistics$IID, het$IID)]
write.table(sample_statistics,
            "sample_statistics.txt", 
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE)

# Update your dataset_statistics dataframe
dataset_statistics = bind_rows(dataset_statistics,
                                  heterozygosity_sample)
# Save as a file
write.table(dataset_statistics,
            "dataset_statistics.txt", 
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE)

print(dataset_statistics)
```

## 7. Check sex {#section-7}

**Preprocess sex chromosomes**

A. Verify that *prefix.fam* has a 'sex' value on the 5 column.

If the values are 0 or -9, you can assign 'sex value' in plink with the --update-sex command. You can also use R to edit the .fam by running the following chunk

```{r edit-fam-2}
# Load your dataframes
fam = read.delim(paste0(prefix,".fam"),
                 sep=' ',
                 header=FALSE)

# If your .fam doesn't have the sex already specified, it will add the SEX from sample_data
if(sum(fam$V5) > length(fam$V5)*0.25) {
  fam$V5 = sample_data$SEX[match(fam$V2, sample_data$IID, nomatch = 0)]
  write.table(fam,
              paste0(prefix,".fam"), 
              col.names = FALSE, 
              row.names = FALSE, 
              quote = FALSE)
} else {print("fam already contains sex values") }
```

B. Split Pseudo-Autosomic region (PAR) of X

It is necessary to remove the PAR region of X for plink's sex check pipeline. *awk* will check if there there are variants in the PAR region of your X chomosome (assumes hg38).

If there are nor variants in the PAR region, then it will append '.split-x' to ensure compatibility further on. Make sure you are using the correct genome alignment for this. e.g. ReDLat files are aligned to hg38. Set up `genome_alignment` variable accordingly You can find more information on this step in [plink documentation](https://www.cog-genomics.org/plink/1.9/data#split_x)

```{bash split-PAR}
if [[ $genome_alignment != 'hg38' ]]; then
 echo 'change start and end bp positions of X PAR in the awk command below'
 elif [[ $(awk '$1 == "23" && ($4 <= 2781479 || $4 >= 155701383) { count++ } END { print count }' "${prefix}.bim") -gt 0 ]]; then
  $plink --bfile ${prefix} --split-x $genome_alignment --make-bed --out ${prefix}.split-x 
else
  mv ${prefix}.bed ${prefix}.split-x.bed
  mv ${prefix}.bim ${prefix}.split-x.bim
  mv ${prefix}.fam ${prefix}.split-x.fam
  echo "plink files have been renamed" ${prefix}.split-x "as there is no PAR in chromosome 23"
fi
```

**Check sex in X chromosome**

```{bash check-sex-X, echo=FALSE}
if [ $(awk -F' ' '{sum+=$5;} END {print sum;}' ${prefix}.split-x.fam) != '0' ]; then
  # Retain X chromosome and calculate r2 between the variants in a given window
  # This step requires that the .bim file has variant IDs in the second column
  # it will generate two files: plink.prune.in and plink.prune.out
  $plink --bfile ${prefix}.split-x --chr 23 --indep-pairwise 50 5 0.2
  
  # Extract unlinked variants in x
  $plink --bfile ${prefix}.split-x --extract plink.prune.in --make-bed --out ${prefix}.split-x.LD
  
  # Calculate heterozygosity F statistic of X chromosome
  $plink --bfile ${prefix}.split-x.LD  --check-sex 0.35 0.7  --out ${prefix}.split-x.chrX
  
  rm plink.*
  rm ${prefix}.split-x.LD.*
else
  echo "Please edit the .fam with the sex of the samples"
fi
```

Plot sex in X chromosome

```{r plot-sex-X}

if (!require("ggbeeswarm", quietly = TRUE))
install.packages("ggbeeswarm")
library(ggbeeswarm)

sex_X = read.delim(file = paste0(prefix,".split-x.chrX.sexcheck"), 
                    header = TRUE, 
                    sep = '')

# Plot these coefficients comparing males vs. females
# Roughly you expect female to have an F coefficient < 0.2-0.3 and males have an F coefficient > 0.7-0.8
# If there are individuals for which the fam has sex = 0 they will be plotted in an additional category called "no data"
if (any(fam$V5 %in% "0")) {
  ggplot(sex_X, 
       aes(x = factor(PEDSEX,
                      labels = c("Not disclosed", "Male", "Female")),
           y = F,
           color = PEDSEX)) +
    geom_quasirandom(alpha = 0.7,
                   size = 1.5) + 
    labs(title = "Chromosomal sex assignement in samples based in X chromosome",
         x = "Disclosed sex",
         y = "F coefficient X chromosome") +
  theme_minimal() +
  theme(legend.position = "none")
} else {
  ggplot(sex_X, 
       aes(x = factor(PEDSEX,
                      labels = c("Male", "Female")),
           y = F,
           color = PEDSEX)) +
    geom_quasirandom(alpha = 0.7,
                   size = 1.5) + 
    labs(title = "Chromosomal sex assignement in samples based in X chromosome",
         x = "Disclosed sex",
         y = "F coefficient X chromosome") +
  theme_minimal() +
  theme(legend.position = "none")
}
```
![F coefficient in X](https://github.com/acostauribe/genetic-data-QC/blob/apr17/redlat_result/7_Xcheck.png?raw=true)

**Check sex according to Y chromosome**

I recommend doing an initial check on the distribution of male vs female Y counts and then you can modify the detection thresholds after.
Detection count thresholds can be modified by changing the following flag: 
``--check-sex y-only [female max Y obs] [male min Y obs]``

```{bash check-sex-Y}
if [[ $(awk '$1 == "24" { count++ } END { print count }' ${prefix}.split-x.bim) != "" ]]; then
  if [ $(awk -F' ' '{sum+=$5;}END{print sum;}' ${prefix}.split-x.fam) != "0" ]; then
    $plink --bfile ${prefix}.split-x --check-sex y-only 2500 4000 --out ${prefix}.split-x.chrY
  else
    echo "Please edit the .fam with the sex of the samples"
  fi
else
  echo "Dataset does not contain Y chromosome. Y chromosome sex check was skipped"
fi
```

Plot sex in Y chromosome

```{r plot-sex-Y}
if (!require("ggbeeswarm", quietly = TRUE))
install.packages("ggbeeswarm")
library(ggbeeswarm)

# File path to search
Y_check_file = paste0(prefix,".split-x.chrY.sexcheck")

# Check if the file exists
if (file.exists(Y_check_file)) {
  sex_Y = read.delim(file = paste0(prefix,".split-x.chrY.sexcheck"),  
                    header = TRUE, 
                    sep = '')
  } else { 
    print("chrY.sexcheck file does not exist. Y chromosome sex check was skipped") 
}

if (exists("sex_Y")) {
# Plot the variant call count in chromosome Y
# If there are individuals for which the fam has sex = 0 they will be plotted in an additional category called "no data"
  if (any(fam$V5 %in% "0")) {
    ggplot(sex_Y, 
         aes(x = factor(PEDSEX,
                        labels = c("No Data", "Male", "Female")),
             y = YCOUNT,
             color = PEDSEX)) +
      geom_quasirandom(alpha = 0.7,
                     size = 1.5) + 
      labs(title = "Chromosomal sex assignement in samples based in Y chromosome",
           x = "Disclosed sex",
           y = "Y chromosome variant count") +
    theme_minimal() +
    theme(legend.position = "none")
    
    } else {
    ggplot(sex_Y, 
         aes(x = factor(PEDSEX,
                        labels = c("Male", "Female")),
             y = YCOUNT,
             color = PEDSEX)) +
      geom_quasirandom(alpha = 0.7,
                     size = 1.5) + 
      labs(title = "Chromosomal sex assignement in samples based in Y chromosome",
           x = "Disclosed sex",
           y = "Y chromosome variant count") +
    theme_minimal() +
    theme(legend.position = "none")
    }
}
```
![Variant Count in Y](https://github.com/acostauribe/genetic-data-QC/blob/apr17/redlat_result/7_Ycheck.png?raw=true)

**Identify individuals that fail sex-check**

```{r identify-sex-check-fails}
# Identify individuals that fail in each chromosome test
fail_x = sex_X[sex_X$STATUS == 'PROBLEM',]
if (exists("sex_Y")) {
  fail_y = sex_Y[sex_Y$STATUS == 'PROBLEM',]
}

# Merge 
if (exists("sex_Y")) {
  sex_fail = bind_rows(fail_x,
                       fail_y)
  } else {
  sex_fail = bind_rows(fail_x)
}

sex_fail_id = unique(select(sex_fail, FID, IID))

write.table(sex_fail_id, file = "sex_fail_samples.txt", 
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
print(paste(nrow(sex_fail_id), "samples have missmatched sex. IDs are written in sex_fail_samples.txt"))
```

These samples will be removed after identifying duplicates.

**Update your Sample Metrics dataframe with sex check**

Notice that this chunk expects genetic data from x chromosome. Y chromosome data is optional.

```{r update-sample-metrics-3}
# Add the disclosed sex to samples
sample_statistics$PEDSEX = sex_X$PEDSEX[match(sample_statistics$IID, sex_X$IID)]
# Add the calculated X chromosome sex
sample_statistics$X_SEX = sex_X$SNPSEX[match(sample_statistics$IID, sex_X$IID)]
# Add the calculated X chromosome sex F statistic
sample_statistics$X_SEX_F = sex_X$F[match(sample_statistics$IID, sex_X$IID)]

if (exists("sex_Y")) {
# Add the calculated Y chromosome sex
sample_statistics$Y_SEX = sex_Y$SNPSEX[match(sample_statistics$IID, sex_Y$IID)]
# Add the Y chromosome snp counts
sample_statistics$Y_SEX_count = sex_Y$YCOUNT[match(sample_statistics$IID, sex_Y$IID)]
}

# Save it as a file
write.table(sample_statistics,
            "sample_statistics.txt", 
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE)
```


## 8. Identify duplicate samples {#section-8}

Notice that KING requires to add the `.bed` at the end of the input

```{bash id-duplicates}
$king -b ${prefix}.split-x.bed --duplicate --rplot --prefix ${prefix}.king
```

`--duplicate` identifies duplicates or MZ twins and generates the file *\${prefix}.con*. Each line represents two samples that have heterozygote concordance rate \> 80%. These samples need to be carefully addressed to determine if they are *known* biological duplicates (e.g. same patient sequenced twice), or if they are supposed to be *different* samples. For the *known* biological duplicates, you should retain the genome with the best metrics. If the genomes are from *different* samples, both genomes should be removed (you cannot tell which one does the genome belong to)

> If you have known biological replicates I suggest to go over the duplicate list manually and decide which genome to exclude from each pair. If you are not expecting any duplicates, eliminate all of them

Make a list of duplicate samples:

```{r list-duplicate-samples}
king_file=paste0(prefix,".king.con")
duplicates_df = read.delim(king_file, header = TRUE)
duplicate_1 = select(duplicates_df, FID1, ID1)
colnames(duplicate_1) = c("FID", "IID")
duplicate_2 = select(duplicates_df, FID2, ID2)
colnames(duplicate_2) = c("FID", "IID")
duplicates_id = unique(bind_rows(duplicate_1, duplicate_2))

write.table(duplicates_id, file = "duplicated_samples.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)
print(paste(nrow(duplicates_id), "samples are found to be duplicates according to the KING algorythm. IDs are written in duplicated_samples.txt"))
```

Generate a list of individuals that are heterozygosity outliers, fail sex checks and/or duplicate samples

```{r write-het-sex-dup-removals}
het_sex_dup = unique(bind_rows(het_outlier_id,
                        sex_fail_id,
                        duplicates_id))

if(data_type=="EXOME" || data_type=="GENOME"){
write.table(select(het_sex_dup, IID),
            "het_sex_dup.irem", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)
  } else {
  write.table(het_sex_dup,
            "het_sex_dup.irem",
            col.names = FALSE,
            row.names = FALSE, 
            quote = FALSE)
  }
print(paste(nrow(het_sex_dup), "samples are found to be heterozygosity outliers, sex-fails and/or duplicates. IDs are written in het_sex_dup.irem and will be removed from dataset in next step"))
```

**Remove heterozygosity outliers, sex-fails and duplicates**

Note: For Plink datasets we will merge X PAR into X and turn haploid calls into missing calls before removing het_sex_dup.txt
```{bash remove-het-sex-dup}
if [[ "$data_type" == "GENOME" || "$data_type" == "EXOME" ]]; then
  vcftools --gzvcf ${prefix}.vcf.gz --remove het_sex_dup.irem --recode --recode-INFO-all --out ${prefix}.het.sex.dup
  mv  ${prefix}.het.sex.dup.recode.vcf ${prefix}.het.sex.dup.vcf
  bgzip ${prefix}.het.sex.dup.vcf
else
  $plink --bfile ${prefix}.split-x --merge-x --make-bed --out ${prefix}.merged-x
  $plink --bfile ${prefix}.merged-x --set-hh-missing --make-bed --out ${prefix}.merged-x.hh
  $plink --bfile ${prefix}.merged-x.hh --remove het_sex_dup.irem --make-bed --out  ${prefix}.het.sex.dup
  rm ${prefix}.merged-x.*
fi 
```

Oraganize all the generated files

```{bash organize-files-4}
if [ "$tidy_up" == "TRUE" ]; then
  mkdir Heterozygosity_stats
  mv *.het ./Heterozygosity_stats
  mv *.het.log ./Heterozygosity_stats
  #mv heterozygosity_outliers.txt ./Heterozygosity_stats
  
  mkdir Sex_check
  mv ${prefix}.split-x.chrX.* ./Sex_check
  if [ -f ${prefix}.split-x.chrY.sexcheck ]; then
    mv ${prefix}.split-x.chrY.* ./Sex_check
  fi
  #mv sex_fail_samples.txt ./Sex_check
  rm *.hh
  
  if [[ "$data_type" == "GENOME" || "$data_type" == "EXOME" ]]; then
  rm *.nosex
  fi
  
  mkdir Duplicate_analysis
  mv ${prefix}.king* ./Duplicate_analysis
  #mv duplicated_samples.txt ./Duplicate_analysis
  
  mv ${prefix}.split-x.bed ./File_evolution
  mv ${prefix}.split-x.bim ./File_evolution
  mv ${prefix}.split-x.fam ./File_evolution 
  if [ -f ${prefix}.split-x.log ]; then
   mv ${prefix}.split-x.log ./File_evolution
  fi
  
  if [ -f ${prefix}.bed ]; then
    mv ${prefix}.bed ./File_evolution
    mv ${prefix}.bim ./File_evolution
    mv ${prefix}.fam ./File_evolution 
    mv ${prefix}.log ./File_evolution
  fi
  
  mv het_sex_dup.irem ./File_evolution
  
  if [ -f ${prefix}.vcf.gz ]; then
    mv ${prefix}.vcf.gz ./File_evolution
  fi
  
  mkdir Missingness_stats
  mv *.miss.imiss ./Missingness_stats
  mv *.miss.lmiss ./Missingness_stats
  mv *.miss.log ./Missingness_stats
  
fi
```

Update prefix

```{r update-prefix-4}
prefix=(paste0(prefix,".het.sex.dup"))
Sys.setenv(prefix=prefix)
```

## 9. Calculate relatedness {#section-9}

This step is especially important if you are working with samples that you know have some relatedness between them, or if you are sampling from a small community. However it is optional and can be turned on/off in 2. Customize the quality control process 'check_relatedness'.

This process is done with the software [king](https://www.kingrelatedness.com/), which uses plink files as an input.

```{bash import-plink}
# If you have a vcf file, import it into plink 
if [[ "$check_relatedness" == 'FALSE' ]]; then
  echo "check_relatedness is set to FALSE. No kinship check will be performed"
elif [[ "$data_type" == "GENOME" || "$data_type" == "EXOME" ]]; then
  $plink --vcf ${prefix}.vcf.gz --keep-allele-order --double-id --vcf-half-call m --set-missing-var-ids @:#$1,$2 --make-bed --out ${prefix}
else
  echo "Data already in plink format"
fi
```

**Check Relatedness**

In this step you verify that disclosed relatedness among individuals matches their genetic relatedness. Additionally you search for the presence of cryptic relatedness between samples.

If the *.fam.* file does not have a Family Identifier per genome (FID), you need to incorporate that information. You can use Plink's [--update-ids](https://www.cog-genomics.org/plink/1.9/data#update_indiv).If you know that in your data set you have parent-offspring samples, its also useful to add this information to the fam with [--update-parents](https://www.cog-genomics.org/plink/1.9/data#update_indiv)

You can also use R to edit the *.fam* by running the following chunk

```{r edit-fam-3}
if(check_relatedness=="TRUE"){
  # Load your dataframes
  fam = read.delim(paste0(prefix,".fam"),
                  sep=' ',
                  header=FALSE)

  if(sum(fam$V5) > length(fam$V5)*0.25) {
    # Use match to replace the existing values in the .fam for those in sample_data
    fam$V5 = sample_data$SEX[match(fam$V2, sample_data$IID)]
    # Replace the NA for 0
    fam$V5[is.na(fam$V5)] = 0
    write.table(fam,
                paste0(prefix,".fam"), 
                col.names = FALSE, 
                row.names = FALSE, 
                quote = FALSE)
  }
  # Use match to replace the existing values in the .fam for those in sample_data
  fam$V1 = sample_data$FID[match(fam$V2, sample_data$IID, nomatch = 0)]
  fam$V3 = sample_data$PID[match(fam$V2, sample_data$IID, nomatch = 0)]
  fam$V4 = sample_data$MID[match(fam$V2, sample_data$IID, nomatch = 0)]

  write.table(fam,
              paste0(prefix,".fam"), 
              col.names = FALSE, 
              row.names = FALSE, 
              quote = FALSE)
  
  print(paste(".fam file will be annotated with:", sample_data_file))
  
} else {
  print("check_relatedness is set to FALSE. No kinship check will be performed")
}
```

After editing the *.fam* you can proceed to check relatedness:

```{bash id-relatedness}
if [ "$check_relatedness" == "TRUE" ]; then
  $king -b ${prefix}.bed --related --rplot --degree 3 --cluster --prefix ${prefix}.king
else
  echo "check_relatedness is set to FALSE. No kinship check will be performed"
fi
```

The output of this command produces (among many others) two important files: 

* *{prefix}.king.kin* Relatedness in disclosed relationships
* *{prefix}.king.kin0* Inferred relatedness in 'unrelated' individuals. 

Each row provides information for one pair of individuals. See [king documentation](https://www.kingrelatedness.com/manual.shtml#RELATED)

If you get `FATAL ERROR - Please correct problems with pedigree structure` it usually means that an individual who is, for example, listed as male but shows up as the mother of another individual or vice versa)

> ⚠️ HUMAN INPUT NEEDED: **YOU** need to manually analyze these two files and generate a table including the samples that have problematic relatedness and have to be excluded (*problematic_relatedness.txt*) If your data type is GENOME or EXOME: Each row should represent the individual ID of each sample. If your data type is an ARRAY: Each row is `Family ID, Individual ID`. The given IDs must match those in the *.fam*.

```{bash remove-problematic-relatedness}
# If you have provided a problematic_relatedness.txt file, it will remove the listed individuals.
if [ -f "problematic_relatedness.txt" ]; then
  if [ "$data_type" == "GENOME" ] || [ "$data_type" == "EXOME" ]; then
    $vcftools --gzvcf ${prefix}.vcf.gz --remove problematic_relatedness.txt --recode --recode-INFO-all --out ${prefix}.rel
    mv  ${prefix}.rel.recode.vcf ${prefix}.rel
    bgzip ${prefix}.rel.vcf
  else
    $plink --bfile ${prefix} --remove problematic_relatedness.txt --make-bed --out ${prefix}.rel
  fi
else
  echo "No 'problematic_relatedness.txt' file. No sample will be removed."
fi
```

Tidy up

```{bash organize-files-5}
if [ "$tidy_up" == "TRUE" ] && [ "$check_relatedness" == "TRUE" ]; then
  mkdir Relatedness_stats
  mv ${prefix}.king* ./Relatedness_stats
  
  if [ -f ${prefix}.rel.vcf.gz ]; then
    mv ${prefix}.vcf.gz ./File_evolution
  fi
  
  if [ -f ${prefix}.rel.bed ] || [ -f ${prefix}.rel.vcf.gz ]; then
    mv ${prefix}.bed ./File_evolution
    mv ${prefix}.bim ./File_evolution
    mv ${prefix}.fam ./File_evolution
    mv ${prefix}.log ./File_evolution
  fi
  
fi
```

Update system prefix if you removed individuals

```{r update-prefix-5}
if (file.exists("problematic_relatedness.txt")) {
prefix=(paste0(prefix,".rel"))
Sys.setenv(prefix=prefix)
}
```

**OPTIONAL: Check for Mendelian errors**

If you know you have trios, parent-offspring duos or multigenerational families in your data set, you can use these to your benefit and check for Mendelian inconsistencies. Plink offers different commands to [identify](https://www.cog-genomics.org/plink/1.9/basic_stats#mendel) and [remove](https://www.cog-genomics.org/plink/1.9/data#set_me_missing) these inconsistencies according to the individuals in your data set.

For plink formatted data we can use:

```{bash mendel-errors, eval=FALSE, include=FALSE}
if [ "$data_type" == "ARRAY" ]; then
  $plink --bfile ${prefix} --set-me-missing --mendel-duos --make-bed --out ${prefix}.me
fi
```

If you are working with a VCF there are [plugins](https://samtools.github.io/bcftools/howtos/plugin.mendelian.html) for bcftools or [programs](https://gatk.broadinstitute.org/hc/en-us/articles/360037594831-FindMendelianViolations-Picard-) that do it, however they usually require trios (both parents + offspring). Please refer to the source pages for its usage

```{bash organize-files-6}
if [ "$tidy_up" == "TRUE" ] && [ -f ${prefix}.me.bed ]; then
  mv ${prefix}.bed ./File_evolution
  mv ${prefix}.bim ./File_evolution
  mv ${prefix}.fam ./File_evolution
  mv ${prefix}.log ./File_evolution
fi
```

Update prefix
```{r update-prefix-6}
if (file.exists(paste0(prefix,".me.bed"))) {
  prefix=(paste0(prefix,".me"))
  Sys.setenv(prefix=prefix)
}
```

## 10. Remove variants & Individuals and with missingness ≥5% {#section-10}

After refining our data set we finally remove the variants and the individuals that are missing more than 5% of the data. 

You can change your threshold of inclusion with the ``--max-missing`` flag for vcf files [0 = all data allowed, 1 = no missing data allowed] or the ``--geno`` flag for plink files [1 = all data allowed, 0 = no missing data allowed].

```{bash final-missingness}
if [ "$data_type" == "EXOME" ] || [ "$data_type" == "GENOME" ]; then
  # Identify individuals with high missingness
  $vcftools --gzvcf ${prefix}.vcf.gz --missing-indv --out ${prefix}.mind-2
  awk '{if ($5 > 0.05) {print $1} }' ${prefix}.mind-2.imiss > ${prefix}.mind-2.irem
  awk '{if ($5 <= 0.05) {print $1} }' ${prefix}.mind-2.imiss > ${prefix}.mind-2.ikeep
  
  # Remove variants and sites missing more than 5% of data
  $vcftools --gzvcf ${prefix}.vcf.gz --keep ${prefix}.mind-2.ikeep --max-missing 0.9 --recode --recode-INFO-all --out ${prefix}.qc
  mv ${prefix}.qc.recode.vcf ${prefix}.qc.vcf
  bgzip ${prefix}.qc.vcf
  
else 
  # Remove variants and samples missing more than 5% of data
  $plink --bfile ${prefix} --geno 0.05 --mind 0.05 --make-bed --out ${prefix}.qc
fi
```

Tidy up

```{bash organize-files-7}
if [ "$tidy_up" == "TRUE" ]; then
  mkdir Post-QC_stats
  
  if [ "$data_type" == "EXOME" ] || [ "$data_type" == "GENOME" ]; then 
  mv *.mind-2.ikeep ./Missingness_stats
  mv *.mind-2.irem ./Missingness_stats
  mv *.mind-2.imiss ./Missingness_stats
  fi
  
  mv ${prefix}.* ./File_evolution
  mv ./File_evolution/${prefix}.qc* ./

fi
```

And now you have a clean data set!

```{r update-prefix-7}
prefix=(paste0(prefix,".qc"))
Sys.setenv(prefix=prefix)
```

If you feel like plotting your results, you can redo the plots from the 1st part:

```{bash get-postQC-stats}
if [ "$data_type" == "EXOME" ] || [ "$data_type" == "GENOME" ]
then
    # Data description
    bcftools stats ${prefix}.vcf.gz > ${prefix}.post-qc.vchk.txt
    
    # Individual missingness
    vcftools --gzvcf ${prefix}.vcf.gz --missing-indv --out ${prefix}.post-qc

    # Individual depth
    vcftools --gzvcf ${prefix}.vcf.gz --depth --out ${prefix}.post-qc

    # Site missingness
    vcftools --gzvcf ${prefix}.vcf.gz --missing-site --out ${prefix}.post-qc
    sed -i 's/chr//g' ${prefix}.post-qc.lmiss 
    
    # Site depth
    vcftools --gzvcf ${prefix}.vcf.gz --site-mean-depth --out ${prefix}.post-qc
    sed -i 's/chr//g' ${prefix}.post-qc.ldepth.mean
    
    # VQRS filtering
    vcftools --gzvcf ${prefix}.vcf.gz --FILTER-summary --out ${prefix}.post-qc
else
    # Use plink for Site and individual missingness in SNP Array. 
    $plink --bfile ${prefix} --missing --out ${prefix}.post-qc
fi
```

Plot missingess per variant and per individual *Before* and *After* the quality control process

**Individual based metrics**

```{r plot-postQC-stats-ind}
# Individual Missingness 
imiss_final = read.delim(file = paste0(prefix,".post-qc.imiss"), 
                    header = TRUE, 
                    sep = '')
# Generate basic summary statistics
imiss_final_F = describe(imiss_final$F_MISS)
rownames(imiss_final_F) = c("final_sample_missingness")
print(imiss_final_F)

# Plot Missingness rate per sample
imiss_hist_final = hist(imiss$F_MISS,
                         xlab="Missingness rate", 
                         ylab="Samples", 
                         main="Missingness rate per sample post QC",
                         col="gold1",
                         breaks=20)

# Compare with initial dataset
imiss_box_final = boxplot(imiss$F_MISS, imiss_final$F_MISS,
                         xlab="Missingness rate", 
                         ylab="Samples", 
                         main="Missingness rate per sample preQC vs postQC", 
                         col = c("cyan3", "gold1"),
                         names = c("preQC", "postQC"))

# Depth per Individual (Only for Exome or Genome data)
if(data_type=="EXOME" || data_type=="GENOME"){
  # Read depth calculations
  idepth_final = read.delim((paste0(prefix,".post-qc.idepth")), header = T, sep = "")
  # Generate basic summary statistics
  idepth_final_mean = describe(idepth_final$MEAN_DEPTH)
  rownames(idepth_final_mean) = c("final_sample_depth")
  print(idepth_final_mean)

  # Plot depth per sample
  idepth_hist_final = hist(idepth_final$MEAN_DEPTH,
                             xlab="Mean Depth",
                             ylab="Samples", 
                             main="Mean Depth per sample postQC", 
                             col="gold1",
                             breaks=50)
  idepth_box_final = boxplot(idepth$MEAN_DEPTH, idepth_final$MEAN_DEPTH,
                         ylab="Mean Depth", 
                         main="Mean depth per sample preQC vs postQC", 
                         col = c("cyan3", "gold1"),
                         names = c("preQC", "postQC"))
}
```
![Mean Depths Per Sample](https://github.com/acostauribe/genetic-data-QC/blob/apr17/redlat_result/10_MeanDepthPerSample.png?raw=true)
![Mean Depths comparison Pre-QC vs Post-QC](https://github.com/acostauribe/genetic-data-QC/blob/apr17/redlat_result/10_MeanSample_prevspost.png?raw=true)
![Missingness Per Sample](https://github.com/acostauribe/genetic-data-QC/blob/apr17/redlat_result/10_MissPerSample.png?raw=true)
![Missingness comparison Pre-QC vs Post-QC](https://github.com/acostauribe/genetic-data-QC/blob/apr17/redlat_result/10_MissSample_prevspost.png?raw=true)

Update the Sample Metrics file:

```{r  update-sample-metrics-4}
if(data_type=="EXOME" || data_type=="GENOME"){
  sample_statistics$Final_n_sites = idepth_final$N_SITES[match(sample_statistics$IID, idepth_final$INDV)]
  sample_statistics$Final_mean_depth = idepth_final$MEAN_DEPTH[match(sample_statistics$IID, idepth_final$INDV)]
  sample_statistics$Final_missingness = imiss_final$F_MISS[match(sample_statistics$IID, imiss_final$INDV)]
}

if(data_type=="ARRAY"){
  sample_statistics$Final_n_sites = imiss_final$N_GENO[match(sample_statistics$IID, imiss_final$IID)]
  sample_statistics$Final_missingness = imiss_final$F_MISS[match(sample_statistics$IID, imiss_final$IID)]
}

write.table(sample_statistics,
            "sample_statistics.txt", 
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE)
```

**Site based metrics**

```{r plot-postQC-stats-site}
# Site Missingness
lmiss_final = read.delim(file = paste0(prefix,".post-qc.lmiss"), 
                    header = TRUE, 
                    sep = '')
# Generate basic summary statistics
lmiss_final_F= describe(lmiss_final$F_MISS)
rownames(lmiss_final_F) = c("final_site_missingness")
print(lmiss_final_F)

# Plot Missingness rate site per chromosome
sorted_chr = unique(lmiss_final$CHR) #organize the chromosomes in ascending order
lmiss_final$CHR = factor(lmiss_final$CHR , levels=sorted_chr)
  
lmiss_box_final = boxplot(F_MISS ~ CHR,
                    data=lmiss_final,
                    xlab="Chromosome",
                    ylab="Missingness rate",
                    cex.axis = 0.5,
                    main="Missingness rate per site post-QC",
                    col="gold1")

# Depth per site - chromosome (Only for Exome or Genome data)
if(data_type=="EXOME" || data_type=="GENOME"){
  # Read depth calculations
  ldepth_final = read.delim((paste0(prefix,".post-qc.ldepth.mean")), header = T, sep = "")

  # Generate basic summary statistics
  ldepth_final_mean = describe(ldepth_final$MEAN_DEPTH)
  rownames(ldepth_final_mean) = c("final_site_depth")
  print(ldepth_final_mean)
  # Plot depth per site
  
  sorted_chr = unique(ldepth_final$CHR) #organize the chromosomes in ascending order
  ldepth_final$CHR = factor(ldepth_final$CHR , levels=sorted_chr)
  ldepth_box_final = boxplot(MEAN_DEPTH ~ CHR,
                       data=ldepth_final,
                       xlab="Chromosome",
                       ylab="Mean Depth per site",
                       cex.axis = 0.5,
                       main="Mean Depth per site post-QC",
                       col="gold1")
}
```
![Mean Depths Per Site](https://github.com/acostauribe/genetic-data-QC/blob/apr17/redlat_result/10_MeanDepthPerSite.png?raw=true)
![Missingness Per Site](https://github.com/acostauribe/genetic-data-QC/blob/apr17/redlat_result/10_MissPerSite.png?raw=true)

Update your dataset_statsitics dataframe

```{r update-dataset_statsitics}
dataset_statistics = bind_rows(dataset_statistics,
                                 imiss_final_F,
                                 lmiss_final_F)

if(data_type=="EXOME" || data_type=="GENOME"){
  dataset_statistics = bind_rows(dataset_statistics,
                                 idepth_final_mean,
                                 ldepth_final_mean)
}
# Save as a file (optional)
write.table(dataset_statistics,
            "dataset_statistics.txt", 
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE)

print(dataset_statistics)
```

Organize all the generated files

```{bash  organize-files-8}
if [ "$tidy_up" == "TRUE" ]; then
  mv ${prefix}.post-qc.* ./Post-QC_stats
  if [ -f "${prefix}.irem" ]; then
  mv ${prefix}.irem ./File_evolution
  fi
fi
```

**And done!! (pat yourself in the back)**
