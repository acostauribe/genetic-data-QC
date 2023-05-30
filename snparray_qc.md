# Quality Control PIPELINE for the RedLat SNP array

## Developed by Juliana Acosta-Uribe 2023

For this pipeline you will need the following:\

**Data:**\
1. A plink formatted `.bed`, `.fam`, `.bim` data set.\
2. The reported sex of the samples, either as an annotation in column 5 of the *.fam* or as a separate text file *sex_info.txt* with 3 columns `Individual_ID, Individual_ID, SEX` (1=male, 2=female, 0=no data) \
3. Information about the individuals that belong to the same family, either as an the same Family ID in column 1 of the *.fam* or as a separate text file *family_info.txt* with4 columns `Individual_ID, Individual_ID, Family_ID, Individual_ID`\
4. If the data set contains parent-offspring samples, the ID of the parents should be stated in columns 3 and 4 of the *.fam*, or given in a separate file *parent_info.txt* with 4 columns:
`Family_ID, Individual_ID, Paternal_ID, Maternal_ID`\

**Tools:**\
[R](https://www.r-project.org/),
[RStudio](https://posit.co/download/rstudio-desktop/),
[plink](https://www.cog-genomics.org/plink2/),
[king](https://www.kingrelatedness.com/)

This pipeline is designed to be run as an [R
markdown](https://rmarkdown.rstudio.com/lesson-1.html) file in R Studio.
This way you can incorporate bash and R commands into a single process.

**Contents Table**

[0. Set up R Markdown]\
[1. Filter for missingness]\
[2. Calculate heterozygosity]\
[3. Check sex]\
[4. Summarize all sample metrics]\
[5. Calculate relatedness]\
[6. Remove variants & Individuals and with missingness ≥5%]

## 0. Set up R Markdown
```{r markdown-setup}
# Install requires packages:
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
#install.packages("igraph")
library(igraph)

# Set your working directory:
knitr::opts_chunk$set(root.dir = "/home/acostauribe/Genetic-Sequencing_raw-data/HudsonAlpha_J.Nicholas.Cochran/redlat/SNP_array/2022.11.18-SA-Batch3/2022.11.18-SA-Batch3-Preimputation", tidy=TRUE)
setwd("/home/acostauribe/Genetic-Sequencing_raw-data/HudsonAlpha_J.Nicholas.Cochran/redlat/SNP_array/2022.11.18-SA-Batch3/2022.11.18-SA-Batch3-Preimputation")

# Set up R variables:
PREFIX='redlat_array'

# Set up environment variables and path to software:
Sys.setenv(PREFIX='redlat_array')
Sys.setenv(PLINK='/home/acostauribe/bin/plink')
Sys.setenv(KING='/home/acostauribe/bin/king')
```

**ALWAYS** Take a look at the files before beginning:\
* Do the Individual IDs (column 2 in *.fam*) match what you were expecting?  \
* Have the families been given a Family ID (column 1 in *.fam*)? \
* Do the Individuals have their sex assigned (column 5 in *.fam*? \
* Do the variants in the *.bim* have an identifier (column 2 in *.bim*)? \
Some analyses will require this information, and we may have to incorporate it to the *.fam*/*.bim* if its not already there.

## 1. Filter for missingness

```{bash filter-missingness}
# A. Calculate missingness rates. Per individual ends in .imiss and per variant ends in .lmiss
$PLINK --bfile ${PREFIX} --missing --out ${PREFIX}.pre-qc

# B. Remove individuals missing more than 5% of data
$PLINK --bfile ${PREFIX} --mind 0.05 --make-bed --out ${PREFIX}.mind

# C. Remove variants missing more than 10% of data
$PLINK --bfile ${PREFIX}.mind --geno 0.1 --make-bed --out ${PREFIX}.mind.geno

# D. Calculate missingness rates after this filter
$PLINK --bfile ${PREFIX}.mind.geno --missing --out ${PREFIX}.mind.geno
```

**NOTICE** that individuals missing more than 5% of the data get removed in this step "X people removed due to missing genotype data (--mind)",l plink writes the ids to PREFIX.mind.irem

**Plot missingess per individual:**

```{r missingess-plots}
# A. Load data into R
## Before QC
imiss = read.delim(file = paste0(PREFIX,".pre-qc.imiss"), 
                    header = TRUE, 
                    sep = '')
## After initial QC
imiss_qc = read.delim(file = paste0(PREFIX,".mind.geno.imiss"), 
                    header = TRUE, 
                    sep = '')

# B. Generate basic summary statistics
## Before QC
imiss_F_MISS = describe(imiss$F_MISS)
print(imiss_F_MISS)
## After initial QC
imiss_qc_F_MISS = describe(imiss_qc$F_MISS)
print(imiss_qc_F_MISS)

# C. Plot Missingness rate per sample
## Before QC
hist(imiss$F_MISS,
     xlab="Missingness rate", #label of x
     ylab="Samples", # label of y
     main="Missingness rate per sample preQC", # main label
     col="lavender",
     breaks=20)
## After initial QC
hist(imiss_qc$F_MISS,
     xlab="Missingness rate",
     ylab="Samples", 
     main="Missingness rate per sample. \nIndividual missing <5% and \nvariants miss <10%", 
     col="hotpink",
     breaks=20)
## Compare data-sets
boxplot(imiss$F_MISS, imiss_qc$F_MISS,
     xlab="Missingness rate",
     ylab="Samples", 
     main="Missingness rate per sample", 
     col = c("lavender", "hotpink"),
     names=c("all", "iMiss <5%, varMiss <10%"))
```

## 2. Calculate heterozygosity

Calculate heterozygosity per sample before QC `${PREFIX}` and after the --mind --geno filtering `${PREFIX}.mind.geno` This step requires that the *.bim* file has variant IDs in the second column. If no variants have been assigned, you could do a preliminary step using [--set-missing-var-ids](https://www.cog-genomics.org/plink/1.9/data#set_missing_var_ids)

```{bash check-het-outliers}
for file in ${PREFIX} ${PREFIX}.mind.geno
do
     # A. Retain variants with MAF > 10% and individuals with low missing %
     $PLINK --bfile ${file} --maf 0.1 --make-bed --out ${file}.maf
     # B. Calculate LD
     #--indep-pairwise <window size>['kb'] <step size (variant ct)> <r^2 threshold>
     $PLINK --bfile ${file}.maf --indep-pairwise 50 10 0.2 
     # C. Retain variants not in LD (independent markers)
     $PLINK --bfile ${file}.maf --extract plink.prune.in --make-bed --out ${file}.maf.ld
     # D. Check heterozygosity
     $PLINK --bfile ${file}.maf.ld --het --out ${file}.het
done
```

This step uses a *sliding window* approach to identify variants in linkage disequilibrium. There are many options to modify the behavior or
this approach in [plink's docummentation](https://www.cog-genomics.org/plink/1.9/ld#indep)

**Plot heterozygosity**

```{r plot-het-outliers}
# A. Load data into R
## i. Before QC
het = read.delim(file = paste0(PREFIX,".het.het"), 
                    header = TRUE, 
                    sep = '')

## ii. After initial QC
het_qc = read.delim(file = paste0(PREFIX,".mind.geno.het.het"), 
                    header = TRUE, 
                    sep = '')

# B. Generate basic summary statistics
## i. Before QC
heterozygosity_sample_raw = describe(het$F)
print(heterozygosity_sample_raw)
## Calculate limits for excluding samples [3 standard deviations]
### Low threshold
heterozygosity_low_limit_raw = mean(het$F)-(3*(sd(het$F)))
print(heterozygosity_low_limit_raw)
### High threshold
heterozygosity_high_limit_raw = mean(het$F)+(3*(sd(het$F)))
print(heterozygosity_high_limit_raw)

## ii. After initial QC
heterozygosity_sample_qc = describe(het_qc$F)
print(heterozygosity_sample_qc)
## Calculate limits for excluding samples [3 standard deviations]
### Low threshold
heterozygosity_low_limit_qc = mean(het_qc$F)-(3*(sd(het_qc$F)))
print(heterozygosity_low_limit_qc)
### High threshold
heterozygosity_high_limit_qc = mean(het_qc$F)+(3*(sd(het_qc$F)))
print(heterozygosity_high_limit_qc)

# C. Plot heterozygosity per sample
## i. Before QC
hist(het$F,  
     freq=TRUE, 
     xlab="Heterozygosity F coefficient",  
     ylab="Samples", 
     main="Heterozygosity rate per sample - Raw dataset ",
     col="lavender",
     breaks=100)
abline(v = (heterozygosity_low_limit_raw), col="red")
abline(v = (heterozygosity_high_limit_raw), col="red")
abline(v = (mean(het$F)), col="blue") 
legend("topleft",
       c("+/-3 SD","mean"),
       col=c("red","blue"),
       pch=16)
### Individuals whose heterozygosity deviated more than 3 SD from the mean should be identified
het_outlier_low = filter(het, F<heterozygosity_low_limit_raw)
het_outlier_high = filter(het, F>heterozygosity_high_limit_raw)
het_outlier_both = bind_rows(het_outlier_low,
                             het_outlier_high)
write.table(het_outlier_both, 
            "samples_heterozygosity_outliers_raw dataset.txt",
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')

## ii. After initial QC
hist(het_qc$F,  
     freq=TRUE, 
     xlab="Heterozygosity F coefficient",  
     ylab="Samples", 
     main="Heterozygosity rate per sample - \nIndividual miss <5% and var miss <10% ",
     col="hotpink",
     breaks=50)
abline(v = (heterozygosity_low_limit_qc), col="red")
abline(v = (heterozygosity_high_limit_qc), col="red")
abline(v = (mean(het_qc$F)), col="blue") 
legend("topleft",
       c("+/-3 SD","mean"),
       col=c("red","blue"),
       pch=16)
### Individuals whose heterozygosity deviated more than 3 SD from the mean should be identified
het_outlier_low_qc = filter(het_qc, F<heterozygosity_low_limit_qc)
het_outlier_high_qc = filter(het_qc, F>heterozygosity_high_limit_qc)
het_outlier_both_qc = bind_rows(het_outlier_low_qc,
                             het_outlier_high_qc)

## Compare data-sets
boxplot(het$F, het_qc$F,
     xlab="Dataset",
     main="Heterozygosity rate", 
     col = c("lavender","paleturquoise3"),
     names=c("raw data", "iMiss <5%, varMiss <10%"))

# D. Identify samples that 'fail' heterozygosity testing
# Generate a table with the heterozygosity outliers
write.table(het_outlier_both, 
            "samples_heterozygosity_outliers.txt",
            col.names = TRUE, 
            row.names = TRUE, 
            quote = FALSE,
            sep = '\t')
```

It is useful to visualize heterozygosity F statistic vs. missingness per sample

```{r heterozygosity-vs-missing}
# Before QC
plot(imiss$F_MISS, het$F,
     xlab="Missingness",  
     ylab="heterozygosity", 
     main="Heterozygosity rate per sample - Data before QC")
abline(h = (heterozygosity_low_limit_raw), col="red")
abline(h = (heterozygosity_high_limit_raw), col="red")
abline(h = (mean(het$F)), col="blue") 

# After initial QC
plot(imiss_qc$F_MISS, het_qc$F,
     xlab="Missingness",  
     ylab="heterozygosity", 
     main="Heterozygosity rate per sample - Data post QC")
abline(h = (heterozygosity_low_limit_qc), col="red")
abline(h = (heterozygosity_high_limit_qc), col="red")
abline(h = (mean(het_qc$F)), col="blue") 
```

## 3. Check sex

**Preprocess sex chromosomes**

A. Verify that *${PREFIX}.mind.geno.fam* has a 'sex' value on the 5 column.\
If the values are 0 or -9, you can assign 'sex value' in plink: *sex_info.txt* has the sex information in three columns.  
Each row is `Family_ID/Individual_ID, Individual_ID, SEX` (1=male, 2=female, 0=no data). The given IDs must match those in the *.fam*

```{bash assign-sex}
# This command will look for a file called sex_info.txt, if it doesn't find it, it will append a .sex to your bed fam bim to ensure compatibility with further steps
if [[ $(find sex_info.txt) == sex_info.txt ]]
then
$PLINK --bfile ${PREFIX}.mind.geno --update-sex sex_info.txt --make-bed --out ${PREFIX}.mind.geno.sex
else
$PLINK --bfile ${PREFIX}.mind.geno --make-bed --out ${PREFIX}.mind.geno.sex
fi
```

B. Split Pseudo-Autosomic region (PAR) of X\
If the files have already been processed to split the PAR, this command won't work.\
Make sure you are using the correct genome alignment for this. e.g. ReDLat files are aligned to hg38. Set up `GENOME_ALIGN` variable accordingly

```{bash split-PAR}
GENOME_ALIGN='hg38'

# This command will check if there is already a XY region in ${PREFIX}.mind.geno.sex.bim, if it doesn't find one it will split the PAR as chromosome 25. If there is already a chromosome 25, then it will append .split-x to ensure compatibility further on
if [[ $(awk '$1 == "25" { count++ } END { print count }' ${PREFIX}.mind.geno.sex.bim) == "" ]]
then 
$PLINK --bfile ${PREFIX}.mind.geno.sex --split-x $GENOME_ALIGN --make-bed --out ${PREFIX}.mind.geno.sex.split-x
else
$PLINK --bfile ${PREFIX}.mind.geno.sex --make-bed --out ${PREFIX}.mind.geno.sex.split-x
fi

```

You can find more information on this step in [plink docummentation](https://www.cog-genomics.org/plink/1.9/data#split_x)

**Check sex in X chromosome**

```{bash check-sex-X}
## a. Retain X chromosome and calculate r2 between the variants in a given window
## This step requires that the .bim file has variant IDs in the second column
$PLINK --bfile ${PREFIX}.mind.geno.sex.split-x --chr 23 --indep-pairwise 50 5 0.2
## This command will generate two files: plink.prune.in y plink.prune.out

## b. Extract unlinked variants in x
$PLINK --bfile ${PREFIX}.mind.geno.sex.split-x --extract plink.prune.in --make-bed --out ${PREFIX}.mind.geno.sex.split-x.LD

## c. Calculate heterozygosity F statistic of X chromosome
$PLINK --bfile ${PREFIX}.mind.geno.sex.split-x.LD  --check-sex 0.3 0.7 --out ${PREFIX}.mind.geno.sex.split-x.LD.chrX
```

Plot sex in X chromosome

```{r plot-sex-X}
#install.packages("ggbeeswarm")
library(ggbeeswarm)

sex_X = read.delim(file = paste0(PREFIX,".mind.geno.sex.split-x.LD.chrX.sexcheck"), 
                    header = TRUE, 
                    sep = '')

## Plot these coefficients comparing males vs. females
## Roughly you expect Female to have an F coefficient < 0.2-0.3 and males have an F coefficient > 0.7-0.8
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
ggsave("Chromosomal X sex assignement in samples.png", width = 8, height = 5) 
```

**Check sex according to Y chromosome**

```{bash check-sex-Y}
$PLINK --bfile ${PREFIX}.mind.geno.sex.split-x --check-sex y-only --out ${PREFIX}.mind.geno.sex.split-x.chrY
```

Plot sex in Y chromosome

```{r plot-sex-Y, echo=TRUE}
#install.packages("ggbeeswarm")
library(ggbeeswarm)

sex_Y = read.delim(file = paste0(PREFIX,".mind.geno.sex.split-x.chrY.sexcheck"), 
                    header = TRUE, 
                    sep = '')

## Plot these coefficients comparing males vs. females
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
ggsave("Chromosomal Y sex assignement in samples.png", width = 8, height = 5) 
```

## 4. Summarize all sample metrics

Generate a data-frame with all th statistics you have calculated so far

```{r generate_sample_metrics}
# A. Start with the file where you calculated 'missingness' in all samples
sample_metrics = imiss
# B. Add the heterozygosity F statistic
sample_metrics$F_het = het$F[match(sample_metrics$IID, het$IID)]
# C. Add missingness rate to samples after filtering with --mind and --geno
sample_metrics$geno_0.1 = imiss_qc$F_MISS[match(sample_metrics$IID, imiss_qc$IID)]
# D. Add the heterozygosity F statistic to samples after filtering with --mind and --geno
sample_metrics$F_hetqc = het_qc$F[match(sample_metrics$IID, het_qc$IID)]
# E. Add the disclosed sex to samples
sample_metrics$PEDSEX = sex_X$PEDSEX[match(sample_metrics$IID, sex_X$IID)]
# F. Add the calculated X chromosome sex
sample_metrics$X_SEX = sex_X$SNPSEX[match(sample_metrics$IID, sex_X$IID)]
# G. Add the calculated X chromosome sex F statistic
sample_metrics$X_SEX_F = sex_X$F[match(sample_metrics$IID, sex_X$IID)]
# H. Add the calculated Y chromosome sex
sample_metrics$Y_SEX = sex_Y$SNPSEX[match(sample_metrics$IID, sex_Y$IID)]
# I. Add the Y chromosome snp counts
sample_metrics$Y_SEX_count = sex_Y$YCOUNT[match(sample_metrics$IID, sex_Y$IID)]
# J. Save it as a file
write.table(sample_metrics,
            "sample_metrics.unimputed-dataset.txt", 
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE,
            sep = '\t')
```

**Remove samples that failed heterozygosity and sex check**

> ⚠️ HUMAN INPUT NEEDED: **You** will have to generate a file (*samples_het_sex_exclude.txt*) with the samples that are heterozygosity outliers or that disclosed sex doesn't match genetic sex. Each row is `Family_ID/Individual_ID, Individual_ID`. The given IDs must match those in the *.fam*

```{bash remove-het-sex-fails}
$PLINK --bfile ${PREFIX}.mind.geno.sex --remove samples_het_sex_exclude.txt --set-hh-missing --make-bed --out ${PREFIX}.mind.geno.sex.het
```

## 5. Calculate relatedness

This is specially important if you are working with samples that you know have some relatedness between them, or if you are sampling from a small community.

**Identify duplicated samples**

Notice that KING requires to add the `.bed` at the end of the input

```{bash id-duplicates}
$KING -b ${PREFIX}.mind.geno.sex.het.bed --duplicate --rplot --prefix ${PREFIX}
```

`--duplicate` identifies duplicates or MZ twins and generates the file *${PREFIX}.con*\
Each line represents two samples that have heterozygote concordance rate > 80%. These samples need to be carefully addressed to determine if they are *known* biological duplicates (e.g. same patient sequenced twice), or if they are supposed to be *different* samples. For the *known* biological duplicates, you should retain the genome with the best metrics. If the genomes are from *different* samples, both genomes should be removed (you cannot tell which one does the genome belong to)

> ⚠️ HUMAN INPUT NEEDED: **You** will have to generate a file (*duplicates.txt*) with the samples that will be excluded. Each row is `Family_ID/Individual_ID, Individual_ID`. The given IDs must match those in the *.fam*

Remove duplicate samples:

```{bash remove-duplicates}
$PLINK --bfile ${PREFIX}.mind.geno.sex.het --remove duplicates.txt --make-bed --out ${PREFIX}.mind.geno.sex.het.dup
```

**Check Relatedness**

In this step you verify that disclosed relatedness among individuals matches their genetic relatedness, additionally you search for the presence of cryptic relatedness between samples.

If the *.fam.* file does not have a Family Identifier per genome (FID), you need to incorporate that information. [--update-ids](https://www.cog-genomics.org/plink/1.9/data#update_indiv) expects a file `family_info.txt` with the following four columns: `Old family ID`, `Old within-family ID`, `New family ID`, `Old/New within-family ID`.

```{bash update-familyID}
# This command will look for a file called family_info.txt if it doesn't find it, it will append a .fid to your bed fam bim to ensure compatibility with further steps
if [[ $(find family_info.txt) == family_info.txt ]]
then
$PLINK --bfile ${PREFIX}.mind.geno.sex.het.dup --update-ids family_info.txt --make-bed --out ${PREFIX}.mind.geno.sex.het.dup.fid
else
$PLINK --bfile ${PREFIX}.mind.geno.sex.het.dup --make-bed --out ${PREFIX}.mind.geno.sex.het.dup.fid
fi
```

If you know that in your data set you have parent-offspring samples, its also useful to add this information to the fam. [--update-parents](https://www.cog-genomics.org/plink/1.9/data#update_indiv) expects a file `parent_info.txt` with the following four columns: `Old family ID`, `Old within-family ID`, `New paternal within-family ID` and `New maternal within-family ID`.

The `within-family ID` given to the parents needs to match the same `within-family ID` given to that individual for their genome.\
A trio .fam you look like this:\
FAM1 | SON | DAD | MOM | 1 | 2\
FAM1 | MOM | 0 | 0 | 2 | 1\
FAM1 | DAD | 0 | 0 | 1 | 1\

```{bash update-parentalID}
# This command will look for a file called parent_info.txt if it doesn't find it, it will append a .pid to your bed fam bim to ensure compatibility with further steps
if [[ $(find parent_info.txt == parent_info.txt) ]]
then
$PLINK --bfile ${PREFIX}.mind.geno.sex.het.dup.fid --update-parents parent_info.txt --make-bed --out ${PREFIX}.mind.geno.sex.het.dup.fid.pid
else
$PLINK --bfile ${PREFIX}.mind.geno.sex.het.dup.fid --make-bed --out ${PREFIX}.mind.geno.sex.het.dup.fid.pid
fi
```

After editing the *.fam* you can proceed to check relatedness:

```{bash id-relatedness}
$KING -b ${PREFIX}.mind.geno.sex.het.dup.fid.pid.bed --related --rplot --degree 4 --cluster --prefix ${PREFIX}
```

The output of this command produces two files: *${PREFIX}.kin* Relatedness in reported relationships *${PREFIX}.kin0*
Inferred relatedness in 'unrelated' individuals. Each row above provides information for one pair of individuals. See
[king documentation](https://www.kingrelatedness.com/manual.shtml#RELATED)

if you get
`FATAL ERROR - Please correct problems with pedigree structure` it usually means that an individual who is, for example, listed as male but shows up as the mother of another individual 9or viceversa)

> ⚠️ HUMAN INPUT NEEDED: **YOU** need to manually analyze these two files and generate a table including the samples that have problematic relatedness and have to be excluded (*problematic_relatedness.txt*)Each row is `Family_ID, Individual_ID`. The given IDs must match those in the *.fam*

```{bash remove-problematic-relatedness}
$PLINK --bfile ${PREFIX}.mind.geno.sex.het.dup.fid.pid --remove problematic_relatedness.txt --make-bed --out ${PREFIX}.mind.geno.sex.het.dup.fid.pid.rel
```

**OPTIONAL: Check for Mendelian errors**

If you know you have trios, parent-offspring duos or multigenerational families in your data set, you can use these to your benefit and check for Mendelian inconsistencies. Plink offers different commands to [identify](https://www.cog-genomics.org/plink/1.9/basic_stats#mendel) and [remove](https://www.cog-genomics.org/plink/1.9/data#set_me_missing) this inconsistencies according to the individuals in your data set.

Since in the ReDLat cohort we have parent-offspring duos, we can use:
```{bash mendel-errors}
$PLINK --bfile ${PREFIX}.mind.geno.sex.het.dup.fid.pid --set-me-missing --mendel-duos --make-bed --out ${PREFIX}.mind.geno.sex.het.dup.fid.pid.rel.me
```

## 6. Remove variants & Individuals and with missingness ≥5%

After refining our data set we finally remove the variants and the individuals that are missing more than 5% of the data
```{bash final-missingness}
# This command will look for a file ${PREFIX}.mind.geno.sex.het.dup.fid.pid.rel.me.bed if it doesn't find it, it will us ${PREFIX}.mind.geno.sex.het.dup.fid.pid.rel.bed instead
if [[ $(find ${PREFIX}.mind.geno.sex.het.dup.fid.pid.rel.me.bed == ${PREFIX}.mind.geno.sex.het.dup.fid.pid.rel.me.bed ]]
then
$PLINK --bfile ${PREFIX}.mind.geno.sex.het.dup.fid.pid.rel.me --mind 0.05 --geno 0.05 --make-bed --out ${PREFIX}.qc
else
$PLINK --bfile ${PREFIX}.mind.geno.sex.het.dup.fid.pid.rel --mind 0.05 --geno 0.05 --make-bed --out ${PREFIX}.qc
fi
```

And now you have a clean data set!

If you feel like plotting your results, you can redo the plots from the 1st part:
```{bash qc-plots}
$PLINK --bfile ${PREFIX}.qc --missing --out ${PREFIX}.qc.stats
```

Plot missingess per variant and per individual *Before* and *After* the quality control process
```{r missingess-plots-qc}
# A. Load data into R
## After Finalizing QC
imiss_final = read.delim(file = paste0(PREFIX,".qc.stats.imiss"), 
                    header = TRUE, 
                    sep = '')

# B. Generate basic summary statistics
## After Finalizing QC
imiss_final_F_MISS = describe(imiss_final$F_MISS)
print(imiss_final_F_MISS)

# C. Plot Missingness rate per sample
## After Finalizing QC
hist(imiss_final$F_MISS,
     xlab="Missingness rate",
     ylab="Samples", 
     main="Missingness rate per sample after quality control", 
     col="aquamarine",
     breaks=20)
## Compare data-sets
boxplot(imiss$F_MISS, imiss_qc$F_MISS, imiss_final$F_MISS,
     xlab="Missingness rate",
     ylab="Samples", 
     main="Missingness rate per sample", 
     col = c("lavender", "hotpink", "aquamarine"),
     names=c("all", "iMiss <5%, varMiss <10%", "final dataset"))
```

Update the Sample Metrics file:
```{r update_sample_metrics}
# A. Start with the file where you calculated 'missingness' in all samples
sample_metrics$Final_MISS = imiss_final$F_MISS[match(sample_metrics$IID, imiss_final$IID)]
write.table(sample_metrics,
            "sample_metrics.unimputed-dataset.txt", 
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE,
            sep = '\t')
```

**And done!! (pat yourself in the back)**
