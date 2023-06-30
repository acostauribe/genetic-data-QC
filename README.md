# Quality control pipeline for the ReD-Lat genomic data
#### Developed by Juliana Acosta-Uribe for the ReD-Lat Consortium 2023

The [redlat_qc.rmd](redlat_qc.rmd) pipeline is designed to be run as an [R markdown](https://rmarkdown.rstudio.com/lesson-1.html) file in R Studio. This way you can run it in a step-by-step mode. 

You could also run it directly from the r command line if you already have the `sample_data.txt` and the `problematic_relatedness.txt` files in your workspace.
```         
library(rmarkdown)
render("path/to/your/file.Rmd")
```

For this pipeline you will need the following:

### Data:

A plink formatted .*bed*, .*fam*, .*bim* set, or a bgzipped *.vcf.gz* file.

> **Always Take a look at the files before beginning:**
> If you start with a plink dataset: 
> - Do the Individual IDs (column 2 in *.fam*) match what you were expecting? 
> - Have the families been given a Family ID (column 1 in *.fam*)?  
> - Do the Individuals have their sex assigned (column 5 in *.fam*?  
> - Do the variants in the *.bim* have an identifier (column 2 in *.bim*)? Some analyses will require this information, and we may have to incorporate it to the *.fam*/*.bim* if its not already there. 
> Make sure your files are properly aligned and the alleles are being called from the correct strand. INDELs should be [left aligned and normalized](https://samtools.github.io/bcftools/bcftools.html#norm).

If you are starting with a *.vcf*, or your *.fam* does not have the sex /family of the samples already specified please provide an additional file `sample_data.txt` with a header as follows:

**IID** ID of the sample as it is in the plink IID or in the VCF\
**SEX** should be specified (1=male, 2=female, 0=no data)\
**FID** Family ID of the sample\
**PID** Paternal ID of the sample, if the individual is also in the dataset, the PID should be identical to their fathers IID\
**MID** Maternal ID of the sample, if the individual is also in the dataset, the PID should be identical to their mothers IID\
**PHENO** Optional, if you are interested in any downstream analyses involving the phenotype. (1=control, 2=case, -9=missing)

⚠️ The **FID** given to the parents needs to match the same **FID** given to that individual for their genome.

A trio would look like this (order of the columns does not matter)

| FID  | IID | PID | MID | SEX | PHENO |
|------|-----|-----|-----|-----|-------|
| FID1 | SON | DAD | MOM | 1   | 2     |
| FID1 | MOM | 0   | 0   | 2   | 1     |
| FID1 | DAD | 0   | 0   | 1   | 1     |

R is expecting a tab delimited file. If your file is delimited by spaces you can fix it with the following bash command `sed -i 's/ /\t/g'  sample_data.txt`


### Tools:

-[R](https://www.r-project.org/)\
-[RStudio](https://posit.co/download/rstudio-desktop/)\
-[plink](https://www.cog-genomics.org/plink2/)\
-[king](https://www.kingrelatedness.com/)

If you are working with Exome or Genome data, you will also need:

-[vcftools](https://vcftools.github.io/man_latest.html)\
-[bcftools](https://samtools.github.io/bcftools/bcftools.html)


### Workflow:

[1. Set up your environment](redlat_qc.rmd/#section-1)\
[2. Customize the quality control process](#section-2)\
[3. Start Quality Control process](#section-3)\
[4. Genotype Quality control](#section-4)\
[5. Filter for missingness](#section-5)\
[6. Calculate heterozygosity](#section-6)\
[7. Check sex](#section-7)\
[8. Identify duplicates](#section-8)\
[9. Calculate relatedness](#section-9)\
[10. Remove variants & Individuals and with missingness ≥5%](#section-10)
