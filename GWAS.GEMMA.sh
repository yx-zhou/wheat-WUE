### To perform GWAS analysis based on GEMMA

## Part1. Data preparation
# Convert vcf files to bed/bim/fam files
plink --vcf genotype.vcf --allow-extra-chr --make-bed --out genotype.data
# Kinship calculation
gemma -bfile genotype.data -gk 1 -o kinship

## Part2. GEMMA performing
gemma -bfile genotype.data -k kinship.cXX.txt -c GAPIT.PCA.txt -n 1 -lmm 1 -o gwas.out.lmm -miss 1.0 -notsnp -r2 1.0 -hwe 0

## Part3. Manhattan and QQ-plot visualization
library(CMplot)
CMplot(gwas.out.lmm, col=mycol,
       plot.type = "m",
       LOG10=TRUE, pch=19,
       threshold =c(3.5e-4),
       threshold.col=c('red'),
       threshold.lty = c(1),
       threshold.lwd = c(2),
       multracks = FALSE, 
       multraits= TRUE,
       amplify = FALSE, 
       file="jpg", dpi = 300,
       file.output = TRUE,
       file.name= "multtrait.Manhattan")
	   
CMplot(gwas, col=mycol,
       plot.type = "q",
       LOG10=TRUE, pch=19,
       threshold =c(3.5e-4),
       threshold.col=c('red'),
       threshold.lty = c(1),
       threshold.lwd = c(2),
       multracks = FALSE,
       multraits= TRUE, 
       amplify = FALSE, 
       file="jpg", dpi = 300,
       file.output = TRUE,
       file.name= "multtrait.qqplot")

## Part4. Phenotypic variation explained (PVE) calculation
library(tidyverse)
data=read_tsv("gwas.out.lmm.assoc.txt")
data$pve = (2*(data$beta^2)*data$af*(1-data$af))/((2*(data$beta^2)*data$af*(1-data$af)) + ((data$se^2)*2*(230-data$n_miss)*data$af*(1-data$af)))
