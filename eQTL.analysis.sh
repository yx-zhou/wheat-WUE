### eQTL Identification by MartixEQTL and GEMMA

## Part1. MartixEQTL perforimg
library(MatrixEQTL)
#Load genotype data
snps = SlicedData$new()
snps$fileDelimiter = "\t"
snps$fileOmitCharacters = "NA"
snps$fileSkipRows = 1
snps$fileSkipColumns = 1
snps$fileSliceSize = 2000
snps$LoadFile("./geno/genotype.012.txt") 
# Load gene expression data
gene = SlicedData$new()
gene$fileDelimiter = "\t"
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile("./cpm/cpm.norm.txt") 
# Load covariates
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"
cvrt$fileOmitCharacters = "NA"
cvrt$fileSkipRows = 1
cvrt$fileSkipColumns = 1
cvrt$LoadFile("./covariates/PCA.txt") 

output_file_name = tempfile()
meq = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = "./eQTL/MatrixEQTL.res.txt", 
  pvOutputThreshold =1e-6,  
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE,
  pvalue.hist = "qqplot")
  
pdf(file = "./eQTL/MatrixEQTL.qqplot.pdf", height=8,width=8)
plot(meq, pch = 16, cex = 0.7)
dev.off()

## Part2. GEMMA perforimg
gemma -bfile genotype.data -k kinship.cXX.txt -c PCA.txt -n 1 -lmm 1 -o eQTL.gemma.out.lmm -miss 1.0 -notsnp -r2 1.0 -hwe 0
