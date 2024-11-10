### To perform co-expression analysis based on WGCNA R package

library(WGCNA)
setwd(outdir)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
ALLOW_WGCNA_THREADS = 48

outdir <- "./WGCNA/"
cutoff <- 0.5
maxBlockSize <- 30000
mergeCutHeight <- 0.25
preser <- as.numeric(1)
notfilter <- as.numeric(1)

## Step1. Data input, cleaning and pre-processing
# Expression data loading
dataExpr <- read.table("./CPM.txt",header=T, row.names=1, sep="\t")
sp_num <- ncol(dataExpr)
datExpr0 <- as.matrix(dataExpr)
if(notfilter==1){
  FltData <- datExpr0
}else{
  FltData <- varFilter(datExpr0, var.func = IQR, var.cutoff = cutoff, filterByQuantile = TRUE)
}
datExpr0 <- as.data.frame(t(FltData))
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
if(!gsg$allOK){
  if (sum(!gsg$goodGenes)>0){  
    write.table(names(datExpr0)[!gsg$goodGenes], file=paste(outdir, "Genes.txt", sep="/"), row.names=F, col.names=F, quote=F)
  }  
  if (sum(!gsg$goodSamples)>0){
    write.table(names(datExpr0)[!gsg$goodSamples], file=paste(outdir, "Samples.txt", sep="/"), row.names=F, col.names=F, quote=F)
  }
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree <- flashClust(dist(datExpr0,method='manhattan'), method = "average")
pdf(file=paste(outdir,"sampleClustering.pdf", sep="/"), width = 12, height = 9)
sizeGrWindow(12, 9)
par(cex = 1)
par(mar = c(2, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.5)
cutHeights <- 400000
abline(h=cutHeights, col="red")
dev.off()

clust <- cutreeStatic(sampleTree, cutHeight=cutHeights, minSize=10)
table(clust)
keepSamples <- (clust==1)
datExpr <- datExpr0[keepSamples, ]
write.table(t(datExpr), file=paste(outdir, "InputData.deal.txt", sep="/"), quote=F, sep="\t")
save(datExpr, file=paste(outdir, "dataInput.RData", sep="/"))

# Phenotype data loading
dataTraits <- read.table("./pheno.txt", header=T, sep="\t")
match_trait <- match(rownames(datExpr), dataTraits$acc)
rownames(dataTraits) <- dataTraits$acc
Traits <- dataTraits[match_trait, -1]

## Step2. One-step network construction and module detection
# sft choosing
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <-  pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
r2 <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
threshold <- 0.9
powerEstimate <- 0
for(i in 1:length(powers)){
  if (r2[i] >= threshold){
    powerEstimate <- powers[i]
    break
  }else if(powerEstimate > 15 || powerEstimate == 0){
    if(nSamples < 20){
      powerEstimate <- 10
    }else if(nSamples >= 20 && nSamples < 30){
      powerEstimate <- 9	
    }else if(nSamples >= 30 && nSamples < 40){
      powerEstimate <- 8
    }else if(nSamples >= 40 && nSamples < 60){
      powerEstimate <- 7
    }else if(nSamples >= 60){
      powerEstimate <- 6
    }
  }
}

for (i in 1:length(powers)){
  if (powerEstimate == powers[i]){
    break
  }
}
r2line <- round(r2[i] * 100 + 0.5)/100
pdf(file=paste(outdir, "soft-thresholding.power.pdf", sep="/"), width = 9, height = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
     type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=r2line, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", 
     type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()
powerEstimate
# One-step network construction
net <- blockwiseModules(datExpr, power=6, TOMType="unsigned", 
                        saveTOMs = TRUE, minModuleSize = 50,
                        reassignThreshold = 0, mergeCutHeight = mergeCutHeight,
                        saveTOMFileBase="blockTOM", verbose = 3,
                        numericLabels = TRUE, pamRespectsDendro = FALSE)
table(net$colors) 
      

## Step3. The module assignment and module eigengene information
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
mergedColorsgene=cbind(as.data.frame(mergedColors),as.data.frame(colnames(datExpr)))
mergedColorsgeneME=cbind(mergedColorsgene,as.data.frame(net$colors))
modulegenenumber=mergedColorsgene %>% aggregate(colnames(datExpr)~mergedColors,FUN=length)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree,file = "networkConstruction-auto.RData")
      
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
modTraitCor <- cor(MEs, Traits, use = "p")
modTraitP <- corPvalueStudent(modTraitCor, nSamples)

## Step4. Exporting to Cytoscape visualization (weight>0.2)
TOM = TOMsimilarityFromExpr(datExpr, power = 3)
module ="red"
probes = names(datExpr)
inModule = (moduleColors == module)
modProbes = probes[inModule]
modGenes = modProbes
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("02.CytoscapeInput-edges-weight-th0.05-", module , ".txt", sep=""),
                               nodeFile = paste("02.CytoscapeInput-nodes-weight-th0.05-", module, ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.05,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
