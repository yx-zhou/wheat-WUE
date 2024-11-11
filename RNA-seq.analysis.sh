### RNA-seq processing, DEGs identification and motif scanning

## Part1. RNA-seq processing
fastp -i {samplename}_1.fq.gz -I {samplename}_2.fq.gz -o clean.{samplename}_R1.fq.gz -O clean.{samplename}_R2.fq.gz --detect_adapter_for_pe -w 4 --compression 9 -h log/{samplename}.html -j log/{samplename}.json
hisat2 --dta -x ChineseSpring.hisat2.index -p 4 -1 clean.{samplename}_R1.fq.gz -2 clean.{samplename}_R2.fq.gz -S samfiles/{samplename}.sam
samtools view -@ 4 -bS samfiles/{samplename}.sam | samtools sort -@ 4 - -o bamfiles/{samplename}.bam
featureCounts -T 4 -t exon -p -P -B -C -g gene_id -a IWGSC_v1.1_HC_20170706.gtf -o featurecount/{samplename}.feacount.txt bamfiles/{samplename}.bam

## Part2. DEGs identification
library(tidyverse)
library(DESeq2)
library(stringr)
allcount <- read_delim("all.rawcount.txt")
colnames(allcount)=gsub("^[0-9]+_","",colnames(allcount))
setwd("compare/") #directory of files contain compare group

# The "comparegroup" files include the sample names for each group of comparison, which correspond to the sample names in the "allcount" table
# for example: call DEGs of CK-WT-Leaf and Drought-WT-Leaf, the "comparegroup" file:
# CFL1
# CFL2
# CFL3
# DFL1
# DFL2
# DFL3

for (comparegroup in list.files()) {
  group <- read.table(comparegroup)
  count <- allcount %>% select(Geneid, group$V1)
  count <- count %>% select(-1)
  rownames(count) <- allcount$Geneid
  group$V1 <- str_sub(group$V1, 1, 2)
  group$V1 <- factor(group$V1,levels=unique(group$V1))
  dds <- DESeqDataSetFromMatrix(countData = count, colData = group, design = ~ V1)
  deg <- DESeq(dds)
  res <- results(deg)
  result <- as.data.frame(res)
  up=subset(result,log2FoldChange > 1 & padj<0.05)
  down=subset(result,log2FoldChange < -1 & padj<0.05)
  print("deseq done")
  deg_file_name <- paste0(comparegroup, ".deg.txt")
  up_file_name <- paste0(comparegroup, ".up.txt")
  down_file_name <- paste0(comparegroup, ".down.txt")
  write.table(result, file = deg_file_name, sep = "\t", quote = FALSE,row.names = T,col.names = T) # write all DESeq2 result
  write.table(rownames(up), file = up_file_name, sep = "\t", quote = FALSE,row.names = F,col.names = F) # write up regulated DEGs 
  write.table(rownames(down), file = down_file_name, sep = "\t", quote = FALSE,row.names = F,col.names = F)  # write down regulated DEGs 
}
sessionInfo()

## Part3. motif scanning
grep -f deg.txt genepro2k.bed > deg.bed
bedtools getfasta -nameOnly -fi 161010_Chinese_Spring_v1.0_pseudomolecules.fasta -bed deg.bed -fo deg.fa
fimo --o deg mybmotif.meme deg.fa
cd deg
awk '{print$3}' fimo.tsv | sort | uniq > deg.target.txt
