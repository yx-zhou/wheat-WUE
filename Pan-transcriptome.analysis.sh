### RNA-seq processing, DEGs identification and traits bias

## Part1. RNA-seq processing
fastp -i ./rawdata/sample1_R1.fq.gz -o ./cleandata/clean.sample1_R1.fq.gz -w 4 --compression 9 -h ./cleandata/sample1.html -j ./cleandata/sample1.json &&
echo ">>>fastp done....." &&
cutadapt -j 4 -a "A{10}" -o ./cleandata/sample1.noPolyA_R1.fq.gz ./cleandata/clean.sample1_R1.fq.gz &&
echo ">>>cut polyA.done....." &&
rm ./cleandata/clean.sample1_R1.fq.gz &&
hisat2 --dta -x ~/genome/ChineseSpring/hisat2_index/ChineseSpring.hisat2 -p 4 -U ./cleandata/sample1.noPolyA_R1.fq.gz -S ./hisat2/sample1.sam &&
echo ">>>hisat2 done....." &&
samtools view -@ 4 -bS ./hisat2/sample1.sam | samtools sort -@ 4 - -o ./hisat2/sample1.bam &&
rm ./hisat2/sample1.sam &&
samtools flagstat -@ 4 ./hisat2/sample1.bam > ./hisat2/sample1.stat.txt &&
echo ">>>samtools done....." &&
featureCounts -T 4 -t exon  -g gene_id -a ~/genome/ChineseSpring/Triticum_aestivum.IWGSC.only.gtf -o ./featureCount/sample1.feacount.txt ./hisat2/sample1.bam &&
echo "done"

## Part2. DEG identification
# Method1. edgeR
library(edgeR)
exprSet=read.table(Counts.txt, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
group_list <- factor(c(rep("CK",1), rep("Treat",1)))
y <- DGEList(counts=exprSet, group = group_list)
y <- calcNormFactors(y)
bcv <- 0.1
et <- exactTest(y, dispersion = bcv ^ 2)
write.table(et$table,file = "edgeR.res.txt",sep = "\t",quote = F,row.names = T)

# Method2. GFOLD
gfold diff -s1 Control.txt -s2  Treat.txt -o ControlVSTreat

## Part3. Traits bias
# Category defining
library(rdist)
data=read_tsv("Traids.txt")
data$Balanced=cdist(data[5:7],t(c(1/3,1/3,1/3)),metric="euclidean",p=2)
data$Asupressed=cdist(data[5:7],t(c(0,0.5,0.5)),metric="euclidean",p=2)
data$Bsupressed=cdist(data[5:7],t(c(0.5,0,0.5)),metric="euclidean",p=2)
data$Dsupressed=cdist(data[5:7],t(c(0.5,0.5,0)),metric="euclidean",p=2)
data$Adominant=cdist(data[5:7],t(c(1,0,0)),metric="euclidean",p=2)
data$Bdominant=cdist(data[5:7],t(c(0,1,0)),metric="euclidean",p=2)
data$Ddominant=cdist(data[5:7],t(c(0,0,1)),metric="euclidean",p=2)
minname.fun=function(data){a=which.min(data);return(a)} 
data$class=apply(data[,8:14],1,minname.fun) 
data=data %>% mutate(group=case_when(class == 1 ~ "Balanced", class == 2 ~ "Asuppressed",class == 3 ~ "Bsuppressed", class == 4 ~ "Dsuppressed",class == 5 ~ "Adominant", class == 6 ~ "Bdominant", class == 7 ~ "Ddominant"))
write_tsv(data,"Traids.category.txt")

# Visualization by ternary
library(ggtern)
p=ggtern(data,aes(Aperave,Bperave,Dperave,color=group)) +
  geom_point(aes(alpha=0.8)) +
  theme_showarrows() +
  scale_color_manual(values=mycol)+
  xlab("A") + ylab("B") + zlab("D")
ggsave("Traids.pdf",p,width =8, height =8)
