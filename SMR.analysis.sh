### SMR software to implement the SMR & HEIDI methods

## Part1. Format file preparation
#--bfile (SNP)
#--gwas-summary (GWAS)
#--beqtl-summary (eQTL)

## Part2. SMR performing
smr-1.3.1 --bfile ./bfile/SNP.out --gwas-summary ./gwassummary/gemma.out --beqtl-summary ./eqtlsummary/eQTL.out --peqtl-smr 1e-6 --heidi-min-m 2 --out ./res/smr.cis.out --thread-num 20 --cis-wind 10000

smr-1.3.1 --bfile ./bfile/SNP.out --gwas-summary ./gwassummary/gemma.out --beqtl-summary ./eqtlsummary/eQTL.out --peqtl-smr 1e-6 --heidi-min-m 2 --out ./res/smr.trans.out --thread-num 20 --trans --trans-wind 5000

padj=p.adjust(smr.cis.out$p_SMR,method="BH")
