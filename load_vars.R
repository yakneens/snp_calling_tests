library(VariantAnnotation)
library(Rsamtools)
library(GenomicAlignments)
# f = readVcf("~/data/snp_calling_tests/freebayes/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.freebayes.vcf.gz", "hs37d5")
# g = readVcf("~/data/snp_calling_tests/hc/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.gatk4.vcf.gz", "hs37d5")
# s = readVcf("~/data/snp_calling_tests/samtools/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.samtools.vcf.gz", "hs37d5")
# p = readVcf("~/data/snp_calling_tests/platypus/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.platypus.vcf.gz", "hs37d5")
# 
# #Genotype likelihoods
# 
# i = 7
# geno(f)$GL[[i]]
# geno(p)$GL[[i]]
# geno(g)$PL[[i]]
# geno(s)$PL[[i]]

my_snp = GRanges(Rle("20"), IRanges(10000758, 10000758))
param = ScanBamParam(which=my_snp, what=scanBamWhat())
bam = scanBam("~/data/snp_calling_tests/inputs/CEUTrio.HiSeq.WGS.b37.NA12878.20.10000758.20.bam", param=param)[[1]]
bam_dt = as.data.table(bam)

pileup_fun <- function(x){
  return(x)
}

aparam = ApplyPileupsParam(which=my_snp, what=c("seq", "qual"))
pfiles = PileupFiles("~/data/snp_calling_tests/inputs/CEUTrio.HiSeq.WGS.b37.NA12878.20.10000758.05.bam")
res = applyPileups(files = pfiles, pileup_fun, param=aparam)

pparam = PileupParam()
x = pileup("~/data/snp_calling_tests/inputs/CEUTrio.HiSeq.WGS.b37.NA12878.20.10000758.20.bam", scanBamParam = param)

aln = readGAlignments(file = "~/data/snp_calling_tests/inputs/CEUTrio.HiSeq.WGS.b37.NA12878.20.10000758.20.bam")
