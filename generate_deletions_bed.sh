#!/bin/bash
export GATK=/Users/siakhnin/tools/gatk4/gatk/build/libs/gatk.jar
export PICARD=/Users/siakhnin/tools/picard/dist/picard.jar

DATA_PREF=/Users/siakhnin/data
mkdir ${DATA_PREF}/reports

SAMPLE=${DATA_PREF}/mnist_na12878_chrom20_100kb.rheos.kafka.no_dels.vcf
GATK_VCF=${DATA_PREF}/mnist_na12878_chrom20.100kb.gatk4.filtered.vcf
FB_VCF=${DATA_PREF}/mnist_na12878_chrom20.100kb.freebayes.minimap.filtered.vcf
GIAB_VCF=${DATA_PREF}/na12878_giab_highconf.9999999-10100000.filtered.vcf

java -jar ${PICARD} GenotypeConcordance \
CALL_VCF=${SAMPLE} \
TRUTH_VCF=${FB_VCF} \
O=${DATA_PREF}/reports/freebayes_to_rheos \
TRUTH_SAMPLE=NA12878 CALL_SAMPLE=NA12878

java -jar ${PICARD} GenotypeConcordance \
CALL_VCF=${SAMPLE} \
TRUTH_VCF=${GIAB_VCF} \
O=${DATA_PREF}/reports/giab_to_rheos \
TRUTH_SAMPLE=HG001 CALL_SAMPLE=NA12878

java -jar ${PICARD} GenotypeConcordance \
CALL_VCF=${SAMPLE} \
TRUTH_VCF=${GATK_VCF} \
O=${DATA_PREF}/reports/gatk4_go_rheos \
TRUTH_SAMPLE=NA12878 CALL_SAMPLE=NA12878

vcftools --vcf \
${GIAB_VCF} \
--diff ${SAMPLE} \
--diff-site \
--out ${DATA_PREF}/reports/giab_to_rheos

vcftools --vcf \
${FB_VCF} \
--diff ${SAMPLE} \
--diff-site \
--out ${DATA_PREF}/reports/freebayes_to_rheos

vcftools --vcf \
${GATK_VCF} \
--diff ${SAMPLE} \
--diff-site \
--out ${DATA_PREF}/reports/gatk4_to_rheos
