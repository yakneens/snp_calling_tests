#!/bin/bash

DATA_PREF=/Users/siakhnin/data
SOURCE_VCF=${DATA_PREF}/giab/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz
REGION="20:9999999-13000000"
FLANK=5
RESULT_FILENAME=${DATA_PREF}/na12878_giab_highconf.9999999-13000000.deletions.5bp_flank_sorted_merged.bed

GET_INDELS_CMD="bcftools view -r ${REGION} -v indels -M2 ${SOURCE_VCF} "
GET_DELS_BED_CMD="vcf2bed --deletions - "
FLANK_CMD="bedtools flank -b 5 -g ${DATA_PREF}/reference/genome.bed.ref "
SORT_CMD="bedtools sort "
MERGE_CMD="bedtools merge "

FULL_CMD="${GET_INDELS_CMD} | \
 ${GET_DELS_BED_CMD} | \
 ${FLANK_CMD} | \
 ${SORT_CMD} | \
 ${MERGE_CMD} > \
 ${RESULT_FILENAME}"
printf "${FULL_CMD}"
eval ${FULL_CMD}