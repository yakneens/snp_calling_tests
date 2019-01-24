#!/bin/bash

DATA_PREF=/Users/siakhnin/data
INPUT_BAM=${DATA_PREF}/giab/RMNISTHS_30xdownsample.chr20.bam
RESULT_FILENAME=${DATA_PREF}/giab/RMNISTHS_30xdownsample_50000000-63025520.name_sorted.fastq
REGION="20:50000000-63025520"

VIEW_CMD="samtools view -b ${INPUT_BAM} ${REGION} "
SORT_CMD="samtools sort -n "
FASTQ_CMD="samtools fastq - "


FULL_CMD="${VIEW_CMD} | \
 ${SORT_CMD} | \
 ${FASTQ_CMD} > \
 ${RESULT_FILENAME}"
printf "${FULL_CMD}"
eval ${FULL_CMD}