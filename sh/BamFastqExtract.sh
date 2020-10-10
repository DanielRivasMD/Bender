#!/bin/bash

################################################################################

# load modules
module load bioinfo-tools
module load samtools
module load BEDTools

################################################################################

# make temporary directory
if [[ "$SNIC_TMP" != "/scratch" ]]
then
  TMP_DIR="${SNIC_TMP}/${SLURMD_NODENAME}_${SLURM_JOBID}/"
  mkdir ${TMP_DIR}
fi
[[ $VERBOSE == "true" ]] && echo "directory created" || echo -n ""

cp ${ADDRESS}${INDIVIDUAL}.bam ${TMP_DIR}
[[ $VERBOSE == "true" ]] && echo "file copied" || echo -n ""

################################################################################

# extract fastq from bam
samtools view -u -f 1 -F 12 ${ADDRESS}${INDIVIDUAL}.bam > ${TMP_DIR}/${INDIVIDUAL}_map_map.bam
# R1 unmapped, R2 mapped
samtools view -u -f 4 -F 264 ${ADDRESS}${INDIVIDUAL}.bam > ${TMP_DIR}/${INDIVIDUAL}_unmap_map.bam
# R1 mapped, R2 unmapped
samtools view -u -f 8 -F 260 ${ADDRESS}${INDIVIDUAL}.bam > ${TMP_DIR}/${INDIVIDUAL}_map_unmap.bam
# R1 & R2 unmapped
samtools view -u -f 12 -F 256 ${ADDRESS}${INDIVIDUAL}.bam > ${TMP_DIR}/${INDIVIDUAL}_unmap_unmap.bam
[[ $VERBOSE == "true" ]] && echo "extracted" || echo -n ""

################################################################################

samtools merge -u ${TMP_DIR}/${INDIVIDUAL}_unmapped.bam ${TMP_DIR}/${INDIVIDUAL}_unmap_map.bam ${TMP_DIR}/${INDIVIDUAL}_map_unmap.bam ${TMP_DIR}/${INDIVIDUAL}_unmap_unmap.bam
[[ $VERBOSE == "true" ]] && echo "merged" || echo -n ""

################################################################################

samtools sort -n -T ${TMP_DIR}/${INDIVIDUAL} ${TMP_DIR}/${INDIVIDUAL}_map_map.bam -o ${TMP_DIR}/${INDIVIDUAL}_mapped.sort.bam
samtools sort -n -T ${TMP_DIR}/${INDIVIDUAL} ${TMP_DIR}/${INDIVIDUAL}_unmapped.bam -o ${TMP_DIR}/${INDIVIDUAL}_unmapped.sort.bam
[[ $VERBOSE == "true" ]] && echo "sorted" || echo -n ""

bamToFastq -i ${TMP_DIR}/${INDIVIDUAL}_mapped.sort.bam -fq ${TMP_DIR}/${INDIVIDUAL}_mapped_1.fastq -fq2 ${TMP_DIR}/${INDIVIDUAL}_mapped_2.fastq
bamToFastq -i ${TMP_DIR}/${INDIVIDUAL}_unmapped.sort.bam -fq ${TMP_DIR}/${INDIVIDUAL}_unmapped_1.fastq -fq2 ${TMP_DIR}/${INDIVIDUAL}_unmapped_2.fastq

################################################################################

cat ${TMP_DIR}/${INDIVIDUAL}_mapped_1.fastq ${TMP_DIR}/${INDIVIDUAL}_unmapped_1.fastq > ${TMP_DIR}/${INDIVIDUAL}_1.fastq
cat ${TMP_DIR}/${INDIVIDUAL}_mapped_2.fastq ${TMP_DIR}/${INDIVIDUAL}_unmapped_2.fastq > ${TMP_DIR}/${INDIVIDUAL}_2.fastq
[[ $VERBOSE == "true" ]] && echo "concatenated" || echo -n ""

################################################################################

# verify fastq extraction
echo -e "\nflagstats${lineBreaks}"
samtools flagstat ${ADDRESS}${INDIVIDUAL}.bam
echo -e "${lineBreaks}"

echo -e "mapped: $( samtools view -c ${TMP_DIR}/${INDIVIDUAL}_mapped.sort.bam )"
echo -e "unmapped: $( samtools view -c ${TMP_DIR}/${INDIVIDUAL}_unmapped.sort.bam )"

[[ $VERBOSE == "true" ]] && echo "verified" || echo -n ""

################################################################################

# compress
gzip ${TMP_DIR}/${INDIVIDUAL}_[12].fastq
[[ $VERBOSE == "true" ]] && echo "compressed" || echo -n ""

################################################################################

# retrieving to storage
mv ${TMP_DIR}/${INDIVIDUAL}_[12].fastq.gz ${STORAGE_folder}
[[ $VERBOSE == "true" ]] && echo "moved" || echo -n ""

################################################################################
