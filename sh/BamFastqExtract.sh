#!/bin/bash

################################################################################

bamFile=$1
bamDir=$2
verbose=$3
storageDir=$4

################################################################################

# load modules
module load bioinfo-tools
module load samtools
module load BEDTools

################################################################################

# make temporary directory
if [[ "$SNIC_TMP" != "/scratch" ]]
then
  tmpDir="${SNIC_TMP}/${SLURMD_NODENAME}_${SLURM_JOBID}/"
  mkdir ${tmpDir}
fi
[[ $verbose == "true" ]] && echo "directory created" || echo -n ""

cp ${bamDir}${bamFile}.bam ${tmpDir}
[[ $verbose == "true" ]] && echo "file copied" || echo -n ""

################################################################################

# extract fastq from bam
samtools view -u -f 1 -F 12 ${bamDir}${bamFile}.bam > ${tmpDir}/${bamFile}_map_map.bam
# R1 unmapped, R2 mapped
samtools view -u -f 4 -F 264 ${bamDir}${bamFile}.bam > ${tmpDir}/${bamFile}_unmap_map.bam
# R1 mapped, R2 unmapped
samtools view -u -f 8 -F 260 ${bamDir}${bamFile}.bam > ${tmpDir}/${bamFile}_map_unmap.bam
# R1 & R2 unmapped
samtools view -u -f 12 -F 256 ${bamDir}${bamFile}.bam > ${tmpDir}/${bamFile}_unmap_unmap.bam
[[ $verbose == "true" ]] && echo "extracted" || echo -n ""

################################################################################

samtools merge -u ${tmpDir}/${bamFile}_unmapped.bam ${tmpDir}/${bamFile}_unmap_map.bam ${tmpDir}/${bamFile}_map_unmap.bam ${tmpDir}/${bamFile}_unmap_unmap.bam
[[ $verbose == "true" ]] && echo "merged" || echo -n ""

################################################################################

samtools sort -n -T ${tmpDir}/${bamFile} ${tmpDir}/${bamFile}_map_map.bam -o ${tmpDir}/${bamFile}_mapped.sort.bam
samtools sort -n -T ${tmpDir}/${bamFile} ${tmpDir}/${bamFile}_unmapped.bam -o ${tmpDir}/${bamFile}_unmapped.sort.bam
[[ $verbose == "true" ]] && echo "sorted" || echo -n ""

bamToFastq -i ${tmpDir}/${bamFile}_mapped.sort.bam -fq ${tmpDir}/${bamFile}_mapped_1.fastq -fq2 ${tmpDir}/${bamFile}_mapped_2.fastq
bamToFastq -i ${tmpDir}/${bamFile}_unmapped.sort.bam -fq ${tmpDir}/${bamFile}_unmapped_1.fastq -fq2 ${tmpDir}/${bamFile}_unmapped_2.fastq

################################################################################

cat ${tmpDir}/${bamFile}_mapped_1.fastq ${tmpDir}/${bamFile}_unmapped_1.fastq > ${tmpDir}/${bamFile}_1.fastq
cat ${tmpDir}/${bamFile}_mapped_2.fastq ${tmpDir}/${bamFile}_unmapped_2.fastq > ${tmpDir}/${bamFile}_2.fastq
[[ $verbose == "true" ]] && echo "concatenated" || echo -n ""

################################################################################

# verify fastq extraction
echo -e "\nflagstats${lineBreaks}"
samtools flagstat ${bamDir}${bamFile}.bam
echo -e "${lineBreaks}"

echo -e "mapped: $( samtools view -c ${tmpDir}/${bamFile}_mapped.sort.bam )"
echo -e "unmapped: $( samtools view -c ${tmpDir}/${bamFile}_unmapped.sort.bam )"

[[ $verbose == "true" ]] && echo "verified" || echo -n ""

################################################################################

# compress
gzip ${tmpDir}/${bamFile}_[12].fastq
[[ $verbose == "true" ]] && echo "compressed" || echo -n ""

################################################################################

# retrieving to storage
mv ${tmpDir}/${bamFile}_[12].fastq.gz ${storageDir}
[[ $verbose == "true" ]] && echo "moved" || echo -n ""

################################################################################
