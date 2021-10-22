#!/bin/bash

################################################################################

fastaFile=$1
fastaDir=$2
verbose=$3
storageDir=$4

echo $fastaFile
echo $fastaDir
echo $verbose
echo $storageDir

# ################################################################################

# # load modules
# module load bioinfo-tools
# module load samtools
# module load bwa

# if [[ "$fastQ" == "true" ]]
# then
#   module load Fastx
# fi

# ################################################################################

# # make temporary directory
# if [[ "$SNIC_TMP" != "/scratch" ]]
# then
#   tmpDir="${SNIC_TMP}/${SLURMD_NODENAME}_${SLURM_JOBID}/"
#   mkdir ${tmpDir}
# fi
# [[ $verbose == "true" ]] && echo "directory created" || echo -n ""

# cp ${fastaDir}${fastaFile}.fa ${tmpDir}
# [[ $verbose == "true" ]] && echo "file copied" || echo -n ""

# ################################################################################

# if [[ "$fastaQ" == "true" ]]
# then

#   # TODO: review fasta_quality_filter
#   # fastq_quality_filter \
#   #   -v \
#   #   -q 30 \
#   #   -p 80 \
#   #   -i ${tmpDir}${fastaFile} \
#   #   -o ${tmpDir}${fastaFile}

# fi

# ################################################################################

# # alignment to reference genome (paired-end)
# bwa mem \
#   -t 40 \
#   ${tmpDir}${reference} \
#   ${tmpDir}${} | \
# samtools view \
#   -F 2048 \
#   -u \
#   -t \
#   ${tmpDir}${reference}.fai - | \
# samtools sort \
#   -l 9 \
#   -T ${tmpDir}${} \
#   - \
#   -o ${tmpDir}${}.sorted.RemoveSupplementary.bam

# ################################################################################



# # extract fastq from bam
# samtools view -u -f 1 -F 12 ${fastaDir}${fastaFile}.bam > ${tmpDir}/${fastaFile}_map_map.bam
# # R1 unmapped, R2 mapped
# samtools view -u -f 4 -F 264 ${fastaDir}${fastaFile}.bam > ${tmpDir}/${fastaFile}_unmap_map.bam
# # R1 mapped, R2 unmapped
# samtools view -u -f 8 -F 260 ${fastaDir}${fastaFile}.bam > ${tmpDir}/${fastaFile}_map_unmap.bam
# # R1 & R2 unmapped
# samtools view -u -f 12 -F 256 ${fastaDir}${fastaFile}.bam > ${tmpDir}/${fastaFile}_unmap_unmap.bam
# [[ $verbose == "true" ]] && echo "extracted" || echo -n ""

# ################################################################################

# samtools merge -u ${tmpDir}/${fastaFile}_unmapped.bam ${tmpDir}/${fastaFile}_unmap_map.bam ${tmpDir}/${fastaFile}_map_unmap.bam ${tmpDir}/${fastaFile}_unmap_unmap.bam
# [[ $verbose == "true" ]] && echo "merged" || echo -n ""

# ################################################################################

# samtools sort -n -T ${tmpDir}/${fastaFile} ${tmpDir}/${fastaFile}_map_map.bam -o ${tmpDir}/${fastaFile}_mapped.sort.bam
# samtools sort -n -T ${tmpDir}/${fastaFile} ${tmpDir}/${fastaFile}_unmapped.bam -o ${tmpDir}/${fastaFile}_unmapped.sort.bam
# [[ $verbose == "true" ]] && echo "sorted" || echo -n ""

# bamToFastq -i ${tmpDir}/${fastaFile}_mapped.sort.bam -fq ${tmpDir}/${fastaFile}_mapped_1.fastq -fq2 ${tmpDir}/${fastaFile}_mapped_2.fastq
# bamToFastq -i ${tmpDir}/${fastaFile}_unmapped.sort.bam -fq ${tmpDir}/${fastaFile}_unmapped_1.fastq -fq2 ${tmpDir}/${fastaFile}_unmapped_2.fastq

# ################################################################################

# cat ${tmpDir}/${fastaFile}_mapped_1.fastq ${tmpDir}/${fastaFile}_unmapped_1.fastq > ${tmpDir}/${fastaFile}_1.fastq
# cat ${tmpDir}/${fastaFile}_mapped_2.fastq ${tmpDir}/${fastaFile}_unmapped_2.fastq > ${tmpDir}/${fastaFile}_2.fastq
# [[ $verbose == "true" ]] && echo "concatenated" || echo -n ""

# ################################################################################

# # verify fastq extraction
# echo -e "\nflagstats${lineBreaks}"
# samtools flagstat ${fastaDir}${fastaFile}.bam
# echo -e "${lineBreaks}"

# echo -e "mapped: $( samtools view -c ${tmpDir}/${fastaFile}_mapped.sort.bam )"
# echo -e "unmapped: $( samtools view -c ${tmpDir}/${fastaFile}_unmapped.sort.bam )"

# [[ $verbose == "true" ]] && echo "verified" || echo -n ""

# ################################################################################

# # compress
# gzip ${tmpDir}/${fastaFile}_[12].fastq
# [[ $verbose == "true" ]] && echo "compressed" || echo -n ""

# ################################################################################

# # retrieving to storage
# mv ${tmpDir}/${fastaFile}_[12].fastq.gz ${storageDir}
# [[ $verbose == "true" ]] && echo "moved" || echo -n ""

# ################################################################################
