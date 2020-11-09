#!/bin/bash

################################################################################

file=$1
projDir=$2
verbose=$3

################################################################################

# load modules
module load bioinfo-tools
module load samtools
module load htseq/0.12.4

################################################################################

# declare variables
align="alignments/"
measure="measurements/"
annot="Annotation/"
# projDir="/crex/proj/snic2020-15-35/private/RNAExpressionProfile/"

################################################################################

# project directory
cd ${projDir}

################################################################################

# create alignment indexes
for alignmentFile in $(ls ${csct}${align}*bam)
do
  fl=${alignmentFile/*\//}
  f=${fl/.bam/}
  [[ $verbose == "true" ]] && echo "Indexing ${f} file..." || echo -n ""
  samtools index ${alignmentFile} &
done
wait

################################################################################

# count reads
for alignmentFile in $(ls ${csct}${align}*bam)
do
  fl=${alignmentFile/*\//}
  f=${fl/.bam/}
  [[ $verbose == "true" ]] && echo "Counting ${f} file..." || echo -n ""
  htseq-count \
    -n 16 \
    ${csct}${align}${f}.bam \
    ${annot}gencode.v35.annotation.gff3.gz \
    > ${csct}${measure}${f}.out \
    2> log/${f}.err &
done

cd - > /dev/null

################################################################################
