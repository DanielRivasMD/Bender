#!/bin/bash

################################################################################

file=$1
verbose=$2

################################################################################

# load modules
module load bioinfo-tools
module load samtools
module load htseq/0.12.4

################################################################################

# declare directory structure variables
align="alignments/"
measure="measurements/"
annot="Annotation/"

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
    -a 20 \
    ${csct}${align}${f}.bam \
    ${annot}gencode.v35.annotation.gff3.gz \
    > ${csct}${measure}${f}.out \
    2> log/${f}.err &
done

################################################################################
