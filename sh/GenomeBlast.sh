#!/bin/bash

################################################################################

# load modules
module load bioinfo-tools
module load blast

################################################################################

# arguments
genomeSuffix=$1
genomeBlast=${genomeSuffix}.out
queryInput=$2
outDir=$3

################################################################################

# set output directory
if [[ ! -d ${outDir} ]]
then
  mkdir -vp ${outDir}
fi

################################################################################

# make database from nucleotide sequence
if [[ ! -f ${genomeSuffix} ]]
then
  echo "Build BLAST database => ${genomeSuffix}"
  makeblastdb -in ${genomeSuffix}.fasta -dbtype nucl -out ${genomeSuffix}
fi

################################################################################

# use protein sequence as query to blast six reading frames
echo "Blasting database..."
tblastn -query ${queryInput}.fa -db ${genomeSuffix} -out ${genomeBlast}

################################################################################
