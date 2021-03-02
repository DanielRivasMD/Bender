#!/bin/bash

################################################################################

genome=$1
genomeDir=$2
library=$3
libraryDir=$4
dbDir=$5
dbDir=${outDir}/${genome}

################################################################################

# set output directory
if [[ ! -d ${dbDir} ]]
then
  mkdir -vp ${dbDir}
fi

################################################################################

# make database from nucleotide sequence
if [[ ! -f ${dbDir}/${genome} ]]
then
  echo "Build BLAST database => ${genome}"
  makeblastdb -in ${genomeDir}/${genome}.fasta -dbtype nucl -out ${dbDir}/${genome}
fi

################################################################################

# use protein sequence as query to blast six reading frames
echo "Blasting database..."
tblastn -query ${libraryDir}/${library} -db ${dbDir}/${genome} -out ${dbDir}/${genome}.out

################################################################################

# purge filtered file
if [[ -f ${dbDir}/${genome}.filtered.out ]]
then
  echo "Cleaning filtered..."
  ${dbDir}/${genome}.filtered.out
fi

################################################################################
