#!/bin/bash

################################################################################

preffixGenome=$1
genomeDir=$2
library=$3
libraryDir=$4
outDir=$5
dbDir=${outDir}db/

################################################################################

# set output directory
if [[ ! -d ${dbDir} ]]
then
  mkdir -vp ${dbDir}
fi

################################################################################

# make database from nucleotide sequence
if [[ ! -f ${dbDir}${preffixGenome} ]]
then
  echo "Build BLAST database => ${preffixGenome}"
  makeblastdb -in ${genomeDir}${preffixGenome}.fasta -dbtype nucl -out ${dbDir}${preffixGenome}
fi

################################################################################

# use protein sequence as query to blast six reading frames
echo "Blasting database..."
tblastn -query ${libraryDir}${library} -db ${dbDir}${preffixGenome} -out ${outDir}${preffixGenome}.out

################################################################################
