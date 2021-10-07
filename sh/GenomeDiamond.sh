#!/bin/bash

################################################################################

species=$1
assembly=$2
assemblyDir=$3
library=$4
libraryDir=$5
outDir=$6
dbDir=${outDir}/${species}

################################################################################

# set output directory
if [[ ! -d ${dbDir} ]]
then
  mkdir -vp ${dbDir}
fi

################################################################################

# make database from syncytin protein sequence
if [[ ! -f ${libraryDir}/${library}.dmnd ]]
then
  echo "Build Diamond database => ${species}"
  diamond makedb \
    --in ${libraryDir}/${library}.fasta \
    --db ${libraryDir}/${library}.dmnd
fi

################################################################################

# use assembly as query with six reading frames
echo "Searching diamond database..."
gzip --decompress --stdout ${assemblyDir}/${assembly}.fasta.gz | \
  diamond blastx \
    --db ${libraryDir}/${library}.dmnd \
    --query - \
    --frameshift 15 \
    --block-size 5 \
    --index-chunks 1 \
    --out ${dbDir}/${species}.tsv

################################################################################

