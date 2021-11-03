#!/bin/bash
set -euo pipefail

################################################################################

species=$1
assemblyT=$2
inDir=$3
libraryT=$4
libraryDir=$5
outDir=$6

################################################################################

# make database from syncytin protein sequence
if [[ ! -f ${libraryDir}/${libraryT}.dmnd ]]
then
  echo "Build Diamond database => ${species}"
  diamond makedb \
    --in ${libraryDir}/${libraryT}.fasta \
    --db ${libraryDir}/${libraryT}.dmnd
fi

################################################################################

# use assembly as query with six reading frames
echo "Searching diamond database..."
diamond blastx \
  --db ${libraryDir}/${libraryT}.dmnd \
  --query ${inDir}/${assemblyT}.fasta \
  --frameshift 15 \
  --block-size 6 \
  --index-chunks 1 \
  --out ${outDir}/${species}.tsv

################################################################################
