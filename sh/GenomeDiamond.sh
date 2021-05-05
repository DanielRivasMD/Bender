#!/bin/bash

################################################################################

species=$1
echo "species = ${species}"
assembly=$2
echo "assembly = ${assembly}"
assemblyDir=$3
echo "assemblyDir = ${assemblyDir}"
library=$4
echo "library = ${library}"
libraryDir=$5
echo "libraryDir = ${libraryDir}"
outDir=$6
echo "outDir = ${outDir}"
dbDir=${outDir}/${species}
echo "dbDir = ${dbDir}"


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
  diamond makedb --in ${libraryDir}/${library}.fasta --db ${libraryDir}/${library}.dmnd
fi

################################################################################

# use assembly as query with six reading frames
echo "Searching diamond database..."
gzip --decompress --stdout ${assemblyDir}/${assembly}.fasta.gz | diamond blastx --db ${libraryDir}/${library}.dmnd --query - --frameshift 15 --block-size 5 --index-chunks 1 --out ${dbDir}/${species}.tsv

################################################################################

