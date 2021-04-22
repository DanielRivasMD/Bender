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
if [[ ! -f ${libraryDir}/${library} ]]
then
  echo "Build BLAST database => ${species}"
  diamond makedb -in ${libraryDir}/${library}.fasta -d ${libraryDir}/${library}
fi

################################################################################

# use assembly as query with six reading frames
echo "Blasting database..."
gzip -dc ${assemblyDir}/${assembly} | diamond blastx -d ${libraryDir}/${library} -q - -o ${dbDir}/${species}.tsv

################################################################################

