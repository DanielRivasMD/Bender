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

# make database from nucleotide sequence
if [[ ! -f ${outDir}/${species} ]]
then
  echo "Build BLAST database => ${species}"
  gzip -dc ${assemblyDir}/${assemblyT} | \
    makeblastdb \
      -in - \
      -dbtype nucl \
      -out ${outDir}/${species} \
      -title ${species}
fi

################################################################################

# use protein sequence as query to blast six reading frames
echo "Blasting database..."
tblastn \
  -query ${libraryDir}/${libraryT} \
  -db ${outDir}/${species} \
  -out ${outDir}/${species}.out

################################################################################
