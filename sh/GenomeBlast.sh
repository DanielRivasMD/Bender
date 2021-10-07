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

# make database from nucleotide sequence
if [[ ! -f ${dbDir}/${species} ]]
then
  echo "Build BLAST database => ${species}"
  gzip -dc ${assemblyDir}/${assembly} | \
    makeblastdb \
      -in - \
      -dbtype nucl \
      -out ${dbDir}/${species} \
      -title ${species}
fi

################################################################################

# use protein sequence as query to blast six reading frames
echo "Blasting database..."
tblastn \
  -query ${libraryDir}/${library} \
  -db ${dbDir}/${species} \
  -out ${dbDir}/${species}.out

################################################################################
