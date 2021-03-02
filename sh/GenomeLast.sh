#!/bin/bash

################################################################################

genome=$1
genomeDir=$2
library=$3
libraryDir=$4
outDir=$5
dbDir=${outDir}/${genome}

################################################################################

# Improved DNA-versus-Protein Homology Search for Protein Fossils
# Yin Yao and Martin C. Frith

# makeblastdb -in RepeatPeps.lib -dbtype prot -out DB
# blastx -query chr21.fa -db DB -evalue 0.1 -outfmt 7 > out

# lastdb -q -c -R01 myDB RepeatPeps.lib
# last-train --codon -X1 myDB chr21.fa > train.out
# lastal -p train.out -D1e9 -m100 -K1 myDB chr21.fa > out

# Option -q appends * to each protein; -R01 lowercases simple sequence with tantan; -c requests masking of lowercase; -X1 treats matches to unknown residues (which are frequent in these proteins) as neutral instead of disfavored; -D1e9 sets the significance threshold to 1 random hit per 109 basepairs; -m100 sets m = 100; -K1 omits alignments whose DNA range lies in that of a higher-scoring alignment.

################################################################################

# set output directory
if [[ ! -d ${dbDir} ]]
then
  mkdir -vp ${dbDir}
fi

################################################################################

# TODO: add comments
# if [[ ! -f ${dbDir}/${genome} ]]
# then
  echo "Build LAST database => ${genome}"
  lastdb -q -c -R01 ${dbDir}/${genome} ${genomeDir}/${genome}.fasta
# fi

################################################################################

echo "Training database..."
last-train --codon -X1 ${dbDir}/${genome} ${libraryDir}/${library} > ${dbDir}/${genome}.train.out

echo "Lasting database..."
lastal -p ${dbDir}/${genome}.train.out -D1e9 -m100 -K1 ${dbDir}/${genome} ${libraryDir}/${library} > ${dbDir}/${genome}.out

################################################################################

# # purge train file
# if [[ -f ${dbDir}/${genome}.filtered.out ]]
# then
#   echo "Cleaning filtered..."
#   ${dbDir}/${genome}.filtered.out
# fi

################################################################################
