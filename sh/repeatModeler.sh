#!/bin/bash

################################################################################

refFile=$1
refDir=$2
verbose=$3
storageDir=$4

################################################################################

# load modules
module load bioinfo-tools
module load RepeatModeler/2.0.1

################################################################################

# make temporary directory
if [[ "$SNIC_TMP" != "/scratch" ]]
then
  tmpDir="${SNIC_TMP}/${SLURMD_NODENAME}_${SLURM_JOBID}/"
  mkdir ${tmpDir}
fi
[[ $verbose == "true" ]] && echo "directory created" || echo -n ""

cp ${refDir}${refFile}.fa ${tmpDir}
[[ $verbose == "true" ]] && echo "file copied" || echo -n ""

################################################################################

# change directory
cd ${tmpDir}

# run repeatModeler
BuildDatabase -name ${refFile}.DB -engine ncbi ${refFile}.fa
[[ $verbose == "true" ]] && echo "database created" || echo -n ""

RepeatModeler -database ${refFile}.DB -engine ncbi -pa 32
[[ $verbose == "true" ]] && echo "library created" || echo -n ""

# storage directory
if [[ ! -d ${storageDir} ]];
then
  mkdir -pv ${storageDir}
fi

################################################################################

# retrieving to storage
mv ${tmpDir}/* ${storageDir}
[[ $verbose == "true" ]] && echo "moved" || echo -n ""

################################################################################
