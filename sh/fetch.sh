#!/bin/bash

################################################################################

inDir=$1
outDir=$2
logFile=$3
maxIt=$4
verbose=$5

################################################################################

# load modules
module load bioinfo-tools
module load sratools/2.9.6-1

################################################################################

# control log
FQDownloadLog="${outDir}${logFile}_${SLURMD_NODENAME}_${SLURM_JOBID}.log"

# initiate errorArray
while read ID
do
  (( idCount++ ))
  errorArray[$idCount]=$ID
done < ${inDir}${logFile}.txt

################################################################################

# errorArray iterartion
while [[ ${#errorArray[@]} != 0 ]] && (( itCount < maxIt ))
do
  # control iterartions as maxIt
  (( itCount++ ))
  # echo -e "${report_break}\tIteration count: ${itCount}\n\tArray length: ${#errorArray[@]}${report_break}" >> $OUT_FILE

  # itirate using the present indexes rather that the length of the array
  for ix in "${!errorArray[@]}"
  do
    # fetch sra file
    fetchCheck=$( prefetch --max-size 35G ${errorArray[$ix]} | tee -a $OUT_FILE 2>> $ERR_FILE )

    if [[ "$fetchCheck" =~ "failed to download" ]]
    then
        echo "${errorArray[$ix]} Fail at ${itCount}" >> $FQDownloadLog
    elif [[ "$fetchCheck" =~ "downloaded successfully" ]] || [[ $fetchCheck =~ "found locally" ]]
    then
      checkCheck=$( vdb-validate ${errorArray[$ix]} 2>&1 | tee -a $OUT_FILE 2>> $ERR_FILE )

      # check download
      if [[ "$checkCheck" =~ "consistent" ]]
      then
        echo "${errorArray[$ix]} Pass at ${itCount}" >> $FQDownloadLog
        unset 'errorArray[$ix]'
      fi
    fi
  done
done

################################################################################
