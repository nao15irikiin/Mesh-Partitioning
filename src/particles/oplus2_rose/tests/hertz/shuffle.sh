#!/bin/bash

if [ -z "${DATADIR}" ]; then
  echo "ERROR: set DATADIR env variable first."
  exit 1
fi

if [ -z "${PROGRAM}" ]; then
  PROGRAM=./op2_cpu_timer
fi
NUMITER=1000

if [ -z "${OUTPUT}" ]; then
  OUTPUT=shuffle.data
fi
rm -f ${OUTPUT}

for seed in {0..1000}; do
  if [ -z "${CUDA_PROFILE}" ]; then
    export CUDA_PROFILE_LOG="cuda_profile.shuffle.$seed.log"
    rm -f ${CUDA_PROFILE_LOG}
  fi
  if [ ${seed} == 0 ]; then
    ${PROGRAM} ${DATADIR}/10000.step -n ${NUMITER} -s ${seed} -v 2> /dev/null | tee -a ${OUTPUT}
  else
    ${PROGRAM} ${DATADIR}/10000.step -n ${NUMITER} -s ${seed}    2> /dev/null | tee -a ${OUTPUT}
  fi
done
