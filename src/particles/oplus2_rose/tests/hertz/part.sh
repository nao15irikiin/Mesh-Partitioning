#!/bin/bash

if [ -z "${DATADIR}" ]; then
  echo "ERROR: set DATADIR env variable first."
  exit 1
fi
PARTDIR=${DATADIR}/partitions

if [ -z "${PROGRAM}" ]; then
  PROGRAM=./op2_cpu_timer
fi
NUMITER=1000

OUTPUT=op2_metis.data
rm -f ${OUTPUT}
FROM=128
TO=512
STEP=32
for p in `seq ${FROM} ${STEP} ${TO}` 240 272; do
  echo -n "${p}, " | tee -a ${OUTPUT};
  ${PROGRAM} ${DATADIR}/10000.step -n ${NUMITER} -p ${PARTDIR}/metis/10000.metis.part.${p} 2> /dev/null | tee -a ${OUTPUT};
done

echo

OUTPUT=op2_rcb.data
rm -f ${OUTPUT}
for p in 5 6 7 8 9; do
  echo -n `head -1 ${PARTDIR}/rcb/10000.part.${p}`", ";
  ${PROGRAM} ${DATADIR}/10000.step -n ${NUMITER} -p ${PARTDIR}/rcb/10000.part.${p}         2> /dev/null | tee -a ${OUTPUT};
done
