#!/bin/bash
set -u

#RUNID_RA_DEC_FILE="refs/RADEC_image_centers.txt"
RUNID_RA_DEC_FILE="$1"

OLDIFS="$IFS" && IFS=$'\n' RUNS=($(tail -n +2 $RUNID_RA_DEC_FILE)) && IFS="$OLDIFS"

export VOUB_AUTOSED='yes'

for LINE in "${RUNS[@]}"
do
    read -a FIELDS <<< $LINE

    export RUN_LABEL="${FIELDS[0]}"
    LOGFILE="RUN_${RUN_LABEL}.log"

    echo "RUNNING FIELD ${FIELDS[0]}" | tee $LOGFILE
    ./bin/vou-blazars ${FIELDS[1]} ${FIELDS[2]} 12  2>&1 | tee -a $LOGFILE
    mv $LOGFILE Results/$RUN_LABEL/.
done

