RUNID_RA_DEC_FILE="refs/RADEC_image_centers.txt"

OLDIFS="$IFS" && IFS=$'\n' RUNS=($(tail -n +2 $RUNID_RA_DEC_FILE)) && IFS="$OLDIFS"

#for LINE in "${RUNS[@]}"
#do 
#    read -a FIELDS <<< $LINE
#    echo "FIELD ${FIELDS[0]}: ${FIELDS[1]} ${FIELDS[2]}"
#done

export VOUB_AUTOSED='yes'

for LINE in "${RUNS[@]}"
do
    read -a FIELDS <<< $LINE
    LOGFILE="RUN_${FIELDS[0]}.log"
    echo "RUNNING FIELD ${FIELDS[0]}" | tee $LOGFILE
    ./bin/vou-blazars ${FIELDS[1]} ${FIELDS[2]} 12  2>&1 | tee -a $LOGFILE
done

