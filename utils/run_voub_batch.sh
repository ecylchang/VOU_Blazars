#!/bin/bash
set -u

# This script run 'vou-blazars' in batch mode.
# ('batch' means a series of jobs, non-interactively)
#
# The script will read a list of coordinates and labels, 
# and will use a default/common radius value.
#
# The file with [labels, ra, dec] columns is given as argument.
# The columns of the file are space-separated.
# The first (header) line of the file is ignored.
#
# An example of such input file follows:
# (notice the first line with column names)
#
#==============================================
#
#   RUNID RA DEC
#   134_332_0 332.82070 -0.3791089
#   149_27_1 27.176890 1.4641790
#   153_349_-1 349.90720 -1.2616670
#   15_40_0 40.799470 0.7736686
#   177_50_-1 50.896390 -1.1488420
#   193_328_0 328.25590 -0.6672633
#   19_46_-1 46.046900 -1.0988520
#   205_5_0 5.4733280 0.1109831
#   218_341_0 341.10090 -0.1489777
#   246_358_0 358.99800 -0.3452113
#   280_30_0 30.340980 0.5499375
#   380_313_0 313.87080 -0.3678985
#   501_303_0 303.59210 -0.7525020
#   506_335_1 335.81130 1.0877520
#   597_322_0 322.31080 0.6661310
#   69_332_0 332.77840 -0.1591273
#
#==============================================

# Input file (label, ra, dec):
#
RUNID_RA_DEC_FILE="$1"

# Field radius (in arcmin)
#
RADIUS='12'

# Variable $VOUB_AUTOSED bypass vou-blazars second phase user-interaction
#
export VOUB_AUTOSED='yes'


OLDIFS="$IFS" && IFS=$'\n' RUNS=($(tail -n +2 $RUNID_RA_DEC_FILE)) && IFS="$OLDIFS"

for LINE in "${RUNS[@]}"
do
    read -a FIELDS <<< $LINE

    export RUN_LABEL="${FIELDS[0]}"
    LOGFILE="RUN_${RUN_LABEL}.log"

    echo "RUNNING FIELD ${FIELDS[0]}" | tee $LOGFILE
    ./bin/vou-blazars ${FIELDS[1]} ${FIELDS[2]} $RADIUS  2>&1 | tee -a $LOGFILE
    mv $LOGFILE Results/$RUN_LABEL/.
done

