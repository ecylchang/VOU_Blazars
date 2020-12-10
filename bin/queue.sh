#!/usr/bin/env bash
set -u

help(){
  echo ""
  echo " Usage: $(basename $0) [-n ] -f <file-commands>"
  echo ""
  echo " Options:"
  echo " -h : this help message"
  echo " -n : number of jobs to run simultaneously [default: $NPROCS]"
  echo " -f : file with the list of commands to run by the queue"
  echo " -q : quiet run"
  echo ""
}
CURDIR=$(cd `dirname $BASH_SOURCE`; pwd)

VERBOSE="1"

# Number of processes
#
NPROCS=1

# Amount of time to refresh the queue (in seconds)
#
SLEEP=0.5

TMPDIR="${PWD}/tmp"
[ -d $TMPDIR ] || mkdir -p $TMPDIR

# Empty init variable
LIST=''

# GetOptions..
#
while getopts ":hqn:f:d:" opt
do
  case $opt in
    h) help; exit 0;;
    q) VERBOSE="0";;
    n) NPROCS="$OPTARG";;
    f) LIST="$OPTARG";;
    \?) echo "ERROR: Wrong option $OPTARG ";;
    :) echo "ERROR: Missing value for $OPTARG ";;
  esac
done


if [ ! -f "$LIST" -o -z "$LIST" ]
then
  echo "ERROR: file '$LIST' does not exist. Finishing."
  help
  exit 2
fi



#===========================
# Check PID..
#
check_process(){
  kill -0 $1 2>/dev/null
  echo $?
}
#===========================

# =====================================================================
# Read pipeline file names of each selected Halo..
#
CNT=0
LINES=()
# for LI in `cat $LIST`
# do
#   CNT=$((CNT+1))
#   LINES[$CNT]="${LI[*]}"
# done
while IFS= read -r LINE
do
  [[ -z "$LINE" ]] && continue
  CNT=$((CNT+1))
  LINES[$CNT]="$LINE"
done < $LIST

echo
#if $CNT
[ "$VERBOSE" = "1" ] && echo "Searching $CNT catalogs"

# Check if any halo was found. If not, finish the run..
#
[ "$CNT" -eq "0" ] && { echo "Empty list of observation?"; exit; }

# Start the queue of jobs..
#
PIDs=()
CNTs=()
NJOBS=0

declare -i nncc=0
while [ "$NJOBS" -le "$CNT" ]
do
  # IF queue is not full && we have not reached to top yet..
  if [ ${#PIDs[*]} -lt $NPROCS -a "$NJOBS" -lt "$CNT" ]
  then
    NJOBS=$((NJOBS+1))

    # each entry in '-f input_list' is expeted to be a json filename
    LINE_COMMAND=${LINES[$NJOBS]}
    1>/dev/null $LINE_COMMAND &

    PID=$!
    PIDs[$PID]=$PID
    CNTs[$PID]=$NJOBS
  else
    sleep $SLEEP
  fi

  # Check the status of each process in queue
  for PID in ${PIDs[*]};
  do
    # If process PID is finished,
    if [ "$(check_process $PID)" -ne "0" ]; then
      _cnt=${CNTs[$PID]}
      _file=${LINES[$_cnt]}
      # get the result status of PID and then remove from queue
      wait $PID
      PSTS=$?
      cats=`echo $_file | awk '{print $5}'`
      ref=`grep ^${cats}[^BGMNW] ${CURDIR}/catrefs.txt | cut -d'!' -f1 | cut -d' ' -f2- `
      
#      echo $cats $PSTS
      if [[ $PSTS -eq 0 ]]; then
         nncc=$nncc+1
         echo "( $nncc / $NJOBS ) \033[32;1m $cats : SUCCESS\033[0m -- ${ref}"
      elif [[ $PSTS -eq 10 ]]; then
         nncc=$nncc+1
         1>&2 echo "( $nncc / $NJOBS ) $cats : NO SOURCES FOUND"
      elif [ $cats == MAGIC -o $cats == VERITAS ]; then
         nncc=$nncc+1
         1>&2 echo "( $nncc / $NJOBS ) $cats : NO SOURCES FOUND"
      else
         nncc=$nncc+1
         1>&2 echo "( $nncc / $NJOBS ) \033[31;1m  $cats : SEARCH FAILED\033[0m "
         echo $_file >> tmp/voerror.txt
      fi
      unset PIDs[$PID]
      unset CNTs[$PID]
    fi
  done

  # BREAK when we reach the top and there are no jobs running
  [ "$NJOBS" -eq "$CNT" -a "${#PIDs[*]}" -eq "0" ] && break

done
