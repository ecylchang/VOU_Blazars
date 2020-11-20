#!/bin/bash

# Directory where the fortran binaries are:
#
HERE=$(cd `dirname $BASH_SOURCE`; pwd) ##point to bin file and fortran file
BINF="${HERE}/fort"

handle_ps () {
    [ `which ps2eps 2> /dev/null` ] || return 0
    FILE_PS=$1
    if [ $FILE_PS == sed.ps -o $FILE_PS == LC.ps -o $FILE_PS == LC_fermi.ps ]; then
       ps2eps -B -q $FILE_PS -R +
    else
       ps2eps -B -q $FILE_PS
    fi
    if [ `which open 2> /dev/null` ]; then
      open ${FILE_PS%.ps}.eps
    else
      gv ${FILE_PS%.ps}.eps &
    fi
    rm -f $FILE_PS
}

help() {
   echo ""
   echo " Usage: $(basename $0) { --ra <degrees> --dec <degrees> --FOV <arcmin> }"
   echo ""
   echo " ARGUMENTS:"
   echo "  --ra     : Right Ascension (in DEGREES)"
   echo "  --dec    : Declination (in DEGREES)"
   echo "  --FOV   : Searchin Radius (in ARC-MINUTES) around RA,DEC to search for observations"
   echo ""
   echo " OPTIONS:"
   echo " --mode    : Running mode"
   echo "       Options are 'f' find candidate mode (default): finding interesing candidates within a specified region"
   echo "                   's' SED mode: obtaining SED for a specified source with given R.A. Dec."
   echo "                   'l' Light curve mode: obtaining light curve for a specified source with given R.A. Dec."
   echo ""
   echo " --nh      : nH column density (in cm^2). Default is 5.e20 cm^2."
   echo "             (If the user has installed Heasoft and did not specify the nh, it will use the value calculated by Heasoft)"
   echo ""
   echo " First and larger error region, circle (--radius) or elliptical (--major --minor --angle)"
   echo " --radius  : Error circle radius (in ARC-MINUTES). Default is 0"
   echo " --major   : Error elliptical major axis (in ARC-MINUTES). Default is 0"
   echo " --minor   : Error elliptical minor axis (in ARC-MINUTES). Default is 0"
   echo " --angle   : Position angle of the error elliptical (in DEGREES), north-east on sky. Default is 0"
   echo ""
   echo " Second and smaller error region, circle (--radius2) or elliptical (--major2 --minor2 --angle2)"
   echo " --radius2 : Second error circle radius (in ARC-MINUTES). Default is 0"
   echo " --major2  : Second error elliptical major axis (in ARC-MINUTES). Default is 0"
   echo " --minor2  : Second error elliptical minor axis (in ARC-MINUTES). Default is 0"
   echo " --angle2  : Position angle of the second error elliptical (in DEGREES), north-east on sky. Default is 0"
   echo ""
}



##############################################################
# Initial set up, read the parameters input
##############################################################

#source activate eada2
#check if docker and heasoft are installed
#dockerthere=$(which docker 2> /dev/null)
nhthere=`which nh`

#####checking the input parameter
while [[ $# -gt 0 ]]
do
   case $1 in
      -h|--help)
         help; exit 0;;
      --ra)
         posra=$2; shift;;
      --dec)
         posdec=$2; shift;;
      --FOV)
         sfov=$2; shift;;
      --nh)
         nhval=$2; shift;;
      --radius)
         r1=$2; shift;;
      --major)
         emaj1=$2; shift;;
      --minor)
         emin1=$2; shift;;
      --angle)
         posa1=$2; shift;;
      --radius2)
         r2=$2; shift;;
      --major2)
         emaj2=$2; shift;;
      --minor2)
         emin2=$2; shift;;
      --angle2)
         posa2=$2; shift;;
      --mode)
         runmode=$2; shift;;
      --*)
         echo "$0: error - unrecognized option $1" 1>&2
         help;exit 1;;
      -?)
         echo "$0: error - unrecognized option $1" 1>&2
         help;exit 1;;
      *)
         break;;
   esac
   shift
done



#read the parameters input
echo Running VOU-Blazars V1.92
echo

#read the XRT Deepsky name and create the file to store the results
raxrt=`echo $posra | sed 's/\./_/g'`
decxrt=`echo $posdec | sed 's/\./_/g' | sed 's/-/m/g'`
rrxrt=`echo $sfov | sed 's/\./_/g'`

# If a variable 'RUN_LABEL' is defined in the environment,
# use it to make the Result/. output directory
#
if [ ! -z "$RUN_LABEL" ]; then
   xrtnm="$RUN_LABEL"
else
   xrtnm=${raxrt}"_"${decxrt}"_"${rrxrt}
fi
if [ -z $pid ]; then
   unset pidnm
   test ! -d Results && mkdir Results
   test ! -d Results/$xrtnm && mkdir Results/$xrtnm
   test ! -d Results/SEDtool && mkdir Results/SEDtool
else
   pidnm=$pid"_"
fi

#adjust the input ra dec from user
nht1=`echo $posra | grep '\.'`
if [ -z $nht1 ]; then
   ranh=$posra".0"
else
   nht1=`echo $posra | cut -d'.' -f2`
   if [ -z $nht1 ]; then
      ranh=$posra"0"
   else
      ranh=$posra
   fi
fi
nht2=`echo $posdec | grep '\.'`
if [ -z $nht2 ]; then
   decnh=$posdec".0"
else
   nht2=`echo $posdec | cut -d'.' -f2`
   if [ -z $nht2 ]; then
      decnh=$posdec"0"
   else
      decnh=$posdec
   fi
fi

#read the parameters input from configuration file
#read the nh, or set the default nh value
if [ -z $nhval ]; then
   if [ -z $nhthere ]; then
      nhval=`grep 'NHVAL' ${HERE}/config_vou.txt | awk '{print $2}'`
   else
#   echo $ranh $decnh
      nh equinox=2000 ra=$ranh dec=$decnh > ${pidnm}nhvalue.txt
      nhval=`tail -1 ${pidnm}nhvalue.txt |  awk '{print $7}'`
      rm -f ${pidnm}nhvalue.txt
   fi
fi
[ -z $r1 ] && r1=`grep 'RADIUS1' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $emaj1 ] && emaj1=`grep 'MAJOR1' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $emin1 ] && emin1=`grep 'MINOR1' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $posa1 ] && posa1=`grep 'POSANG1' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $r2 ] && r2=`grep 'RADIUS2' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $emaj2 ] && emaj2=`grep 'MAJOR2' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $emin2 ] && emin2=`grep 'MINOR2' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $posa2 ] && posa2=`grep 'POSANG2' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $runmode ] && runmode=`grep 'MODE' ${HERE}/config_vou.txt | awk '{print $2}'`
if [ $runmode != f ]; then
   sfov=1.
   r1=0
   emaj1=0
   emin1=0
   posa1=0
   r2=0
   emaj2=0
   emin2=0
   posa2=0
fi

#set $ranh $decnh $sfov $nhval $r1 $emaj1 $emin1 $posa1 $r2 $emaj2 $emin2 $posa2 $runmode



#check if the SED mode need plotting or not
if [ $runmode == -s ]; then #last parameter
   runmode=s
   plotsed=N
else
   plotsed=Y
fi



echo $ranh $decnh $sfov $nhval $r1 $emaj1 $emin1 $posa1 $r2 $emaj2 $emin2 $posa2 $runmode

#analyzing the XRT Deepsky data
#if [ -n $dockerthere ]; then
#docker run --name caldb chbrandt/heasoft_caldb:swift
#alias swift_deepsky='docker run --rm -it --volumes-from caldb -v $PWD/work:/work chbrandt/swift_deepsky:latest'
#swift_deepsky --ra $1 --dec $2 --radius $3
#echo finish searching swift deep data
#fi





##############################################################
# FIRST PHASE
# aim: find the candidates from available radio and X-ray sources
##############################################################

#FIRST phase data retrieving, runing the conesearch 1st phase
rm -f ${pidnm}*.1.csv
echo conesearch --db ${HERE}/cats1.ini --catalog NVSS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}nvss.1.csv > ${pidnm}vosearch.txt
echo conesearch --db ${HERE}/cats1.ini --catalog FIRST --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}first.1.csv >> ${pidnm}vosearch.txt
echo conesearch --db ${HERE}/cats1.ini --catalog SUMSS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}sumss.1.csv >> ${pidnm}vosearch.txt
echo conesearch --db ${HERE}/cats1.ini --catalog VLASSQL --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}vlassql.1.csv >> ${pidnm}vosearch.txt
echo conesearch --db ${HERE}/cats1.ini --catalog 2SXPS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}2sxps.1.csv >> ${pidnm}vosearch.txt
echo conesearch --db ${HERE}/cats1.ini --catalog SDS82 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}sds82.1.csv >> ${pidnm}vosearch.txt
echo conesearch --db ${HERE}/cats1.ini --catalog 1OUSX --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}1ousx.1.csv >> ${pidnm}vosearch.txt
echo conesearch --db ${HERE}/cats1.ini --catalog RASS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}rass.1.csv >> ${pidnm}vosearch.txt
echo conesearch --db ${HERE}/cats1.ini --catalog XMMSL2 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}xmmsl2.1.csv >> ${pidnm}vosearch.txt
echo conesearch --db ${HERE}/cats1.ini --catalog 4XMM-DR9 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}4xmmdr9.1.csv >> ${pidnm}vosearch.txt
echo conesearch --db ${HERE}/cats1.ini --catalog BMW --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}bmw.1.csv >> ${pidnm}vosearch.txt
echo conesearch --db ${HERE}/cats1.ini --catalog WGACAT --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}wgacat.1.csv >> ${pidnm}vosearch.txt
echo conesearch --db ${HERE}/cats1.ini --catalog IPC2E --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}ipc.1.csv >> ${pidnm}vosearch.txt
echo conesearch --db ${HERE}/cats1.ini --catalog IPCSL --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}ipcsl.1.csv >> ${pidnm}vosearch.txt
echo conesearch --db ${HERE}/cats1.ini --catalog Chandra-CSC2 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}chandracsc2.1.csv >> ${pidnm}vosearch.txt
#echo conesearch --db ${HERE}/cats1.ini --catalog MAXISSC --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}maxissc.1.csv >> ${pidnm}vosearch.txt
#echo conesearch --db ${HERE}/cats1.ini --catalog MAXIGSC --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}maxigsc.1.csv >> ${pidnm}vosearch.txt
if [ $runmode == f ]; then #skip the catalogs that we don't plot its data on SED
   echo conesearch --db ${HERE}/cats1.ini --catalog ZWCLUSTERS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}zw.1.csv >> ${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats1.ini --catalog PSZ2 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}psz2.1.csv >> ${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats1.ini --catalog ABELL --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}abell.1.csv >> ${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats1.ini --catalog MCXC --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}mcxc.1.csv >> ${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats1.ini --catalog 5BZCat --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}5bzcat.1.csv >> ${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats1.ini --catalog SDSSWHL --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}whl.1.csv >> ${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats1.ini --catalog SWXCS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}swxcs.1.csv >> ${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats1.ini --catalog 3HSP --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}3hsp.1.csv >> ${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats2.ini --catalog CRATES --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}crates.1.csv >> ${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats1.ini --catalog BROS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}bros.1.csv >> ${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats1.ini --catalog PULSAR --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}pulsar.1.csv >> ${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats1.ini --catalog F2PSR --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}f2psr.1.csv >> ${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats2.ini --catalog 4FGL-DR2 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}4fgldr2.1.csv >> ${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats2.ini --catalog 3FHL --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}3fhl.1.csv >> ${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats2.ini --catalog 3FGL --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}3fgl.1.csv >> ${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats2.ini --catalog 2BIGB --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}2bigb.1.csv >> ${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats2.ini --catalog MST9Y --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}mst9y.1.csv >> ${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats2.ini --catalog FermiMeV --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}fmev.1.csv >> ${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats2.ini --catalog 2AGILE --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}2agile.1.csv >> ${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats1.ini --catalog FermiGRB --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}fgrb.1.csv >> ${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats1.ini --catalog MilliQuas --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}mquas.1.csv >> ${pidnm}vosearch.txt
#   echo conesearch --db ${HERE}/cats1.ini --catalog 4LAC --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --columns default -o ${pidnm}4lac.1.csv >> ${pidnm}vosearch.txt
fi

#check if the phase 1 data is already there
ncat1=`cat ${pidnm}vosearch.txt | wc -l`
rm -f voerror.txt
sh ${HERE}/queue.sh -f ${pidnm}vosearch.txt -n $ncat1 #2> voerror #!2>&1
#   cat voerror | grep 'ERROR'
#see if there is an error when searching through VO
declare -i irunvo=1
if [ -s voerror.txt ]; then
   voerror=yes
else
   voerror=no
fi
#echo $voerror
until [ $voerror == no -o $irunvo -ge 3 ]
do
   runvo=y
#   [ -z $pid ] && read -p "Found error when searhcing through VO. Type 'n' to skip and continue; otherwise, rerun the VO conesearch" runvo
   [ -z $runvo ] && runvo=y
   cat voerror.txt > vosearch.txt
   rm -f voerror.txt
   irunvo=${irunvo}+1
   echo
   echo VO search returned errors for one or more catalogues. Run conesearch again $irunvo
   [ $runvo != n -a $runvo != N ] && sh ${HERE}/queue.sh -f ${pidnm}vosearch.txt -n $ncat1
   if [ -s voerror.txt ]; then
      voerror=yes
   else
      voerror=no
   fi
done
echo
echo finish conesearch phase 1
rm -f ${pidnm}vosearch.txt

#read the data
ls ${pidnm}*.1.csv > ${pidnm}catlist1.txt
#copy the XRTDEEP result and add to the catlist1
rm -f ${pidnm}xrtdeep.csv
if [ -s work/$xrtnm/table_flux_detections.csv ]; then
   cat work/$xrtnm/table_flux_detections.csv | sed 's/;/,/g' | sed 's/:/ /g' > ${pidnm}xrtdeep.csv
fi
[ -s ${pidnm}xrtdeep.csv ] && echo ${pidnm}xrtdeep.csv >> ${pidnm}catlist1.txt
rm -f ${pidnm}output1.csv
echo
[ -s ${pidnm}catlist1.txt ] && ${BINF}/readcat ${pidnm}catlist1.txt ${pidnm}output1.csv $ranh $decnh $sfov $nhval $r1 $emaj1 $emin1 $posa1

#running the PHASAE 1
rm -f ${pidnm}*temp.txt
if [ $runmode == s -o $runmode == l ]; then
   ${BINF}/find_candidates1 ${pidnm}output1 ${pidnm}find_out_temp.txt ${pidnm}RX_temp.txt ${pidnm}Sed_temp.txt ${pidnm}no_matched_temp.txt ${BINF} 0 > ${pidnm}phase1
   echo
   cat ${pidnm}phase1
else
   ${BINF}/find_candidates1 ${pidnm}output1 ${pidnm}find_out_temp.txt ${pidnm}RX_temp.txt ${pidnm}Sed_temp.txt ${pidnm}no_matched_temp.txt ${BINF} 1 > ${pidnm}phase1
   echo
   cat ${pidnm}phase1

###############################PLOTTING###############################
#plot the result, candidates map and radio-X-ray source map
   rm -f ${pidnm}candidates.*ps
   rm -f ${pidnm}RX_map.*ps
   sort -n -k 3 ${pidnm}RX_temp.txt > ${pidnm}RX_sorted.txt
   rm -f ${pidnm}RX_temp.txt
   ${BINF}/gnomo_plot_types ${pidnm}RX_sorted.txt,${pidnm}candidates_posix.txt,${pidnm}RX_map.ps/vcps, $ranh $decnh $r1 $sfov $ranh $decnh $emaj1 $emin1 $posa1 $r2 $emaj2 $emin2 $posa2 $ranh $decnh
   ${BINF}/gnomo_plot_types ${pidnm}find_out_temp.txt,${pidnm}candidates_posix.txt,${pidnm}candidates.ps/vcps, $ranh $decnh $r1 $sfov $ranh $decnh $emaj1 $emin1 $posa1 $r2 $emaj2 $emin2 $posa2 $ranh $decnh
   handle_ps ${pidnm}RX_map.ps
   handle_ps ${pidnm}candidates.ps
#   ps2eps -B -q RX_map.ps
#   open RX_map.eps
#   ps2eps -B -q candidates.ps
#  open candidates.eps
#   rm -f candidates.ps
#   rm -f RX_map.ps
##############################################################
fi





##############################################################
# Intermediate PHASE
# aim: find more available candiates from single radio/ single X-ray sources.
##############################################################

rm -f ${pidnm}Intermediate_out.txt
rm -f ${pidnm}output_int.csv
if [ -s ${pidnm}no_matched_temp.txt ]; then
#check the error circle
   if [ $r1 != 0. -a $emaj1 == 0. ]; then
      radint=$r1
   elif [ $r1 == 0. -a $emaj1 != 0. ]; then
      radint=$emaj1
   else
      radint=1.
      echo
      echo set default error circle area to 1 arcmin in Intermediate phase
   fi
   radcrit=`echo $radint | cut -d'.' -f 1`
   [ -z $radcrit ] && radcrit=1
   [ $radcrit -lt 1 ] && radcrit=1
#   echo 'radcrit' $radcrit
   if [ -z $pid ]; then
      process=y
   else
      process=n
   fi

#   rcut=0
#   declare -i radint
   if [ $radcrit -ge 20 ]; then
      echo
      echo Warning! Large searching radius in Intermediate phase for conesearch.
      [ -z $pid ] && read -p "Type 'y' to continue! Otherwise, stop processing Intermediate phase." process #rcut
      [ -z $process ] && process=n
      [ $process != y -a $process != Y ] && echo Not searching sources in Intermediate phase.
   fi

#running conesearch for Intermediate phase
   if [ $process == y -o $process == Y ]; then
      rm -f *.i.csv
      ln=`cat ${pidnm}no_matched_temp.txt | wc -l`
      declare -i cc
      cc=0
      for (( jj=1; jj<=$ln; jj=jj+1 ))
      do
         raint=`head -$jj ${pidnm}no_matched_temp.txt | tail -1 |  awk '{print $5}'`
         decint=`head -$jj ${pidnm}no_matched_temp.txt | tail -1 |  awk '{print $6}'`
         typeint=`head -$jj ${pidnm}no_matched_temp.txt | tail -1 |  awk '{print $10}'`
         #declare -i typeint
         if [ $typeint -gt 50 -o $typeint -lt 0 ]; then
            cc=$cc+1
         fi
      done
      echo conesearch --db ${HERE}/cats2.ini --catalog GALEX --ra $ranh --dec $decnh --radius $radint --runit arcmin --columns default -o ${pidnm}galex.i.csv > ${pidnm}vosearch.txt
#      echo conesearch --db ${HERE}/cats2.ini --catalog Fermi8YL --ra $ranh --dec $decnh --radius $radint --runit arcmin --columns default -o fermi8yr.i.csv >> vosearch.txt
      echo conesearch --db ${HERE}/cats2.ini --catalog PMN --ra $ranh --dec $decnh --radius $radint --runit arcmin --columns default -o ${pidnm}pmn.i.csv >> ${pidnm}vosearch.txt
      echo conesearch --db ${HERE}/cats2.ini --catalog GB6 --ra $ranh --dec $decnh --radius $radint --runit arcmin --columns default -o ${pidnm}gb6.i.csv >> ${pidnm}vosearch.txt
      echo Searching further data in intermediate phase for $cc source
      ncatint=`cat ${pidnm}vosearch.txt | wc -l`
      rm -f voerror.txt
      sh ${HERE}/queue.sh -f ${pidnm}vosearch.txt -n $ncatint #1> voerror 2>&1
#      cat voerror | grep 'ERROR'

# check if there is error during the searching
      declare -i irunvo=1
      if [ -s voerror.txt ]; then
         voerror=yes
      else
         voerror=no
      fi
#      echo $voerror
      until [ $voerror == no -o $irunvo -ge 3 ]
      do
         runvo=y
      #   [ -z $pid ] && read -p "Found error when searhcing through VO. Type 'n' to skip and continue; otherwise, rerun the VO conesearch" runvo
         [ -z $runvo ] && runvo=y
         cat voerror.txt > vosearch.txt
         rm -f voerror.txt
         irunvo=${irunvo}+1
         echo
         echo VO search returned errors for one or more catalogues. Run conesearch again $irunvo
         [ $runvo != n -a $runvo != N ] && sh ${HERE}/queue.sh -f ${pidnm}vosearch.txt -n $ncat1
         if [ -s voerror.txt ]; then
            voerror=yes
         else
            voerror=no
         fi
      done
      echo
      echo finish conesearch intermediate phase
      rm -f ${pidnm}vosearch.txt

#read the Intermediate phase data
      ls ${pidnm}*.i.csv > ${pidnm}catlist_int.txt
      rm -f ${pidnm}output_int.csv
      echo
      [ -s ${pidnm}catlist_int.txt ] && ${BINF}/readcat ${pidnm}catlist_int.txt ${pidnm}output_int.csv $ranh $decnh $sfov $nhval $r1 $emaj1 $emin1 $posa1

#run the Intermediate phase
      echo
      ${BINF}/find_candidates_int ${pidnm}output_int ${pidnm}no_matched_temp.txt ${pidnm}find_out_temp.txt ${pidnm}Intermediate_out.txt
   fi
fi

#transfer the data from Intermediate phase to 2 phase
if [ -s ${pidnm}Intermediate_out.txt ]; then
   nint=`tail -1 ${pidnm}Intermediate_out.txt | awk '{print $1}' `
#   declare -i nint
   head -$nint ${pidnm}Intermediate_out.txt >> ${pidnm}find_out_temp.txt
   nnllint=`cat ${pidnm}Intermediate_out.txt | wc -l`
   declare -i nsedl=$nnllint-$nint
   declare -i nsedu=$nnllint-$nint-1
   tail -$nsedl ${pidnm}Intermediate_out.txt | head -$nsedu >> ${pidnm}Sed_temp.txt

###############################PLOTTING###############################
#plot the candidates again with Intermediate phase finish ???
   rm -f ${pidnm}candidates.*ps
   ${BINF}/gnomo_plot_types ${pidnm}find_out_temp.txt,${pidnm}candidates_posix.txt,${pidnm}candidates.ps/vcps, $ranh $decnh $r1 $sfov $ranh $decnh $emaj1 $emin1 $posa1 $r2 $emaj2 $emin2 $posa2 $ranh $decnh
   handle_ps ${pidnm}candidates.ps
#    ps2eps -B -q candidates.ps
#    open candidates.eps
#    rm -f candidates.ps
######################################################################
fi

#save the phase 1 and phase intermediate results in the folder
[ $runmode == f -a -d Results/$xrtnm ] && cp candidates.*ps Results/$xrtnm/.
[ $runmode == f -a -d Results/$xrtnm ] && cp RX_map.*ps Results/$xrtnm/.
[ -f output1.csv -a -d Results/$xrtnm ] && cp output1.csv Results/$xrtnm/output1.csv
[ -f output_int.csv -a -d Results/$xrtnm ] && cp output_int.csv Results/$xrtnm/output_int.csv

[ -d Results/$xrtnm ] && echo \\section{$ranh $decnh} > texcand.txt
if [ $runmode == f -a -d Results/$xrtnm ]; then
   echo \\begin{figure}[h!] >> texcand.txt
   echo \\centering >> texcand.txt
   echo \\includegraphics[width=0.49\\linewidth]{Results/$xrtnm/candidates.eps} >> texcand.txt
   echo \\includegraphics[width=0.49\\linewidth]{Results/$xrtnm/RX_map.eps} >> texcand.txt
   echo \\caption{Left: Candidates Map. Right: Radio-X-ray map.} >> texcand.txt
   echo \\end{figure} >> texcand.txt
fi





#######################################################
# Second PHASE
# aim: Plot the SED and error circle map for a specified source
#######################################################

#begin phase 2 setup, leave the tool when no candidates
rm -f ${pidnm}*.2.csv
rm -f texsed*
source=1

if [ ! -s ${pidnm}find_out_temp.txt ]; then
   source=nocand
   ln=0
fi

#read the number interested
until [ $source == sed -o $source == q -o $source == lcurve ]
do
   if [ $runmode == f -a $source != nocand ]; then
      echo
      if [ -z $VOUB_AUTOSED ]; then
#        cat ${pidnm}phase1
        if [ -z $pid ]; then
           read -p "Please enter the source number and zoom in area (in arcsec) intered: (Type 'q' to quit) (Type 'candlist' to show all candidates)" source zzinput zoomin showsed
        else
           source=q
        fi
      else
        runiter=`ls ${pidnm}*.$source.2.csv 2>/dev/null | wc -l`
        if [ $VOUB_AUTOSED == all -a ${runiter} -eq 0 ]; then
           [ -f ${pidnm}Sed_temp.txt ] && `grep "matched source" ${pidnm}Sed_temp.txt | awk '{print $1}' | uniq > ${pidnm}sources.list`
#           cat sources.list
        elif [ $VOUB_AUTOSED == part -a ${runiter} -eq 0 ]; then
           cat phase1 | tail -2 | head -1 |awk '{print $7 "\n" $8 "\n" $9 "\n" $10 "\n" $11}' > sources.list
        fi
        #cat sources.list
        #[ `wc -l ${pidnm}sources.list | cut -d' ' -f1` == 0 ] && { source='q'; continue; }
        declare -i sourcenb
        sourcenb=`cat ${pidnm}sources.list | wc -l`
        sourcenb=${sourcenb}-1
        if [ $sourcenb -ge 0 ]; then
           source=`head -1 ${pidnm}sources.list`
        #echo $sourcenb
           cat ${pidnm}sources.list | tail -${sourcenb} > ${pidnm}sources.tmp
           mv ${pidnm}sources.tmp ${pidnm}sources.list
        else
           source=q
        fi
      #elif [ $VOUB_AUTOSED=part ]; then
      fi
      ln=`cat ${pidnm}find_out_temp.txt | wc -l`
   elif [ $runmode == s -o $runmode == l ]; then
      ln=1
      [ $runmode == s ] && source=sed
      [ $runmode == l ] && source=lcurve
#      [ ! -f ${pidnm}find_out_temp.txt ] && echo $1 $decnh 99 > ${pidnm}find_out_temp.txt
      declare -i zoomintt #deal with non interger input
      declare -i zoomindd
      declare -i zoomin
      nht3=`echo $sfov | grep '\.'`
      if [ $nht3 ]; then
         nht31=`echo $sfov | cut -d'.' -f1`
         nht32=`echo $sfov | cut -d'.' -f2`
         if [ -z $nht32 ]; then
            zoomintt=$nht31
            zoomindd=0
         else
            zoomintt=${nht31}
            zoomindd=${nht32:0:1}
            zoomindd=${zoomindd}*6
            if [ ${nht32:1:2} ]; then
               if [ ${nht32:1:2} -lt 3 ]; then
                  zoomindd=${zoomindd}+1
               elif [ ${nht32:1:2} -lt 5 ]; then
                  zoomindd=${zoomindd}+2
               elif [ ${nht32:1:2} -eq 5 ]; then
                  zoomindd=${zoomindd}+3
               elif [ ${nht32:1:2} -lt 8 ]; then
                  zoomindd=${zoomindd}+4
               else
                  zoomindd=${zoomindd}+5
               fi
            fi
         fi
      else
         zoomintt=$sfov
         zoomindd=0
      fi
      zoomin=${zoomintt}*60+${zoomindd}
      ####echo testtesttest $source $ln
   fi

#check if need to output the SED file or not and deal with no input
   [ -z $source ] && source=0 #no input
   [ -z $showsed ] && showsed=nosed
   if [ $zoomin ]; then
       if [ $zoomin == osed ]; then
         showsed=osed
         zzinput1=`echo $zzinput | cut -d'.' -f1`
         if [ $zzinput1 -gt 5 ]; then
            zoomin=$zzinput
            zzinput=0
         else
            zoomin=60
         fi
      fi
   else
      if [ $zzinput ]; then
         if [ $zzinput == osed ]; then
            showsed=osed
            zzinput=0
         else
            zzinput1=`echo $zzinput | cut -d'.' -f1`
            if [ $zzinput1 -gt 5 ]; then
               showsed=nosed
               zoomin=$zzinput
               zzinput=0
            fi
         fi
      fi
   fi
#   echo snr=$source redshift=$zzinput area=$zoomin sed=$showsed

#read the find_out list and the second phase ra dec
   declare -i nn
   nn=0
   [ -z $ln ] && ln=0
   for (( ii=1; ii<=$ln; ii=ii+1 ))
   do
      if [ $runmode == f ]; then
         rar=`head -$ii ${pidnm}find_out_temp.txt | tail -1 |  awk '{print $1}'`
         decr=`head -$ii ${pidnm}find_out_temp.txt | tail -1 |  awk '{print $2}'`
         typer=`head -$ii ${pidnm}find_out_temp.txt | tail -1 |  awk '{print $3}'`
      else
         rar=$ranh
         decr=$decnh
         typer=99
#         echo the mode $source
      fi

#check if the source is BZ/WHSP already has matched, and echo the conessearch
# 99 for sed moed, -9999 for pulsar. -1111 for Gamma, -2222 for GRB, no number
      if  [ $typer -gt 10000 -o $typer -lt -40000 -o $typer -eq 99 -o $typer -eq -9999 ]; then
         nn=$nn+1
         if [ $nn = $source -o $source == sed -o $source == lcurve ]; then
            echo conesearch --db ${HERE}/cats2.ini --catalog WISH352 --ra $rar --dec $decr --radius 15 --runit arcsec --columns default -o ${pidnm}wish352.$nn.2.csv > ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog TGSS150 --ra $rar --dec $decr --radius 15 --runit arcsec --columns default -o ${pidnm}tgss150.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog GLEAM --ra $rar --dec $decr --radius 15 --runit arcsec --columns default -o ${pidnm}gleam.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog LoTSS --ra $rar --dec $decr --radius 15 --runit arcsec --columns default -o ${pidnm}lotss.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog VLSSR --ra $rar --dec $decr --radius 15 --runit arcsec --columns default -o ${pidnm}vlssr.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog PMN --ra $rar --dec $decr --radius 30 --runit arcsec --columns default -o ${pidnm}pmn.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog GB6 --ra $rar --dec $decr --radius 30 --runit arcsec --columns default -o ${pidnm}gb6.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog GB87 --ra $rar --dec $decr --radius 30 --runit arcsec --columns default -o ${pidnm}gb87.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog AT20G --ra $rar --dec $decr --radius 15 --runit arcsec --columns default -o ${pidnm}at20.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog ATPMN --ra $rar --dec $decr --radius 15 --runit arcsec --columns default -o ${pidnm}atpmn.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog CRATES --ra $rar --dec $decr --radius 15 --runit arcsec --columns default -o ${pidnm}crates.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog NORTH20 --ra $rar --dec $decr --radius 2 --runit arcmin --columns default -o ${pidnm}north20.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog PCCS44 --ra $rar --dec $decr --radius 3 --runit arcmin --columns default -o ${pidnm}pccs44.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog PCCS70 --ra $rar --dec $decr --radius 3 --runit arcmin --columns default -o ${pidnm}pccs70.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog PCCS100 --ra $rar --dec $decr --radius 3 --runit arcmin --columns default -o ${pidnm}pccs100.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog PCCS143 --ra $rar --dec $decr --radius 3 --runit arcmin --columns default -o ${pidnm}pccs143.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog PCCS217 --ra $rar --dec $decr --radius 3 --runit arcmin --columns default -o ${pidnm}pccs217.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog PCCS353 --ra $rar --dec $decr --radius 3 --runit arcmin --columns default -o ${pidnm}pccs353.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog PCNT --ra $rar --dec $decr --radius 3 --runit arcmin --columns default -o ${pidnm}pcnt.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog ALMA --ra $rar --dec $decr --radius 15 --runit arcsec --columns default -o ${pidnm}alma.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog WISE --ra $rar --dec $decr --radius 10 --runit arcsec --columns default -o ${pidnm}wise.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog WISEME --ra $rar --dec $decr --radius 10 --runit arcsec --columns default -o ${pidnm}wiseme.$nn.2.csv >> ${pidnm}vosearch.txt
#            echo conesearch --db ${HERE}/cats2.ini --catalog NEOWISE --ra $rar --dec $decr --radius 10 --runit arcsec --columns default -o ${pidnm}neowise.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog 2MASS --ra $rar --dec $decr --radius 10 --runit arcsec --columns default -o ${pidnm}2mass.$nn.2.csv >> ${pidnm}vosearch.txt
#            echo conesearch --db ${HERE}/cats2.ini --catalog SPIRE250 --ra $rar --dec $decr --radius 15 --runit arcsec --columns default -o ${pidnm}spire250.$nn.2.csv >> ${pidnm}vosearch.txt
#            echo conesearch --db ${HERE}/cats2.ini --catalog SPIRE350 --ra $rar --dec $decr --radius 15 --runit arcsec --columns default -o ${pidnm}spire350.$nn.2.csv >> ${pidnm}vosearch.txt
#            echo conesearch --db ${HERE}/cats2.ini --catalog SPIRE500 --ra $rar --dec $decr --radius 15 --runit arcsec --columns default -o ${pidnm}spire500.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog H-ATLAS-DR1 --ra $rar --dec $decr --radius 15 --runit arcsec --columns default -o ${pidnm}hatlas1.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog H-ATLAS-DR2NGP --ra $rar --dec $decr --radius 15 --runit arcsec --columns default -o ${pidnm}hatlas2ngp.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog H-ATLAS-DR2SGP --ra $rar --dec $decr --radius 15 --runit arcsec --columns default -o ${pidnm}hatlas2sgp.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog AKARIBSC --ra $rar --dec $decr --radius 20 --runit arcsec --columns default -o ${pidnm}akaribsc.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog IRAS-PSC --ra $rar --dec $decr --radius 1 --runit arcmin --columns default -o ${pidnm}iraspsc.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog SDSS --ra $rar --dec $decr --radius 10 --runit arcsec --columns default -o ${pidnm}sdss.$nn.2.csv >> ${pidnm}vosearch.txt
            #echo conesearch --db ${HERE}/cats2.ini --catalog USNO --ra $rar --dec $decr --radius 10 --runit arcsec --columns default -o ${pidnm}usno.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog HSTGSC --ra $rar --dec $decr --radius 10 --runit arcsec --columns default -o ${pidnm}hst.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog PanSTARRS --ra $rar --dec $decr --radius 10 --runit arcsec --columns default -o ${pidnm}panstarrs.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog GAIA --ra $rar --dec $decr --radius 10 --runit arcsec --columns default -o ${pidnm}gaia.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog GALEX --ra $rar --dec $decr --radius 15 --runit arcsec --columns default -o ${pidnm}galex.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog XMMOM --ra $rar --dec $decr --radius 15 --runit arcsec --columns default -o ${pidnm}xmmom.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog UVOT --ra $rar  --dec $decr --radius 15 --runit arcsec --columns default -o ${pidnm}uvot.$nn.2.csv >> ${pidnm}vosearch.txt
#            echo conesearch --db ${HERE}/cats2.ini --catalog CMA --ra $rar  --dec $decr --radius 30 --runit arcsec --columns default -o ${pidnm}cma.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog XRTSPEC --ra $rar --dec $decr --radius 15 --runit arcsec --columns default -o ${pidnm}xrtspec.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog OULC --ra $rar --dec $decr --radius 15 --runit arcsec --columns default -o ${pidnm}oulc.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats1.ini --catalog OUSXB --ra $rar --dec $decr --radius 15 --runit arcmin --columns default -o ${pidnm}ousxb.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats1.ini --catalog OUSXG --ra $rar --dec $decr --radius 15 --runit arcmin --columns default -o ${pidnm}ousxg.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog BEPPOSAX --ra $rar --dec $decr --radius 30 --runit arcsec --columns default -o ${pidnm}bepposax.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog BAT105m --ra $rar --dec $decr --radius 10 --runit arcmin --columns default -o ${pidnm}bat105.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog 3FGL --ra $rar  --dec $decr --radius 20 --runit arcmin --columns default -o ${pidnm}3fgl.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog 2FHL --ra $rar  --dec $decr --radius 20 --runit arcmin --columns default -o ${pidnm}2fhl.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog 4FGL-DR2 --ra $rar --dec $decr --radius 20 --runit arcmin --columns default -o ${pidnm}4fgldr2.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog 3FHL --ra $rar  --dec $decr --radius 20 --runit arcmin --columns default -o ${pidnm}3fhl.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog 2BIGB --ra $rar  --dec $decr --radius 10 --runit arcmin --columns default -o ${pidnm}2bigb.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog 2AGILE --ra $rar  --dec $decr --radius 50 --runit arcmin --columns default -o ${pidnm}2agile.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog FermiMeV --ra $rar  --dec $decr --radius 30 --runit arcmin --columns default -o ${pidnm}fmev.$nn.2.csv >> ${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog FMonLC --ra $rar  --dec $decr --radius 30 --runit arcmin --columns default -o ${pidnm}fmonlc.$nn.2.csv >> ${pidnm}vosearch.txt
            echo specsearch --db ${HERE}/cats2.ini --service MAGIC --ra $rar  --dec $decr --radius 10 --runit arcmin --columns default -o ${pidnm}magictt.$nn.2.csv >> ${pidnm}vosearch.txt
            echo specsearch --db ${HERE}/cats2.ini --service VERITAS --ra $rar  --dec $decr --radius 10 --runit arcmin --columns default -o ${pidnm}veritas.$nn.2.csv >> ${pidnm}vosearch.txt
            echo $rar $decr $typer $nn
            racand=$rar
            deccand=$decr
         fi
      fi
   done

#check if the number input are not valid, or the user want to quit
   if [ ! -f ${pidnm}vosearch.txt -a ${source} != q -a ${source} != candlist -a ${source} != nocand  -a ${source} != sed -a ${source} != a -a ${source} != lcurve -a ${source} != n ]; then
      source=nosource
   fi
   if [ $source == nosource ]; then
      echo SOURCE NUMBER NOT FOUND!
   elif [ $source == candlist ]; then
      echo
      cat phase1
   elif [ $source == nocand ]; then
      echo Candidates NOT Found!!!
      source=q
   elif [ $source == q ]; then
      echo
   elif [ $source == n ]; then
      cat Sed.txt | grep -v = > sed4nupeak.txt
      python ~/app/nu_peak.py --dec $decr  --sed_path sed4nupeak.txt
   elif [ $source == a ]; then
      rm -rf ${pidnm}vou-aladin-cand.html
      ${BINF}/aladin_interface ${pidnm}output1.csv ${pidnm}find_out_temp.txt ${pidnm}candidates_posix.txt ${HERE}/aladin_script.js ${pidnm}vou-aladin-cand.html $sfov
#      mv vou-aladin.html ~/
#      open ${pidnm}vou-aladin-cand.html
#check if the phase 2 data are already there, if no, retrieving PHASE 2 data
   else
      ncat2=`cat ${pidnm}vosearch.txt | wc -l`
      catthere=`ls ${pidnm}*.$source.2.csv 2>/dev/null | wc -l`
      if [ $catthere != 0 ]; then
         echo
         echo data already downloaded!
         rm -f ${pidnm}vosearch.txt
      else
         rm -f ${pidnm}voerror.txt
         sh ${HERE}/queue.sh -f ${pidnm}vosearch.txt -n $ncat2 #1> voerror 2>&1
#         cat voerror | grep 'ERROR'
         declare -i irunvo=1
         if [ -s ${pidnm}voerror.txt ]; then
            voerror=yes
         else
            voerror=no
         fi
#      echo $voerror
         until [ $voerror == no -o $irunvo -ge 3 ]
         do
            runvo=y
#   [ -z $pid ] && read -p "Found error when searhcing through VO. Type 'n' to skip and continue; otherwise, rerun the VO conesearch" runvo
            [ -z $runvo ] && runvo=y
            cat ${pidnm}voerror.txt > ${pidnm}vosearch.txt
            rm -f ${pidnm}voerror.txt
            irunvo=${irunvo}+1
            echo
            echo VO search returned errors for one or more catalogues. Run conesearch again $irunvo
            [ $runvo != n -a $runvo != N ] && sh ${HERE}/queue.sh -f ${pidnm}vosearch.txt -n $ncat1
            if [ -s ${pidnm}voerror.txt ]; then
               voerror=yes
            else
               voerror=no
            fi
         done
         rm -f ${pidnm}vosearch.txt
         echo
      fi
      echo finish conesearch phase 2
      echo ra , dec = $racand , $deccand

#the nh value in phase 2
      if [ ! -z $nhthere ]; then
         nh equinox=2000 ra=$racand dec=$deccand > ${pidnm}nhvalue.txt
         nhval=`tail -1 ${pidnm}nhvalue.txt |  awk '{print $7}'`
         rm -f ${pidnm}nhvalue.txt
      fi

#read the phase 2 data
      if [ -f ${pidnm}magictt.*.2.csv ]; then
         if [ $source == sed -o $source == lcurve ]; then
            magic=1
         else
            magic=$source
         fi
         cat ${pidnm}magictt.$magic.2.csv | grep -v '#' > ${pidnm}magic.$magic.2.csv #deal with MAGIC extension
         rm -f ${pidnm}magictt.$magic.2.csv
      fi
#### check the extra Fermi LC
      if [ -s ftaptlc.csv ]; then
         if [ $source == sed -o $source == lcurve ]; then
            ftaptlc=1
         else
            ftaptlc=$source
         fi
         cp ftaptlc.csv ftaptlc.$ftaptlc.2.csv
      fi
      ls ${pidnm}*.2.csv > ${pidnm}catlist2.txt
      rm -f ${pidnm}output2.csv
      echo
      [ -s ${pidnm}catlist2.txt ] && ${BINF}/readcat ${pidnm}catlist2.txt ${pidnm}output2.csv $racand $deccand $sfov $nhval 0. 0. 0. 0.
      if [ ! -s ${pidnm}output2.csv ]; then
         echo RA= $rar Dec= $decr radius= $sfov > output2.csv
         echo nH= 0.0 Error circle/elliptical= 0.  0.  0.  0.  0. >> output2.csv
      fi

#running the PHASAE 2
      rm -f ${pidnm}sed.txt
      rm -f ${pidnm}error_map.txt
      [ -z $zzinput ] && zzinput=0
      ${BINF}/find_candidates2 ${pidnm}output2 ${pidnm}find_out_temp.txt ${pidnm}Sed_temp.txt ${pidnm}error_map.txt ${pidnm}Sed.txt ${HERE}/catrefs.txt ${BINF} $zzinput $source > ${pidnm}phase2
      echo
      cat ${pidnm}phase2
      #open sed.txt

###############################PLOTTING###############################
#plot the PHASE 2 error map and set default area 1 arcmin
      if [ -z $zoomin ]; then
         zoomin=60.
         echo set default second phase zoomin area to 1.0 arcmin.
      fi
      rm -f ${pidnm}LC.txt
      rm -f ${pidnm}error_map.*ps
      rm -f ${pidnm}sed.*ps
      rm -f ${pidnm}LC.*ps
      rm -f ${pidnm}LC_fermi.*ps
      if [ $plotsed != N ]; then
         if [ $runmode == f ]; then
            ${BINF}/gnomo_plot_types ${pidnm}error_map.txt,${pidnm}candidates_posix.txt,${pidnm}error_map.ps/vcps, $racand $deccand 0. $zoomin $racand $deccand 0 0 0 0 0 0 0 $racand $deccand
            handle_ps ${pidnm}error_map.ps
         fi
      # ps2eps -B -q error_map.ps
      # open error_map.eps

#plot the SED
         echo
         if [ $runmode != l ]; then
            ${BINF}/plot_sed ${pidnm}Sed.txt ${pidnm}sed.ps/cps $source
            handle_ps ${pidnm}sed.ps
         fi
         echo
         if [ $runmode != s ]; then
            ${BINF}/plot_lc ${pidnm}Sed.txt ${pidnm}Lc.ps/cps ${pidnm}LC.txt ${pidnm}LC_fermi.ps/cps $source
            [ -f ${pidnm}LC.ps ] && handle_ps ${pidnm}LC.ps
            [ -f ${pidnm}LC_fermi.ps ] && handle_ps ${pidnm}LC_fermi.ps
         fi
         rm -rf ${pidnm}vou-aladin-error.html
         ${BINF}/aladin_error_map ${pidnm}phase2 ${pidnm}error_map.txt ${HERE}/aladin_script.js ${pidnm}vou-aladin-error.html
#         mv vou-aladin.html ~/
#         open ${pidnm}vou-aladin-error.html
      fi
      # ps2eps -B -q sed.ps -R +
      # rm -f error_map.ps
      # rm -f sed.ps
      # open sed.eps
#      sort -n -k 5 ${pidnm}LC.txt > ${pidnm}LC_sorted.txt
      [ $showsed == 'osed' ] && open ${pidnm}Sed.txt
##############################################################

#for SED builder tool and future catalog
      rm -f ${pidnm}Out4SedTool.txt
      rased=`echo $racand | sed 's/\./_/g'`
      decsed=`echo $deccand | sed 's/\./_/g' | sed 's/-/m/g'`
      ${BINF}/convert_sed ${pidnm}Sed.txt ${pidnm}Out4SedTool.txt ${pidnm}Sed.csv
      [ -d Results/SEDtool -a -f Out4SedTool.txt ] && cp Out4SedTool.txt Results/SEDtool/$rased"_"$decsed"_"sedtool.txt
      if [ $runmode != l -a $plotsed != N ]; then
         rm -f ${pidnm}PySED.png
         python ${HERE}/VOU-plotSED.py --xaxis f --infile ${pidnm}Sed.csv --outfile ${pidnm}PySED.png --title 'Source nr. '$source' RA='$racand' Dec='$deccand --upperl 'yes'
         open ${pidnm}PySED.png
      fi
##############################################################

#save the phase 2 results in the folder
      [ -d Results/$xrtnm -a -f sed.*ps ] && cp sed.*ps Results/$xrtnm/$source"_"sed.eps
      [ -d Results/$xrtnm -a -f error_map.*ps ] && cp error_map.*ps Results/$xrtnm/$source"_"error_map.eps
      [ -d Results/$xrtnm -a -f LC.*ps ] && cp LC.*ps Results/$xrtnm/$source"_"LC.eps
      [ -d Results/$xrtnm -a -f LC_fermi.*ps ] && cp LC_fermi.*ps Results/$xrtnm/$source"_"LC_fermi.eps
      [ -d Results/$xrtnm -a -f Sed.txt ] && cp Sed.txt Results/$xrtnm/$source"_"Sed.txt
      [ -d Results/$xrtnm -a -f output2.csv ] && cp output2.csv Results/$xrtnm/$source"_"output2.csv

#save the results as a tex pdf
#      cp sed.eps fortex/$source"_"sed.eps
#      cp error_map.eps fortex/$source"_"error_map.eps
#      cp LC.eps fortex/$source"_"LC.eps
#      cp candidates.eps fortex/.
#      cp RX_map.eps fortex/.
      
#copy figure text to tex file
      [ $runmode == f -a -d Results/$xrtnm ] && echo \\subsection{Candidate No. $source}  > texsed_$source.txt
      [ $runmode == f -a -d Results/$xrtnm ] && echo Candidate position: $racand , $deccand >> texsed_$source.txt
      [ $runmode == f -a -d Results/$xrtnm ] && echo \\begin{figure}[h!] >> texsed_$source.txt
      [ $runmode == f -a -d Results/$xrtnm ] && echo \\centering >> texsed_$source.txt
      [ -d Results/$xrtnm -a -f sed.*ps ] && echo \\includegraphics[width=0.3\\linewidth]{Results/$xrtnm/$source"_"error_map.eps} >> texsed_$source.txt
      [ -d Results/$xrtnm -a -f error_map.*ps ] && echo \\includegraphics[width=0.69\\linewidth]{Results/$xrtnm/$source"_"sed.eps} >> texsed_$source.txt
      [ $runmode == f -a -d Results/$xrtnm ] && echo \\caption{Left: Error Circle Map. Right: SED.} >> texsed_$source.txt
      [ $runmode == f -a -d Results/$xrtnm ] && echo \\end{figure} >> texsed_$source.txt
      if [ -f LC.*ps -a -d Results/$xrtnm ]; then
         echo \\begin{figure}[h!] >> texsed_$source.txt
         echo \\centering >> texsed_$source.txt
         echo \\includegraphics[width=0.95\\linewidth]{Results/$xrtnm/$source"_"LC.eps} >> texsed_$source.txt
         echo \\caption{Light Curve.} >> texsed_$source.txt
         echo \\end{figure} >> texsed_$source.txt
      fi
      if [ -f LC_fermi.*ps -a -d Results/$xrtnm ]; then
         echo \\begin{figure}[h!] >> texsed_$source.txt
         echo \\centering >> texsed_$source.txt
         echo \\includegraphics[width=0.95\\linewidth]{Results/$xrtnm/$source"_"LC_fermi.eps} >> texsed_$source.txt
         echo \\caption{Fermi gamma-ray Light Curve.} >> texsed_$source.txt
         echo \\end{figure} >> texsed_$source.txt
      fi

#PID number for online version
      if [ $pid ]; then
         echo PID number = $pid
         [ -f ${pidnm}sed.eps ] && mv ${pidnm}sed.eps $pid"_"$source"_"sed.eps
         [ -f ${pidnm}error_map.eps ] && mv ${pidnm}error_map.eps $pid"_"$source"_"error_map.eps
         [ -f ${pidnm}LC.eps ] && mv ${pidnm}LC.eps $pid"_"$source"_"LC.eps
         [ -f ${pidnm}LC_fermi.eps ] && mv ${pidnm}LC_fermi.eps $pid"_"$source"_"LC_fermi.eps
         [ -f ${pidnm}error_map.txt ] && mv ${pidnm}error_map.txt $pid"_"$source"_"error_map.txt
         [ -f ${pidnm}Sed.txt ] && mv ${pidnm}Sed.txt $pid"_"$source"_"sed.txt
         [ -f ${pidnm}Out4SedTool.txt ] && mv ${pidnm}Out4SedTool.txt $pid"_"$source"_"Out4SedTool.txt
      fi
   fi
done

# Make a tex pdf file
if [ -d Results/$xrtnm ]; then
   cat ${HERE}/headtex.txt > vou-blazars.tex
   cat texcand.txt 2>/dev/null >> vou-blazars.tex
   cat texsed_*.txt 2>/dev/null >> vou-blazars.tex
   echo \\end{document} >> vou-blazars.tex
#   pdflatex vou-blazars.tex
#   open vou-blazars.pdf
fi

#remove the vo files from various catalogs
rm -rf eada_files
rm -f ${pidnm}*.1.csv
rm -f ${pidnm}*.i.csv
rm -f ${pidnm}*.2.csv
rm -f ${pidnm}xrtdeep.csv
rm -f Results/$xrtnm/*.pdf
rm -f vou-blazars.aux
rm -f vou-blazars.log
unset zoomin zzinput
unset posra posdec sfov nhval r1 major1 minor1 posa1 r2 major2 minor2 posa2 runmode
unset VOUB_AUTOSED