#!/bin/bash

HERE=$(cd `dirname $BASH_SOURCE`; pwd) ##point to bin file and fortran file
# Directory where the fortran binaries are:
BINF="${HERE}/fort"
#
SITE='standard'
PESCARA='online'
#PESCARA='offline'
#SITE='alternative'
if [ $SITE == 'alternative' ]; then
   cp ${HERE}/cats1.ini.alternative ${HERE}/cats1.ini
   cp ${HERE}/cats2.ini.alternative ${HERE}/cats2.ini
else
   cp ${HERE}/cats1.ini.orig.hybrid ${HERE}/cats1.ini
   cp ${HERE}/cats2.ini.orig.hybrid ${HERE}/cats2.ini
fi

handle_ps () {
    [ `which ps2eps 2> /dev/null` ] || return 0
    FILE_PS=$1
    if [ $plotsed != N ]; then
    if [ `which open 2> /dev/null` ]; then
       echo ""
    else
      gv ${FILE_PS%.ps}.eps &
    fi
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


# If no arguments given, print Help and exit.
[ "${#@}" -eq 0 ] && { help; exit 0; }

##############################################################
# Initial set up, read the parameters input
##############################################################

checkvo=fine

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
      --auto)
         VOUB_AUTOSED=$2; shift;;
      --plot)
         plotsed=$2; shift;;
      --allcats)
         allcatalog=$2; shift;;
      --light)
        light=$2; shift;;
      --legend)
         plotlab=$2; shift;;
      --PID)
         pid=$2; shift;;
      --IDPATH)
         oupath=$2; shift;;
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
echo Running VOU-Blazars V2.23Hybrid
echo

if (( $(echo "$posdec < -90" | bc -l) || $(echo "$posdec > 90" | bc -l) )); then
   echo "Declination outside the allowed range (-90/+90)"
   exit 1
fi

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
else
   pidnm=${oupath}/${pid}"_"
fi
test ! -d tmp && mkdir tmp
test ! -d tmp/${oupath} && mkdir tmp/${oupath}

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
   nhval=`python3 ${HERE}/nh.py --ra $ranh --dec $decnh `
fi
echo "Nh from heasarc "$nhval
[ -z $r1 ] && r1=`grep 'RADIUS1' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $emaj1 ] && emaj1=`grep 'MAJOR1' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $emin1 ] && emin1=`grep 'MINOR1' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $posa1 ] && posa1=`grep 'POSANG1' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $r2 ] && r2=`grep 'RADIUS2' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $emaj2 ] && emaj2=`grep 'MAJOR2' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $emin2 ] && emin2=`grep 'MINOR2' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $posa2 ] && posa2=`grep 'POSANG2' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $runmode ] && runmode=`grep 'MODE' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $plotsed ] && plotsed=`grep 'PLOTSED' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $allcatalog ] && allcatalog=`grep 'ALLCATS' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $light ] && light=`grep 'LIGHT' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $plotlab ] && plotlab=`grep 'LEGEND' ${HERE}/config_vou.txt | awk '{print $2}'`
if [ $runmode != f ]; then
   sfov=1.
   r1=0.
   emaj1=0.
   emin1=0.
   posa1=0.
   r2=0.
   emaj2=0.
   emin2=0.
   posa2=0.
fi

#set $ranh $decnh $sfov $nhval $r1 $emaj1 $emin1 $posa1 $r2 $emaj2 $emin2 $posa2 $runmode
#check if the SED mode need plotting or not

if [ $runmode == -s ]; then #last parameter
   runmode=s
   #plotsed=N
fi

##############################################################
# FIRST PHASE
# aim: find the candidates from available radio and X-ray sources
##############################################################

#FIRST phase data retrieving, runing the conesearch/conesearch2.py 1st phase
rm -f tmp/${pidnm}*.1.csv
#rm -f tmp/${pidnm}*_out.txt
rm -f tmp/${pidnm}vosearch.txt
rm -f tmp/gammacandidates.csv
rm -f tmp/catalog_error.txt
if [ $light != y -a $light != Y ]; then
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog NVSS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}nvss.1.csv > tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog FIRST --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}first.1.csv >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog SUMSS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}sumss.1.csv >> tmp/${pidnm}vosearch.txt
if [ $runmode != f ]; then 
      echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog VLASSQL --ra $ranh --dec $decnh --radius 4 --runit arcsec --o tmp/${pidnm}vlassql.1.csv >> tmp/${pidnm}vosearch.txt
else
      echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog VLASSQL --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}vlassql.1.csv >> tmp/${pidnm}vosearch.txt
fi
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog RACS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}racs.1.csv >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog 2SXPS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}2sxps.1.csv >> tmp/${pidnm}vosearch.txt
#if [ $PESCARA == 'online' ]; then
#   echo conesearch --db ${HERE}/cats1.ini --catalog 1OUSX --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o tmp/${pidnm}1ousx.1.csv >> tmp/${pidnm}vosearch.txt
#fi
   echo python3.10 ${HERE}/conesearch2.py  --db ${HERE}/cats1.ini --catalog RASS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}rass.1.csv >> tmp/${pidnm}vosearch.txt
if [ $SITE == 'alternative' ]; then
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog 4XMM-DR12 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}4xmmdr12.1.csv >> tmp/${pidnm}vosearch.txt
else
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog 4XMM-DR13 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}4xmmdr13.1.csv >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog BMW --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}bmw.1.csv >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog XMMSL2 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}xmmsl2.1.csv >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog IPC2E --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}ipc.1.csv >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py  --db ${HERE}/cats1.ini --catalog IPCSL --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}ipcsl.1.csv >> tmp/${pidnm}vosearch.txt
fi
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog WGACAT --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}wgacat.1.csv >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog Chandra-CSC2 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}chandracsc2.1.csv >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog eRASS1 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}erass1.1.csv  >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog eRASS1-S --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}erass1-s.1.csv  >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog eFEDS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}efeds.1.csv >> tmp/${pidnm}vosearch.txt
fi
if [ $runmode == f ]; then #skip the catalogs that we don't plot its data on SED
   if [ $light != y -a $light != Y ]; then
     echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog ZWCLUSTERS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}zw.1.csv >> tmp/${pidnm}vosearch.txt
     echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog PSZ2 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}psz2.1.csv >> tmp/${pidnm}vosearch.txt
     echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog MCXC --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}mcxc.1.csv >> tmp/${pidnm}vosearch.txt
#     echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog SDSSWHL --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}whl.1.csv >> tmp/${pidnm}vosearch.txt
      if [[ $SITE == 'standard' ]]; then
        echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog ABELL --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}abell.1.csv >> tmp/${pidnm}vosearch.txt
         echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog SWXCS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}swxcs.1.csv >> tmp/${pidnm}vosearch.txt
      fi
     echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog F2PSR --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}f2psr.1.csv >> tmp/${pidnm}vosearch.txt
     echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog MilliQuas --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}mquas.1.csv >> tmp/${pidnm}vosearch.txt
     echo conesearch --db ${HERE}/cats1.ini --catalog F357cat --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --columns default -o tmp/${pidnm}f357cat.1.csv >> tmp/${pidnm}vosearch.txt
   fi
     echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog CRATES --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}crates.1.csv >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog 5BZCat --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}5bzcat.1.csv >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog PULSAR --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}pulsar.1.csv >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog 3HSP --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}3hsp.1.csv >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog 4LAC-DR3 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}4lacdr3.1.csv >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog 4FGL-DR3 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}4fgldr3.1.csv >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog 3FHL --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}3fhl.1.csv >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog 2BIGB --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}2bigb.1.csv >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog 4FGL-DR4 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}4fgldr4.1.csv >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog MST12Y --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}mst12y.1.csv >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog 1FLE --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}fmev.1.csv >> tmp/${pidnm}vosearch.txt
#   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog FermiMeV --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}oldfmev.1.csv >> tmp/${pidnm}vosearch.txt
   echo conesearch --db ${HERE}/cats2.ini --catalog 2AGILE --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --columns default -o tmp/${pidnm}2agile.1.csv >> tmp/${pidnm}vosearch.txt
fi

rm -f tmp/${pidnm}voerror.txt
bash ${HERE}/queue.sh -f tmp/${pidnm}vosearch.txt -p $pidnm #2> voerror #!2>&1
#see if there is an error when searching through VO
declare -i irunvo=1
if [ -s tmp/${pidnm}voerror.txt ]; then
   voerror=yes
else
   voerror=no
fi
#until [ $voerror == no -o $irunvo -ge 3 ]
until [ $voerror == no -o $irunvo -ge 2 ]
do
   runvo=y
   [ -z $runvo ] && runvo=y
   cat tmp/${pidnm}voerror.txt > tmp/${pidnm}vosearch.txt
   rm -f tmp/${pidnm}voerror.txt
   irunvo=${irunvo}+1
   echo
   echo VO search returned errors or timed out for one or more catalogues. Running conesearch again $irunvo
   [ $runvo != n -a $runvo != N ] && bash ${HERE}/queue.sh -f tmp/${pidnm}vosearch.txt -p $pidnm
   if [ -s tmp/${pidnm}voerror.txt ]; then
      voerror=yes
   else
      voerror=no
   fi
done

if [ -s tmp/${pidnm}erass1.1.csv ]; then
  if [ -s tmp/${pidnm}erass1-s.1.csv ]; then
     cat tmp/${pidnm}erass1-s.1.csv | grep -v IAUName >> tmp/erass1.1.csv
  fi
else
  if [ -s tmp/${pidnm}erass1-s.1.csv ]; then
     cp tmp/${pidnm}erass1-s.1.csv  tmp/erass1.1.csv
  fi 
fi
if [ $voerror == yes ]; then
    echo "Warning some catalogs were not downloaded " > tmp/catalog_error.txt
    cat tmp/voerror.txt | awk -F'--catalog ' '{print $2}' | awk '{print $1}' >> tmp/catalog_error.txt
fi

echo  phase 1 completed
noOfCats=`ls tmp/${pidnm}*.1.csv 2>/dev/null | wc -l`
if [ $noOfCats == 0 ] && [ $runmode != s ]; then
  echo "There are no blazar candidates in this field"
  exit 0;
fi
if [ -s tmp/${pidnm}voerror.txt ]; then
   checkvo=check
fi
#read the data
ls tmp/${pidnm}*.1.csv > tmp/${pidnm}catlist1.txt

#copy the XRTDEEP result and add to the catlist1
rm -f tmp/${pidnm}xrtdeep.csv
if [ -s work/$xrtnm/table_flux_detections.csv ]; then
   cat work/$xrtnm/table_flux_detections.csv | sed 's/;/,/g' | sed 's/:/ /g' > tmp/${pidnm}xrtdeep.csv
fi
[ -s tmp/${pidnm}xrtdeep.csv ] && echo tmp/${pidnm}xrtdeep.csv >> tmp/${pidnm}catlist1.txt
rm -f tmp/${pidnm}output1.csv
#convert format from new site to heasarc output format
if [ $SITE == 'alternative' ]; then
   python3 ${HERE}/vizier2heasarc.py tmp/${pidnm}crates.1.csv crates
   python3 ${HERE}/vizier2heasarc.py tmp/${pidnm}pulsar.1.csv pulsar
   python3 ${HERE}/vizier2heasarc.py tmp/${pidnm}nvss.1.csv   nvss
   python3 ${HERE}/vizier2heasarc.py tmp/${pidnm}3fhl.1.csv   3fhl
fi
echo
[ -s tmp/${pidnm}catlist1.txt ] && ${BINF}/readcat tmp/${pidnm}catlist1.txt tmp/${pidnm}output1.csv $ranh $decnh $sfov $nhval $r1 $emaj1 $emin1 $posa1
#running the PHASE 1
rm -f tmp/${pidnm}*temp.txt
if [ $runmode == s -o $runmode == l ]; then
   ${BINF}/find_candidates1 tmp/${pidnm}output1 tmp/${pidnm}find_out_temp.txt tmp/${pidnm}RX_temp.txt tmp/${pidnm}Sed_temp.txt tmp/${pidnm}no_matched_temp.txt ${BINF} 1 | cat > tmp/${pidnm}phase1
   echo
   #cat tmp/${pidnm}phase1
else
   ${BINF}/find_candidates1 tmp/${pidnm}output1 tmp/${pidnm}find_out_temp.txt tmp/${pidnm}RX_temp.txt tmp/${pidnm}Sed_temp.txt tmp/${pidnm}no_matched_temp.txt ${BINF} 1 | cat > tmp/${pidnm}phase1
   echo
   cat tmp/${pidnm}phase1

###############################PLOTTING###############################
#plot the result, candidates map and radio-X-ray source map
   rm -f tmp/${pidnm}candidates.*ps
   rm -f tmp/${pidnm}RX_map.*ps
   sort -n -k 3 tmp/${pidnm}RX_temp.txt > tmp/${pidnm}RX_sorted.txt
   rm -f tmp/${pidnm}RX_temp.txt
   ${BINF}/gnomo_plot_types tmp/${pidnm}RX_sorted.txt,tmp/${pidnm}candidates_posix.txt,tmp/${pidnm}RX_map.ps/vcps,${BINF}, $ranh $decnh $r1 $sfov $ranh $decnh $emaj1 $emin1 $posa1 $r2 $emaj2 $emin2 $posa2 $ranh $decnh
   ${BINF}/gnomo_plot_types tmp/${pidnm}find_out_temp.txt,tmp/${pidnm}candidates_posix.txt,tmp/${pidnm}candidates.ps/vcps,${BINF}, $ranh $decnh $r1 $sfov $ranh $decnh $emaj1 $emin1 $posa1 $r2 $emaj2 $emin2 $posa2 $ranh $decnh
   if [ $plotlab == Y -o $plotlab == y ]; then
      open ${HERE}/../data/legend_cand.png
   fi
##############################################################
fi
##############################################################
# Intermediate PHASE
# aim: find more available candiates from single radio/ single X-ray sources.
##############################################################

rm -f tmp/${pidnm}Intermediate_out.txt
rm -f tmp/${pidnm}output_int.csv
if [ -s tmp/${pidnm}no_matched_temp.txt ]; then
#check the error circle
   #echo checking intermediate $ra $emaj1
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
   
   process=y
   [ ${pid} ] && process=n
   [ ${VOUB_AUTOSED} ] && process=n
   if [ $radcrit -ge 20 ]; then
      if [ -z ${VOUB_AUTOSED} ]; then
         echo
         echo Warning! Large searching radius in Intermediate phase for conesearch.
         [ -z $pid ] && read -p "Type 'y' to continue! Otherwise, stop processing Intermediate phase." process #rcut
         [ -z $process ] && process=n
         [ $process != y -a $process != Y ] && echo Not searching sources in Intermediate phase.
      fi
   fi

#running conesearch for Intermediate phase
   if [ $process == y -o $process == Y ]; then
      rm -f tmp/${pidnm}*.i.csv
      ln=`cat tmp/${pidnm}no_matched_temp.txt | wc -l`
      declare -i cc
      cc=0
      for (( jj=1; jj<=$ln; jj=jj+1 ))
      do
         raint=`head -$jj tmp/${pidnm}no_matched_temp.txt | tail -1 |  awk '{print $5}'`
         decint=`head -$jj tmp/${pidnm}no_matched_temp.txt | tail -1 |  awk '{print $6}'`
         typeint=`head -$jj tmp/${pidnm}no_matched_temp.txt | tail -1 |  awk '{print $10}'`
         #declare -i typeint
         if [ $typeint -gt 50 -o $typeint -lt 0 ]; then
            cc=$cc+1
         fi
      done
      echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog GALEX --ra $ranh --dec $decnh --radius $radint --runit arcmin --o tmp/${pidnm}galex.i.csv > tmp/${pidnm}vosearch.txt
      echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog PMN --ra $ranh --dec $decnh --radius $radint --runit arcmin --o tmp/${pidnm}pmn.i.csv >> tmp/${pidnm}vosearch.txt
      echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog GB6 --ra $ranh --dec $decnh --radius $radint --runit arcmin --o tmp/${pidnm}gb6.i.csv >> tmp/${pidnm}vosearch.txt
      echo Searching further data in intermediate phase for $cc source
      rm -f tmp/${pidnm}voerror.txt
      bash ${HERE}/queue.sh -f tmp/${pidnm}vosearch.txt -p $pidnm #1> voerror 2>&1

# check if there is error during the searching
      declare -i irunvo=1
      if [ -s tmp/${pidnm}voerror.txt ]; then
         voerror=yes
      else
         voerror=no
      fi
      #until [ $voerror == no -o $irunvo -ge 3 ]
      until [ $voerror == no -o $irunvo -ge 2 ]
      do
         runvo=y
         [ -z $runvo ] && runvo=y
         cat tmp/${pidnm}voerror.txt > tmp/${pidnm}vosearch.txt
         rm -f tmp/${pidnm}voerror.txt
         irunvo=${irunvo}+1
         echo
         echo VO search returned errors or timed out for one or more catalogues. Running conesearch again $irunvo
         [ $runvo != n -a $runvo != N ] && bash ${HERE}/queue.sh -f tmp/${pidnm}vosearch.txt -p $pidnm
         if [ -s tmp/${pidnm}voerror.txt ]; then
            voerror=yes
         else
            voerror=no
         fi
      done
      if [ $voerror == yes ]; then
         if [ ! -f  tmp/catalog_error.txt ]; then
            echo "Warning some catalogs were not downloaded " > tmp/catalog_error.txt
         fi
            cat tmp/voerror.txt | awk -F'--catalog ' '{print $2}' | awk '{print $1}' >> tmp/catalog_error.txt
      fi
      echo intermediate phase completed
      if [ -s tmp/${pidnm}voerror.txt ]; then
         checkvo=check
      fi
      rm -f tmp/${pidnm}vosearch.txt

#read the Intermediate phase data
      ls tmp/${pidnm}*.i.csv > tmp/${pidnm}catlist_int.txt
      rm -f tmp/${pidnm}output_int.csv
      echo
      [ -s tmp/${pidnm}catlist_int.txt ] && ${BINF}/readcat tmp/${pidnm}catlist_int.txt tmp/${pidnm}output_int.csv $ranh $decnh $sfov $nhval $r1 $emaj1 $emin1 $posa1

#run the Intermediate phase
      ${BINF}/find_candidates_int tmp/${pidnm}output_int tmp/${pidnm}no_matched_temp.txt tmp/${pidnm}find_out_temp.txt tmp/${pidnm}Intermediate_out.txt ${BINF}
   fi
fi
#transfer the data from Intermediate phase to 2 phase
if [ -s tmp/${pidnm}Intermediate_out.txt ]; then
   nint=`tail -1 tmp/${pidnm}Intermediate_out.txt | awk '{print $1}' `
   if [ $nint != 'No' ]; then 
      head -$nint tmp/${pidnm}Intermediate_out.txt >> tmp/${pidnm}find_out_temp.txt
      nnllint=`cat tmp/${pidnm}Intermediate_out.txt | wc -l`
      #echo "cat tmp/Intermediate_out.txt"
      #cat tmp/Intermediate_out.txt
      declare -i nsedl=$nnllint-$nint
      declare -i nsedu=$nnllint-$nint-1
      tail -$nsedl tmp/${pidnm}Intermediate_out.txt | head -$nsedu >> tmp/${pidnm}Sed_temp.txt
      #echo "cat tmp/Sed_temp.txt"
      #cat tmp/Sed_temp.txt
   fi
###############################PLOTTING###############################
#plot the candidates again with Intermediate phase finish ???
#   ${BINF}/gnomo_plot_types tmp/${pidnm}find_out_temp.txt,tmp/${pidnm}candidates_posix.txt,tmp/${pidnm}candidates.ps/vcps,${BINF}, $ranh $decnh $r1 $sfov $ranh $decnh $emaj1 $emin1 $posa1 $r2 $emaj2 $emin2 $posa2 $ranh $decnh
######################################################################
fi

###############################  HTML table file  ###############################
if [ $runmode == f ]; then
   echo "Number R.A. Dec. Type Name" > tmp/${pidnm}candidates.csv
   cat tmp/${pidnm}find_out_temp.txt | awk ' $3>10000 {print $1, $2, int($3/10000), $4 $5}  $3<-40000 {print $1, $2, int($3/10000), $4 $5} $3==-9999 {print $1, $2, int($3/10000), $4 $5} ' | awk '{print NR, $1, $2, "type"$3, $4 $5}' >> tmp/${pidnm}candidates.csv
   cat tmp/${pidnm}candidates.csv | sed 's/type1/HBL/g' | sed 's/type2/IBL/g' | sed 's/type3/LBL/g' | sed 's/type4/jetted-AGN/g' | sed 's/type5/Unknown/g' | sed 's/type-5/3HSP/g' | sed 's/type-6/5BZCat/g' | sed 's/type-7/CRATES/g' | sed 's/type0/Pulsar/g' > tmp/${pidnm}candidates.csv
#echo "candidates.txt"
#cat tmp/candidates.csv
#echo "find_out_temp.txt"
#cat tmp/find_out_temp.txt
   
   if [ $pidnm ]; then
      echo "<table>" > tmp/${pidnm}candidates.html
      ln=`cat tmp/${pidnm}candidates.csv | wc -l`
      for (( ii=1; ii<=${ln}; ii=ii+1 ))
      do
         echo "<tr>" >> tmp/${pidnm}candidates.html
         for (( jj=1; jj<=4; jj=jj+1 ))
         do
            htmlpar=`head -$ii tmp/${pidnm}candidates.csv | tail -1 |  awk '{print $'$jj'}'`
            if [ $ii == 1 ]; then
               echo "<th> " $htmlpar " </th>" >> tmp/${pidnm}candidates.html
            else
               echo "<td> " $htmlpar " </td>" >> tmp/${pidnm}candidates.html
            fi
         done
         echo "</tr>" >> tmp/${pidnm}candidates.html
      done
      echo "</table>" >> tmp/${pidnm}candidates.html
   fi
   
   sed -i '' 's/\ /\,/g' tmp/${pidnm}candidates.csv
   
fi

#######################################################
# Second PHASE
# aim: Plot the SED and error circle map for a specified source
#######################################################

#begin phase 2 setup, leave the tool when no candidates
rm -f tmp/${pidnm}*.2.csv
rm -f tmp/${pidnm}*.2.txt
rm -f tmp/${pidnm}texsed*
source=1

if [ ! -s tmp/${pidnm}find_out_temp.txt ]; then
   source=nocand
   ln=0
fi
if [ $runmode != s ]; then
#######################################################
# python plotting of the candidates map
#######################################################
   rrad=0
   rm -f candidates.png
   if [ $emaj1  != "0." ]; then 
     python3.10 ${HERE}/skymap.py --ra $posra --dec $posdec --FOV $sfov --major $emaj1  --minor $emin1 --angle $posa1 --major2 $emaj2  --minor2 $emin2 --angle2 $posa2 --infile_candidates tmp/candidates.csv --infile_multi tmp/find_out_temp.txt --infile_local IHBL-catalog --out candidates.png 
#     echo "python3.10 ${HERE}/skymap.py --ra $posra --dec $posdec --FOV $sfov --major $emaj1  --minor $emin1 --angle $posa1 --major2 $emaj2  --minor2 $emin2 --angle2 $posa2 --infile_candidates tmp/candidates.csv --infile_multi tmp/find_out_temp.txt --infile_local IHBL-catalog --out candidates.png "
   elif [ $r1 != "0." ]; then
     python3.10 ${HERE}/skymap.py --ra $posra --dec $posdec --FOV $sfov --radius $r1 --infile_candidates tmp/candidates.csv --infile_multi tmp/find_out_temp.txt --infile_local IHBL-catalog --out candidates.png 
   else  
     python3.10 ${HERE}/skymap.py --ra $posra --dec $posdec --FOV $sfov --infile_candidates tmp/candidates.csv --infile_multi tmp/find_out_temp.txt --infile_local IHBL-catalog --out candidates.png 
   fi
   #open candidates.png
fi
######
#read the number interested
until [ $source == sed -o $source == q -o $source == lcurve ]
do
   if [ $runmode == f -a $source != nocand ]; then
      echo
      if [ -z $VOUB_AUTOSED ]; then
        if [ -z $pid ]; then
           read -p "Please enter the source number and zoom in area (in arcsec) intered: (Type 'q' to quit) (Type 'candlist' to show all candidates)" source zzinput zoomin showsed
        else
           source=q
        fi
      else
        runiter=`ls tmp/${pidnm}*.$source.2.csv 2>/dev/null | wc -l`
        if [ ${runiter} -eq 0 ]; then
           [ -f tmp/${pidnm}Sed_temp.txt ] && `grep "matched source" tmp/${pidnm}Sed_temp.txt | awk '{print $1}' | uniq > tmp/${pidnm}sources.list`
           declare -i sourcenb
           sourcenb=`cat tmp/${pidnm}sources.list | wc -l`
           #sourcenb=${sourcenb}-1
           if [ $VOUB_AUTOSED == part -a ${sourcenb} -gt 5 ]; then
             python ${HERE}/select_cand.py --RA_cent $ranh --Dec_cent $decnh --quiet --input_list "tmp/"{$pidnm}"find_out_temp.txt" --input_file "tmp/"{$pidnm}"Sed_temp.txt" --output_list "tmp/"$pidnm"sources.list"
           fi
        fi
        declare -i sourcenb
        sourcenb=`cat tmp/${pidnm}sources.list | wc -l`
        sourcenb=${sourcenb}-1
        if [ $sourcenb -ge 0 ]; then
           source=`head -1 tmp/${pidnm}sources.list`
           cat tmp/${pidnm}sources.list | tail -${sourcenb} > tmp/${pidnm}sources.tmp
           mv tmp/${pidnm}sources.tmp tmp/${pidnm}sources.list
        else
           source=q
        fi
      fi
      ln=`cat tmp/${pidnm}find_out_temp.txt | wc -l`
   elif [ $runmode == s -o $runmode == l ]; then
      ln=1
      [ $runmode == s ] && source=sed
      [ $runmode == l ] && source=lcurve
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

#read the find_out list and the second phase ra dec
   declare -i nn
   nn=0
   [ -z $ln ] && ln=0
   for (( ii=1; ii<=$ln; ii=ii+1 ))
   do
      if [ $runmode == f ]; then
         rar=`head -$ii tmp/${pidnm}find_out_temp.txt | tail -1 |  awk '{print $1}'`
         decr=`head -$ii tmp/${pidnm}find_out_temp.txt | tail -1 |  awk '{print $2}'`
         typer=`head -$ii tmp/${pidnm}find_out_temp.txt | tail -1 |  awk '{print $3}'`
      else
         rar=$ranh
         decr=$decnh
         typer=99
      fi

#check if the source is BZ/WHSP already has matched, and echo the conessearch
# 99 for sed moed, -9999 for pulsar. -1111 for Gamma, -2222 for GRB, no number
      if  [ $typer -gt 10000 -o $typer -lt -40000 -o $typer -eq 99 -o $typer -eq -9999 ]; then
         nn=$nn+1
         if [ $nn = $source -o $source == sed -o $source == lcurve ]; then
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog WISH352 --ra $rar --dec $decr --radius 20 --runit arcsec --o tmp/${pidnm}wish352.$nn.2.csv > tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog TGSS150 --ra $rar --dec $decr --radius 1.5 --runit arcmin --o tmp/${pidnm}tgss150.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog LoTSS --ra $rar --dec $decr --radius 30 --runit arcsec --o tmp/${pidnm}lotss.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog VLSSR --ra $rar --dec $decr --radius 30 --runit arcsec --o tmp/${pidnm}vlssr.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog PMN --ra $rar --dec $decr --radius 45 --runit arcsec --o tmp/${pidnm}pmn.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog GB6 --ra $rar --dec $decr --radius 45 --runit arcsec --o tmp/${pidnm}gb6.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo  python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog GB87 --ra $rar --dec $decr --radius 45 --runit arcsec --o tmp/${pidnm}gb87.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py  --db ${HERE}/cats2.ini --catalog AT20G --ra $rar --dec $decr --radius 30 --runit arcsec --o tmp/${pidnm}at20.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog ATPMN --ra $rar --dec $decr --radius 15 --runit arcsec --o tmp/${pidnm}atpmn.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog CRATES --ra $rar --dec $decr --radius 5 --runit arcsec --o tmp/${pidnm}crates.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog NORTH20 --ra $rar --dec $decr --radius 1.5 --runit arcmin --o tmp/${pidnm}north20.$nn.2.csv >> tmp/${pidnm}vosearch.txt
#            echo conesearch --db ${HERE}/cats2.ini --catalog F357det --ra $rar --dec $decr --radius 10 --runit arcsec --columns default -o tmp/${pidnm}f357det.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog PCCS30 --ra $rar --dec $decr --radius 10 --runit arcmin --o tmp/${pidnm}pccs30.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog PCCS44 --ra $rar --dec $decr --radius 8 --runit arcmin  --o tmp/${pidnm}pccs44.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog PCCS70 --ra $rar --dec $decr --radius 6 --runit arcmin  --o tmp/${pidnm}pccs70.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog PCCS100 --ra $rar --dec $decr --radius 6 --runit arcmin  --o tmp/${pidnm}pccs100.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog PCCS143 --ra $rar --dec $decr --radius 6 --runit arcmin  --o tmp/${pidnm}pccs143.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog PCCS217 --ra $rar --dec $decr --radius 6 --runit arcmin  --o tmp/${pidnm}pccs217.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog PCCS353 --ra $rar --dec $decr --radius 6 --runit arcmin --o tmp/${pidnm}pccs353.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog PCCS545 --ra $rar --dec $decr --radius 7 --runit arcmin  --o tmp/${pidnm}pccs545.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog PCCS857 --ra $rar --dec $decr --radius 8 --runit arcmin  --o tmp/${pidnm}pccs857.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog PCNT --ra $rar --dec $decr --radius 5 --runit arcmin --o tmp/${pidnm}pcnt.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog ALMA --ra $rar --dec $decr --radius 15 --runit arcsec --o tmp/${pidnm}alma.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog WISE --ra $rar --dec $decr --radius 5 --runit arcsec --o tmp/${pidnm}wise.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog UNWISE --ra $rar --dec $decr --radius 5 --runit arcsec --o tmp/${pidnm}unwise.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog 2MASS --ra $rar --dec $decr --radius 5 --runit arcsec --o tmp/${pidnm}2mass.$nn.2.csv >> tmp/${pidnm}vosearch.txt
#            echo conesearch --db ${HERE}/cats2.ini --catalog SPIRE250 --ra $rar --dec $decr --radius 5 --runit arcsec --columns default -o tmp/${pidnm}spire250.$nn.2.csv >> tmp/${pidnm}vosearch.txt
#            echo conesearch --db ${HERE}/cats2.ini --catalog SPIRE350 --ra $rar --dec $decr --radius 5 --runit arcsec --columns default -o tmp/${pidnm}spire350.$nn.2.csv >> tmp/${pidnm}vosearch.txt
#            echo conesearch --db ${HERE}/cats2.ini --catalog SPIRE500 --ra $rar --dec $decr --radius 5 --runit arcsec --columns default -o tmp/${pidnm}spire500.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog H-ATLAS-DR1 --ra $rar --dec $decr --radius 10 --runit arcsec --o tmp/${pidnm}hatlas1.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog H-ATLAS-DR2NGP --ra $rar --dec $decr --radius 10 --runit arcsec --o tmp/${pidnm}hatlas2ngp.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog H-ATLAS-DR2SGP --ra $rar --dec $decr --radius 10 --runit arcsec --o tmp/${pidnm}hatlas2sgp.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog AKARIBSC --ra $rar --dec $decr --radius 10 --runit arcsec --o tmp/${pidnm}akaribsc.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog IRAS-PSC --ra $rar --dec $decr --radius 1 --runit arcmin --o tmp/${pidnm}iraspsc.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog SDSS --ra $rar --dec $decr --radius 2 --runit arcsec --o tmp/${pidnm}sdss.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            #echo conesearch --db ${HERE}/cats2.ini --catalog USNO --ra $rar --dec $decr --radius 10 --runit arcsec --columns default -o tmp/${pidnm}usno.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog HSTGSC --ra $rar --dec $decr --radius 10 --runit arcsec --o tmp/${pidnm}hst.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog PanSTARRS --ra $rar --dec $decr --radius 10 --runit arcsec --o tmp/${pidnm}panstarrs.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog GAIA2 --ra $rar --dec $decr --radius 10 --runit arcsec --o tmp/${pidnm}gaia2.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog GALEX --ra $rar --dec $decr --radius 5 --runit arcsec --o tmp/${pidnm}galex.$nn.2.csv >> tmp/${pidnm}vosearch.txt
#            echo conesearch --db ${HERE}/cats2.ini --catalog CMA --ra $rar  --dec $decr --radius 1 --runit arcmin --columns default -o tmp/${pidnm}cma.$nn.2.csv >> tmp/${pidnm}vosearch.txt
if [ $PESCARA == 'online' ]; then
            echo specsearch --db ${HERE}/cats2.ini --service MAGIC --ra $rar  --dec $decr --radius 5 --runit arcmin --columns default -o tmp/${pidnm}magictt.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo specsearch --db ${HERE}/cats2.ini --service MAGIC --ra $rar  --dec $decr --radius 5 --runit arcmin --columns default -o tmp/${pidnm}magictt.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo specsearch --db ${HERE}/cats2.ini --service VERITAS --ra $rar  --dec $decr --radius 5 --runit arcmin --columns default -o tmp/${pidnm}veritas.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog OULC --ra $rar --dec $decr --radius 15 --runit arcsec --columns default -o tmp/${pidnm}oulc.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog 2AGILE --ra $rar  --dec $decr --radius 30 --runit arcmin --columns default -o tmp/${pidnm}2agile.$nn.2.csv >> tmp/${pidnm}vosearch.txt
#            echo conesearch --db ${HERE}/cats2.ini --catalog SMARTS --ra $rar --dec $decr --radius 20 --runit arcsec --columns default -o tmp/${pidnm}smarts.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog XRTSPEC --ra $rar --dec $decr --radius 1 --runit arcmin --columns default -o tmp/${pidnm}xrtspec.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats1.ini --catalog OUSXB --ra $rar --dec $decr --radius 1 --runit arcmin --columns default -o tmp/${pidnm}ousxb.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats1.ini --catalog OUSXG --ra $rar --dec $decr --radius 1 --runit arcmin --columns default -o tmp/${pidnm}ousxg.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog BEPPOSAX --ra $rar --dec $decr --radius 2 --runit arcmin --columns default -o tmp/${pidnm}bepposax.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo conesearch --db ${HERE}/cats2.ini --catalog NuBlazar --ra $rar --dec $decr --radius 20 --runit arcsec --columns default -o tmp/${pidnm}nublazar.$nn.2.csv >> tmp/${pidnm}vosearch.txt
fi
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog BAT105m --ra $rar --dec $decr --radius 5 --runit arcmin --o tmp/${pidnm}bat105.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog 2FHL --ra $rar  --dec $decr --radius 10 --runit arcmin --o tmp/${pidnm}2fhl.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog 4FGL-DR3 --ra $rar --dec $decr --radius 20  --runit arcmin --o tmp/${pidnm}4fgldr3.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog 4FGL-DR4 --ra $rar --dec $decr --radius 20 --runit arcmin --o tmp/${pidnm}4fgldr4.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog 3FHL --ra $rar  --dec $decr --radius 7 --runit arcmin --o tmp/${pidnm}3fhl.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog 2BIGB --ra $rar  --dec $decr --radius 10 --runit arcmin --o tmp/${pidnm}2bigb.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog SPECFINDV3 --ra $rar  --dec $decr --radius 20 --runit arcsec --o tmp/${pidnm}specfindv3.$nn.2.txt >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog EPRS --ra $rar  --dec $decr --radius 1 --runit arcmin --o tmp/${pidnm}eprs.$nn.2.txt >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog RATAN-600 --ra $rar  --dec $decr --radius 20 --runit arcsec --o tmp/${pidnm}ratan600.$nn.2.txt >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog MM-MONITORING4 --ra $rar  --dec $decr --radius 20 --runit arcsec --o tmp/${pidnm}mmmonitoring4.$nn.2.txt >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog PACO --ra $rar  --dec $decr --radius 20 --runit arcsec --o tmp/${pidnm}paco.$nn.2.txt >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog PACOPCCS --ra $rar  --dec $decr --radius 20 --runit arcsec --o tmp/${pidnm}pacopccs.$nn.2.txt >> tmp/${pidnm}vosearch.txt
           echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog SPTSZ --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}sptsz.$nn.2.txt >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog GLEAMV2 --ra $rar  --dec $decr --radius 20 --runit arcsec --o tmp/${pidnm}gleamv2.$nn.2.txt >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog 1FLE --ra $rar  --dec $decr --radius 20 --runit arcmin --o tmp/${pidnm}fmev.$nn.2.csv >> tmp/${pidnm}vosearch.txt
#   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog FermiMeV --ra $ranh  --dec $decnh --radius 20  --runit arcmin --o tmp/${pidnm}oldfmev.1.csv >> tmp/${pidnm}vosearch.txt
            echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog WISEME --ra $rar --dec $decr --radius 5 --runit arcsec --o tmp/${pidnm}wiseme.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            if [ $allcatalog == y -o $allcatalog == Y ]; then
               echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog NEOWISE --ra $rar --dec $decr --radius 10 --runit arcsec --o tmp/${pidnm}neowise.$nn.2.csv --timeout 40 >> tmp/${pidnm}vosearch.txt
            fi
            if [ $SITE == 'standard' ]; then
               echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog KUEHR --ra $rar --dec $decr --radius 1 --runit arcmin --o tmp/${pidnm}kuehr.$nn.2.csv >> tmp/${pidnm}vosearch.txt
#               echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog GLEAM --ra $rar --dec $decr --radius 50 --runit arcsec --o tmp/${pidnm}gleam.$nn.2.csv >> tmp/${pidnm}vosearch.txt
               echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog XMMOM --ra $rar --dec $decr --radius 5 --runit arcsec --o tmp/${pidnm}xmmom.$nn.2.csv >> tmp/${pidnm}vosearch.txt
               echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog UVOT --ra $rar  --dec $decr --radius 5 --runit arcsec --o tmp/${pidnm}uvot.$nn.2.csv >> tmp/${pidnm}vosearch.txt
               echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog FMonLC --ra $rar  --dec $decr --radius 5 --runit arcmin --o tmp/${pidnm}fmonlc.$nn.2.csv >> tmp/${pidnm}vosearch.txt
            fi
            racand=$rar
            deccand=$decr
         fi
      fi
   done
#check if the number input are not valid, or the user want to quit
   if [ ! -f tmp/${pidnm}vosearch.txt -a ${source} != q -a ${source} != candlist -a ${source} != nocand  -a ${source} != sed -a ${source} != a -a ${source} != lcurve -a ${source} != n ]; then
      source=nosource
   fi
   if [ $source == nosource ]; then
      echo SOURCE NUMBER NOT FOUND!
   elif [ $source == candlist ]; then
      echo
      cat tmp/${pidnm}phase1
   elif [ $source == nocand ]; then
      echo Candidates NOT Found!!!
      source=q
   elif [ $source == q ]; then
      echo
   elif [ $source == n ]; then
      cat tmp/${pidnm}Sed.txt | grep -v = > sed4nupeak.txt
      python ~/app/nu_peak.py --dec $decr  --sed_path sed4nupeak.txt
   elif [ $source == a ]; then
      rm -rf tmp/${pidnm}vou-aladin-cand.html
      ${BINF}/aladin_interface tmp/${pidnm}output1.csv tmp/${pidnm}find_out_temp.txt tmp/${pidnm}candidates_posix.txt ${HERE}/aladin_script.js tmp/${pidnm}vou-aladin-cand.html $sfov
      open tmp/${pidnm}vou-aladin-cand.html
#check if the phase 2 data are already there, if no, retrieving PHASE 2 data
   else
      catthere=`ls tmp/${pidnm}*.$source.2.csv 2>/dev/null | wc -l`
      if [ $catthere != 0 ]; then
         echo
         echo data already downloaded!
         rm -f tmp/${pidnm}vosearch.txt
      else
         rm -f tmp/${pidnm}voerror.txt
         bash ${HERE}/queue.sh -f tmp/${pidnm}vosearch.txt -p $pidnm #1> voerror 2>&1
         declare -i irunvo=1
         if [ -s tmp/${pidnm}voerror.txt ]; then
            voerror=yes
         else
            voerror=no
         fi
         until [ $voerror == no -o $irunvo -ge 3 ]
#         until [ $voerror == no -o $irunvo -ge 1 ]
         do
            runvo=y
            [ -z $runvo ] && runvo=y
            cat tmp/${pidnm}voerror.txt > tmp/${pidnm}vosearch.txt
            rm -f tmp/${pidnm}voerror.txt
            irunvo=${irunvo}+1
            echo
            echo VO search returned errors or timed out for one or more catalogues. Running conesearch again $irunvo
            [ $runvo != n -a $runvo != N ] && bash ${HERE}/queue.sh -f tmp/${pidnm}vosearch.txt -p $pidnm
            if [ -s tmp/${pidnm}voerror.txt ]; then
               voerror=yes
            else
               voerror=no
            fi
         done
         if [ $voerror == yes ]; then
             if [ ! -f  tmp/catalog_error.txt ]; then
                echo "Warning some catalogs were not downloaded " > tmp/catalog_error.txt
             fi
             cat tmp/voerror.txt | awk -F'--catalog ' '{print $2}' | awk '{print $1}' >> tmp/catalog_error.txt
         fi
#New convert format from new site to heasarc output format
         if [ $source == 'sed' ]; then
            nnn=1
         else 
            nnn=$source
         fi 
         if [ $SITE == 'alternative' ]; then
            python3 ${HERE}/vizier2heasarc.py tmp/${pidnm}gb6.$nnn.2.csv gb6
            python3 ${HERE}/vizier2heasarc.py tmp/${pidnm}pccs44.$nnn.2.csv pcc
            python3 ${HERE}/vizier2heasarc.py tmp/${pidnm}pccs70.$nnn.2.csv pcc
            python3 ${HERE}/vizier2heasarc.py tmp/${pidnm}pccs100.$nnn.2.csv pcc
            python3 ${HERE}/vizier2heasarc.py tmp/${pidnm}pccs143.$nnn.2.csv pcc
            python3 ${HERE}/vizier2heasarc.py tmp/${pidnm}pccs217.$nnn.2.csv pcc
            python3 ${HERE}/vizier2heasarc.py tmp/${pidnm}pccs353.$nnn.2.csv pcc
            python3 ${HERE}/vizier2heasarc.py tmp/${pidnm}crates.$nnn.2.csv crates
            python3 ${HERE}/vizier2heasarc.py tmp/${pidnm}atpmn.$nnn.2.csv atpmn
            python3 ${HERE}/vizier2heasarc.py tmp/${pidnm}tgss150.$nnn.2.csv tgss150
         fi
# panstarrs is not working from vizier anymore so the format conversion need to be here
         #python3 ${HERE}/vizier2heasarc.py tmp/${pidnm}panstarrs.$nnn.2.csv panstarrs
         if [ -s tmp/${pidnm}voerror.txt ]; then
            checkvo=check
         fi
         rm -f tmp/${pidnm}vosearch.txt
         echo
      fi
      echo phase 2 completed

#the nh value in phase 2
      if [ ! -z $nhthere ]; then
         nh equinox=2000 ra=$racand dec=$deccand > tmp/${pidnm}nhvalue.txt
         nhval=`tail -1 tmp/${pidnm}nhvalue.txt |  awk '{print $7}'`
         rm -f tmp/${pidnm}nhvalue.txt
      fi
#read the phase 2 data
      if [ -f tmp/${pidnm}magictt.*.2.csv ]; then
         if [ $source == sed -o $source == lcurve ]; then
            magic=1
         else
            magic=$source
         fi
         cat tmp/${pidnm}magictt.$magic.2.csv | grep -v '#' > tmp/${pidnm}magic.$magic.2.csv #deal with MAGIC extension
         rm -f tmp/${pidnm}magictt.$magic.2.csv
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
      ls tmp/${pidnm}*.2.csv > tmp/${pidnm}catlist2.txt
      rm -f tmp/${pidnm}output2.csv
      echo
      [ -s tmp/${pidnm}catlist2.txt ] && ${BINF}/readcat tmp/${pidnm}catlist2.txt tmp/${pidnm}output2.csv $racand $deccand $sfov $nhval 0. 0. 0. 0.
      if [ ! -s tmp/${pidnm}output2.csv ]; then
         echo RA= $rar Dec= $decr radius= $sfov > tmp/${pidnm}output2.csv
         echo nH= 0.0 Error circle/elliptical= 0.  0.  0.  0.  0. >> tmp/${pidnm}output2.csv
      fi

#running the PHASE 2
      rm -f tmp/${pidnm}sed.txt
      rm -f tmp/${pidnm}error_map.txt
      [ -z $zzinput ] && zzinput=0

      ${BINF}/find_candidates2 tmp/${pidnm}output2 tmp/${pidnm}find_out_temp.txt tmp/${pidnm}Sed_temp.txt tmp/${pidnm}error_map.txt tmp/${pidnm}Sed.txt ${HERE}/catrefs.txt ${BINF} $zzinput $source | cat > tmp/${pidnm}phase2
      #cat tmp/${pidnm}phase2

###############################PLOTTING###############################
#plot the PHASE 2 error map and set default area 1 arcmin
      if [ -z $zoomin ]; then
         zoomin=60.
         echo set default second phase zoomin area to 1.0 arcmin.
      fi
      rm -f tmp/${pidnm}LC.txt
#         ${BINF}/gnomo_plot_types tmp/${pidnm}error_map.txt,tmp/${pidnm}candidates_posix.txt,tmp/${pidnm}error_map.ps/vcps,${BINF}, $racand $deccand 0. $zoomin $racand $deccand 0 0 0 0 0 0 0 $racand $deccand

#plot the SED
      echo
      if [ $runmode != l ]; then
            handle_ps tmp/${pidnm}sed.ps
      fi
      echo
      if [ $runmode != s ]; then
         ${BINF}/plot_lc tmp/${pidnm}Sed.txt tmp/${pidnm}Lc.ps/cps tmp/${pidnm}LC.txt tmp/${pidnm}LC_fermi.ps/cps ${BINF} $source
         [ -f tmp/${pidnm}LC.ps ] && handle_ps tmp/${pidnm}LC.ps
         [ -f tmp/${pidnm}LC_fermi.ps ] && handle_ps tmp/${pidnm}LC_fermi.ps
      fi
      rm -rf tmp/${pidnm}vou-aladin-error.html
      #echo python3.10 ${HERE}/aladin_error_map.py --ra $racand --dec $deccand --input_file_error_circles tmp/${pidnm}error_map.txt --out tmp/vou-aladin-error-py.html
      python3.10 ${HERE}/aladin_error_map.py --ra $racand --dec $deccand --input_file_error_circles tmp/${pidnm}error_map.txt --out tmp/vou-aladin-error-py.html
      if [ $plotsed != N ]; then
#plotting aladin map here
         echo " "
         #open tmp/${pidnm}vou-aladin-error-py.html
      fi
      [ $showsed == 'osed' ] && open tmp/${pidnm}Sed.txt
##############################################################
# Write data from additional catalogs via cat2sed.py and add to Sed.txt
      if [ -f tmp/specfindv3.1.2.txt ]; then
         python3.10 ${HERE}/cat2sed.py --input_file tmp/specfindv3.1.2.txt --output tmp/4sed.csv --refs_file  ${HERE}/catrefs.txt
         cat tmp/4sed.csv >> tmp/Sed.txt
      fi
      if [ -f tmp/mmmonitoring4.1.2.txt ]; then
         python3.10 ${HERE}/cat2sed.py --input_file tmp/mmmonitoring4.1.2.txt --output tmp/4sed.csv --refs_file  ${HERE}/catrefs.txt
         cat tmp/4sed.csv >> tmp/Sed.txt
      fi
      if [ -f tmp/eprs.1.2.txt ]; then
         python3.10 ${HERE}/cat2sed.py --input_file tmp/eprs.1.2.txt --output tmp/4sed.csv --refs_file  ${HERE}/catrefs.txt
         cat tmp/4sed.csv >> tmp/Sed.txt
      fi
      if [ -f tmp/gleamv2.1.2.txt ]; then
         python3.10 ${HERE}/cat2sed.py --input_file tmp/gleamv2.1.2.txt --output tmp/4sed.csv --refs_file  ${HERE}/catrefs.txt
         cat tmp/4sed.csv >> tmp/Sed.txt
      fi
      if [ -f tmp/ratan600.1.2.txt ]; then
         python3.10 ${HERE}/cat2sed.py --input_file tmp/ratan600.1.2.txt --output tmp/4sed.csv --refs_file  ${HERE}/catrefs.txt
         cat tmp/4sed.csv >> tmp/Sed.txt
      fi
      if [ -f tmp/paco.1.2.txt ]; then
         python3.10 ${HERE}/cat2sed.py --input_file tmp/paco.1.2.txt --output tmp/4sed.csv --refs_file  ${HERE}/catrefs.txt
         cat tmp/4sed.csv >> tmp/Sed.txt
      fi
      if [ -f tmp/pacopccs.1.2.txt ]; then
         python3.10 ${HERE}/cat2sed.py --input_file tmp/pacopccs.1.2.txt --output tmp/4sed.csv --refs_file  ${HERE}/catrefs.txt
         cat tmp/4sed.csv >> tmp/Sed.txt
      fi
      if [ -f tmp/sptsz.1.2.txt ]; then
         python3.10 ${HERE}/cat2sed.py --input_file tmp/sptsz.1.2.txt --output tmp/4sed.csv --refs_file  ${HERE}/catrefs.txt
         cat tmp/4sed.csv >> tmp/Sed.txt
      fi
      rm -f matching_rows.csv
      echo "Searching for Swift SUFST data"
      python3.10 ${HERE}/getSUFST.py $racand $deccand 10
      cat matching_rows.csv >> tmp/Sed.txt
      echo "Searching for SMARTS data"
      python3.10 ${HERE}/getSMARTS.py $racand $deccand 10
      cat matching_rows.csv >> tmp/Sed.txt
      echo "Searching for Fermi LC data"
      python3.10 ${HERE}/getFermi.py $racand $deccand 10
      cat Fermi_matching_rows.csv >> tmp/Sed.txt

#for SED builder tool and future catalog
      rm -f tmp/${pidnm}Out4SedTool.txt
      rased=`echo $racand | sed 's/\./_/g'`
      decsed=`echo $deccand | sed 's/\./_/g' | sed 's/-/m/g'`
      python3.10 ${HERE}/convert_sed.py tmp/${pidnm}Sed.txt tmp/${pidnm}Out4SedTool.txt tmp/${pidnm}Sed.csv
      if [ $runmode != l -a $plotsed != N ]; then
         rm -f tmp/${pidnm}PySED.png
         python3.10 ${HERE}/All-SED.py --xaxis f --infile1 tmp/${pidnm}Sed.csv --outfile tmp/${pidnm}PySED.png --title 'Source nr. '$source' RA='$racand' Dec='$deccand --upperlimits 'no' --erosita 'yes' --Show 'no' --fermi 'yes'
#Plotting SED here
         #open tmp/${pidnm}PySED.png
      fi
##############################################################

#PID number for online version
      if [ $pid ]; then
         ###echo PID number = $pid
         [ -f tmp/${pidnm}sed.eps ] && mv tmp/${pidnm}sed.eps tmp/${pidnm}${source}"_"sed.eps
         [ -f tmp/${pidnm}error_map.eps ] && mv tmp/${pidnm}error_map.eps tmp/${pidnm}${source}"_"error_map.eps
         [ -f tmp/${pidnm}LC.eps ] && mv tmp/${pidnm}LC.eps tmp/${pidnm}${source}"_"LC.eps
         [ -f tmp/${pidnm}LC_fermi.eps ] && mv tmp/${pidnm}LC_fermi.eps tmp/${pidnm}${source}"_"LC_fermi.eps
         [ -f tmp/${pidnm}error_map.txt ] && mv tmp/${pidnm}error_map.txt tmp/${pidnm}${source}"_"error_map.txt
         [ -f tmp/${pidnm}Sed.txt ] && mv tmp/${pidnm}Sed.txt tmp/${pidnm}${source}"_"sed.txt
         [ -f tmp/${pidnm}Out4SedTool.txt ] && mv tmp/${pidnm}Out4SedTool.txt tmp/${pidnm}${source}"_"Out4SedTool.txt
      fi
   fi
done

if [ -f tmp/catalog_error.txt ]; then 
   python3.10 ${HERE}/printMissedCatalogs.py tmp/catalog_error.txt ${HERE}/cats1.ini ${HERE}/cats2.ini
else
   echo "All catalogs have been queried regularly"
fi
#remove the vo files from various catalogs
rm -rf eada_files
#rm -f tmp/${pidnm}*.1.csv
rm -f tmp/${pidnm}*.i.csv
#rm -f tmp/${pidnm}*.2.csv
rm -f tmp/${pidnm}xrtdeep.csv
###rm -f Results/$xrtnm/*.pdf
rm -f vou-blazars.aux
rm -f vou-blazars.log
unset zoomin zzinput
unset posra posdec sfov nhval r1 emaj1 emin1 posa1 r2 emaj2 emin2 posa2 runmode
unset VOUB_AUTOSED plotsed allcatalog
