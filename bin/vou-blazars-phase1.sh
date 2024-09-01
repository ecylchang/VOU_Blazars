#!/bin/bash

HERE=$(cd `dirname $BASH_SOURCE`; pwd) ##point to bin file and fortran file
# Directory where the fortran binaries are:
BINF="${HERE}/fort"
#
SITE='standard'
MinFermiSearchRadius='10'
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
    if [ $ploterrorregion != N ]; then
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
   echo "  --FOV   : Searching Radius (in ARC-MINUTES) around RA,DEC to search for observations"
   echo ""
   echo " OPTIONS:"
   echo " --mode    : Running mode"
   echo "       Options are 'f' find candidate mode (default): finding interesing candidates within a specified region"
   echo "                   's' SED mode: obtaining SED for a specified source with given R.A. Dec."
   echo "                   'l' Light curve mode: obtaining light curve for a specified source with given R.A. Dec."
   echo ""
   echo "             (If the user has installed Heasoft and did not specify the nh, it will use the value calculated by Heasoft)"
   echo " --nh      : nH column density (in cm^2). Default is 5.e20 cm^2."
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
   echo ""
}


# If no arguments given, print Help and exit.
[ "${#@}" -eq 0 ] && { help; exit 0; }

##############################################################
# Initial set up, read the parameters input
##############################################################

#nhthere=`which nh`
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
         ploterrorregion=$2; shift;;
      --allcats)
         allcatalog=$2; shift;;
      --light)
        light=$2; shift;;
      --optical)
        optical=$2; shift;;
      --radio)
        radio=$2; shift;;
      --xray)
        xray=$2; shift;;
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

echo Running VOU-BlazarsHybrid Phase1 V1.27

echo
if [ "$optical" == "" ]; then
   #optical='Y'
   optical='NI'
fi
optical=$(echo "$optical" | tr '[:lower:]' '[:upper:]')
light=$(echo "$light" | tr '[:lower:]' '[:upper:]')

ifov="${sfov%%.*}"
if [[ "$ifov" -gt 11 && ( "$optical" == "Y" || "$optical" == "NI" ) ]]; then
   echo "================================================================================"
   echo "The requested FOV is larger than 12.0 arcminutes"
   echo "Flag to enable the retrieval of optical sources will be set to 'NO' to avoid dealing with too many objects"
   echo "================================================================================"
   optical="NO"
fi

if (( $(echo "$posdec < -90" | bc -l) || $(echo "$posdec > 90" | bc -l) )); then
   echo "Declination outside the allowed range (-90/+90)"
   exit 1
fi
#read the input parameters 

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

if [ -z $nhval ]; then
   nhval=`python3.10 ${HERE}/nh.py --ra $ranh --dec $decnh `
   echo "nH value from HEASARC "$nhval
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
[ -z $ploterrorregion ] && ploterrorregion=`grep 'PLOTERRORREGION' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $allcatalog ] && allcatalog=`grep 'ALLCATS' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $light ] && light=`grep 'LIGHT' ${HERE}/config_vou.txt | awk '{print $2}'`
[ -z $plotlab ] && plotlab=`grep 'LEGEND' ${HERE}/config_vou.txt | awk '{print $2}'`
runmode='m'

##############################################################
# FIRST PHASE
# aim: find the candidates from available radio and X-ray sources
##############################################################

rm -f tmp/${pidnm}*.1.csv
rm -f tmp/${pidnm}*.i.csv
rm -f tmp/${pidnm}*.2.csv
rm -f tmp/${pidnm}*_out.csv
rm -f tmp/${pidnm}vosearch.txt
rm -f tmp/gammacandidates.csv
rm -f tmp/Sed_temp.txt 
rm -f tmp/candidates_int.txt
rm -f tmp/candidates.csv
rm -f tmp/phase1_candidates*.csv
rm -f tmp/phase1_find_out.txt
rm -f tmp/find_out_temp.txt
rm -f tmp/output_int.csv
rm -f tmp/output1.csv
rm -f tmp/output2.csv
rm -f tmp/hstgsc_out.txt
rm -f tmp/catlist_int.txt
rm -f tmp/catalog_error.txt
rm -f tmp/slopes.csv
rm -f tmp/viz-wise.csv
rm -f *.npy
if [ $light != y -a $light != Y ]; then
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog NVSS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}nvss.1.csv  > tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog FIRST --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}first.1.csv  >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog VLASSQL --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}vlassql.1.csv  >> tmp/${pidnm}vosearch.txt

if (( $(echo "$decnh < -39" | bc -l) )); then
      echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog SUMSS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}sumss.1.csv >> tmp/${pidnm}vosearch.txt
      echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog RACS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}racs.1.csv  >> tmp/${pidnm}vosearch.txt
fi 

   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog 2SXPS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}2sxps.1.csv  >> tmp/${pidnm}vosearch.txt
#   echo conesearch --db ${HERE}/cats1.ini --catalog 1OUSX --ra $ranh --dec $decnh --radius $sfov --runit arcmin --columns default -o tmp/${pidnm}1ousx.1.csv >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py  --db ${HERE}/cats1.ini --catalog RASS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}rass.1.csv  >> tmp/${pidnm}vosearch.txt


if [ $SITE == 'alternative' ]; then
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog 4XMM-DR12 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}4xmmdr13.1.csv >> tmp/${pidnm}vosearch.txt
else
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog 4XMM-DR13 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}4xmmdr13.1.csv  >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog BMW --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}bmw.1.csv  >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog XMMSL2 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}xmmsl2.1.csv  >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog IPC2E --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}ipc.1.csv  >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py  --db ${HERE}/cats1.ini --catalog IPCSL --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}ipcsl.1.csv  >> tmp/${pidnm}vosearch.txt
fi
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog WGACAT --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}wgacat.1.csv  >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog Chandra-CSC2 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}chandracsc2.1.csv  >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog eRASS1 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}erass1.1.csv  >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog eRASS1-S --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}erass1-s.1.csv  >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog eFEDS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}efeds.1.csv  >> tmp/${pidnm}vosearch.txt
fi
if [ $light != y -a $light != Y ]; then
     echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog ZWCLUSTERS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}zw.1.csv  >> tmp/${pidnm}vosearch.txt
     echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog PSZ2 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}psz2.1.csv  >> tmp/${pidnm}vosearch.txt
     echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog MCXC --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}mcxc.1.csv  >> tmp/${pidnm}vosearch.txt
      if [[ $SITE == 'standard' ]]; then
        echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog ABELL --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}abell.1.csv  >> tmp/${pidnm}vosearch.txt
         echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog SWXCS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}swxcs.1.csv  >> tmp/${pidnm}vosearch.txt
      fi
     echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog F2PSR --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}f2psr.1.csv  >> tmp/${pidnm}vosearch.txt
#     echo conesearch --db ${HERE}/cats1.ini --catalog F357cat --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --columns default -o tmp/${pidnm}f357cat.1.csv >> tmp/${pidnm}vosearch.txt
  if [ "$optical" == "Y" ] || [ "$optical" == "NI" ]; then
     echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog MilliQuas --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}mquas.1.csv  >> tmp/${pidnm}vosearch.txt
     echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog WISE-VIZIER --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}viz-wise.csv >> tmp/${pidnm}vosearch.txt
     echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog UNWISE --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}unwise.i.csv >> tmp/${pidnm}vosearch.txt
     echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog HSTGSC --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/hstgsc_out.txt  >> tmp/${pidnm}vosearch.txt 
     echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog GAIA2 --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}gaia2.i.csv >> tmp/${pidnm}vosearch.txt
    if [ "$ifov" -le 12 ]; then
#       echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog PanSTARRS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/panstarrs.i.csv --timeout 30 
       echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog PanSTARRS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/panstarrs.i.csv --timeout 30 >> tmp/${pidnm}vosearch.txt
    fi
    if [ "$ifov" -le 5 ]; then
       echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog SDSS --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}sdss_out.txt >> tmp/${pidnm}vosearch.txt
    fi
  fi
fi
echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog CRATES --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}crates.1.csv  >> tmp/${pidnm}vosearch.txt
echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog 5BZCat --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}5bzcat.1.csv  >> tmp/${pidnm}vosearch.txt
echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog 6DF --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}6df.1.csv  >> tmp/${pidnm}vosearch.txt
#echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog HD --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}hd.c.csv  >> tmp/${pidnm}vosearch.txt
echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog SAO --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}sao.c.csv  >> tmp/${pidnm}vosearch.txt
echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog PULSAR --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}pulsar.1.csv  >> tmp/${pidnm}vosearch.txt
echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog 3HSP --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}3hsp.1.csv  >> tmp/${pidnm}vosearch.txt
echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog 3FHL --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}3fhl.1.csv  >> tmp/${pidnm}vosearch.txt
if [ "$ifov" -le "$MinFermiSearchRadius" ]; then
  sfovFermi="$MinFermiSearchRadius"
else
  sfovFermi="$sfov"
fi
echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog 4LAC-DR3 --ra $ranh --dec $decnh --radius $sfovFermi --runit arcmin --o tmp/${pidnm}4lacdr3.1.csv  >> tmp/${pidnm}vosearch.txt
echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog 4FGL-DR3 --ra $ranh --dec $decnh --radius $sfovFermi --runit arcmin --o tmp/${pidnm}4fgldr3.1.csv  >> tmp/${pidnm}vosearch.txt
echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog 4FGL-DR4 --ra $ranh --dec $decnh --radius $sfovFermi --runit arcmin --o tmp/${pidnm}4fgldr4.1.csv  >> tmp/${pidnm}vosearch.txt
echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog 2AGILE-VIZIER --ra $ranh --dec $decnh --radius $sfovFermi --runit arcmin --o tmp/${pidnm}2agile-vizier.1.csv  >> tmp/${pidnm}vosearch.txt
echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats2.ini --catalog 1FLE --ra $ranh  --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}fmev.1.csv >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog CVCAT --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}cvcat_out.csv  >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog XRBCAT --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}xrbcat_out.csv  >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog SNRGREEN --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}snrgreen_out.csv  >> tmp/${pidnm}vosearch.txt
#   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog MWMC --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}mwmc_out.csv  
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog MWMC --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}mwmc_out.csv  >> tmp/${pidnm}vosearch.txt
   echo python3.10 ${HERE}/conesearch2.py --db ${HERE}/cats1.ini --catalog MWSC --ra $ranh --dec $decnh --radius $sfov --runit arcmin --o tmp/${pidnm}mwsc_out.csv  >> tmp/${pidnm}vosearch.txt
#fi

rm -f tmp/${pidnm}voerror.txt
bash ${HERE}/queue.sh -f tmp/${pidnm}vosearch.txt -p $pidnm #2> voerror #!2>&1
declare -i irunvo=1
if [ -s tmp/${pidnm}voerror.txt ]; then
   voerror=yes
else
   voerror=no
fi
#until [ $voerror == no -o $irunvo -ge 4 ]
until [ $voerror == no -o $irunvo -ge 3 ]
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
   echo "Warning some catalogs were not downloaded - main list" > tmp/catalog_error.txt
   cat tmp/voerror.txt | awk -F'--catalog ' '{print $2}' | awk '{print $1}' >> tmp/catalog_error.txt
fi
echo

if [ -s tmp/${pidnm}erass1.1.csv ]; then
  if [ -s tmp/${pidnm}erass1-s.1.csv ]; then
     cat tmp/${pidnm}erass1-s.1.csv | grep -v IAUName >> tmp/erass1.1.csv
  fi
else
  if [ -s tmp/${pidnm}erass1-s.1.csv ]; then
     cp tmp/${pidnm}erass1-s.1.csv  tmp/erass1.1.csv
  fi
fi

echo  phase 1 completed
noOfCats=`ls tmp/${pidnm}*.1.csv 2>/dev/null | wc -l`
if [ $noOfCats == 0 ] && [ $runmode != s ]; then
  echo "There are no blazar candidates in this field"
#  exit 0;
fi
if [ -s tmp/${pidnm}voerror.txt ]; then
   checkvo=check
fi
#rm -f tmp/${pidnm}vosearch.txt
#read the data ensuring that racs catalog is the last one
ls tmp/${pidnm}*.1.csv  > tmp/${pidnm}catlist1.txt

if [ -s tmp/${pidnm}hstgsc_out.txt ]; then
   ls tmp/${pidnm}hstgsc_out.txt >> tmp/${pidnm}catlist1.txt
fi
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
if [ -s tmp/vlassql.1.csv ]; then
   grep -v "0.0000,0.0000,0.0000,0.000,0.000" tmp/vlassql.1.csv > tmp/1.txt
   mv tmp/1.txt tmp/vlassql.1.csv
fi
echo
[ -s tmp/${pidnm}catlist1.txt ] && ${BINF}/readcat tmp/${pidnm}catlist1.txt tmp/${pidnm}output1.csv $ranh $decnh $sfov $nhval $r1 $emaj1 $emin1 $posa1
#[ -s tmp/${pidnm}catlist1.txt ] && ${BINF}/readcat tmp/${pidnm}catlist1.txt tmp/${pidnm}output1.csv $ranh $decnh $sfov 3.e20 $r1 $emaj1 $emin1 $posa1
rm -f tmp/${pidnm}*temp.txt
cp tmp/output1.csv tmp/output1.csv.save
python3.10 ${HERE}/remove_old_unnecessary_xraypoints.py 
if [ $runmode == s -o $runmode == l ]; then
   ${BINF}/find_candidates1 tmp/${pidnm}output1 tmp/${pidnm}find_out_temp.txt tmp/${pidnm}RX_temp.txt tmp/${pidnm}Sed_temp.txt tmp/${pidnm}no_matched_temp.txt ${BINF} 0 | cat > tmp/${pidnm}phase1
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
#if [ -s tmp/${pidnm}no_matched_temp.txt ]; then
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
      #if [ -z ${VOUB_AUTOSED} ]; then
      #   echo
      #   echo Warning! Large searching radius in Intermediate phase for conesearch.
      #   [ -z $pid ] && read -p "Type 'y' to continue! Otherwise, stop processing Intermediate phase." process #rcut
      #   [ -z $process ] && process=n
      #   [ $process != y -a $process != Y ] && echo Not searching sources in Intermediate phase.
      #fi
      process=n
   fi

#running conesearch for Intermediate phase
   if [ $process == y -o $process == Y ]; then
      #rm -f tmp/${pidnm}*.i.csv
      ln=`cat tmp/${pidnm}no_matched_temp.txt | wc -l`
      #if [ $ln -gt 1000 ]; then
      if [ $ln -gt 3000 ]; then
         echo "Too many candidates, stopping intermediate phase"
         exit
      fi
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
      until [ $voerror == no -o $irunvo -ge 4 ]
#      until [ $voerror == no -o $irunvo -ge 3 ]
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
          echo "Warning some catalogs were not downloaded - second list " >> tmp/catalog_error.txt
          cat tmp/voerror.txt | awk -F'--catalog ' '{print $2}' | awk '{print $1}' >> tmp/catalog_error.txt
      fi
      echo
      echo intermediate phase completed
      if [ -s tmp/${pidnm}voerror.txt ]; then
         checkvo=check
      fi
      rm -f tmp/${pidnm}vosearch.txt

#read the Intermediate phase data
      shopt -s nullglob
      files=(tmp/*.i.csv)
      if [ ${#files[@]} -gt 0 ]; then
          ls "${files[@]}" > tmp/${pidnm}catlist_int.txt
          #ls "${files[@]}" 
      fi
      #ls tmp/${pidnm}*.i.csv > tmp/${pidnm}catlist_int.txt
      if [ -s tmp/hstgsc_out.txt ]; then
         echo tmp/hstgsc_out.txt >> tmp/${pidnm}catlist_int.txt
      fi
      if [ -s tmp/${pidnm}first.1.csv ]; then
         ls tmp/${pidnm}first.1.csv >> tmp/${pidnm}catlist_int.txt
      fi
      if [ -s tmp/${pidnm}sumss.1.csv ]; then
         ls tmp/${pidnm}sumss.1.csv >> tmp/${pidnm}catlist_int.txt
      fi
      if [ -s tmp/${pidnm}nvss.1.csv ]; then
         ls tmp/${pidnm}nvss.1.csv >> tmp/${pidnm}catlist_int.txt
      fi
      if [ -s tmp/${pidnm}vlassql.1.csv ]; then
         ls tmp/${pidnm}vlassql.1.csv >> tmp/${pidnm}catlist_int.txt
      fi
      if [ -s tmp/${pidnm}racs.1.csv ]; then
         ls tmp/${pidnm}racs.1.csv >> tmp/${pidnm}catlist_int.txt
      fi
      rm -f tmp/${pidnm}output_int.csv
      echo
      #[ -s tmp/${pidnm}catlist_int.txt ] && ${BINF}/readcat tmp/${pidnm}catlist_int.txt tmp/${pidnm}output_int.csv $ranh $decnh $sfov 3.e20 $r1 $emaj1 $emin1 $posa1
      [ -s tmp/${pidnm}catlist_int.txt ] && ${BINF}/readcat tmp/${pidnm}catlist_int.txt tmp/${pidnm}output_int.csv $ranh $decnh $sfov $nhval $r1 $emaj1 $emin1 $posa1

#run the Intermediate phase
      #cat tmp/Sed_temp.txt
      python3.10 ${HERE}/get_radio_xray_matches.py tmp/Sed_temp.txt >> tmp/no_matched_temp.txt
      echo
      ${BINF}/find_candidates_int tmp/${pidnm}output_int tmp/${pidnm}no_matched_temp.txt tmp/${pidnm}find_out_temp.txt tmp/${pidnm}Intermediate_out.txt ${BINF}
   fi
if [ -s tmp/${pidnm}Intermediate_out.txt ]; then
   nint=`tail -1 tmp/${pidnm}Intermediate_out.txt | awk '{print $1}' `
   head -$nint tmp/${pidnm}Intermediate_out.txt >> tmp/${pidnm}find_out_temp.txt
   cat tmp/${pidnm}Intermediate_out.txt  >> tmp/${pidnm}Sed_temp.txt
fi
#cat tmp/${pidnm}find_out_temp.txt
echo "Number,R.A.,Dec.,Type Name" > tmp/${pidnm}candidates.csv
cat tmp/${pidnm}find_out_temp.txt | awk ' $3>10000 {print $1, $2, int($3/10000), $4 $5}  $3<-40000 {print $1, $2, int($3/10000), $4 $5} $3==-9999 {print $1, $2, int($3/10000), $4 $5} ' | awk '{print NR, $1, $2, "type"$3, $4 $5}' >> tmp/${pidnm}candidates.csv
#echo "---"
#cat tmp/${pidnm}candidates.csv
#echo "---"
cat tmp/${pidnm}candidates.csv | sed 's/type1/HBL/g' | sed 's/type2/IBL/g' | sed 's/type3/LBL/g' | sed 's/type4/radio-AGN/g' | sed 's/type5/Unknown/g' | sed 's/type-5/3HSP/g' | sed 's/type-6/5BZCat/g' | sed 's/type-7/CRATES/g' | sed 's/type0/Pulsar/g' > tmp/${pidnm}candidates.csv 

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
   
#sed -i '' 's/\ /\,/g' tmp/${pidnm}candidates.csv
sed 's/ /,/g' tmp/${pidnm}candidates.csv > tmp/candidates_comma.csv
mv tmp/candidates_comma.csv tmp/candidates.csv

#######################################################
if [ "$optical" == "Y" ] || [ "$optical" == "NI" ]; then
    # add optical sources and merge the results
    # option NI is for the case where optical data is used to find counterparts but it is not plotted
    python3.10 "${HERE}/complete_phase1.py" "$ranh" "$decnh" "$optical"
    if [ -s tmp/${pidnm}phase1_find_out.txt ]; then
       if [ -f "tmp/${pidnm}find_out_temp.txt" ]; then
          cat "tmp/phase1_find_out.txt" >> "tmp/${pidnm}find_out_temp.txt"
       else
          cat "tmp/phase1_find_out.txt" > "tmp/${pidnm}find_out_temp.txt"
       fi
    fi
    if [ -f "tmp/find_out_temp.txt" ]; then
       python3.10 "${HERE}/remove_duplicates.py"  tmp/find_out_temp.txt
    fi
    if [ -f "tmp/phase1_candidates.csv" ]; then
       python3.10 "${HERE}/remove_duplicates.py"  tmp/phase1_candidates.csv
    fi
    python3.10 "${HERE}/remove_no_opt_cntp_sources.py"
else
    if [ -s tmp/candidates.csv ]; then
        cat "tmp/candidates.csv" > "tmp/phase1_candidates.csv"
    fi
fi
if [ "$radio" == "y" ] || [ "$radio" == "Y" ]; then
   if [ -f "tmp/${pidnm}RX_sorted.txt" ]; then
        cat "tmp/${pidnm}RX_sorted.txt" | grep "\-9" | sed 's/$/ test/' >> "tmp/${pidnm}find_out_temp.txt"
   else 
        echo "File with radio sources not found"
   fi 
fi
if [ "$xray" == "y" ] || [ "$xray" == "Y" ]; then
   if [ -f "tmp/${pidnm}RX_sorted.txt" ]; then
        cat "tmp/${pidnm}RX_sorted.txt" | grep "\-8" | sed 's/$/ test/' >> "tmp/${pidnm}find_out_temp.txt"
   else 
        echo "File with radio sources not found"
   fi 
fi
#######################################################
# python plotting of the candidates map
#######################################################
rrad=0
rm -f candidates.png
python3.10 "${HERE}/complete_known_sources.py"  
if [ ! -e tmp/phase1_candidates.csv ]; then
    echo "Number,R.A.,Dec.,Type,Name" > tmp/phase1_candidates.csv
fi
python3.10 "${HERE}/remove_nearby_duplicates.py" tmp/phase1_candidates.csv $sfov
python3.10 "${HERE}/add_catalog_id.py"
cat tmp/phase1_candidates_final.csv
plot_opt='y'
if [ -e tmp/phase1_candidates.csv ] && [ -e tmp/find_out_temp.txt ]; then
   if [ $emaj1  != "0." ]; then 
        #echo python3.10 ${HERE}/skymap.py --ra $posra --dec $posdec --FOV $sfov --major $emaj1  --minor $emin1 --angle $posa1 --major2 $emaj2  --minor2 $emin2 --angle2 $posa2 --infile_candidates tmp/phase1_candidates_clean.csv --infile_multi tmp/find_out_temp.txt --infile_local IHBL-catalog --out candidates.png --plot_star $plot_opt
        python3.10 ${HERE}/skymap.py --ra $posra --dec $posdec --FOV $sfov --major $emaj1  --minor $emin1 --angle $posa1 --major2 $emaj2  --minor2 $emin2 --angle2 $posa2 --infile_candidates tmp/phase1_candidates_clean.csv --infile_multi tmp/find_out_temp.txt --infile_local IHBL-catalog --out candidates.png --plot_star $plot_opt
   elif [ $r1 != "0." ]; then
        python3.10 ${HERE}/skymap.py --ra $posra --dec $posdec --FOV $sfov --radius $r1 --infile_candidates tmp/phase1_candidates_clean.csv --infile_multi tmp/find_out_temp.txt --infile_local IHBL-catalog --out candidates.png --plot_star $plot_opt
   else  
        python3.10 ${HERE}/skymap.py --ra $posra --dec $posdec --FOV $sfov --infile_candidates tmp/phase1_candidates_clean.csv --infile_multi tmp/find_out_temp.txt --infile_local IHBL-catalog --out candidates.png --plot_star $plot_opt
   fi
else
   python3.10 ${HERE}/skymap.py --ra $posra --dec $posdec --FOV $sfov --radius $r1 --infile_candidates tmp/phase1_candidates_clean.csv --infile_multi tmp/find_out_temp.txt --infile_local IHBL-catalog --out candidates.png 
fi
if [ $emaj1  != "0." ]; then 
   python3.10 ${HERE}/aladin_error_region.py --ra $posra  --dec $posdec --FOV $sfov --maj_axis $emaj1 --min_axis $emin1 --angle $posa1 --out tmp/vou-aladin-error-region.html --plotErrorEllipse $ploterrorregion
elif [ $r1 != "0." ]; then
   python3.10 ${HERE}/aladin_error_region.py --ra $posra  --dec $posdec --FOV $sfov --radius $r1 --out tmp/vou-aladin-error-region.html --plotErrorEllipse $ploterrorregion
fi
open tmp/vou-aladin-error-region.html 
open candidates.png
