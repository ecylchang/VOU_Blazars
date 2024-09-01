import os
import csv
import argparse
import re
import math

path ='tmp/'
input_file_4fgl = path+'4fgldr4.1.csv'
input_file_4lac = path+'4lacdr3.1.csv'
input_file_efeds = path+'efeds.1.csv'
input_file_erass1 = path+'erass1.1.csv'
input_file_2sxps = path+'2sxps.1.csv'
input_file_2agile = path+'2agile-vizier.1.csv'
input_file_ousx = path+'1ousx.1.csv'
input_file_4xmmdr13 = path+'4xmmdr13.1.csv'
input_file_chandra = path+'chandracsc2.1.csv'
input_file_xmmslew = path+'xmmsl2.1.csv'
input_file_ipc = path+'ipc.1.csv'
input_file_rass = path+'rass.1.csv'
input_file_wga = path+'wgacat.1.csv'
input_file_nvss = path+'nvss.1.csv'
input_file_racs = path+'racs.1.csv'
input_file_first = path+'first.1.csv'
input_file_vlass = path+'vlassql.1.csv'
input_file_icecube = '/Users/paologiommi/cats4topcat/IceCubeEllipses.txt'

input_file_candidates = path+'phase1_candidates_final.csv'

parser = argparse.ArgumentParser(description='FirmamentoAladin')
parser.add_argument('--ra', type=float, help='R.A. (degrees d/f=0.0)', default='0.0')
parser.add_argument('--dec', type=float, help='Dec. (degrees d/f=0.0)', default='0.0')
parser.add_argument('--radius', type=float, help='radius (arcmin d/f=0.0)', default='0.0')
parser.add_argument('--FOV', type=float, help='radius of searched region (arcmin d/f=0.0)', default='0.0')
parser.add_argument('--maj_axis', type=float, help='error ellipse major axis (arcmin d/f=0.0)', default='0.0')
parser.add_argument('--min_axis', type=float, help='error ellipse minor axis (arcmin d/f=0.0)', default='0.0')
parser.add_argument('--angle', type=float, help='error ellipse ange (arcmin d/f=0.0)', default='0.0')
parser.add_argument('--plotErrorEllipse', type=str, help='plot error ellipse? (d/f=no)', default='no')
parser.add_argument('--output_file', type=str, help='output file d/f=vou-aladin-py.html)', default='vou-aladin-py.html')

def sky_distance(ra1, dec1, ra2, dec2):
    # Convert coordinates to radians
    ra1 = math.radians(ra1)
    dec1 = math.radians(dec1)
    ra2 = math.radians(ra2)
    dec2 = math.radians(dec2)

    # Haversine formula for sky distance
    d_ra = ra2 - ra1
    d_dec = dec2 - dec1
    a = math.sin(d_dec/2)**2 + math.cos(dec1) * math.cos(dec2) * math.sin(d_ra/2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a)) * 57.32
    #57.32 converts radians to deg
    return c

def ra_to_deg(ra_str):
    # Extract the hours, minutes, and seconds from the RA string using regular expressions
    match = re.search(r"(\d+) (\d+) (\d+)", ra_str)
    if not match:
        raise ValueError("Invalid RA string format")
    hours = int(match.group(1))
    minutes = int(match.group(2))
    seconds = int(match.group(3))
    
    # Convert the hours, minutes, and seconds to decimal degrees
    deg = 15.0 * (hours + minutes/60.0 + seconds/3600.0)
    return deg

def dec_to_deg(dec_str):
    # Extract the degrees, minutes, and seconds from the Dec string using regular expressions
    match = re.search(r"([+-]?)(\d+) (\d+) (\d+)", dec_str)
    if not match:
        raise ValueError("Invalid Dec string format")
    sign = -1.0 if match.group(1) == "-" else 1.0
    degrees = int(match.group(2))
    minutes = int(match.group(3))
    seconds = int(match.group(4))
    
    # Convert the degrees, minutes, and seconds to decimal degrees
    deg = sign * (degrees + minutes/60.0 + seconds/3600.0)
    return deg


ok = True
a1 = "'"
a2 = "}"
a0 = "{"

args = parser.parse_args()
ra = args.ra
dec = args.dec 
ra_center = ra
dec_center = dec
plot_error_ellipse = args.plotErrorEllipse
fov = args.FOV/60
ell_maj = args.maj_axis/60
ell_min = args.min_axis/60
ell_angle = -args.angle-90
field_radius = args.radius/60
output_file =args.output_file
sra = (format(ra,"10.4f"))
sdec = (format(dec,"10.4f"))
open(output_file, 'w').close()
with open(output_file, 'w') as of:
   of.write('<HTML style="background-color:#ffffff">\n')
   of.write('<style>\n')
   of.write('sm {\n')
   of.write('  font-size: 10px;\n')
   of.write(' }\n')
   of.write('.title {\n')
   of.write('  font-size: 13px;\n')
   of.write(' }\n')
   of.write('body {\n')
   of.write(' font-family: Tahoma, sans-serif;\n')
   of.write(' margin:0;\n')
   of.write(' padding:0;\n')
   of.write('}\n')
   of.write('</style>\n')
   of.write('<HEAD>\n')
   of.write('   <script type="text/javascript" src="https://code.jquery.com/jquery-1.9.1.min.js" charset="utf-8"></script>\n')
   of.write('   <link rel="stylesheet" href="https://aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.css" >\n')
   #of.write('   <link rel="stylesheet" href="https://aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.css" >\n')
   of.write('   <script type="text/javascript">var jqMenu = jQuery.noConflict();</script>\n')
   of.write('   <script type="text/javascript">\n')
   of.write('   var hipsDir=null;</script>\n')
   of.write('</HEAD>\n')

   of.write('&emsp;&emsp;&emsp;&emsp;&emsp;\n')
   #of.write('<script type="text/javascript"  src="https://aladin.cds.unistra.fr/AladinLite/api/v3/latest/aladin.js" charset="utf-8"></script>\n')
   of.write('<script type="text/javascript" src="https://aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.js" charset="utf-8"></script>\n')
   of.write("<legend class='title'>SED error regions map for<sub><sub><sub>\n")
   of.write('<img class="logoImg" src="https://vouwebsite.web.app/blackLogo.png" alt="LOGO" width="100" height="30"></sub></sub></sub>&emsp;Based on "Aladin Lite" developed at CDS, Strasbourg Observatory, France<br></legend>\n')
   of.write('<table><tr valign="top"><td>\n')
   of.write('<div id="aladin-lite-div" style="width:80vw; height:100vh;"></div>\n')
   of.write('</td><td style="background-color: #ededed"><sm> <div style="width:18vw;">\n')
   of.write('<b><center>Radio/IR surveys</center></b><br>\n')
   if dec < -40: 
      of.write('<input id="SUMSS" type="radio"   name="survey" value="SUMSS"><label for="SUMSS">SUMSS</label><br>\n')
      if dec > -53: 
         of.write('<input id="TGSSADR" type="radio"   name="survey" value="TGSSADR"><label for="TGSSADR">TGSSADR (150 MHz)</label><br>\n')
   else:
      of.write('<input id="NVSS" type="radio"   name="survey" value="NVSS"><label for="NVSS">NVSS (1.4 GHz)</label><br>\n' )
      of.write('--&nbsp;&nbsp;&nbsp;<b>VLASS</b> (3GHz)&nbsp;&nbsp;&nbsp;--<br>\n')
      of.write('<input id="VLASS_Epoch1" type="radio"   name="survey" value="VLASS_Epoch1"><label for="VLASS_Epoch1">Epoch 1</label>\n' )
      of.write('<input id="VLASS2.2" type="radio"   name="survey" value="VLASS2.2"><label for="VLASS2.2">2.2</label>\n' )
      of.write('<input id="VLASS3.1" type="radio"   name="survey" value="VLASS3.1"><label for="VLASS3.1">3.1</label><br>\n' )
#      of.write('<input id="VLASSall" type="radio"   name="survey" value="VLASSall"><label for="VLASSall">All epochs</label><br>\n' )
      of.write('<input id="VCSS1" type="radio"   name="survey" value="VCSS1"><label for="VCSS1">VCSS1 (340 MHz)</label><br>\n' )
      of.write('<input id="TGSSADR" type="radio"   name="survey" value="TGSSADR"><label for="TGSSADR">TGSSADR (150 MHz)</label><br>\n')
   if dec < 8: 
      of.write('<input id="EMU" type="radio"   name="survey" value="EMU"><label for="EMU">EMU 943 MHz</label><br>\n')
   if dec < 49: 
      of.write('--&nbsp;&nbsp;&nbsp;<b>RACS</b>&nbsp;&nbsp;&nbsp;--<br>\n')
   if dec < 40: 
      of.write('<input id="RACS-low" type="radio"   name="survey" value="RACS-low"><label for="RACS-low">Low-887 MHz</label>\n')
   if dec < 49: 
      of.write('<input id="RACS-mid" type="radio"   name="survey" value="RACS-mid"><label for="RACS-mid">Mid-1367 MHz</label><br>\n')
   of.write('<b><center>--</center></b>\n')
#   of.write('<br><b><center>IR surveys</center></b><br>\n')
   of.write('<input id="allwise" type="radio"   name="survey" value="P/allWISE/color"><label for="allwise">AllWISE</label><br>\n')
   of.write('<input id="unwise"  type="radio"   name="survey" value="unWISE"><label for="unWISE">unWISE color</label><br>\n')

   of.write('<input id="2MASS" type="radio"   name="survey" value="P/2MASS/color"><label for="2MASS">2MASS</label><br>\n')
   of.write('<input id="Spitzer/IRAC134" type="radio"   name="survey" value="Spitzer/IRAC134"><label for="Spitzer/IRAC134">Spitzer/IRAC134</label><br>\n')
   of.write('<input id="Herschel/color" type="radio" name="survey" value="Herschel/color"><label for="Herschel/color">Herschel/color</label><br>\n')
   of.write('<br><b><center>OpticalUV surveys</center></b><br>\n')

   if dec > -30:
      of.write('<input id="DSS2" type="radio"   name="survey" value="P/DSS2/color"><label for="DSS2">DSS2</label><br>\n')
      of.write('<input id="PANSTARRS" type="radio"   name="survey" value="P/PanSTARRS/DR1/color-z-zg-g"  checked="checked"><label for="P/PanSTARRS/DR1/color">PanSTARRS-DR1</label><br>\n')
   else:
      of.write('<input id="DSS2" type="radio"   name="survey" value="P/DSS2/color" checked="checked"><label for="DSS2">DSS2</label><br>\n')
      of.write('<input id="DES-DR2" type="radio"   name="survey" value="DES/DR2/CDS_P_DES-DR2_ColorIRG"><label for="DES-DR2">DES-DR2</label><br>\n')
   if dec > -8:
      of.write('<input id="SDSS9" type="radio"   name="survey" value="P/SDSS9/color"><label for="SDSS9">SDSS9</label><br>\n')
   if ( (dec < 35 ) & ( dec > -10 ) ) & ( ( ra > 100 ) & ( ra < 280 ) ): 
      of.write('<input id="DECaLS/DR5/color"   type="radio" name="survey" value="DECaLS/DR5/color"><label for="DECaLS-DR5">DECaLS-DR5</label><br>\n')
   if ( dec > 30 ) & ( ( ra > 100 ) & (ra < 280 ) ):
      of.write('<input id="BASS/DR3/image"   type="radio" name="survey" value="BASS/DR3/image"><label for="BASS/DR3">BASS-DR3</label><br>\n')
   if dec > -28: 
      of.write('<input id="ZTF/DR7/color"   type="radio" name="survey" value="ZTF/DR7/color"><label for="ZTF/DR7/color">ZTF-DR7</label><br>\n')
   if ( dec < 2 ) & ( dec > -30 ): 
      of.write('<input id="MAMA color"   type="radio" name="survey" value="MAMA color"><label for="MAMA color">MAMA</label><br>\n')
   if dec < 0 : 
      of.write('<input id="SkyMapper" type="radio"   name="survey" value="SkyMapper"><label for="Skymapper">SkyMapper-DR1</label><br>\n')

#   of.write("<br><b><center>UV surveys</center></b><br>\n")
   of.write('<input id="GALEXGR6_7" type="radio"  name="survey" value="GALEXGR6_7"><label for="GALEXGR6_7">GALEXGR6_7</label><br>\n')
   of.write('--&nbsp;&nbsp;&nbsp;<b>Swift-UVOT</b>&nbsp;&nbsp;&nbsp;--<br>\n')
   of.write('<input id="Swift-UVOT-M2" type="radio"  name="survey" value="Swift-UVOT-M2"><label for="Swift-UVOT-M2">M2</label>\n')
   of.write('<input id="Swift-UVOT-W1" type="radio"  name="survey" value="Swift-UVOT-W1"><label for="Swift-UVOT-W1">W1</label>\n')
   of.write('<input id="Swift-UVOT-W2" type="radio"  name="survey" value="Swift-UVOT-W2"><label for="Swift-UVOT-W2">W2</label><br>\n')

   of.write('<br><b><center>X-ray/&gamma;-ray surveys</center></b><br>\n ')
   of.write('--&nbsp;&nbsp;&nbsp;<b>eRASS1</b>&nbsp;&nbsp;&nbsp;--<br>\n')
   of.write('<input id="eRosita" type="radio"  name="survey" value="eRASS1_RGB_Rate_c010"><label for="eROSITA-RGB">Color</label>\n')
   of.write('<input id="eRosita" type="radio"  name="survey" value="eRASS1_024_Rate_c010"><label for="eROSITA-024">0.2-2.3 keV</label>\n')
   of.write('<input id="eRosita" type="radio"  name="survey" value="eRASS1_023_Rate_c010"><label for="eROSITA-023">2.3-5 keV</label><br>\n')
   of.write('<input id="Swift-XRT" type="radio"  name="survey" value="Swift/XRT/OpenUniverse color"><label for="Swift-XRT">Swift-XRT</label><br>\n')
   of.write('<input id="XMM-PN" type="radio"  name="survey" value="XMM/PN/color"><label for="XMM/PN/color">XMM-PN</label><br>\n')
   of.write('<input id="Chandra" type="radio"  name="survey" value="Chandra color"><label for="Chandra color">Chandra</label><br>\n')
   of.write('<input id="RASS" type="radio"  name="survey" value="P/RASS"><label for="P/RASS">RASS</label><br>\n')
#   of.write('<input type="radio" name="survey" id="erass1_rgb" value="eRASS1"><label for="erass1_rgb">eRASS1</label><br>\n')
   of.write('<input id="Swift-BAT" type="radio"  name="survey" value="Swift-BAT"><label for="Swift-BAT">Swift-BAT(14-20keV)</label><br>\n')
   of.write('<input id="MAXI-GSC" type="radio"  name="survey" value="MAXI-GSC"><label for="MAXI-GSC">MAXI-GSC</label><br>\n')
#   of.write('<br><b><center>&gamma;-ray surveys</center></b><br> \n')
   of.write('<input id="Fermi-LAT" type="radio"  name="survey" value="P/Fermi/color"><label for="P/Fermi/color">Fermi-LAT</label><br>\n')
   of.write('</div></sm></td></tr><tr><td> <center> \n')
   of.write('<input id="CONSTELLAIONS3" type="radio"  name="survey" value="CONSTELLATIONS3"><label for="Constellations">Constellations</label>\n')
   of.write('Projection <select name="projection"><option>aitoff</option><option selected>sin</option></select>\n')
   of.write('</center></td> </tr>\n')
   of.write('</table>\n')
   of.write('<script type="text/javascript">\n')

   string = "var aladin = $.aladin(\"#aladin-lite-div\", {survey: 'P/DSS/colored', showSimbadPointerControl: true,fov:0.2,target:'"+sra.strip()+" "+sdec.strip()+"'});\n"
   of.write(string)

   of.write('hipsDir = "https://openuniverse.asi.it/SwiftXRT4OU/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("hipsDir = hipsDir+'/outSwiftHiPS_ALLupAug2020-RGB/'\n")
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('Swift/XRT/OpenUniverse color','Swift/XRT/OpenUniverse color',hipsDir,'equatorial',7,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")

   of.write('hipsDir = "https://erosita.mpe.mpg.de/dr1/erodat/static/hips/eRASS1_RGB_Rate_c010/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('eRASS1_RGB_Rate_c010','eRASS1_RGB_Rate_c010',hipsDir,'equatorial',6,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")

   of.write('hipsDir = "https://erosita.mpe.mpg.de/dr1/erodat/static/hips/eRASS1_024_Rate_c010/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('eRASS1_024_Rate_c010','eRASS1_024_Rate_c010',hipsDir,'equatorial',6,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")

   of.write('hipsDir = "https://erosita.mpe.mpg.de/dr1/erodat/static/hips/eRASS1_023_Rate_c010/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('eRASS1_023_Rate_c010','eRASS1_023_Rate_c010',hipsDir,'equatorial',6,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")

   of.write('hipsDir = "https://alaskybis.cds.unistra.fr/NVSS/intensity/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('NVSS','NVSS',hipsDir,'equatorial',5,{imgFormat:'jpg'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('cubehelix')\n")

   of.write('hipsDir = "https://www.aoc.nrao.edu/~jmarvil/emuHiPS/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('EMU','EMU',hipsDir,'equatorial',7,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")

   of.write('hipsDir = "https://archive-new.nrao.edu/vlass/HiPS/VLASS2.2/Quicklook/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('VLASS2.2','VLASS2.2',hipsDir,'equatorial',9,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('rainbow')\n")

   of.write('hipsDir = "https://archive-new.nrao.edu/vlass/HiPS/VLASS3.1/Quicklook/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('VLASS3.1','VLASS3.1',hipsDir,'equatorial',9,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('rainbow')\n")

#   of.write('hipsDir = "https://archive-new.nrao.edu/vlass/HiPS/All_VLASS/Quicklook";\n')
#   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
#   of.write("aladin.setImageSurvey(aladin.createImageSurvey('VLASSall','VLASSall',hipsDir,'equatorial',9,{imgFormat:'png'}))\n")
#   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('eosb')\n")

   of.write('hipsDir = "http://archive-new.nrao.edu/vlass/HiPS/VLASS_Epoch1/Quicklook/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('VLASS_Epoch1','VLASS_Epoch1',hipsDir,'equatorial',9,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('rainbow')\n")

   of.write('hipsDir = "https://archive-new.nrao.edu/vlass/HiPS/VCSS1/Images_vcss1_HiPS/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('VCSS1','VCSS1',hipsDir,'equatorial',8,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('rainbow')\n")

   of.write('hipsDir = "https://casda.csiro.au/hips/RACS/low/I/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('RACS-low','RACS-low',hipsDir,'equatorial',8,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('rainbow')\n")

   of.write('hipsDir = "https://casda.csiro.au/hips/RACS/mid/I/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('RACS-mid','RACS-mid',hipsDir,'equatorial',8,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('rainbow')\n")

   of.write('hipsDir = "https://alaskybis.cds.unistra.fr/DECaLS/DR5/color/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('DECaLS/DR5/color','DECaLS/DR5/color',hipsDir,'equatorial',11,{imgFormat:'jpg'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")
   of.write('hipsDir = "http://skymapper.anu.edu.au/CDS_P_skymapper-color/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('SkyMapper','SkyMapper',hipsDir,'equatorial',9,{imgFormat:'jpg'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")

   of.write('hipsDir = "https://alaskybis.cds.unistra.fr/RASS/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('P/RASS','RASS',hipsDir,'equatorial',4,{imgFormat:'jpg'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('rainbow')\n")

#   of.write('hipsDir = "https://catsweb.oas.inaf.it/RASS_P/rgb/";\n')
#   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
#   of.write("aladin.setImageSurvey(aladin.createImageSurvey('eRASS1','erass1_rgb',hipsDir,'equatorial',5,{imgFormat:'jpg'}))\n")
#   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('rainbow')\n")
#
   of.write('hipsDir = "https://erosita.mpe.mpg.de/ExposureHiPS/eRASS3/946/eROSITA_exposure_021/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('eROSITA exposure','eROSITA exposure',hipsDir,'equatorial',5,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('rainbow')\n")
   of.write('hipsDir = "https://alaskybis.cds.unistra.fr/SUMSS/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('SUMSS','SUMSS',hipsDir,'equatorial',6,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")

   of.write('hipsDir = "https://xcatdb.unistra.fr/PNColor/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('XMM/PN/color','XMM/PN/color',hipsDir,'equatorial',7,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")
   of.write('hipsDir = "https://alaskybis.cds.unistra.fr/CHANDRA/cxc.harvard.edu_P_cda_hips_allsky_rgb/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('Chandra color','Chandra color',hipsDir,'equatorial',11,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")
   of.write('hipsDir = "https://data.darts.isas.jaxa.jp/pub/HiPS/MAXI-GSC_Image/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('MAXI-GSC','MAXI-GSC',hipsDir,'equatorial',3,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")
   of.write('hipsDir = "https://data.darts.isas.jaxa.jp/pub/HiPS/SWIFT-BAT_Image/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('Swift-BAT image','Swift-BAT image',hipsDir,'equatorial',6,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")
   of.write('hipsDir = "http://cade.irap.omp.eu/documents/Ancillary/4Aladin/BAT_14_20/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('Swift-BAT','Swift-BAT',hipsDir,'galactic',3,{imgFormat:'jpg'}))\n")

   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('cubehelix')\n")
   of.write('hipsDir = "https://alaskybis.cds.unistra.fr/JAXA/JAXA_P_CONSTELLATIONS3/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('CONSTELLATIONS3','CONSTELLATIONS3',hipsDir,'equatorial',6,{imgFormat:'png'}))\n")

   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")
   of.write('hipsDir = "https://alaskybis.cds.unistra.fr/JAXA/JAXA_P_CONSTELLATIONS6/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('CONSTELLATIONS6','CONSTELLATIONS6',hipsDir,'equatorial',6,{imgFormat:'png'}))\n")

   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")
   of.write('hipsDir = "https://alaskybis.cds.unistra.fr/ZTF/DR7/CDS_P_ZTF_DR7_g/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('ZTF/DR7/color','ZTF/DR7/color',hipsDir,'equatorial',9,{imgFormat:'png'}))\n")

   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")
   of.write('hipsDir = "https://alaskybis.cds.unistra.fr/China-VO/China-VO_P_BASS_DR3_image/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('BASS/DR3/image','BASS/DR3/image',hipsDir,'equatorial',10,{imgFormat:'png'}))\n")

   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")
   of.write('hipsDir = "https://alaskybis.cds.unistra.fr/ESAC/ESAVO_P_HST_WFPC/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('HST/WFPC','HST/WFPC',hipsDir,'equatorial',13,{imgFormat:'jpg'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")
   of.write('hipsDir = "https://alaskybis.cds.unistra.fr/MAMA/CDS_P_MAMA_color_posse_srcj/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('MAMA color','MAMA color',hipsDir,'equatorial',10,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")
   of.write('hipsDir = "https://alaskybis.cds.unistra.fr/ESAC/ESAVO_P_Herschel_color/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('Herschel/color','Herschel/color',hipsDir,'equatorial',8,{imgFormat:'png'}))\n")

   of.write('hipsDir = "https://alasky.cds.unistra.fr/unWISE/color-W2-W1W2-W1/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('unWISE','unWISE',hipsDir,'equatorial',8,{imgFormat:'jpg'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('rainbow')\n")

   of.write('hipsDir = "https://alaskybis.cds.unistra.fr/ESAC/ESAVO_P_Spitzer_IRAC134-RGB-bright/";\n')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('Spitzer/IRAC134','Spitzer/IRAC134',hipsDir,'equatorial',10,{imgFormat:'png'}))\n")

   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")
   of.write('hipsDir = "http://hips.canfar.net/CDS_P_HST_wideV/"; ')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('HST_wideV','HST_wideV',hipsDir,'equatorial',13,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")
   of.write('hipsDir = "https://tgssadr.strw.leidenuniv.nl/hips/"; ')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('TGSSADR','TGSSADR',hipsDir,'equatorial',7,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")

   of.write('hipsDir = "http://alasky.cds.unistra.fr/GALEX/GALEXGR6_7_color/"; ')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('GALEXGR6_7','GALEXGR6_7',hipsDir,'equatorial',9,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('rainbow')\n")

   of.write('hipsDir = "http://alasky.cds.unistra.fr/DES/DR2/CDS_P_DES-DR2_ColorIRG/"; ')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('DES/DR2/CDS_P_DES-DR2_ColorIRG','DES-DR2',hipsDir,'equatorial',11,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")

   of.write('hipsDir = "https://skyview.gsfc.nasa.gov/surveys/uvot/cnt/UVM2/"; ')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('Swift-UVOT-M2','Swift-UVOT-M2',hipsDir,'equatorial',9,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")
   of.write('hipsDir = "https://skyview.gsfc.nasa.gov/surveys/uvot/cnt/UVW1/"; ')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('Swift-UVOT-W1','Swift-UVOT-W1',hipsDir,'equatorial',9,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")

   of.write('hipsDir = "https://skyview.gsfc.nasa.gov/surveys/uvot/cnt/UVW2/"; ')
   of.write('hipsDir = hipsDir.substring(0,hipsDir.lastIndexOf("/",hipsDir.length));\n')
   of.write("aladin.setImageSurvey(aladin.createImageSurvey('Swift-UVOT-W2','Swift-UVOT-W2',hipsDir,'equatorial',9,{imgFormat:'png'}))\n")
   of.write("aladin.getBaseImageLayer(hipsDir).getColorMap(hipsDir).update('native')\n")
   if dec > -30:
       of.write("aladin.setBaseImageLayer('P/PanSTARRS/DR1/color-z-zg-g')\n")
       of.write("aladin.setBaseImageLayer('P/PanSTARRS/DR1/color-z-zg-g')\n")
   else: 
       of.write("aladin.setBaseImageLayer('P/DSS/colored')\n")
   of.write('\n')

   of.write('function drawRotatedEllipse(centerRA, centerDec, angle, majorAxis, minorAxis, numSegments, color) {\n')
   of.write(' var angleStep = 2 * Math.PI / numSegments;\n')
   of.write(' var angleRad = (angle * Math.PI) / 180;\n')
   of.write(' var DecRad = (centerDec * Math.PI) / 180;\n')
   of.write(' var ellipsePoints = [];\n')
   of.write(' for (var i = 0; i < numSegments + 1; i++) {\n')
   of.write('     var angleCurrent = i * angleStep;\n')
   of.write('     var x =\n')
   of.write('         centerRA +\n')
   of.write('         (majorAxis * Math.cos(angleCurrent) * Math.cos(angleRad) -\n')
   of.write('         minorAxis * Math.sin(angleCurrent) * Math.sin(angleRad)) / Math.cos(DecRad) ;\n')
   of.write('     var y =\n')
   of.write('         centerDec +\n')
   of.write('         majorAxis * Math.cos(angleCurrent) * Math.sin(angleRad) +\n')
   of.write('         minorAxis * Math.sin(angleCurrent) * Math.cos(angleRad);\n')
   of.write('     ellipsePoints.push([x, y]);\n')
   of.write(' }\n')
   of.write(' overlay.add(A.polyline(ellipsePoints, { color: color, lineWidth: 1 }));\n')
   of.write('}\n')

   of.write('var cat = A.catalog();\n')
   of.write('aladin.addCatalog(cat);\n')
   of.write("aladin.on('objectHovered', function(object) {\n")
   of.write('    var msg;\n')
   of.write('    if (object) {\n')
   of.write("        msg = 'Object ' + object.data.name + ' R.A. = ' + object.ra + ', Dec. = ' + object.dec;\n")
   of.write('    }\n')
   of.write('    else {\n')
   of.write("        msg = ' Move the mouse (or click) on a source for details  ';\n")
   of.write('    }\n')
   of.write("    $('#infoDiv').html(msg);\n")
   of.write('});\n')
   of.write('// define function triggered when an object is clicked\n')
   of.write('var objClicked;\n')
   of.write("aladin.on('objectClicked', function(object) {\n")
   of.write('    var msg;\n')
   of.write('    if (object) {\n')
   of.write('        objClicked = object;\n')
   of.write('      object.select();\n')
   of.write("        msg = 'You clicked on object ' + object.data.name + ' located at R.A. = ' + object.ra + ', Dec. = ' + object.dec;\n")
   of.write('    }\n')
   of.write('    else {\n')
   of.write('        objClicked.deselect(); \n')
   of.write("        msg = 'You clicked in void';\n")
   of.write('    }\n')
   of.write("    $('#infoDiv').html(msg);\n")
   of.write('});\n')
   of.write("var overlay = A.graphicOverlay({color: 'cyan', lineWidth: 3});\n")
   of.write('aladin.addOverlay(overlay);\n')

   of.write(' var drawFunction = function(source, canvasCtx, viewParams) {\n')
   of.write('     canvasCtx.beginPath();\n')
   of.write("     canvasCtx.arc(source.x, source.y, source.data['size'] * 2, 0, 2 * Math.PI, false);\n")
   of.write('     canvasCtx.closePath();\n')
   of.write("     canvasCtx.strokeStyle = '#c38';\n")
   of.write('     canvasCtx.lineWidth = 3;\n')
   of.write('     canvasCtx.globalAlpha = 0.7,\n')
   of.write('     canvasCtx.stroke();\n')
   of.write("     var fov = Math.max(viewParams['fov'][0], viewParams['fov'][1]);\n")
   of.write('     // object name is displayed only if fov<10°\n')
   of.write('     if (fov>10) {\n')
   of.write('         return;\n')
   of.write('     }\n')
   of.write('     canvasCtx.globalAlpha = 0.9;\n')
   of.write('     canvasCtx.globalAlpha = 1;\n')
   of.write('     var xShift = 20;\n')
   of.write("     canvasCtx.font = '15px Arial'\n")
   of.write("     canvasCtx.fillStyle = '#eee';\n")
   of.write("     canvasCtx.fillText(source.data['name'], source.x + xShift, source.y -4);\n")
   of.write('     // object type is displayed only if fov<2°\n')
   of.write('     if (fov>10) {\n')
   of.write('         return;\n')
   of.write('     }\n')
   of.write("     canvasCtx.font = '12px Arial'\n")
   of.write("     canvasCtx.fillStyle = '#abc';\n")
   of.write("     canvasCtx.fillText(source.data['otype'], source.x + 2 + xShift, source.y + 10);\n")
   of.write(' };\n')
 
   with open(input_file_candidates, 'r') as candidates:
       csv_reader_candidates = csv.reader(candidates)
       header = next(csv_reader_candidates) 
       i = 0
       candidates =''
       of.write('// create catalog layer with custom draw function\n')
       of.write("var cat = A.catalog({}name: 'Candidates', shape: drawFunction{});\n".format(a0,a2))
       for row in csv_reader_candidates:
          i +=1
          number = int(row[0])       
          ra_cand = float(row[1])
          dec_cand = float(row[2])
          type = row[3]
          if len(row) < 6:
             name = row[4]
          elif row[5] == '---':
             name = row[4]
          else:
             name = row[5]
          nnam = 'A'+str(i)
          namenodots = name.replace(".", "")
          of.write("var {} = A.source({:.5f},{:.5f}, {} name: '{}- {}',size: 4.5,otype: '{}' {});\n".format(nnam,ra_cand, dec_cand, a0, number, name, type, a2))
          if i == 1:
             candidates = nnam
          else:
             candidates = candidates+','+nnam
       of.write('cat.addSources([{}]);\n'.format(candidates))
       of.write('aladin.addCatalog(cat);\n')
       of.write('\n\n')
   if os.path.exists(input_file_icecube):
      zero = 0
      of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: 'red' {} ));\n".format(zero, zero, zero, a0, a2))
      max_distance = 8
      with open(input_file_icecube, 'rt') as input_icecube:
          next(input_icecube)
          for line in input_icecube:
              # Split the line into columns using any whitespace characters
              columns = line.split()
              name = columns[0]
              ra_icecube = float(columns[1])
              dec_icecube = float(columns[2])
              maj= float(columns[3])/60
              min= float(columns[4])/60
              angle = -float(columns[5])-90
              if abs(dec_icecube - dec_center) < max_distance:
                  dist = sky_distance(ra_icecube, dec_icecube, ra_center, dec_center)
                  if dist < max_distance:
                     maj = maj*1.3
                     min = min*1.3
                     of.write("drawRotatedEllipse({:.5f},{:.5f},{:.2f},{:.5f},{:.5f},100,'#ffffe6');\n".format(ra_icecube,dec_icecube,angle,maj,min))


   if os.path.exists(input_file_2agile):
      with open(input_file_2agile, 'r') as inf:
        csv_reader_2agile = csv.reader(inf)
        header = next(csv_reader_2agile) 
        for row in csv_reader_2agile:
          ra2agile = float(row[1])
          dec2agile = float(row[2])
          maj = float(row[3])       
          min = float(row[4])      
          angle = -float(row[5])-90
          of.write("drawRotatedEllipse({:.5f},{:.5f},{:.2f},{:.5f},{:.5f},100,'#99ffcc');\n".format(ra2agile,dec2agile,angle,maj,min))

   if os.path.exists(input_file_efeds):
      with open(input_file_efeds, 'r') as inf:
        csv_reader = csv.reader(inf)
        header = next(csv_reader) 
        for row in csv_reader:
          ra = float(row[1])       
          dec = float(row[2])       
          radius = float(row[3]) / 3600
          #print("ra dec radius ",ra,dec,radius)
          #of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: '#1a1aff' {} ));\n".format(ra, dec, radius, a0, a2))
          of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: '#ff9966' {} ));\n".format(ra, dec, radius, a0, a2))

   if os.path.exists(input_file_erass1):
      with open(input_file_erass1, 'r') as inf:
        csv_reader = csv.reader(inf)
        header = next(csv_reader) 
        for row in csv_reader:
          ra = float(row[1])       
          dec = float(row[2])       
          if row[3] != '':
             radius = float(row[3]) / 3600 * 2
          else:
             radius = 3./3600.
          #of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: '#1a1aff' {} ));\n".format(ra, dec, radius, a0, a2))
          of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: '#ff9933' {} ));\n".format(ra, dec, radius, a0, a2))

   if os.path.exists(input_file_rass):
      with open(input_file_rass, 'r') as inf:
        csv_reader = csv.reader(inf)
        header = next(csv_reader) 
        for row in csv_reader:
          ra = float(row[1])       
          dec = float(row[2])       
          radius = 30 / 3600
          of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: '#1a1aff' {} ));\n".format(ra, dec, radius, a0, a2))

   if os.path.exists(input_file_2sxps):
      with open(input_file_2sxps, 'r') as inf:
        csv_reader = csv.reader(inf)
        header = next(csv_reader) 
        for row in csv_reader:
          ra = float(row[1])       
          dec = float(row[2])       
          radius = float(row[3]) / 3600
          of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: '#1a1aff' {} ));\n".format(ra, dec, radius, a0, a2))

   if os.path.exists(input_file_ousx):
      with open(input_file_ousx, 'r') as inf:
        csv_reader = csv.reader(inf)
        header = next(csv_reader) 
        for row in csv_reader:
          ra = float(row[1])       
          dec = float(row[2])       
          radius = 5 / 3600
          of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: '#1a1aff' {} ));\n".format(ra, dec, radius, a0, a2))

   if os.path.exists(input_file_4xmmdr13):
      with open(input_file_4xmmdr13, 'r') as inf:
        csv_reader = csv.reader(inf)
        header = next(csv_reader) 
        for row in csv_reader:
          ra = float(row[1])       
          dec = float(row[2])       
          radius = float(row[15]) 
          if radius < 2:
             radius = 2
          radius /=3600
          of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: '#1a1aff' {} ));\n".format(ra, dec, radius, a0, a2))

   if os.path.exists(input_file_chandra):
      with open(input_file_chandra, 'r') as inf:
        csv_reader = csv.reader(inf)
        header = next(csv_reader) 
        for row in csv_reader:
          ra = float(row[1])       
          dec = float(row[2])       
          radius = float(row[5]) 
          if radius < 2:
             radius = 2
          radius /=3600
          of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: '#1a1aff' {} ));\n".format(ra, dec, radius, a0, a2))

   if os.path.exists(input_file_xmmslew):
      with open(input_file_xmmslew, 'r') as inf:
        csv_reader = csv.reader(inf)
        header = next(csv_reader) 
        for row in csv_reader:
          ra = float(row[1])       
          dec = float(row[2])       
          radius = float(row[9])*1.5
          if radius < 1:
             radius = 1
          radius /=3600
          of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: '#1a1aff' {} ));\n".format(ra, dec, radius, a0, a2))

   if os.path.exists(input_file_ipc):
      with open(input_file_ipc, 'r') as inf:
        csv_reader = csv.reader(inf)
        header = next(csv_reader) 
        for row in csv_reader:
          ra = float(row[1])       
          dec = float(row[2])       
          radius = float(row[5]) / 3600
          of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: '#1a1aff' {} ));\n".format(ra, dec, radius, a0, a2))

   if os.path.exists(input_file_wga):
      with open(input_file_wga, 'r') as inf:
        csv_reader = csv.reader(inf)
        header = next(csv_reader) 
        for row in csv_reader:
          ra = float(row[1])       
          dec = float(row[2])       
          radius = float(row[5]) / 3600
          of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: '#1a1aff' {} ));\n".format(ra, dec, radius, a0, a2))

   if os.path.exists(input_file_first):
      with open(input_file_first, 'r') as inf:
        csv_reader = csv.reader(inf)
        header = next(csv_reader) 
        for row in csv_reader:
          ra = ra_to_deg(row[1])
          dec = dec_to_deg(row[2])
          radius = 1 / 3600 
          of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: 'red' {} ));\n".format(ra, dec, radius, a0, a2))

   if os.path.exists(input_file_vlass):
      with open(input_file_vlass, 'r') as inf:
        csv_reader = csv.reader(inf)
        header = next(csv_reader) 
        for row in csv_reader:
          ra = float(row[1])
          dec = float(row[2])
          radius = 1 / 3600 
          of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: 'red' {} ));\n".format(ra, dec, radius, a0, a2))

   if os.path.exists(input_file_nvss):
      with open(input_file_nvss, 'r') as inf:
        csv_reader = csv.reader(inf)
        header = next(csv_reader) 
        for row in csv_reader:
          ra = float(row[1])
          dec = float(row[2])
          radius = (float(row[3])+float(row[4])) / 2
          if radius < 1:
             radius = 1
          radius /=3600
          of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: 'red' {} ));\n".format(ra, dec, radius, a0, a2))

   if os.path.exists(input_file_4fgl):
      with open(input_file_4fgl, 'r') as inf:
        csv_reader_4fgl = csv.reader(inf)
        header = next(csv_reader_4fgl) 
        for row in csv_reader_4fgl:
          ra4fgl = float(row[1])       
          dec4fgl = float(row[2])       
          maj = 0.
          if row[3] !='':
             maj = float(row[3])       
          min = 0.
          if row[4] !='':
             min = float(row[4])       
          angle = 0
          if row[5] !='':
             angle = -float(row[5])-90
          of.write("drawRotatedEllipse({:.5f},{:.5f},{:.2f},{:.5f},{:.5f},100,'#ff99ff');\n".format(ra4fgl,dec4fgl,angle,maj,min))
          if row[11] != '':
             ra_counterpart = float(row[11])
             dec_counterpart = float(row[12])
             radius = 10/3600
             of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: '#99e6ff' {} ));\n".format(ra_counterpart, dec_counterpart, radius, a0, a2))

   if os.path.exists(input_file_racs):
      with open(input_file_racs, 'r') as inf:
        csv_reader = csv.reader(inf)
        header = next(csv_reader) 
        for row in csv_reader:
          ra = float(row[1])
          dec = float(row[2])
          radius = (float(row[8])+float(row[9])) / 2
          if radius < 1:
             radius = 1
          radius /=3600
          of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: 'red' {} ));\n".format(ra, dec, radius, a0, a2))
   if plot_error_ellipse == 'yes':
       of.write("drawRotatedEllipse({:.5f},{:.5f},{:.2f},{:.5f},{:.5f},100,'#cc3300');\n".format(ra_center,dec_center,ell_angle,ell_maj,ell_min))

   if os.path.exists(input_file_4lac):
      with open(input_file_4lac, 'r') as inf:
        csv_reader_4lac = csv.reader(inf)
        header = next(csv_reader_4lac) 
        for row in csv_reader_4lac:
          name = row[1]
          ra4lac = float(row[1])       
          dec4lac = float(row[2])       
          radius = 10/3600
          of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: '#ffff99' {} ));\n".format(ra4lac, dec4lac, radius, a0, a2))

   of.write('$("input[name=survey]").change(function() {  \n')
   of.write('   aladin.setImageSurvey($(this).val());\n')
   of.write('});\n')
   of.write('$("select[name=projection]").change(function() {\n')
   of.write('   aladin.setProjection($(this).val());\n')
   of.write('});\n')
   of.write('</script>\n')
   of.write('<script type="text/javascript">\n')
   of.write('document.getElementById("hipsBase").innerHTML=hips\n')
   of.write('</script>\n')
   of.write('</HTML>\n')
