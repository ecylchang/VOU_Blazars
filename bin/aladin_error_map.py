import os
import csv
import argparse

path ='tmp/'
input_file_4fgl = path+'4fgldr3.1.2.csv'
parser = argparse.ArgumentParser(description='FirmamentoAladin')
parser.add_argument('--ra', type=float, help='R.A. (degrees d/f=0.0)', default='0.0')
parser.add_argument('--dec', type=float, help='Dec. (degrees d/f=0.0)', default='0.0')
parser.add_argument('--input_file_error_circles', type=str, help='input_file error_circles d/f=error_map.txt)', default=path+'error_map.txt')
parser.add_argument('--output_file', type=str, help='output file d/f=vou-aladin-py.html)', default='vou-aladin-py.html')


ok = True
a1 = "'"
a2 = "}"
a0 = "{"

args = parser.parse_args()
ra = args.ra
dec = args.dec 
input_file_error_circles =args.input_file_error_circles
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
#   of.write('   <script type="text/javascript" src="https://code.jquery.com/jquery-1.12.1.min.js" charset="utf-8"></script>\n')
   of.write('   <script type="text/javascript" src="https://code.jquery.com/jquery-1.9.1.min.js" charset="utf-8"></script>\n')
   of.write('   <link rel="stylesheet" href="https://aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.css" >\n')
   of.write('   <script type="text/javascript">var jqMenu = jQuery.noConflict();</script>\n')
   of.write('   <script type="text/javascript">\n')
   of.write('   var hipsDir=null;</script>\n')
   of.write('</HEAD>\n')
   of.write('&emsp;&emsp;&emsp;&emsp;&emsp;\n')
   of.write('<script type="text/javascript" src="https://aladin.u-strasbg.fr/AladinLite/api/v2/latest/aladin.min.js" charset="utf-8"></script>\n')
   of.write("<legend class='title'>VOU-Blazars - SED error circles map for<sub><sub><sub>\n")
   of.write('<img class="logoImg" src="https://vouwebsite.web.app/blackLogo.png" alt="LOGO" width="100" height="30"></sub></sub></sub>&emsp;Based on "Aladin Lite" developed at CDS, Strasbourg Observatory, France<br></legend>\n')
   of.write('<table><tr valign="top"><td>\n')
   of.write('<div id="aladin-lite-div" style="width:80vw; height:80vh;"></div>\n')
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
   of.write('<input id="DSS2" type="radio"   name="survey" value="P/DSS2/color" checked="checked"><label for="DSS2">DSS2</label><br>\n')

   if dec > -30:
      of.write('<input id="PANSTARRS" type="radio"   name="survey" value="P/PanSTARRS/DR1/color-z-zg-g"><label for="P/PanSTARRS/DR1/color">PanSTARRS-DR1</label><br>\n')
   else:
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
   of.write('</td></tr><tr><td> <center> \n')
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

#   of.write('hipsDir = "https://alaskybis.cds.unistra.fr/SSC/xcatdb_P_XMM_PN_color/";\n')
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

   if os.path.exists(input_file_4fgl):
      with open(input_file_4fgl, 'r') as inf:
        csv_reader_4fgl = csv.reader(inf)
        header = next(csv_reader_4fgl) 
        for row in csv_reader_4fgl:
          ra4fgl = float(row[1])       
          dec4fgl = float(row[2])       
          maj = float(row[3])       
          min = float(row[4])       
          angle = -float(row[5])-90.
          #of.write("drawRotatedEllipse({:.5f},{:.5f},{:.2f},{:.5f},{:.5f},100,'BlueViolet');\n".format(ra4fgl,dec4fgl,angle,maj,min))
          of.write("drawRotatedEllipse({:.5f},{:.5f},{:.2f},{:.5f},{:.5f},100,'#ff99ff');\n".format(ra4fgl,dec4fgl,angle,maj,min))

   if not os.path.isfile(input_file_error_circles):
       print("file " + input_file_error_circles + " not found")
       print("Error circles will not be generated")
   else:    
      with open(input_file_error_circles, 'r') as infil:
          i = 0
          while ok:
              line = infil.readline()
              if not line:
                  break
              line = line.strip()
              values = line.split()
              ra = float(values[0])
              dec = float(values[1])
              ii = int(values[2])
              radius = float(values[3]) / 3600.0
              ia = abs(ii) // 1000
              i += 1
              if ia == 1:
                  of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: 'red' {} ));\n".format(ra, dec, radius, a0, a2))
              elif ia == 5:
                  of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: 'orange' {} ));\n".format(ra, dec, radius, a0, a2))
              elif ia == 6:
                  of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: 'silver' {} ));\n".format(ra, dec, radius, a0, a2))
              elif ia == 7:
                  of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: 'aquamarine' {} ));\n".format(ra, dec, radius, a0, a2))
              elif ia == 8:
                  of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: '#1a1aff' {} ));\n".format(ra, dec, radius, a0, a2))
              #elif ia == 9:
              #    of.write("overlay.add(A.circle({:.5f} , {:.5f} , {:.5f} , {} color: 'BlueViolet' {} ));\n".format(ra, dec, radius, a0, a2))
   of.write("overlay.add(A.circle(0. , 0. , 0. , {} color: 'white' {} ));\n".format(a0, a2))

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
