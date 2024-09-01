import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import os 
from astropy.constants import kpc
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
mpl.rcParams['axes.linewidth']  = 1
plt.style.use('seaborn-v0_8-deep')
#plt.style.use('seaborn-talk')
import argparse

def is_file_empty(file_path):
    return os.path.getsize(file_path) == 0


parser = argparse.ArgumentParser(description='SED plotting tool')
parser.add_argument('--infile1' , type=str, help='Input (SED) csv file name generated with VOU-Blazars (d/f=--)', default='--')
parser.add_argument('--infile2' , type=str, help='Input (SED) file name formatted for SSDC SED tool (e.g. sed_pow_4SED.txt)', default='--')
parser.add_argument('--infile3' , type=str, help='Input (SED) file from Browse (d/f=None)', default='--')
parser.add_argument('--infile4' , type=str, help='Input (SED) file from UVOT or ASASSN (d/f= None )', default='--')
parser.add_argument('--infile5' , type=str, help='Input average SED file for comparison (d/f= None )', default='--')
parser.add_argument('--infile6' , type=str, help='Input NuSTAR SED (d/f= None )', default='--')
parser.add_argument('--filetimeslice' , type=str, help='Input file with time filtered data (d/f= None )', default='--')
parser.add_argument('--fromSSDC' , type=str, help='Input exported SSDC SED tool file (d/f= None )', default='--')
parser.add_argument('--fromVIZIER' , type=str, help='Input exported Vizier SED tool (d/f= None )', default='--')
parser.add_argument('--fromZTF' , type=str, help='Input from ZTF  (d/f= None )', default='--')
parser.add_argument('--fromNED' , type=str, help='Input from NED photometric tool (d/f= None )', default='--')
parser.add_argument('--fromWPEAK' , type=str, help='W-Peak value from neowise_data program (d/f= None )', default='--')
parser.add_argument('--fromBLAST' , type=str, help='Input from BLAST (d/f= None )', default='--')
parser.add_argument('--outfile', type=str, help='Sed output file name (d/f=PySED.png)', default='PySED.png')
parser.add_argument('--dataoutfile', type=str, help='Data output file name (d/f=PySED_data.txt)', default='PySED_data.txt')
parser.add_argument('--xaxis'  , type=str, help='X-axis type (f=Frequency,e=Energy eV,k=Energy KeV,g=Energy GeV, t=Energy TeV, d/f=f)', default='f')
parser.add_argument('--title'  , type=str, help='SED title (centred) ', default=' ')
parser.add_argument('--ltitle'  , type=str, help='SED title (left justified)', default=' ')
parser.add_argument('--rtitle'  , type=str, help='SED title (rigth justified)', default=' ')
parser.add_argument('--redshift'  , type=float, help='Redshift', default='0.0')
parser.add_argument('--scaling_factor'  , type=float, help='Scaling factor for data in infil5', default='1.0')
parser.add_argument('--upperlimits'  , type=str, help='Upper limits', default='yes')

parser.add_argument('--CTAdetectability', type=str, help='Show approx. CTA detectability level(d/f=no)', default='no')
parser.add_argument('--GETemplate' , type=str, help='Add a giant elliptical template? (d/f=no)', default='no')
parser.add_argument('--BBTemplate' , type=str, help='Add a Blue Bump + accretion template? (d/f=no)', default='no')
parser.add_argument('--BBScalFactor' , type=float, help='Blue Bump Scaling factor (d/f=5.e-12)', default='5.e-12')
parser.add_argument('--C279Template' , type=str, help='Add a 3C279 template? (LBL, <nupeak> ~1.e13 d/f=no)', default='no')
parser.add_argument('--C279ScalFactor' , type=float, help='3C279 Scaling factor (d/f=1)', default='1')
parser.add_argument('--OJ287Template' , type=str, help='Add a OJ287 template? (IBL <nupeak> ~7.e13 d/f=no)', default='no')
parser.add_argument('--OJ287ScalFactor' , type=float, help='OJ287 Scaling factor (d/f=1)', default='1')
parser.add_argument('--BLLACTemplate' , type=str, help='Add a BL Lac template? (IBL <nupeak> ~1.e14 d/f=no)', default='no')
parser.add_argument('--BLLACScalFactor' , type=float, help='BL Lac Scaling factor (d/f=1)', default='1')
parser.add_argument('--C371Template' , type=str, help='Add a 3C371 template? (IBL <nupeak> ~2.e14 d/f=no)', default='no')
parser.add_argument('--C371ScalFactor' , type=float, help='3C371 Scaling factor (d/f=1)', default='1')
parser.add_argument('--TXS0506Template' , type=str, help='Add a TXS0506+056 template? (IBL, <nupeak> ~5.e14 d/f=no)', default='no')
parser.add_argument('--TXS0506ScalFactor' , type=float, help='TXS0506+056 Scaling factor (d/f=1)', default='1')
parser.add_argument('--PKS2155Template' , type=str, help='Add a PKS2155-304 template? (HBL, <nupeak> ~5.E15 d/f=no)', default='no')
parser.add_argument('--PKS2155ScalFactor' , type=float, help='PKS2155-304 Scaling factor (d/f=1)', default='1')
parser.add_argument('--MKN421Template' , type=str, help='Add a MKN421 template? (HBL, <nupeak> ~1.e17 d/f=no)', default='no')
parser.add_argument('--MKN421ScalFactor' , type=float, help='MKN421 Scaling factor (d/f=1)', default='1')
parser.add_argument('--ES1959Template' , type=str, help='Add a 1ES1959+650 template? (HBL, <nupeak> ~5.e17 d/f=no)', default='no')
parser.add_argument('--ES1959ScalFactor' , type=float, help='1ES1959+650 Scaling factor (d/f=1)', default='1')
parser.add_argument('--MKN501Template' , type=str, help='Add a MRK501 template? (HBL, <nupeak> ~5.e17 d/f=no)', default='no')
parser.add_argument('--MKN501ScalFactor' , type=float, help='MRK501 Scaling factor (d/f=1)', default='1')
parser.add_argument('--ShowLines' , type=str, help='PLot lines at 1e12, 1e13, 1e14 and 1e15 Hz (d/f=no)', default='no')

parser.add_argument('--eblcorrected'  , type=str, help='Show data ebl-corrected', default='no')
parser.add_argument('--erosita'  , type=str, help='Highlight eROSITA data', default='no')
parser.add_argument('--alma'  , type=str, help='Highlight ALMA data', default='no')
parser.add_argument('--neowise'  , type=str, help='Highlight NEOWISE data', default='no')
parser.add_argument('--smarts'  , type=str, help='Highlight SMARTS data', default='no')
parser.add_argument('--swift'  , type=str, help='Highlight Swift-XRT data', default='no')
parser.add_argument('--fermi'  , type=str, help='Highlight Fermi data', default='no')
parser.add_argument('--xmm'  , type=str, help='Highlight XMM data', default='no')
parser.add_argument('--nustar'  , type=str, help='Highlight NuSTAR data', default='no')
parser.add_argument('--panstarrs'  , type=str, help='Highlight PanSTARRS data', default='no')
parser.add_argument('--asassn'  , type=str, help='Highlight ASAS-SN data', default='no')
parser.add_argument('--ztf'  , type=str, help='Highlight ZTF data', default='no')

args = parser.parse_args()
type = args.xaxis
infile1 = args.infile1
infile2 = args.infile2
infile3 = args.infile3
infile4 = args.infile4
infile5 = args.infile5
infile6 = args.infile6
filetimeslice = args.filetimeslice
fromSSDC = args.fromSSDC
fromVIZIER = args.fromVIZIER
fromNED = args.fromNED
fromWPEAK = args.fromWPEAK
fromBLAST = args.fromBLAST
infile4 = args.fromZTF
outfile = args.outfile
dataoutfile = args.dataoutfile
sed_title = args.title
sed_ltitle = args.ltitle
sed_rtitle = args.rtitle
z = args.redshift
scaling_factor = args.scaling_factor
ulimit = args.upperlimits
EBL = args.eblcorrected
eROSITA = args.erosita
ALMA = args.alma
NEOWISE = args.neowise
XMM = args.xmm
SWIFT = args.swift
FERMI = args.fermi
NUSTAR = args.nustar
PANSTARRS = args.panstarrs
ASASSN = args.asassn
ZTF = args.ztf
SMARTS = args.smarts

ShowCTAdetectability = args.CTAdetectability
Templ_GE = args.GETemplate
Templ_BB = args.BBTemplate
ScalingFactorBB = args.BBScalFactor
Templ_3C279 = args.C279Template
ScalingFactor3C279 = args.C279ScalFactor
Templ_3C371 = args.C371Template
ScalingFactor3C371 = args.C371ScalFactor
Templ_OJ287 = args.OJ287Template
ScalingFactorOJ287 = args.OJ287ScalFactor
Templ_BLLAC = args.BLLACTemplate
ScalingFactorBLLAC = args.BLLACScalFactor
Templ_TXS0506 = args.TXS0506Template
ScalingFactorTXS0506 = args.TXS0506ScalFactor
Templ_PKS2155 = args.PKS2155Template
ScalingFactorPKS2155 = args.PKS2155ScalFactor
Templ_MKN421 = args.MKN421Template
ScalingFactorMKN421 = args.MKN421ScalFactor
Templ_MKN501 = args.MKN501Template
ScalingFactorMKN501 = args.MKN501ScalFactor
Templ_1ES1959 = args.ES1959Template
ScalingFactor1ES1959 = args.ES1959ScalFactor
ShowLines = args.ShowLines

f = open(dataoutfile,"w")

if infile1 != '--':
   data  = pd.read_csv(infile1,  delimiter="," , names=["Frequency", "Flux", "Flux_err", "MJD", "MJD_end", "IsDet", "Cat", "Reference"],  skiprows=1, header=None)
#   data = dat[dat["Flux"] != 'Infinity' ]
if infile2 != '--':
   data2 = pd.read_csv(infile2,  delimiter="|" , names=["ra","dec","Frequency", "Freq_err", "Flux", "Flux_err", "MJD", "MJD_end","IsDet","a","b"],  skiprows=0, header=None)
if infile3 != '--':
   data3 = pd.read_csv(infile3,  delimiter="," , names=["ra","dec","MJD","F5kev","F5kev_err","F05kev","F05kev_err","F15kev","F15kev_err","F3kev","F3kev_err","F45kev","F45kev_err","F1kev","F1kev_err","IsDet"],  skiprows=0, header=None, index_col=False)
#   data3["Freq_1kev"] = data3["ra"] * 0.0 + 2.418*10**(17)
   data3["Freq_1kev"] =  2.418*10**(17)
   data3["Freq_05kev"] = data3["Freq_1kev"] * 0.5 
   data3["Freq_15kev"] = data3["Freq_1kev"] * 1.5 
   data3["Freq_3kev"] = data3["Freq_1kev"] * 3.0 
   data3["Freq_45kev"] = data3["Freq_1kev"] * 4.5 
if infile4 != '--':
   data4  = pd.read_csv(infile4,  delimiter="," , names=["Frequency", "Flux", "Flux_err", "MJD", "MJD_end" , "IsDet" , "b"],  skiprows=1, header=None)
if infile5 != '--':
   data5  = pd.read_csv(infile5,  delim_whitespace=True , names=["LogFreq", "zero", "LogFlux", "zero1"],  skiprows=1, header=None)
   data5["xtype"] = 10**data5["LogFreq"]
   data5["Flux"] = 10**data5["LogFlux"]/scaling_factor
if infile6 != '--':
   data6 = pd.read_csv(infile6,  delimiter="|" , names=["ra","dec","Frequency", "Freq_error", "Flux", "Flux_err", "MJD", "d"],  skiprows=0, header=None)
if fromSSDC != '--':
   data7 = pd.read_csv(fromSSDC, delim_whitespace=True , names=["Frequency", "Ferr", "Flux", "Flux_err", "MJD", "MJD_end"],  skiprows=0, header=None)
if filetimeslice != '--':
   dtime  = pd.read_csv(filetimeslice,  delimiter="," , names=["Frequency", "Flux", "Flux_err", "MJD", "MJD_end", "IsDet", "Cat", "Reference"],  skiprows=1, header=None)
if fromVIZIER != '--' :
    there = os.path.isfile(fromVIZIER)
    empty = is_file_empty(fromVIZIER)
    if not there or empty: 
      fromVIZIER == '--'
    else:
      dataviz = pd.read_csv(fromVIZIER, delimiter="," , names=["Frequency", "Flux", "Flux_err"],  skiprows=1, header=None)
if fromNED != '--':
   dataNED = pd.read_csv(fromNED, delimiter="," , names=["Frequency", "Flux", "Flux_err", "MJD", "MJD_end", "flag", "Cat", "Reference"],  skiprows=1, header=None)
fig,ax1 = plt.subplots(figsize = (12, 8), facecolor = 'white', dpi=500)
fig,ax1 = plt.subplots(figsize = (12, 8), facecolor = 'white', dpi=500)
ax1.set_ylabel(r'E$\cdot$F$_{\rmE}$ [erg cm$^{-2}$ s$^{-1}$]', fontsize=24, fontweight='bold')

if type == 'f':
   if infile1 != '--':
      data["xtype"] = data["Frequency"] 
   if filetimeslice != '--':
      dtime["xtype"] = dtime["Frequency"] 
   if infile2 != '--':
      data2["xtype"] = data2["Frequency"] 
   if infile4 != '--':
      data4["xtype"] = data4["Frequency"] 
   if infile6 != '--':
      data6["xtype"] = data6["Frequency"] 
      dd6 = data6[(data6["Flux"]/data6["Flux_err"] > 2.0) & ( data6["Frequency"] < 1.7e19) ]
   if fromVIZIER != '--':
      dataviz["xtype"] = dataviz["Frequency"] 
   if fromSSDC != '--':
      data7["xtype"] = data7["Frequency"] 
#   ax1.set_xlabel(r'Frequency ($\nu$, observer frame) [Hz]', fontsize=24, fontweight='bold')
   if fromNED != '--':
      dataNED["xtype"] = dataNED["Frequency"] 
      dNED = dataNED[(dataNED["Flux"]/dataNED["Flux_err"] > 2.0)]
   ax1.set_xlabel(r'${\bf Frequency}$ $\nu_{\rm\bf observer~frame}{\bf [Hz]}$', fontsize=24 )
   ax1.set_ylabel(r'$\nu$F$_{\nu}$ [erg cm$^{-2}$ s$^{-1}$]', fontsize=24, fontweight='bold')
elif type == 'e':
   data["xtype"] = data["Frequency"] / (2.418*10**(14))
   data2["xtype"] = data2["Frequency"] / (2.418*10**(14))
   ax1.set_xlabel("Energy [eV]" , fontsize=24, fontweight='bold')
elif type == 'k':
   data["xtype"] = data["Frequency"] / (2.418*10**(17))
   data2["xtype"] = data2["Frequency"] / (2.418*10**(17))
   ax1.set_xlabel("Energy [KeV]" , fontsize=24, fontweight='bold')
elif type == 'g':
   data["xtype"] = data["Frequency"] / (2.418*10**(23))
   data2["xtype"] = data2["Frequency"] / (2.418*10**(23))
   ax1.set_xlabel("Energy [GeV]" , fontsize=24, fontweight='bold')
elif type == 't':
   data["xtype"] = data["Frequency"] / (2.418*10**(26))
   data2["xtype"] = data2["Frequency"] / (2.418*10**(26))
   ax1.set_xlabel("Energy [TeV]" , fontsize=24, fontweight='bold')

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.tick_params(axis='both', which='both', direction='in', length=5 , labelsize=18, top=True, right=True, labeltop=False, labelright=False)

if infile1 != '--':
   dda = data[data["Cat"].str.strip() != 'DEBL']
   debl = data[data["Cat"].str.strip() == 'DEBL']
   deros = data[( (data["Cat"].str.strip()  == 'eROSITA-EDR' ) | (data["Cat"].str.strip()  == 'eRASS1' ) ) & (data["IsDet"] != 'UL')]
   dalma = data[data["Cat"].str.strip()   == 'ALMA']
   dneowise = data[data["Cat"].str.strip()   == 'NEOWISE']
   dsmarts = data[data["Cat"].str.strip() == 'SMARTS']
   dnustar = data[data["Cat"].str.strip() == 'NuBlazar']
   dpanstarrs = data[data["Cat"].str.strip() == 'Pan-STARRS-LC']
   dasas = data[data["Cat"].str.strip() == 'ASAS-SN-LC']
   dztf = data[data["Cat"].str.strip() == 'ZTF-LC']
   dxmm = data[(data["Cat"].str.strip()  == 'XMMSL2') | (data["Cat"].str.strip()  == '4XMM-DR11') | (data["Cat"].str.strip()  == '4XMM-DR13')]
#   dswift = data[(data["Cat"].str.strip()  == '1OUSX') | (data["Cat"].str.strip()  == '2SXPS')  | (data["Cat"].str.strip()  == 'SUFST') | (data["Cat"].str.strip()  == 'XRTSPEC')]
   dswift = data[(data["Cat"].str.strip()  == 'SUFST') & (data["IsDet"] != 'UL')]
   dfermi = data[(data["Cat"].str.strip()  == 'Fermi')]
   dd  = dda[(dda["IsDet"] != 'UL') & (dda["Flux_err"] >= 0.) & ( ( (dda["Frequency"] > 1.4e13) | (dda["Frequency"] < 1.3e13) ) | ( ( (dda["Frequency"] > 1.3e13) & (dda["Frequency"] < 1.4e13) ) & (dda["Flux"] > 1.5e-12) ) ) ]
   dlimit = data[ (data["IsDet"] == 'UL') | (data["Flux_err"] >= data["Flux"]) ]
f.write ("Frequency (Hz), nuFnu flux (erg/cm2/s), nuFnu flux error (erg/cm2/s) , MJD (days) \n")
if fromNED != '--':
#   ax1.errorbar(dNED["xtype"],dNED["Flux"],xerr=None,yerr=dNED["Flux_err"], fmt='o', markeredgecolor='#339933', color = 'none', markersize=8)
   ax1.errorbar(dNED["xtype"],dNED["Flux"],xerr=None,yerr=dNED["Flux_err"], fmt='o',color = '#339933', markersize='4')
   f.write ("#Data from NED \n")
   df_list = [dataNED[['xtype', 'Flux', 'Flux_err']]]
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['xtype'] = df_list['xtype'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if fromVIZIER != '--':
   ax1.errorbar(dataviz["xtype"],dataviz["Flux"],xerr=None,yerr=dataviz["Flux_err"], fmt='o',markeredgecolor='#ff9900', color = '#ff9900', markersize=3,zorder=1)
if fromSSDC != '--':
   ax1.errorbar(data7["xtype"],data7["Flux"],xerr=None,yerr=data7["Flux_err"], fmt='o',markeredgecolor='#58aaee', color = '#58aaee', markersize=6,zorder=1)
   f.write ("#Data from SSDC SED tool \n")
   df_list = [data7[['xtype', 'Flux', 'Flux_err', 'MJD']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['xtype'] = df_list['xtype'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if infile2 != '--':
   dd2 = data2[(data2["Flux"]/data2["Flux_err"] > 1.8) & ( (data2["Frequency"] < 1.9e18) | (data2["Frequency"] > 1.e21) )]
if infile3 != '--':
   dd3  = data3[(data3["F1kev"]/data3["F1kev_err"] > 1.5) & (data3["IsDet"] != 'UL')]
   d305 = dd3[dd3["F05kev_err"] > 0.]
   d315 = dd3[dd3["F15kev_err"] > 0.]
   d330 = dd3[dd3["F3kev_err"] > 0.]
   d345 = dd3[dd3["F45kev_err"] > 0.]
   dlimit3 = data3[data3["IsDet"] == 'UL']
if infile4 != '--':
   d4  = data4[(data4["IsDet"] != 'UL') & (data4["Flux"]/data4["Flux_err"] > 2.0)]
   dlimit4  = data4[(data4["IsDet"] == 'UL')]
if infile1 != '--':
   ax1.errorbar(dd["xtype"],dd["Flux"],xerr=None,yerr=dd["Flux_err"], fmt='o',color = '#0066cc', markeredgecolor='black', markersize='5',zorder=8)
   f.write ("#Data from VOU-BLazars tool \n")
   df_list = [dd[['xtype', 'Flux', 'Flux_err' , 'MJD']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['xtype'] = df_list['xtype'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   if filetimeslice != '--':
      ax1.errorbar(dtime["xtype"],dtime["Flux"],xerr=None,yerr=dtime["Flux_err"], fmt='o',color = '#ff0055', markersize='4',zorder=18)
   if EBL == 'yes':
      ax1.errorbar(debl["xtype"],debl["Flux"],xerr=None, yerr=debl["Flux_err"], fmt='D',color = '#ff4d4d',fillstyle='none', markersize='5') 
   if eROSITA == 'yes':
      ax1.errorbar(deros["xtype"],deros["Flux"],xerr=None, yerr=deros["Flux_err"], fmt='o',color = '#ff2222', markersize='5',zorder=10) 
   if ALMA == 'yes':
      ax1.errorbar(dalma["xtype"],dalma["Flux"],xerr=None, yerr=dalma["Flux_err"], fmt='o',color = '#009933', markersize='5',zorder=10) 
   if NEOWISE == 'yes':
      ax1.errorbar(dneowise["xtype"],dneowise["Flux"],xerr=None, yerr=dneowise["Flux_err"], fmt='o',color = '#ff0066', markersize='5',zorder=10) 
   if SMARTS == 'yes':
      ax1.errorbar(dsmarts["xtype"],dsmarts["Flux"],xerr=None, yerr=dsmarts["Flux_err"], fmt='o',color = '#00cc99', markersize='5',zorder=10) 
   if SWIFT == 'yes':
      ax1.errorbar(dswift["xtype"],dswift["Flux"],xerr=None, yerr=dswift["Flux_err"], fmt='o',color = '#ff00ff', markersize='5',zorder=11) 
   if FERMI == 'yes':
      ax1.errorbar(dfermi["xtype"],dfermi["Flux"],xerr=None, yerr=dfermi["Flux_err"], fmt='o',color = '#ff00ff', markersize='5',zorder=11) 
   if NUSTAR == 'yes':
      ax1.errorbar(dnustar["xtype"],dnustar["Flux"],xerr=None, yerr=dnustar["Flux_err"], fmt='o',color = '#ff9933', markersize='5',zorder=10) 
   if ZTF == 'yes':
      ax1.errorbar(dztf["xtype"],dztf["Flux"],xerr=None, yerr=dztf["Flux_err"], fmt='o',color = '#ff9966', markersize='2',zorder=10) 
   if PANSTARRS == 'yes':
      ax1.errorbar(dpanstarrs["xtype"],dpanstarrs["Flux"],xerr=None, yerr=dpanstarrs["Flux_err"], fmt='o',color = '#339966', markersize='2',zorder=10) 
   if ASASSN == 'yes':
      ax1.errorbar(dasas["xtype"],dasas["Flux"],xerr=None, yerr=dasas["Flux_err"], fmt='o',color = '#66ff66', markersize='2',zorder=10) 
   if XMM == 'yes':
      ax1.errorbar(dxmm["xtype"],dxmm["Flux"],xerr=None, yerr=dxmm["Flux_err"], fmt='o',color = '#ff6699', markersize='5',zorder=10) 
if infile2 != '--':
   ax1.errorbar(dd2["xtype"],dd2["Flux"],xerr=None,yerr=dd2["Flux_err"], fmt='o',color = '#1f7604', markersize='3')
   f.write ("#Data from Swift_xrtproc XRT analysis \n")
   df_list = [dd2[['xtype', 'Flux', 'Flux_err', 'MJD']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['xtype'] = df_list['xtype'].astype(float)
   df_list.to_csv(f,mode = 'a' , index=None, header=None , sep=',',  float_format='%.8E')
if infile3 != '--':
   ax1.errorbar(dd3["Freq_1kev"] ,dd3["F1kev"] ,xerr=None,yerr=dd3["F1kev_err"] , fmt='o', color = '#fbd799', markersize='4')
   ax1.errorbar(d305["Freq_05kev"],d305["F05kev"],xerr=None,yerr=d305["F05kev_err"], fmt='o',color = '#e33e1d', markersize='3')
   ax1.errorbar(d315["Freq_15kev"],d315["F15kev"],xerr=None,yerr=d315["F15kev_err"], fmt='o',color = '#e33e1d', markersize='3')
   ax1.errorbar(d345["Freq_45kev"],d345["F45kev"],xerr=None,yerr=d345["F45kev_err"], fmt='o',color = '#e33e1d', markersize='3')
   f.write ("#1KeV data from Swift_xrtproc useful for lightcurve use \n")
   df_list = [dd3[['Freq_1kev', 'F1kev', 'F1kev_err', 'MJD']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['Freq_1kev'] = df_list['Freq_1kev'].astype(float)
   df_list.to_csv(f, index=None,  sep=',',  header=None ,float_format='%.8E')
   f.write ("#Data from Swift_xrtproc XIMAGE analysis \n")
   df_list = [d305[['Freq_05kev', 'F05kev', 'F05kev_err', 'MJD']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['Freq_05kev'] = df_list['Freq_05kev'].astype(float)
   df_list.to_csv(f, index=None,  sep=',',  header=None ,float_format='%.8E')
   df_list = [d315[['Freq_15kev', 'F15kev', 'F15kev_err', 'MJD']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['Freq_15kev'] = df_list['Freq_15kev'].astype(float)
   df_list.to_csv(f, index=None,  sep=',',  header=None ,float_format='%.8E')
   df_list = [d345[['Freq_45kev', 'F45kev', 'F45kev_err', 'MJD']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['Freq_45kev'] = df_list['Freq_45kev'].astype(float)
   df_list.to_csv(f, index=None,  sep=',',  header=None , float_format='%.8E')
xmin=1.e7
xmax=5.e27
ax1.set_xlim(xmin,xmax)
if infile4 != '--':
   ax1.errorbar(d4["xtype"],d4["Flux"],xerr=None,yerr=d4["Flux_err"], fmt='d',color = '#ff9966', markersize='3')
#   ax1.errorbar(data4["xtype"],data4["Flux"],xerr=None,yerr=data4["Flux_err"], fmt='o',color = '#ffcccc', markersize='5')
   f.write ("#Data from UVOT or ASASSN  \n")
   df_list = [d4[['xtype', 'Flux', 'Flux_err' , 'MJD']]]
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['xtype'] = df_list['xtype'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if infile5 != '--':
   ax1.errorbar(data5["xtype"],data5["Flux"],xerr=None,yerr=None, fmt='o',color = '#ffb84d', markersize='4')
if infile6 != '--':
   ax1.errorbar(dd6["xtype"],dd6["Flux"],xerr=None, yerr=dd6["Flux_err"], fmt='o',color = '#ff00ee', markeredgecolor='black', markersize='7')
   f.write ("#Data from NuSTAR analysis (Middei et al. 2021) \n")
   df_list = [dd6[['xtype', 'Flux', 'Flux_err' , 'MJD']]]
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['xtype'] = df_list['xtype'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if ulimit == 'yes':
   if infile1 != '--':
      ax1.errorbar(dlimit["xtype"],dlimit["Flux"], yerr=0.15*dlimit["Flux"], fmt='o', markersize='0', color='#8b908e', uplims=True) 
      f.write ("#Upper limits from VOU_blazars \n")
      f.write ("#Frequency, nufnu flux (erg/cm2/s), MJD (days) \n")
      df_list = [dlimit[['xtype', 'Flux', 'MJD']]]
      df_list = pd.concat(df_list, ignore_index=True)
      df_list['xtype'] = df_list['xtype'].astype(float)
      df_list.to_csv(f, index=None,  sep=',',  header=None ,float_format='%.8E')
   if infile3 != '--':
      ax1.errorbar(dlimit3["Freq_1kev"],dlimit3["F1kev"],xerr=None,yerr=0.15*dlimit3["F1kev"], fmt='o', markersize='0', color='#888888', uplims=True)
      f.write ("#Upper limits from XRT XIMAGE analysis \n")
      df_list = [dlimit3[['Freq_1kev', 'F1kev', 'MJD']]]
      df_list = pd.concat(df_list, ignore_index=True)
      df_list['Freq_1kev'] = df_list['Freq_1kev'].astype(float)
      df_list.to_csv(f, index=None,  sep=',',  header=None ,float_format='%.8E')
   if infile4 != '--':
      ax1.errorbar(dlimit4["xtype"],dlimit4["Flux"],xerr=None,yerr=0.15*dlimit4["Flux"], fmt='o', markersize='0', color='#888888', uplims=True)
      f.write ("#Upper limits from ASAS-SN data \n")
      df_list = [dlimit4[['xtype', 'Flux', 'MJD']]]
      df_list = pd.concat(df_list, ignore_index=True)
      df_list['xtype'] = df_list['xtype'].astype(float)
      df_list.to_csv(f, index=None,  sep=',',  header=None ,float_format='%.8E')
mnl, mxl = ax1.get_ylim()
f.write ("#Data from template(s) \n")
f.write ("#Frequency, nufnu flux (erg/cm2/s) \n")
if Templ_BB =='yes': 
# alpha_ox = -0.137 Log(nuLnu @1.2e15Hz)+4.704 
#Steffen ApJ 2006, 131,2826
   infile ='/Users/paologiommi/app/Templates/BlueBumpTemplNormalised_4py.txt'
   dataTempl_BB = pd.read_csv(infile, delimiter="\s+", names=["Frequency", "Flux"],  skiprows=1, header=None)
   dataTempl_BB["FreqObserverframe"] = dataTempl_BB["Frequency"] / (1.+z) 
   dataTempl_BB["Flux5000"] = dataTempl_BB["Flux"] * ScalingFactorBB
   if z > 0. :
     distance = cosmo.luminosity_distance(z)*kpc*1e5
     nulnuat1p2e15 = 4.*np.pi*distance.value**2*ScalingFactorBB*1.938
# 1.e15 = 2500 A
# 1.938 is the flux value at nu = 1.2e15 in the BlueBumpTemplNormalised_4py.txt file 
# 0.154 is the flux value at nu = 2.41E17 in the BlueBumpTemplNormalised_4py.txt file
     scalingXflux = 10.**((-0.137*np.log10(nulnuat1p2e15)+4.704+1.0)*2.605)/0.154
     ddUVOTTB = dataTempl_BB [dataTempl_BB["Frequency"] < 1.e16]
     ddXTB = dataTempl_BB [dataTempl_BB["Frequency"] > 1.e16]
     ax1.errorbar(ddUVOTTB["FreqObserverframe"],ddUVOTTB["Flux5000"],xerr=None, yerr=None, fmt='o', color = '#009900', markersize='1') 
     ax1.errorbar(ddXTB["FreqObserverframe"],ddXTB["Flux5000"]*scalingXflux,xerr=None, yerr=None, fmt='o', color = '#990099', markersize='1') 
   ax1.errorbar(dataTempl_BB["FreqObserverframe"],dataTempl_BB["Flux5000"],xerr=None, yerr=None, fmt='o', color = '#009900', markersize='1') 
if Templ_MKN501 =='yes': 
   infile ='/Users/paologiommi/app/Templates/TemplateMRK501SyncLow.csv'
   dataTempl_MKN501  = pd.read_csv(infile, delimiter="," , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl_MKN501["FreqObserverframe"] = 10**dataTempl_MKN501["LogFreq"] * 1.03 / (1.+z) 
   dataTempl_MKN501["FluxScaled"] = 10**dataTempl_MKN501["LogFlux"] * ScalingFactorMKN501
   if ScalingFactorMKN501 < 0. :
     distance = cosmo.luminosity_distance(z)*kpc*1e5
     distance_mkn501 = cosmo.luminosity_distance(0.033)*kpc*1e5
     scal = -(distance.value/distance_mkn501.value)**-2*ScalingFactorMKN501
     dataTempl_MKN501["FluxScaled"] = 10**dataTempl_MKN501["LogFlux"]*scal
   ax1.errorbar(dataTempl_MKN501["FreqObserverframe"],dataTempl_MKN501["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#ff9900', markersize='0', ls='--') 
   f.write ("#MKN501 template -1st part\n")
   df_list = [dataTempl_MKN501[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateMRK501SyncHigh.csv'
   dataTempl_MKN501  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl_MKN501["FreqObserverframe"] = 10**dataTempl_MKN501["LogFreq"] * 1.03 / (1.+z) 
   dataTempl_MKN501["FluxScaled"] = 10**dataTempl_MKN501["LogFlux"] * ScalingFactorMKN501
   if ScalingFactorMKN501 < 0. :
      dataTempl_MKN501["FluxScaled"] = 10**dataTempl_MKN501["LogFlux"] * scal 
   ax1.errorbar(dataTempl_MKN501["FreqObserverframe"],dataTempl_MKN501["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#ff9900', markersize='0', ls='--') 
   f.write ("#MKN501 template -2nd part \n")
   df_list = [dataTempl_MKN501[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateMRK501ICLow.csv'
   dataTempl_MKN501  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl_MKN501["FreqObserverframe"] = 10**dataTempl_MKN501["LogFreq"] * 1.03 / (1.+z) 
   dataTempl_MKN501["FluxScaled"] = 10**dataTempl_MKN501["LogFlux"] * ScalingFactorMKN501
   if ScalingFactorMKN501 < 0. :
      dataTempl_MKN501["FluxScaled"] = 10**dataTempl_MKN501["LogFlux"] * scal
   ax1.errorbar(dataTempl_MKN501["FreqObserverframe"],dataTempl_MKN501["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#ff9900', markersize='0', ls='--') 
   f.write ("#MKN501 template -3rd part \n")
   df_list = [dataTempl_MKN501[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateMRK501ICHigh.csv'
   dataTempl_MKN501  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl_MKN501["FreqObserverframe"] = 10**dataTempl_MKN501["LogFreq"] * 1.03 / (1.+z) 
   dataTempl_MKN501["FluxScaled"] = 10**dataTempl_MKN501["LogFlux"] * ScalingFactorMKN501
   if ScalingFactorMKN501 < 0. :
      dataTempl_MKN501["FluxScaled"] = 10**dataTempl_MKN501["LogFlux"] * scal
   ax1.errorbar(dataTempl_MKN501["FreqObserverframe"],dataTempl_MKN501["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#ff9900', markersize='0', ls='--') 
   f.write ("#MKN501 template -4th part \n")
   df_list = [dataTempl_MKN501[['FreqObserverframe', 'FluxScaled']]] 
   df_list = [dataTempl_MKN501[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if Templ_PKS2155 =='yes': 
   infile ='/Users/paologiommi/app/Templates/TemplatePKS2155SyncLow.csv'
   dataTempl  = pd.read_csv(infile, delimiter="," , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.116 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorPKS2155
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#ff0000', markersize='0', ls='--') 
   f.write ("#PKS2155 template -1st part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplatePKS2155SyncHigh.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.116 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorPKS2155
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#ff0000', markersize='0', ls='--') 
   f.write ("#PKS2155 template -2nd part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplatePKS2155ICLow.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.116 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorPKS2155
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#ff0000', markersize='0', ls='--') 
   f.write ("#PKS2155 template -3rd part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplatePKS2155ICHigh.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.116 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorPKS2155
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#ff0000', markersize='0', ls='--') 
   f.write ("#PKS2155 template -4th part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if Templ_OJ287 =='yes': 
   infile ='/Users/paologiommi/app/Templates/TemplateOJ287SyncLow.csv'
   dataTempl  = pd.read_csv(infile, delimiter="," , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.3056 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorOJ287
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#996633', markersize='0', ls='--') 
   f.write ("#OJ287 template -1st part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateOJ287SyncHigh.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.3056 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorOJ287
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#996633', markersize='0', ls='--') 
   f.write ("#OJ287 template -2nd part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateOJ287ICLow.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.3056 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorOJ287
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#996633', markersize='0', ls='--') 
   f.write ("#OJ287 template -3rd part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateOJ287ICHigh.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.3056 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorOJ287
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#996633', markersize='0', ls='--') 
   f.write ("#OJ287 template -4th part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if Templ_BLLAC =='yes': 
   infile ='/Users/paologiommi/app/Templates/TemplateBLLacSyncLow.csv'
   dataTempl  = pd.read_csv(infile, delimiter="," , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.0686 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorBLLAC
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#996633', markersize='0', ls='--') 
   f.write ("#BL Lac template -1st part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateBLLacSyncHigh.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.0686 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorBLLAC
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#996633', markersize='0', ls='--') 
   f.write ("#BL Lac template -2nd part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateBLLacICLow.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.0686 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorBLLAC
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#996633', markersize='0', ls='--') 
   f.write ("#BL Lac template -3rd part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateBLLacICHigh.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.0686 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorBLLAC
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#996633', markersize='0', ls='--') 
   f.write ("#BL Lac template -4th part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if Templ_1ES1959 =='yes': 
   infile ='/Users/paologiommi/app/Templates/Template1ES1959SyncLow.csv'
   dataTempl  = pd.read_csv(infile, delimiter="," , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.047 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactor1ES1959
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#0099cc', markersize='0', ls='--') 
   f.write ("#1ES1959+650 template -1st part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/Template1ES1959SyncHigh.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.047 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactor1ES1959
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#0099cc', markersize='0', ls='--') 
   f.write ("#1ES1959+650 template -2nd part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/Template1ES1959ICLow.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.047 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactor1ES1959
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#0099cc', markersize='0', ls='--') 
   f.write ("#1ES1959+650 template -3rd part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/Template1ES1959ICHigh.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.047 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactor1ES1959
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#0099cc', markersize='0', ls='--') 
   f.write ("#1ES1959+650 template -4th part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if Templ_MKN421 =='yes': 
   infile ='/Users/paologiommi/app/Templates/TemplateMKN421SyncLow.csv'
   dataTempl  = pd.read_csv(infile, delimiter="," , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.116 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorMKN421
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#006600', markersize='0', ls='--') 
   f.write ("#MKN421 template -1st part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateMKN421SyncHigh.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.116 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorMKN421
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#006600', markersize='0', ls='--') 
   f.write ("#MKN421 template -2nd part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateMKN421ICLow.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.116 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorMKN421
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#006600', markersize='0', ls='--') 
   f.write ("#MKN421 template -3rd part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateMKN421ICHigh.csv'
   dataTempl  = pd.read_csv(infile, delimiter=","  , names=["LogFreq", "LogFlux" ],  skiprows=1, header=None)
   dataTempl["FreqObserverframe"] = 10**dataTempl["LogFreq"] * 1.116 / (1.+z) 
   dataTempl["FluxScaled"] = 10**dataTempl["LogFlux"] * ScalingFactorMKN421
   ax1.errorbar(dataTempl["FreqObserverframe"],dataTempl["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#006600', markersize='0', ls='--') 
   f.write ("#MKN421 template -4th part\n")
   df_list = [dataTempl[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if Templ_TXS0506 =='yes': 
   infile ='/Users/paologiommi/app/Templates/TemplateTXS0506SyncLow.csv'
   dataTempl_TXS0506  = pd.read_csv(infile,  delimiter=",", names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_TXS0506["FreqObserverframe"] = 10**dataTempl_TXS0506["LogFreq"] * 1.3365 / (1.+z) 
   dataTempl_TXS0506["FluxScaled"] = 10**dataTempl_TXS0506["LogFlux"] * ScalingFactorTXS0506
   ax1.errorbar(dataTempl_TXS0506["FreqObserverframe"],dataTempl_TXS0506["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#00e6e6', markersize='0', ls='--') 
   f.write ("#TXS0506-056 template -1st part\n")
   df_list = [dataTempl_TXS0506[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateTXS0506SyncHigh.csv'
   dataTempl_TXS0506  = pd.read_csv(infile,  delimiter=",", names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_TXS0506["FreqObserverframe"] = 10**dataTempl_TXS0506["LogFreq"] * 1.3365 / (1.+z) 
   dataTempl_TXS0506["FluxScaled"] = 10**dataTempl_TXS0506["LogFlux"] * ScalingFactorTXS0506
   ax1.errorbar(dataTempl_TXS0506["FreqObserverframe"],dataTempl_TXS0506["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#00e6e6', markersize='0', ls='--') 
   f.write ("#TXS0506-056 template -2nd part\n")
   df_list = [dataTempl_TXS0506[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateTXS0506ICLow.csv'
   dataTempl_TXS0506  = pd.read_csv(infile,  delimiter=",", names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_TXS0506["FreqObserverframe"] = 10**dataTempl_TXS0506["LogFreq"] * 1.3365 / (1.+z) 
   dataTempl_TXS0506["FluxScaled"] = 10**dataTempl_TXS0506["LogFlux"] * ScalingFactorTXS0506
   ax1.errorbar(dataTempl_TXS0506["FreqObserverframe"],dataTempl_TXS0506["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#00e6e6', markersize='0', ls='--') 
   f.write ("#TXS0506-056 template -3rd part\n")
   df_list = [dataTempl_TXS0506[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/TemplateTXS0506ICHigh.csv'
   dataTempl_TXS0506  = pd.read_csv(infile,  delimiter=",", names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_TXS0506["FreqObserverframe"] = 10**dataTempl_TXS0506["LogFreq"] * 1.3365 / (1.+z) 
   dataTempl_TXS0506["FluxScaled"] = 10**dataTempl_TXS0506["LogFlux"] * ScalingFactorTXS0506
   ax1.errorbar(dataTempl_TXS0506["FreqObserverframe"],dataTempl_TXS0506["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#00e6e6', markersize='0', ls='--') 
   f.write ("#TXS0506-056 template -4th part\n")
   df_list = [dataTempl_TXS0506[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if Templ_3C371 =='yes': 
   infile ='/Users/paologiommi/app/Templates/Template3C371SyncLow.csv'
   dataTempl_3C371  = pd.read_csv(infile,  delimiter=",", names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_3C371["FreqObserverframe"] = 10**dataTempl_3C371["LogFreq"] * 1.3365 / (1.+z) 
   dataTempl_3C371["FluxScaled"] = 10**dataTempl_3C371["LogFlux"] * ScalingFactor3C371
   ax1.errorbar(dataTempl_3C371["FreqObserverframe"],dataTempl_3C371["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#e6ac00', markersize='0', ls='--') 
   f.write ("#3C371 template -1st part\n")
   df_list = [dataTempl_3C371[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/Template3C371SyncHigh.csv'
   dataTempl_3C371  = pd.read_csv(infile,  delimiter=",", names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_3C371["FreqObserverframe"] = 10**dataTempl_3C371["LogFreq"] * 1.051 / (1.+z) 
   dataTempl_3C371["FluxScaled"] = 10**dataTempl_3C371["LogFlux"] * ScalingFactor3C371
   ax1.errorbar(dataTempl_3C371["FreqObserverframe"],dataTempl_3C371["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#e6ac00', markersize='0', ls='--') 
   f.write ("#3C371 template -2nd part\n")
   df_list = [dataTempl_3C371[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/Template3C371ICLow.csv'
   dataTempl_3C371  = pd.read_csv(infile,  delimiter=",", names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_3C371["FreqObserverframe"] = 10**dataTempl_3C371["LogFreq"] * 1.051 / (1.+z) 
   dataTempl_3C371["FluxScaled"] = 10**dataTempl_3C371["LogFlux"] * ScalingFactor3C371
   ax1.errorbar(dataTempl_3C371["FreqObserverframe"],dataTempl_3C371["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#e6ac00', markersize='0', ls='--') 
   f.write ("#3C371 template -3rd part\n")
   df_list = [dataTempl_3C371[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/Template3C371ICHigh.csv'
   dataTempl_3C371  = pd.read_csv(infile,  delimiter=",", names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_3C371["FreqObserverframe"] = 10**dataTempl_3C371["LogFreq"] * 1.051 / (1.+z) 
   dataTempl_3C371["FluxScaled"] = 10**dataTempl_3C371["LogFlux"] * ScalingFactor3C371
   ax1.errorbar(dataTempl_3C371["FreqObserverframe"],dataTempl_3C371["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#e6ac00', markersize='0', ls='--') 
   f.write ("#3C371 template -4th part\n")
   df_list = [dataTempl_3C371[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if Templ_3C279 =='yes': 
   infile ='/Users/paologiommi/app/Templates/Template3C279SyncLow.csv'
   dataTempl_3C279  = pd.read_csv(infile,  delimiter="," , names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_3C279["FreqObserverframe"] = 10**dataTempl_3C279["LogFreq"] * 1.536 / (1.+z) 
   dataTempl_3C279["FluxScaled"] = 10**dataTempl_3C279["LogFlux"] * ScalingFactor3C279
   ax1.errorbar(dataTempl_3C279["FreqObserverframe"],dataTempl_3C279["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#cc3300', markersize='0', ls='--') 
   f.write ("#3C279 template -1st part\n")
   df_list = [dataTempl_3C279[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/Template3C279SyncHigh.csv'
   dataTempl_3C279  = pd.read_csv(infile,  delimiter="," , names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_3C279["FreqObserverframe"] = 10**dataTempl_3C279["LogFreq"] * 1.536 / (1.+z) 
   dataTempl_3C279["FluxScaled"] = 10**dataTempl_3C279["LogFlux"] * ScalingFactor3C279
   ax1.errorbar(dataTempl_3C279["FreqObserverframe"],dataTempl_3C279["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#cc3300', markersize='0', ls='--') 
   f.write ("#3C279 template -2nd part \n")
   df_list = [dataTempl_3C279[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/Template3C279ICLow.csv'
   dataTempl_3C279  = pd.read_csv(infile,  delimiter="," , names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_3C279["FreqObserverframe"] = 10**dataTempl_3C279["LogFreq"] * 1.536 / (1.+z) 
   dataTempl_3C279["FluxScaled"] = 10**dataTempl_3C279["LogFlux"] * ScalingFactor3C279
   ax1.errorbar(dataTempl_3C279["FreqObserverframe"],dataTempl_3C279["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#cc3300', markersize='0', ls='--') 
   f.write ("#3C279 template -3rd part \n")
   df_list = [dataTempl_3C279[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
   infile ='/Users/paologiommi/app/Templates/Template3C279ICHigh.csv'
   dataTempl_3C279  = pd.read_csv(infile,  delimiter="," , names=["LogFreq", "LogFlux"],  skiprows=0, header=None)
   dataTempl_3C279["FreqObserverframe"] = 10**dataTempl_3C279["LogFreq"] * 1.536 / (1.+z) 
   dataTempl_3C279["FluxScaled"] = 10**dataTempl_3C279["LogFlux"] * ScalingFactor3C279
   ax1.errorbar(dataTempl_3C279["FreqObserverframe"],dataTempl_3C279["FluxScaled"],xerr=None, yerr=None, fmt='o', color = '#cc3300', markersize='0', ls='--') 
   f.write ("#3C279 template -4th part \n")
   df_list = [dataTempl_3C279[['FreqObserverframe', 'FluxScaled']]] 
   df_list = pd.concat(df_list, ignore_index=True)
   df_list['FreqObserverframe'] = df_list['FreqObserverframe'].astype(float)
   df_list.to_csv(f, index=None,  sep=',', header=None , float_format='%.8E')
if z != 0.0:
   distance =  cosmo.luminosity_distance(z)*kpc*1e5
   const = 4.*np.pi*distance.value**2
   if Templ_GE =='yes':
      infile ='/Users/paologiommi/app/Templates/GiantEllipticalTemplate_4py.txt'
      dataTempl_GE = pd.read_csv(infile, delimiter="\s+", names=["Frequency", "Luminosity"],  skiprows=1, header=None)
      dataTempl_GE["FreqObserverframe"] = dataTempl_GE["Frequency"] / (1.+z)
      dataTempl_GE["Flux"] = dataTempl_GE["Luminosity"] / const * 1.0
      ax1.errorbar(dataTempl_GE["FreqObserverframe"],dataTempl_GE["Flux"],xerr=None, yerr=None, fmt='o', color = '#ff33cc', markersize='1', ls ='-')
   ax2 = ax1.twinx()
   mn, mx = ax1.get_ylim()
   ax2.set_ylim(mn*const, mx*const)
   ax2.set_yscale('log')
   ax1.tick_params(axis='both', which='both', direction='in', labelsize=18, top=True, right=False, labeltop=False, labelright=False)
   ax2.tick_params(axis='both', which='both', direction='in', labelsize=18, top=True, right=True, labeltop=False )
   if type == 'f':
      ax2.set_ylabel(r'$\nu$L$_{\nu}$ [erg s$^{-1}$]', fontsize=18, fontweight='bold') 
   else:
      ax2.set_ylabel(r'E$\cdot$L$_{\rmE}$ [erg s$^{-1}$]', fontsize=18, fontweight='bold')
mnl, mxl = ax1.get_ylim()
if ShowLines =='yes': 
   ax1.axvline(1.e12,ymin=0,ymax=1,linestyle='--', color = '#aaaaaa',zorder=9)
   ax1.axvline(1.e13,ymin=0,ymax=1,linestyle='--', color = '#aaaaaa',zorder=9)
   ax1.axvline(1.e14,ymin=0,ymax=1,linestyle='--', color = '#aaaaaa',zorder=9)
   ax1.axvline(1.e15,ymin=0,ymax=1,linestyle='--', color = '#aaaaaa',zorder=9) 
   ax1.axvline(1.e16,ymin=0,ymax=1,linestyle='--', color = '#aaaaaa',zorder=9) 
   ax1.axvline(1.e17,ymin=0,ymax=1,linestyle='--', color = '#aaaaaa',zorder=9) 
   ax1.axvline(1.e18,ymin=0,ymax=1,linestyle='--', color = '#aaaaaa',zorder=9) 
   yt = mnl*(mxl/mnl)**0.01
   ax1.text(1.e12,yt,'10$^{12}$',fontsize=10 , color = '#888888')
   ax1.text(1.e14,yt,'10$^{14}$',fontsize=10 , color = '#888888')
   ax1.text(1.e15,yt,'10$^{15}$',fontsize=10 , color = '#888888')
   ax1.text(1.e17,yt,'10$^{17}$',fontsize=10 , color = '#888888')
   ax1.text(1.e18,yt,'10$^{18}$',fontsize=10 , color = '#888888')
nupeak_blast = 1.0
if fromBLAST != '--':
   with open(fromBLAST, 'r') as f:
      contents = f.read().strip()
      numbers = contents.split()
      nupeak_blast = 10**float(numbers[0])
      nupeak_blast_error = float(numbers[1])
yyst = 0.02
if fromWPEAK != '--':
   dataWPEAK = pd.read_csv(fromWPEAK, delimiter="," , names=["nupeak", "nufnuMean", "nufnu_Min", "nufnu_Max"],  skiprows=1, header=None)
   nup = float(dataWPEAK["nupeak"])
   nupeak_Wpeak=10**nup
   nufnumean=10**dataWPEAK["nufnuMean"]
   nufnumin=10**dataWPEAK["nufnu_Min"]
   nufnumax=10**dataWPEAK["nufnu_Max"]
   if (nup > 1):
      yye = yyst+ 0.03
      yy = mnl*(mxl/mnl)**yye
      if ((nupeak_blast/nupeak_Wpeak > 1.5) | (nupeak_blast/nupeak_Wpeak < 0.7) ):
        ax1.text(nupeak_Wpeak*0.6,yy,r'W-Peak $\nu_{\rm peak}$ + $\nu$f$\nu_{\rm peak}$ range',fontsize=9,rotation='vertical',color = '#4ce600',zorder=12)
      else:
        if (nupeak_blast > nupeak_Wpeak):
           ax1.text(nupeak_Wpeak/2.7,yy,r'W-Peak $\nu_{\rm peak}$ + $\nu$f$\nu_{\rm peak}$ range',fontsize=9,rotation='vertical',color = '#4ce600',zorder=12)
        else:
           ax1.text(nupeak_Wpeak*1.1,yy,r'W-Peak $\nu_{\rm peak}$ + $\nu$f$\nu_{\rm peak}$ range',fontsize=9,rotation='vertical',color = '#4ce600',zorder=12)
      ax1.vlines(nupeak_Wpeak,ymin=nufnumin,ymax=nufnumax,linestyle='-',color='#66ff33',linewidth=4,zorder=12) 
      ax1.plot(nupeak_Wpeak,nufnumean,marker="s",markersize=4, markeredgecolor="none",markerfacecolor="#ff6699",zorder=12)
      if ((nup > 14.5) & (nup < 15.5)):
         ax1.plot(nupeak_Wpeak,nufnumax,marker="^",markersize=7, markerfacecolor="#66ff33",zorder=12)
      if (nup >= 15.5):
         ax1.plot(nupeak_Wpeak,nufnumax,marker="^",markersize=7, markeredgecolor="none", markerfacecolor="#66ff33",zorder=12)
         ax1.plot(nupeak_Wpeak,nufnumin,marker="v",markersize=6, markeredgecolor="none", markerfacecolor="#66ff33",zorder=12)
mnl, mxl = ax1.get_ylim()
yyst = 0.02
if fromBLAST != '--':
   ax1.axvline(nupeak_blast,ymin=0.35,ymax=0.98,linestyle=':',color='#b36b00',linewidth=2,zorder=16) 
   yye = yyst + 0.03
   yy = mnl*(mxl/mnl)**yye
   ax1.text(nupeak_blast*0.6,yy,r'BLAST $\nu_{\rm peak}$',fontsize=9,rotation='vertical',color = '#b36b00')
   yye = 0.38
   yy = mnl*(mxl/mnl)**yye
   ax1.plot(nupeak_blast,yy,marker="^",markersize=8, markeredgecolor="none", markerfacecolor="#b36b00",zorder=16)
   ax1.plot(nupeak_blast,yy*2,marker="^",markersize=8, markeredgecolor="none", markerfacecolor="#b36b00",zorder=16)
   yye = 0.98
   yy = mnl*(mxl/mnl)**yye
   ax1.plot(nupeak_blast,yy,marker="v",markersize=8, markeredgecolor="none", markerfacecolor="#b36b00",zorder=16)
   ax1.plot(nupeak_blast,yy/2,marker="v",markersize=8, markeredgecolor="none", markerfacecolor="#b36b00",zorder=16)
#
if ShowCTAdetectability == 'yes':
   if nupeak_blast > 13.0:
      nufnuminCTA=10**-12.1
      nufnuminIACTs=10**-11.1
      nufnuminLHAASO=10**-10.7
   else:
      nufnuminCTA=10**-12.3
      nufnuminIACTs=10**-11.3
   if nupeak_Wpeak > 10.:
      ax1.vlines(nupeak_Wpeak,ymin=nufnuminCTA,ymax=nufnuminIACTs,linestyle=':',color='#999966',linewidth=7,zorder=9,alpha=0.7) 
      ax1.vlines(nupeak_Wpeak,ymin=nufnuminCTA,ymax=nufnuminIACTs,linestyle='-',color='#999966',linewidth=8,zorder=9,alpha=0.4) 
      mnl, mxl = ax1.get_ylim()
      ax1.vlines(nupeak_Wpeak,ymin=nufnuminIACTs,ymax=mxl,linestyle=':',color='#ff6666',linewidth=8,zorder=9,alpha=0.8) 
      ax1.vlines(nupeak_Wpeak,ymin=nufnuminIACTs,ymax=mxl,linestyle='-',color='#ff6666',linewidth=9,zorder=9,alpha=0.6) 
   else:
      ax1.vlines(nupeak_blast,ymin=nufnuminCTA,ymax=nufnuminIACTs,linestyle=':',color='#999966',linewidth=7,zorder=9,alpha=0.7) 
      ax1.vlines(nupeak_blast,ymin=nufnuminCTA,ymax=nufnuminIACTs,linestyle='-',color='#999966',linewidth=8,zorder=9,alpha=0.4) 
      mnl, mxl = ax1.get_ylim()
      ax1.vlines(nupeak_blast,ymin=nufnuminIACTs,ymax=mxl,linestyle=':',color='#ff6666',linewidth=8,zorder=9,alpha=0.8) 
      ax1.vlines(nupeak_blast,ymin=nufnuminIACTs,ymax=mxl,linestyle='-',color='#ff6666',linewidth=9,zorder=9,alpha=0.6) 
   infile ='/Users/paologiommi/app/Templates/CTA-Sensitivity.csv'
   dataCTASens  = pd.read_csv(infile, delimiter="," , names=["En-TeVUnits", "Flux" ],  skiprows=0, header=None)
   dataCTASens["Freq"] = 2.418e26*dataCTASens["En-TeVUnits"] 
   ax1.errorbar(dataCTASens["Freq"],dataCTASens["Flux"],xerr=None, yerr=None, fmt='o', color = '#999966', markersize='1', ls='--')
   str ='CTAO best sensitivity '
   ax1.text(8.0e24,8.e-14,str,fontsize=8 , rotation=-58, color = '#999966')
#
yye = 0.02
if (Templ_GE =='yes') & (z != 0.0) : 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   ax1.text(2.e22,yy,'Giant Elliptical template',fontsize=10 , color = '#ff33cc')
if Templ_BB =='yes': 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   ax1.text(2.e22,yy,'QSO template',fontsize=10 , color = '#009900')
if Templ_MKN501 =='yes': 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   str ='MRK501 template * '+ str(ScalingFactorMKN501)
   ax1.text(2.e22,yy,str,fontsize=10 , color = '#ff9900')
if Templ_PKS2155 =='yes': 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   str ='PKS2155-304 template * '+ str(ScalingFactorPKS2155)
   ax1.text(2.e22,yy,str,fontsize=10 , color = '#ff0000')
if Templ_TXS0506 =='yes': 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   str ='TXS0506+056 template * '+ str(ScalingFactorTXS0506)
   ax1.text(2.e22,yy,str,fontsize=10 , color = '#00e6e6')
if Templ_3C371 =='yes': 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   str ='3C371 template * '+ str(ScalingFactor3C371)
   ax1.text(2.e22,yy,str,fontsize=10 , color = '#e6ac00')
if Templ_3C279 =='yes': 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   str ='3C279 template * '+ str(ScalingFactor3C279)
   ax1.text(2.e22,yy,str,fontsize=10 , color = '#cc3300')
if Templ_MKN421 =='yes': 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   str ='MKN421 template * '+ str(ScalingFactorMKN421)
   ax1.text(2.e22,yy,str,fontsize=10 , color = '#006600')
if Templ_1ES1959 =='yes': 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   str ='1ES1959+650 template * '+ str(ScalingFactor1ES1959)
   ax1.text(2.e22,yy,str,fontsize=10 , color = '#0099cc')
if Templ_OJ287 =='yes': 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   str ='OJ287 template * '+ str(ScalingFactorOJ287)
   ax1.text(2.e22,yy,str,fontsize=10 , color = '#996633')
if Templ_BLLAC =='yes': 
   yye = yye + 0.03
   yy = mnl*(mxl/mnl)**yye
   str ='BL Lac template * '+ str(ScalingFactorBLLAC)
   ax1.text(2.e22,yy,str,fontsize=10 , color = '#996633')

ax2 = ax1.twiny()
ax2.set_xlabel(r'Energy$_{\rm\bf observer ~frame}$ [eV]', fontweight='bold', fontsize=16)
ax2.set_xscale('log')
const = 2.418e14
ax2.set_xlim(xmin/const,xmax/const)
ax2.tick_params(axis='both', which='both', direction='in', labelsize=18, top=True, right=False)

plt.title(sed_title, fontsize=18, fontweight='bold')
plt.title(sed_ltitle, loc='left' , color='#800000', weight='bold')
plt.title(sed_rtitle, loc='right' , color='#800000', weight='bold')
#plt.savefig(outfile, bbox_inches='tight', format='pdf')
plt.savefig(outfile, bbox_inches='tight', format='png')
