import argparse
from math import pi

import numpy as np
import pandas as pd
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse

import ligo.skymap.plot


def convert_to_skycoord(data_df, ra_col_name, dec_col_name):
    """
    Will convert all ra and dec coordinates to astropy.coordinates.SkyCoord format.

    Parameters
    ----------
    data_df : pd.DataFrame
        Data in DataFrame format.
    ra_col_name : str
        Right ascension column name.
    dec_col_name : str
        Declination column name.

    Returns
    -------
    astropy.coordinates.SkyCoord
    """
    return SkyCoord(f"{data_df[ra_col_name]} {data_df[dec_col_name]}", unit='deg')


def plot_ellipse(ax, center, major, minor, angle, facecolor='#d4deea', edgecolor='#1e4b9f', linestyle='dashdot'):
    """
    Parameters
    ----------
    ax : matplotlib.axes._subplots.AstroDegreesZoomAxesSubplot
        Declared axes, on which ellipses will be drawn.
    center : list or tuple
        RA and DEC coordinates in a list
    major : int or float
        Major axis size of ellipse in arcmin.
    minor : int or float
        Minor axis size of ellipse in arcmin.
    angle : int or float
        Angle size of ellipse in arcmin. This parameter will rotate ellipse.
    facecolor : str
        Color of the inner part of the ellipse.
    edgecolor : str
        Color of the edge of the ellipse.
    linestyle : str
        Line Style of the edge of the ellipse.

    Returns
    -------
    matplotlib.axes._subplots.AstroDegreesZoomAxesSubplotcdelta)
    """
    major_axes_deg = Angle('{}'.format(major), unit='arcmin').deg
    minor_axes_deg = Angle('{}'.format(minor), unit='arcmin').deg
    angle_axes_deg = Angle('{}'.format(angle), unit='deg').deg
    cdelta = np.cos(Angle(center[1], unit='deg').rad)
    aa = angle_axes_deg * pi / 180.
    #w = 2 * minor_axes_deg / cdelta
    #h = 2 * major_axes_deg 
    w = 2 * minor_axes_deg * np.sqrt( (np.cos(aa)/cdelta) ** 2 +np.sin(aa)**2 )
    h = 2 * major_axes_deg * np.sqrt( (np.sin(aa)/cdelta) ** 2 +np.cos(aa)**2 )
#    h = 2 * minor_axes_deg  
#    w = 2 * major_axes_deg  / cdelta
#    print (h,w)
    angle_of_ellipse =  -angle_axes_deg 
    ell = Ellipse(xy=center, width=w, height=h, angle=angle_of_ellipse,
                 transform=ax.get_transform('world'), facecolor=facecolor, linestyle=linestyle, edgecolor=edgecolor)
    ax.add_artist(ell)
    return ax


def plot_points(data_to_plot, ax, marker='s', marker_size=2, color='None', marker_edge_width=2,
                marker_face_color='color', numbering=True, zorder=1):
    """
    Will check data and plot errorbar in axes, also will add annotation by argument value.

    Parameters
    ----------
    data_to_plot : pd.DataFrame
        Data which will be plotted.
    ax : matplotlib.axes._subplots.AstroDegreesZoomAxesSubplot
        Declared axes, on which plot will be drawn.
    marker : str
        Type of marker, which will be used.
        Check available values here
        https://matplotlib.org/3.1.0/api/markers_api.html#module-matplotlib.markers
    marker_size : int or float
        Size of the marker.
    color : str
        Color of the marker in hex or string format.
    marker_edge_width : int or float
        Markers edge width parameter.
    marker_face_color : str
        Color for filling inner part of marker.
    numbering : Boolean
        Apply annotation with numbers. 
    zorder : int
        Order of the plot points in z-scale.
         
    Returns
    -------
    matplotlib.axes._subplots.AstroDegreesZoomAxesSubplot
    """
    if marker_face_color == 'color':
        marker_face_color = color
    annotation_x_diff = siz / 20
    annotation_y_diff = -siz / 20
    #annotation_x_diff = siz / 30
    #annotation_y_diff = -siz / 40
    if data_to_plot.size > 0:
        for idx_num in data_to_plot.sky_coordinates.index:
            ax.errorbar(data_to_plot.loc[idx_num, 'sky_coordinates'].ra.deg,
                        data_to_plot.loc[idx_num, 'sky_coordinates'].dec.deg, transform=ax.get_transform('world'),
                        marker=marker, markersize=marker_size, markeredgewidth=marker_edge_width, color=color,
                        markerfacecolor=marker_face_color, zorder=zorder)
            if 'number' in data_to_plot.columns and numbering:
                ax.text(data_to_plot.loc[idx_num, 'sky_coordinates'].ra.deg - annotation_x_diff,
                        data_to_plot.loc[idx_num, 'sky_coordinates'].dec.deg - annotation_y_diff,
                        f"{data_to_plot.loc[idx_num, 'number']}", transform=ax.get_transform('world'),
                        fontsize=annotation_fontsize, color='#555555', weight='bold')

    return ax


# Parse all necessary arguments.
parser = argparse.ArgumentParser(description='Sky map')
parser.add_argument('--ra', type=float, help='R.A. (degrees d/f=0.0)', default='0.0')
parser.add_argument('--dec', type=float, help='Dec. (degrees d/f=0.0)', default='0.0')
parser.add_argument('--FOV', type=float, help='size (degrees d/f=3.0)', default='3.0')
parser.add_argument('--radius', type=float, help='radius (degrees d/f=0)', default='0')
parser.add_argument('--major', type=float, help='Ellipse1 major axis (d/f=0)', default='0.0')
parser.add_argument('--minor', type=float, help='Ellipse1 minor axis  (d/f=0)', default='0.0')
parser.add_argument('--angle', type=float, help='Ellipse1 rotation angle (d/f=90)', default='90')
parser.add_argument('--major2', type=float, help='Ellipse2 major axis (d/f=0)', default='0.0')
parser.add_argument('--minor2', type=float, help='Ellipse2 minor axis  (d/f=0)', default='0.0')
parser.add_argument('--angle2', type=float, help='Ellipse2 rotation angle  (d/f=90)', default='90')
parser.add_argument('--infile_candidates', type=str, help='Candidates file name (d/f=candidates.csv)',
                    default='candidates.csv')
parser.add_argument('--infile_multi', type=str, help='Multi-frequency sources file name (d/f=find_out_temp.txt)',
                    default='find_out_temp.txt')
parser.add_argument('--infile_local', type=str, help='Local filename (d/f=IHBL-catalog, use "none" if not requested )',
                    default='IHBL-catalog')
parser.add_argument('--infile_local2', type=str, help='Local filename (d/f=none)',
                    default='none')
parser.add_argument('--infile_IceCube', type=str, help='Local filename (d/f=IceCubeTracksPositions.txt)',
                    default='IceCubeTracksPositions.txt')
parser.add_argument('--outfile', type=str, help='Output file name (d/f=skymap.png)', default='skymap.png')
parser.add_argument('--title', type=str, help='Output file name (d/f=--)', default='--')
parser.add_argument('--catsonly', type=str, help='Show only candidates from existing catalogues (d/f=--)', default='--')
parser.add_argument('--plot_stars', type=str, help='plot stars from HSTGSC catalog (d/f=no)', default='no')

args = parser.parse_args()
ra = args.ra
dec = args.dec
dra = "{} d".format(ra)
ddec = "{} d".format(dec)
radec = "{} d {} d".format(ra, dec)
if dec > 0:
    radec = "{} d + {} d".format(ra, dec)
radec_sky = SkyCoord(radec)
siz = args.FOV / 60
size = "{} deg".format(siz)
markersize = 9
markersize_gamma = 6
markersize_local = 10
markersize_catalogued = 3
markersize_blue = 2
markersize_hbl = 3
markersize_orange = 2
markersize_cyan = 2
marker_edge_width = 1
marker_edge_width_hbl_ibl_lbl = 0.2
rad = args.radius
annotation_fontsize = 5
title_fontsize = 10
major = args.major
minor = args.minor
angle = args.angle
if rad > 0:
    major = rad
    minor = rad
    angle = 90
major2 = args.major2
minor2 = args.minor2
angle2 = args.angle2
infile_cand = args.infile_candidates
infile_multi = args.infile_multi
infile_local = args.infile_local
infile_local2 = args.infile_local2
infile_IceCube = args.infile_IceCube
plot_stars = args.plot_stars
outfile = args.outfile
title = args.title
catsonly = args.catsonly
# Reading all datafiles.

if infile_local != 'none':
   if dec > 0:
     infile_local = '~/cats4topcat/'+infile_local+'-north.txt'
   else:
     infile_local = '~/cats4topcat/'+infile_local+'-south.txt'
#if infile_local2 != 'none':
#   if dec > 0:
#     infile_local2 = '~/cats4topcat/'+infile_4lac+'-north.txt'
#   else:
#     infile_local2 = '~/cats4topcat/'+infile_4lac+'-south.txt'
#
infile_IceCube = '~/cats4topcat/'+infile_IceCube

data = pd.read_csv(infile_cand, delimiter=",", names=["number", "rra", "ddec", "type", "source_name"], skiprows=1,
                   header=None, index_col=False)
data2 = pd.read_csv(infile_multi, delimiter="\s+", names=["rra2", "ddec2", "code", "source_name"], skiprows=0,
                    header=None, index_col=False)
if infile_local != 'none':
   data3 = pd.read_csv(infile_local, delimiter=",", names=["source_name", "rra3", "ddec3"], skiprows=1, header=None,index_col=False)
if infile_local2 != 'none':
   data4 = pd.read_csv(infile_local2, delimiter=",", names=["source_name", "rra4", "ddec4"], skiprows=1, header=None,index_col=False)
dataIC = pd.read_csv(infile_IceCube, delimiter=",", names=["source_name", "rra", "ddec"], skiprows=1, header=None,
                    index_col=False)
data2['radiocode'] = abs(data2['code'] - (data2['code'] / 100).astype(int) * 100)

# Selecting necessary data.
sizra = siz / np.cos(dec * pi / 180.)
#print(len(data))
if len(data) > 0:
    # Calculating all coordinates for all available data.
    if (ra-sizra) < 0:
       add = 360+ra-sizra
       data = data[(data["ddec"] < (dec + siz)) & (data["ddec"] > (dec - siz)) & 
            (data["rra"] > (ra - sizra)) & ( (data["rra"] < (ra + sizra)) | (data["rra"] > add ) )]
    elif (ra+sizra) > 360:
       add = ra+sizra-360
       data = data[(data["ddec"] < (dec + siz)) & (data["ddec"] > (dec - siz)) &
           ( (data["rra"] > (ra - sizra)) | (data["rra"] < add ) ) & (data["rra"] < (ra + sizra))]
    else:
       data = data[(data["ddec"] < (dec + siz)) & (data["ddec"] > (dec - siz)) & (data["rra"] > (ra - sizra)) &
             (data["rra"] < (ra + sizra))]
    data['sky_coordinates'] = data.apply(convert_to_skycoord, args=('rra', 'ddec'), axis=1)
    duncl = data[(data["type"] == "UNCL")]
    dhbl_orange = data[(data["type"] == "HBL")]
    dibl_cyan = data[(data["type"] == "IBL")]
    dlbl_blue = data[(data["type"] == "LBL")]
    dother_black = data[((data["type"] == "Unknown") | (data["type"] == "radio-AGN") | (data["type"] == "radio-UNCL")) & (data["ddec"] < 90.)]
    dother_gray = data[(data["type"] == "HSTGSC") | (data["type"] == "OptCand") | (data["type"] == "RQ-AGN") ]
    dother_cat = data[(data["type"] == "Catalog")]
    dother_star = data[(data["type"] == "X-Star")]
    #dpulsar_purple = data[(data["type"] == "Pulsar")]
    dxcandidate_x_values = data[(data["type"] == "type-8")]
    drcandidate_r_values = data[(data["type"] == "type-9")]
    dhspcat = data[data["source_name"].astype(str).str.contains('3HSP')]
    dbzbcat = data[data["source_name"].astype(str).str.contains('5BZ')]
    dcrates_red = data[(data["type"] == "CRATES") | (data["type"] == "5BZCat") | (data["type"] == "3HSP")]
else:
    # If the data is missing
    dhbl_orange, dibl_cyan, dlbl_blue, dother_black, dother_gray, dother_cat , dother_star = [pd.DataFrame()] * 7
    dxcandidate_x_values, drcandidate_r_values, dhspcat, dbzbcat, dcrates_red = [pd.DataFrame()] * 5
if len(data2) > 0:
    if (ra-sizra) < 0:
       add = 360+ra-sizra
       data2 = data2[(data2["ddec2"] < (dec + siz)) & (data2["ddec2"] > (dec - siz)) &
            (data2["rra2"] > (ra - sizra)) & ( (data2["rra2"] < (ra + sizra)) | (data2["rra2"] > add ) )]
    elif (ra+sizra) > 360:
       add = ra+sizra-360
       data2 = data2[(data2["ddec2"] < (dec + siz)) & (data2["ddec2"] > (dec - siz)) &
           ( (data2["rra2"] > (ra - sizra)) | (data2["rra2"] < add ) ) & (data2["rra2"] < (ra + sizra))]
    else:
       data2 = data2[(data2["ddec2"] < (dec + siz)) & (data2["ddec2"] > (dec - siz)) & (data2["rra2"] > (ra - sizra)) &
              (data2["rra2"] < (ra + sizra))]
    data2['sky_coordinates'] = data2.apply(convert_to_skycoord, args=('rra2', 'ddec2'), axis=1)
    dpulsar_purple = data2[(data2["code"] == -8888)]
    dothersmall_star = data2[(data2["code"] > 100000) & (data2["code"] <= 103000)]
    dothermedium_star = data2[(data2["code"] > 103000) & (data2["code"] <= 106000)]
    dotherlarge_star = data2[(data2["code"] > 106000) & (data2["code"] < 110000)]
    dothersverysmall_agn = data2[(data2["code"] > 110000) & (data2["code"] < 111001)]
    dothersmall_agn = data2[(data2["code"] >= 111001) & (data2["code"] < 113000)]
    dothermedium_agn = data2[(data2["code"] > 113001) & (data2["code"] < 115000)]
    dotherlarge_agn = data2[(data2["code"] > 115001) & (data2["code"] < 120000)]
    dothersmall_gray = data2[(data2["code"] > 70000) & (data2["code"] <= 73000)]
    dothermedium_gray = data2[(data2["code"] > 73000) & (data2["code"] <= 76000)]
    dotherlarge_gray = data2[(data2["code"] > 76000) & (data2["code"] < 80000)]
    dother_black_small = data2[(data2["code"] > 40000) & (data2["code"] < 43000)]
    dother_black_medium = data2[(data2["code"] > 43100) & (data2["code"] < 46000)]
    dother_black_large = data2[(data2["code"] > 46100) & (data2["code"] < 50000)]
    dhblsmall_orange = data2[(data2["code"] > 10000) & (data2["code"] < 13000)]
    dhblmedium_orange = data2[(data2["code"] > 13000) & (data2["code"] < 16000)]
    dhbllarge_orange = data2[(data2["code"] > 16000) & (data2["code"] < 20000)]
    diblsmall_cyan = data2[(data2["code"] > 20000) & (data2["code"] < 23000)]
    diblmedium_cyan = data2[(data2["code"] > 23000) & (data2["code"] < 26000)]
    dibllarge_cyan = data2[(data2["code"] > 26000) & (data2["code"] < 30000)]
    dlblsmall_blue = data2[(data2["code"] > 30000) & (data2["code"] < 40000) & (data2["radiocode"] <= 30)]
    dlblmedium_blue = data2[(data2["code"] > 30000) & (data2["code"] < 40000) & (data2["radiocode"] > 30) & (data2["radiocode"] <= 60)]
    dlbllarge_blue = data2[(data2["code"] > 30000) & (data2["code"] < 40000) & (data2["radiocode"] > 60)]
    dcrates_small = data2[(data2["code"] < -30000) & (data2["code"] > -30030)]
    dcrates_medium = data2[( (data2["code"] < -30030) & (data2["code"] > -30060) ) | (data2["code"] == -50000)]
    dcrates_large = data2[(data2["code"] < -30060) & (data2["code"] > -40000)]
    dradio_small = data2[(data2["code"] < -90000) & (data2["code"] > -90030)]
    dradio_medium = data2[( (data2["code"] < -90030) & (data2["code"] > -90060) ) | (data2["code"] == -50000)]
    dradio_large = data2[(data2["code"] < -90060) & (data2["code"] > -100000)]
    dxray_small = data2[(data2["code"] < -80000) & (data2["code"] > -83000)]
    dxray_medium = data2[( (data2["code"] < -83000) & (data2["code"] >= -86000) ) | (data2["code"] == -50000)]
    dxray_large = data2[(data2["code"] < -86000) & (data2["code"] > -90000)]
    dgrb = data2[data2["code"] == -2222]
    d4lac = data2[data2["code"] == -4444]
    dmilliquas = data2[data2["code"] == -7000]
    dclusters = data2[data2["code"] == -40000]
    dgamma = data2[( (data2["code"] == -1111) | (data2["code"] == -1234) ) & (data2["ddec2"] < 90.)]
    dstellar_cluster = data2[(data2["code"] == 160000)]
    dcv = data2[(data2["code"] == 120000)]
    dxray_binaries = data2[(data2["code"] == 130000)]
    dsnr = data2[(data2["code"] == 140000)]
    dmolecular_clouds = data2[(data2["code"] == 150000)]
else:
    # If the data is missing
    dhblsmall_orange, dhblmedium_orange, dhbllarge_orange = [pd.DataFrame()] * 3
    diblsmall_cyan, diblmedium_cyan, dibllarge_cyan = [pd.DataFrame()] * 3
    dlblsmall_blue, dlblmedium_blue, dlbllarge_blue = [pd.DataFrame()] * 3
    dcrates_small, dcrates_medium, dcrates_large = [pd.DataFrame()] * 3
    dradio_small, dradio_medium, dradio_large = [pd.DataFrame()] * 3
    dxray_small, dxray_medium, dxray_large = [pd.DataFrame()] * 3
    dother_black_small, dother_black_medium, dother_black_large = [pd.DataFrame()] * 3
    dmilliquas, dgrb, d4lac, dclusters, dgamma, dstellar_cluster, dcv, dxray_binaries, dsnr, dmolecular_clouds = [pd.DataFrame()] * 10


if infile_local != 'none':
   if len(data3) > 0:
      data3['sky_coordinates'] = data3.apply(convert_to_skycoord, args=('rra3', 'ddec3'), axis=1)
      if (ra-sizra) < 0:
         add = 360+ra-sizra
         d3good = data3[(data3["ddec3"] < (dec + siz)) & (data3["ddec3"] > (dec - siz)) &
            (data3["rra3"] > (ra - sizra)) & ( (data3["rra3"] < (ra + sizra)) | (data3["rra3"] > add) )]
      elif (ra+sizra) > 360:
         add = ra+sizra-360
         d3good = data3[(data3["ddec3"] < (dec + siz)) & (data3["ddec3"] > (dec - siz)) &
            ((data3["rra3"] > (ra - sizra)) | (data3["rra3"] < add)) & (data3["rra3"] < (ra + sizra))]
      else:
         d3good = data3[(data3["ddec3"] < (dec + siz)) & (data3["ddec3"] > (dec - siz)) & (data3["rra3"] > (ra - sizra)) &
                     (data3["rra3"] < (ra + sizra))]
    # d3good local_values.
      if len(d3good) > 0:
        toprint = d3good.apply(
            lambda x: f"{x.rra3} +{x.ddec3}  {x.source_name}" if (
                    x.ddec3 > 0) else f"{x.rra3} {x.ddec3}  {x.source_name}",
            axis=1)
        print(toprint)
   else:
       # If the data is missing
       d3good = pd.DataFrame()

if infile_local2 != 'none':
   if len(data4) > 0:
      data4['sky_coordinates'] = data4.apply(convert_to_skycoord, args=('rra4', 'ddec4'), axis=1)
      d4good = data4[(data4["ddec4"] < (dec + siz)) & (data4["ddec4"] > (dec - siz)) & (data4["rra4"] > (ra - sizra)) &
                     (data4["rra4"] < (ra + sizra))]
    # d4good local_values.
      if len(d4good) > 0:
        toprint = d4good.apply(
            lambda x: f"{x.rra4} +{x.ddec4}  {x.source_name}" if (
                    x.ddec4 > 0) else f"{x.rra4} {x.ddec4}  {x.source_name}",
            axis=1)
        print(toprint)
   else:
       # If the data is missing
       d4good = pd.DataFrame()

if len(dataIC) > 0:
    dataIC['sky_coordinates'] = dataIC.apply(convert_to_skycoord, args=('rra', 'ddec'), axis=1)
    dICgood = dataIC[(dataIC["ddec"] < (dec + siz)) & (dataIC["ddec"] > (dec - siz)) & (dataIC["rra"] > (ra - sizra)) &
                   (dataIC["rra"] < (ra + sizra))]
    # dICgood local_values.
    if len(dICgood) > 0:
        toprint = dICgood.apply(
            lambda x: f"{x.rra} +{x.ddec}  {x.source_name}" if (
                    x.ddec > 0) else f"{x.rra} {x.ddec}  {x.source_name}",
            axis=1)
        print(toprint)
else:
    # If the data is missing
    dICgood = pd.DataFrame()

# Correcting title.
if title == '--':
    title = f"Image Centre R.A.={radec_sky.ra.to_string(unit='hour', decimal=False, sep=' ')} " \
            f"DEC={radec_sky.dec.to_string(decimal=False, sep=' ')}"
else:
    title = args.title

# Plotter part.
fig = plt.figure(figsize=(10, 4), facecolor='white', dpi=400)
ax_map = plt.axes(projection='astro degrees zoom',
                  # projection='astro hours zoom',
                  center=radec,
                  radius=size)
ax_map.grid(color='lightgray', linestyle='dashed')
ax_map.set_xlabel('R.A. (degrees)')
ax_map.set_ylabel('Dec. (degrees)')
ax_map.set_title(title, fontsize=title_fontsize)
ax_map.tick_params(axis='both', which='both', direction='in')

# Centre of the ellipses.
ellipse_center = (SkyCoord(radec).ra.deg, SkyCoord(radec).dec.deg)

# Drawing ellipses.
#plot_ellipse(ax_map, ellipse_center, major, minor, angle, facecolor='#d8e5f3', edgecolor='#1e4b9f', linestyle='dashdot')
plot_ellipse(ax_map, ellipse_center, major, minor, angle, facecolor='#ddffff', edgecolor='#1e4b9f', linestyle='dashdot')
plot_ellipse(ax_map, ellipse_center, major2, minor2, angle2, facecolor='#eef5f7', edgecolor='#89b2a1',
             linestyle='dashdot')

# All plots done with the same function.
# Please check and correct all parameters by yourself.
# For color change color param, for annotation number number_for_annotation parameter, and for marker correspondingly
# marker parameter.
ax_map.errorbar(ellipse_center[0], ellipse_center[1], transform=ax_map.get_transform('world'),
                markersize=3, markeredgewidth=0.2, color='red', marker='*', markerfacecolor='none')
# Plotting points
if infile_local != 'none':
   # d3good is local_values
   if len(d3good) > 0:
       plot_points(d3good, ax_map, marker='o', marker_size=markersize_local, color='green',
                   marker_edge_width=1, marker_face_color='None', zorder=13)
if infile_local2 != 'none':
   # d4good is local_values
   if len(d4good) > 0:
       plot_points(d4good, ax_map, marker='o', marker_size=markersize_local*0.85, color='#ff1aff',
                   marker_edge_width=1, marker_face_color='None', zorder=13)
if len(data) > 0:
   plot_numbers = 'False'
   hsp_color = '#ffcc99' 
   bzb_color = '#ffcc99' 
   hsp_size = 3.0
   bzb_size = 2.0
   if catsonly != 'y':
      plot_numbers = 'True'
      hsp_color = '#ffff99' 
      bzb_color = '#ffff99' 
      hsp_size = 2.5
      bzb_size = 1.0
      # UNCL (unclassified) objects.
      plot_points(duncl, ax_map, marker='$?$', marker_size=markersize_hbl*3, color='#5555aa',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl,zorder=11)
      # HBL objects.
      plot_points(dhbl_orange, ax_map, marker='o', marker_size=markersize_hbl, color='#ff9933',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl,zorder=11)
      # IBL objects.
      plot_points(dibl_cyan, ax_map, marker='o', marker_size=markersize_cyan, color='#66ccff',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl,zorder=11)
      # LBL objects.
      plot_points(dlbl_blue, ax_map, marker='o', marker_size=markersize_blue, color='#0033cc',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl,zorder=11)
      # Unknown type AGN or Non-jetted AGN or radio-AGN.
      plot_points(dother_black, ax_map, marker=(12,1,0), marker_size=markersize/2, color='#555566',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl)
      # iX-ray emitting Star  
      plot_points(dother_star, ax_map, marker=(12,1,0), marker_size=markersize/2, color='green',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl)
      if plot_stars == 'y':
          # plot stars from HSTGSC catalog
          plot_points(dother_gray, ax_map, marker=',', marker_size=1, color='gray',
                      marker_edge_width=marker_edge_width_hbl_ibl_lbl)
      # catalogues 
      plot_points(dother_cat, ax_map, marker='.', marker_size=markersize*0.6, color='#1a48ff',
                  marker_edge_width=marker_edge_width, marker_face_color='None')
      # Type 8 x candidate source.
      plot_points(dxcandidate_x_values, ax_map, marker='o', marker_size=markersize*0.6, color='#1a48ff',
                  marker_edge_width=marker_edge_width, marker_face_color='None')
      # Type 9 radio candidate source.
      plot_points(drcandidate_r_values, ax_map, marker='o', marker_size=markersize_blue, color='#e70b5c',
                  marker_edge_width=marker_edge_width)
      # Crates or 5BZB catalog.
      plot_points(dcrates_red, ax_map, marker='o', marker_size=markersize_blue, color='#e70b5c',
                  marker_edge_width=marker_edge_width)
      plot_points(dcrates_red, ax_map, marker='s', marker_size=4, color='#6267ca',
                  marker_edge_width=0.5, marker_face_color='None', numbering=plot_numbers, zorder=15)
   # 3HSP catalog objects.
   plot_points(dhspcat, ax_map, marker='*', marker_size=markersize_catalogued * hsp_size, color=hsp_color,
               marker_edge_width=0.5, marker_face_color=hsp_color, numbering=plot_numbers, zorder=10)
   # 5BZ catalog objects.
   plot_points(dbzbcat, ax_map, marker='D', marker_size=markersize_catalogued * bzb_size, color=bzb_color,
               marker_edge_width=0.5, marker_face_color=bzb_color, numbering=plot_numbers, zorder=15)
if len(data2) > 0:
   if catsonly != 'y':
      if plot_stars == 'y':
          plot_points(dothersmall_gray, ax_map, marker=(4,1,0), marker_size=markersize_orange * 2, color='gray',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl)
          plot_points(dothermedium_gray, ax_map, marker=(4,1,0), marker_size=markersize_orange * 3, color='gray',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl)
          plot_points(dotherlarge_gray, ax_map, marker=(4,1,0), marker_size=markersize_orange * 6, color='gray',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl)
      plot_points(dothersverysmall_agn, ax_map, marker=(8,1,0), marker_size=markersize_orange * 1.5, color='#ff0066',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl, zorder=10)
      plot_points(dothersmall_agn, ax_map, marker=(8,1,0), marker_size=markersize_orange * 3.0, color='#ff0066',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl, zorder=10)
      plot_points(dothermedium_agn, ax_map, marker=(8,1,0), marker_size=markersize_orange * 4, color='#ff0066',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl, zorder=10)
      plot_points(dotherlarge_agn, ax_map, marker=(8,1,0), marker_size=markersize_orange * 6, color='#ff0066',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl, zorder=10)
      plot_points(dothersmall_star, ax_map, marker=(4,1,0), marker_size=markersize_orange * 2, color='#00cc66',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl)
      plot_points(dothermedium_star, ax_map, marker=(4,1,0), marker_size=markersize_orange * 3, color='#00cc66',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl)
      plot_points(dotherlarge_star, ax_map, marker=(4,1,0), marker_size=markersize_orange * 6, color='#00cc66',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl)
      plot_points(dhblsmall_orange, ax_map, marker=(12,1,0), marker_size=markersize_orange * 3, color='#ff9933',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl)
      plot_points(dhblmedium_orange, ax_map, marker=(12,1,0), marker_size=markersize_orange * 4, color='#ff9933',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl)
      plot_points(dhbllarge_orange, ax_map, marker=(12,1,0), marker_size=markersize_orange * 6, color='#ff9933',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl)
      plot_points(diblsmall_cyan, ax_map, marker=(12,1,0), marker_size=markersize_cyan * 3, color='#66ccff',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl)
      plot_points(diblmedium_cyan, ax_map, marker=(12,1,0), marker_size=markersize_cyan * 4, color='#66ccff',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl)
      plot_points(dibllarge_cyan, ax_map, marker=(12,1,0), marker_size=markersize_cyan * 6, color='#66ccff',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl)
      plot_points(dlblsmall_blue, ax_map, marker=(12,1,0), marker_size=markersize_blue * 3, color='#0033cc',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl)
      plot_points(dlblmedium_blue, ax_map, marker=(12,1,0), marker_size=markersize_blue * 4, color='#0033cc',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl)
      plot_points(dlbllarge_blue, ax_map, marker=(12,1,0), marker_size=markersize_blue * 6, color='#0033cc',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl)
      plot_points(dother_black_small, ax_map, marker=(12,1,0), marker_size=markersize*0.5, color='#555566',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl)
      plot_points(dother_black_medium, ax_map, marker=(12,1,0), marker_size=markersize*0.8, color='#555566',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl)
      plot_points(dother_black_large, ax_map, marker=(12,1,0), marker_size=markersize*1.2, color='#555566',
                  marker_edge_width=marker_edge_width_hbl_ibl_lbl)
      plot_points(dcrates_small, ax_map, marker='o', marker_size=markersize_blue * 2, color='#e70b5c',
                  marker_edge_width=marker_edge_width)
      plot_points(dcrates_medium, ax_map, marker='o', marker_size=markersize_blue * 3, color='#e70b5c',
                  marker_edge_width=marker_edge_width)
      plot_points(dcrates_large, ax_map, marker='o', marker_size=markersize_blue * 5, color='#e70b5c',
                  marker_edge_width=marker_edge_width)
      plot_points(dradio_small, ax_map, marker='o', marker_size=markersize_blue * 1.0, color='#ff0066',
                  marker_edge_width=marker_edge_width, marker_face_color='None')
      plot_points(dradio_medium, ax_map, marker='o', marker_size=markersize_blue * 2.0, color='#ff0066',
                  marker_edge_width=marker_edge_width, marker_face_color='None')
      plot_points(dradio_large, ax_map, marker='o', marker_size=markersize_blue * 3.2, color='#ff0066',
                  marker_edge_width=marker_edge_width, marker_face_color='None')
      plot_points(dxray_small, ax_map, marker='x', marker_size=markersize_blue * 1.0, color='#0033cc',
                  marker_edge_width=marker_edge_width)
      plot_points(dxray_medium, ax_map, marker='x', marker_size=markersize_blue * 1.4, color='#0033cc',
                  marker_edge_width=marker_edge_width)
      plot_points(dxray_large, ax_map, marker='x', marker_size=markersize_blue * 2.4, color='#0033cc',
                  marker_edge_width=marker_edge_width)
      # Cluster of galaxies
      plot_points(dclusters, ax_map, marker='o', marker_size=markersize_catalogued * 3.0, color='#993300',
                  marker_edge_width=marker_edge_width, zorder=14, marker_face_color='None')
      plot_points(dclusters, ax_map, marker='$Clu$', marker_size=markersize_catalogued * 4.0, color='#993300',
                  marker_edge_width=marker_edge_width, zorder=15, marker_face_color='None')
      # Milliquas sources
      plot_points(dmilliquas, ax_map, marker='o', marker_size=markersize_catalogued, color='#84e184',
                  marker_edge_width=marker_edge_width, marker_face_color='None')
      # GRB sources
      plot_points(dgrb, ax_map, marker='H', marker_size=markersize_catalogued*3, color='#555555',
                  marker_edge_width=marker_edge_width, marker_face_color='None')
      # 4LAC 
      plot_points(d4lac, ax_map, marker='o', marker_size=markersize_local*0.85, color='#ff1aff',
                  marker_edge_width=marker_edge_width, marker_face_color='None',zorder=13)
      # Stellar clusters
      plot_points(dstellar_cluster, ax_map, marker='o', marker_size=markersize_catalogued*3, color='#ff0066',
                  marker_edge_width=0.5, numbering=False, zorder=20, marker_face_color='None')
      plot_points(dstellar_cluster, ax_map, marker='$SC$', marker_size=markersize_catalogued*3, color='#ff0066',
                  marker_edge_width=0.5, numbering=False, zorder=20, marker_face_color='None')
      # Cataclysmic variabile
      plot_points(dcv, ax_map, marker='o', marker_size=markersize_catalogued*3, color='#993300',
                  marker_edge_width=0.5, numbering=False, zorder=14, marker_face_color='None')
      plot_points(dcv, ax_map, marker='$CV$', marker_size=markersize_catalogued*3, color='#993300',
                  marker_edge_width=0.5, numbering=False, zorder=15, marker_face_color='None')
      # SNR
      plot_points(dsnr, ax_map, marker='o', marker_size=markersize_catalogued*3, color='#ff0066',
                  marker_edge_width=0.5, numbering=False, zorder=20, marker_face_color='None')
      plot_points(dsnr, ax_map, marker='$SNR$', marker_size=markersize_catalogued*3, color='#ff0066',
                  marker_edge_width=0.5, numbering=False, zorder=20, marker_face_color='None')
      # X-ray binaries 
      plot_points(dxray_binaries, ax_map, marker='o', marker_size=markersize_catalogued*3, color='#ff0066',
                  marker_edge_width=0.5, numbering=False, zorder=20, marker_face_color='None')
      plot_points(dxray_binaries, ax_map, marker='$XRB$', marker_size=markersize_catalogued*4, color='#ff0066',
                  marker_edge_width=0.5, numbering=False, zorder=20, marker_face_color='None')
      # Molecular clouds 
      plot_points(dmolecular_clouds, ax_map, marker='o', marker_size=markersize_catalogued*3, color='#ff0066',
                  marker_edge_width=0.5, numbering=False, zorder=20, marker_face_color='None')
      plot_points(dmolecular_clouds, ax_map, marker='$MC$', marker_size=markersize_catalogued*3, color='#ff0066',
                  marker_edge_width=0.5, numbering=False, zorder=20, marker_face_color='None')
   # Pulsars
   plot_points(dpulsar_purple, ax_map, marker='p', marker_size=markersize_catalogued*2, color='#8000ff',
               marker_edge_width=marker_edge_width)
   # Gamma
   plot_points(dgamma, ax_map, marker='^', marker_size=markersize_gamma, color='#7e0de8',
               marker_edge_width=0.5, numbering=False, zorder=15, marker_face_color='None')
   if len(dICgood) > 0:
    plot_points(dICgood, ax_map, marker='*', marker_size=markersize_local*0.7, color='red',
                marker_edge_width=1, marker_face_color='None', zorder=16)

plt.savefig(outfile, bbox_inches='tight', format='png')
