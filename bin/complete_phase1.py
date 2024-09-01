import csv
import math
import argparse
import sys
import os
import numpy as np
from scipy.optimize import curve_fit

if len(sys.argv) > 1:
    ra_center = float(sys.argv[1])
    dec_center = float(sys.argv[2])
    flag =sys.argv[3]
    #gammaray = True
    gammaray = False
else:
    print("No arguments (RA Dec) provided.")
    exit()

path ='tmp/'
sed_file = path+'Sed_temp.txt'

input_file_hstgsc = path+'hstgsc_out.txt'
input_file_panstarrs = path+'panstarrs.i.csv'
input_file_gaia = path+'gaia2.i.csv'
input_file_candidates = path+'candidates.csv'
input_file_erass1 = path+'erass1.1.csv'
input_file_efeds = path+'efeds.1.csv'
input_file_4xmm = path+'4xmmdr13.1.csv'
input_file_4fgl = path+'4fgldr4.1.csv'
input_file_4lac = path+'4lacdr3.1.csv'
input_file_wise = path+'viz-wise.csv'
output_file_candidates = path+'phase1_candidates.csv'
output_file_slopes = path+'slopes.csv'
file_find_out = path+'phase1_find_out.txt'
file_find_out_temp = path+'find_out_temp.txt'
input_file_candidates_int = path+'candidates_int.txt'
input_file_sdss = path+'sdss_out.txt'
temporary_file = path+'temp.txt'

def wisemag2flux(m_band, filter):
    if (filter[:3] != 'ww1' and filter[:3] != 'ww2' and filter[:3] != 'ww3' and filter[:3] != 'ww4'):
        print("wisemag2flux: Filter not supported")
        return
    av = 0
    if filter[:3] == 'ww1':
        lambda_ =34000.
        const=np.log10(309.540)-23.
    elif filter[:3] == 'ww2':
        lambda_ =46000.
        const=np.log10(171.787)-23.
    elif filter[:3] == 'ww3':
        lambda_ =120000.
        const=np.log10(31.674)-23.
    elif filter[:3] == 'ww4':
        lambda_ =220000.
        const=np.log10(8.363)-23.
    else:
        print(' wisemag2flux: Filter not supported  ')
        return
    a_band = 0
    c = 2.9979e10
    a = 1.0
    frequency = c / (lambda_ * 1.e-8)
# Extinction law
    flux = 10. ** (-0.4 * (m_band - a_band) + const) * frequency
    if m_band == 0:
        flux = 0
    return flux,frequency


def log_model_function(log_x, m, c):
    return m * log_x + c

def check_host_galaxy_sdss(ra,dec):
    host_galaxy_image = False 
    if os.path.exists(input_file_sdss):
        for srows in sdss_data:
            ra_sdss = float(srows[1])
            dec_sdss = float(srows[2])
            sdss_class = srows[13]
            min_dists = 5
            dists = sky_distance(ra, dec, ra_sdss, dec_sdss)
            dists *=3600
            if dists < min_dists and sdss_class =='3':
                host_galaxy_image = True
                break 
    if host_galaxy_image:
       return True
    else:
       return False

def check_host_galaxy_panstarrs(ra,dec):
    host_galaxy_image = False 
    if os.path.exists(input_file_panstarrs):
        for row in panstarrs_data:
            ra_panstarrs = float(row[0])
            dec_panstarrs = float(row[1])
            magnitude = float(row[4])
            gmag = float(row[4])
            rmag = float(row[6])
            kron_gmag = float(row[15])
            kron_rmag = float(row[16])
            deltegmag= -99.
            if rmag > 0 and kron_rmag >0:
                deltarmag = rmag - kron_rmag
            else:
                deltarmag = -99
            if gmag > 0 and kron_gmag > 0:
                deltagmag = gmag - kron_gmag
            else:
                deltagmag = -99
            dists = sky_distance(ra, dec, ra_panstarrs, dec_panstarrs)
            dists *=3600
            if dists < 5:
                if deltagmag > 0.05 or deltarmag > 0.05:
                   host_galaxy_image = True
                break 
    else:
       return False
    if host_galaxy_image:
       return True
    else:
       return False

def read_sed_data(file_path,ra,dec):
    # Read data from the file and filter rows containing "E+17"
    # if the row includes code 51 data if from XMM-slew
    # if the row includes code 2 or 52 data if from XMM-DR13
    # if the row includes code 5 or 55 data if from Swift-2SXPS
    # if the row includes code 14 or 64 data if from erosita
    data = []
    data_xmm = []
    data_swift = []
    data_erosita = []
    with open(file_path, 'r') as file:
        dist = 100
        for line in file:
            if "matched" in line:
                c = line.strip().split()
                ra_source = float(c[3])
                dec_source = float(c[4])
                dist = sky_distance(ra, dec, ra_source, dec_source) 
                dist *=3600
                line = next(file)
            maxd = 9
            if dist < maxd:
                check = line[-3:].replace("\n", "")
                if "E+16" in line or"E+17" in line or "E+18" in line:
                    columns = line.strip().split()
                    if float(columns[1]) > 0:
                        data.append([float(columns[0]), float(columns[1]), (float(columns[2])-float(columns[3]))/2])
                if ("E+16" in line or "E+17" in line or "E+18" in line) and (check == " 2" or check == "52"):
                    columns = line.strip().split()
                    if float(columns[1]) > 0:
                        snr = 1
                        if float(columns[2]) > float(columns[3]):
                           snr = float(columns[1])/((float(columns[2])-float(columns[3]))/2)
                        if snr > 2.0:
                           data_xmm.append([float(columns[0]), float(columns[1]), (float(columns[2])-float(columns[3]))/2])
                if ("E+16" in line or "E+17" in line or "E+18" in line) and (check == "14" or check == "64"):
                    columns = line.strip().split()
                    if float(columns[1]) > 0:
                        error = (float(columns[2])-float(columns[3]))/2
                        snr = float(columns[1])/((float(columns[2])-float(columns[3]))/2)
                        if snr > 1.6:
                           data_erosita.append([float(columns[0]), float(columns[1]), (float(columns[2])-float(columns[3]))/2])
                if ("E+16" in line or "E+17" in line or "E+18" in line) and (check == " 5" or check == "55"):
                    columns = line.strip().split()
                    if float(columns[1]) > 0:
                        snr = float(columns[1])/((float(columns[2])-float(columns[3]))/2)
                        if snr > 1.6:
                            data_swift.append([float(columns[0]), float(columns[1]), (float(columns[2])-float(columns[3]))/2])
    if len(data_xmm) > 0:
        first_column = [row[0] for row in data_xmm]
        min_value = min(first_column)
        max_value = max(first_column)
    if len(data_xmm) > 2 and min_value != max_value:
        #print("XMM")
        return np.array(data_xmm)
    elif len(data_swift) > 2:
        #print("Swift")
        return np.array(data_swift)
    elif len(data_erosita) > 2:
        #print("eROSITA")
        return np.array(data_erosita)
    else:
        #print("else")
        return np.array(data)

def least_squares_fit(x, y):
    # Perform least squares linear fit in log space
    log_x = np.log10(x)
    log_y = np.log10(y)
    A = np.vstack([log_x, np.ones_like(log_x)]).T
    m, c = np.linalg.lstsq(A, log_y, rcond=None)[0]
    return m, c

def get_x_slope(sed_file,ra,dec):
    if not os.path.exists(sed_file):
       return -99
    data = read_sed_data(sed_file,ra,dec)
    #print("x-ray data ",data)
    x_slope = -99
    if data.size == 0:
       return x_slope
    if len(data[:,0]) > 0:
        Ma = max(data[:,0])
        Mi = min(data[:,0])
        if Ma > Mi: 
            # Separate data into x (first column), y (second column), and error (statistical uncertainty)
            x = data[:, 0]
            y = data[:, 1]
            error = data[:, 2]

            log_x = np.log(x)
            log_y = np.log(y)
            params, covariance = curve_fit(log_model_function, log_x, log_y, sigma=error/y)
            # Extract the fitted parameters
            m, c = params
            variances = np.diag(covariance)
            statistical_errors = np.sqrt(variances)
            #print("errors ",statistical_errors)

            # Display the slope and intercept of the linear fit
            #print(f"Linear Fit: y = 10^({m:.4f}log10(x) + {c:.4f})")
            x_slope = -m + 1
            #print("x_slope +/- error ",x_slope,statistical_errors[0])
            if statistical_errors[0] < 0.0001 or statistical_errors[0] > 0.3:
               x_slope = -99
    return x_slope

def check_if_incandidates(ra,dec):
    min_dist = 10000.
    phase1_candidates = np.load('candidates_data.npy')
    for row in phase1_candidates:
        ra_cand=float(row[1])
        dec_cand=float(row[2])
        dist = sky_distance(ra, dec, ra_cand, dec_cand)
        dist = dist*3600
        if dist < min_dist:
            min_dist = dist
            type=row[3]
        if min_dist < 2 and (type!='type-8' and type!='type-9' or type=='radio-AGN' or type=='X-Star'):
            return False
        else:
            return True 

# Function to calculate the sky distance between two points
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


if os.path.exists(input_file_sdss):
    sdss_data = []
    with open(input_file_sdss, mode='r') as sdss:
        csv_reader_sdss = csv.reader(sdss)
        header = next(csv_reader_sdss)  # Skip the header row
        for row in csv_reader_sdss:
            sdss_data.append(row)

if os.path.exists(input_file_panstarrs):
    panstarrs_data = []
    with open(input_file_panstarrs, mode='r') as panstarrs:
        csv_reader_panstarrs = csv.reader(panstarrs)
        header = next(csv_reader_panstarrs)  # Skip the header row
        for row in csv_reader_panstarrs:
            panstarrs_data.append(row)

with open(file_find_out_temp, 'r') as fo, open(temporary_file, 'w') as te:
    for line in fo:
        parts = line.split('"')
        if len(parts) > 1:
            parts[1] = parts[1].replace(" ", "")
            parts[0] = parts[0].replace('"', '')
            parts[1] = parts[1].replace('"', '')
        modified_line = ''.join(parts)
        elements = modified_line.split()
        newline = ' '.join(elements)
        te.write(newline+'\n')

# Read the modified file and create a NumPy array
with open(temporary_file, mode='r') as ffo:
    lines = ffo.readlines()
# Convert to NumPy array
data = np.array(lines)
np.save('findout_data.npy', data)

if os.path.exists(input_file_gaia):
    input_file_optical=input_file_gaia
elif os.path.exists(input_file_panstarrs):
    input_file_optical=input_file_panstarrs
elif os.path.exists(input_file_hstgsc):
    input_file_optical=input_file_hstgsc
else:
   print(f"No optical sources file found")
   exit()

#print("input_file_optical ",input_file_optical)

if not os.path.exists(input_file_candidates):
    cand_phase1 = False
else:
   with open(input_file_candidates, 'r') as filec:
      linesc = filec.readlines()
      if len(linesc) > 1:
          cand_phase1 = True
      else:
          cand_phase1 = False 

if not os.path.exists(input_file_candidates_int):
    cand_int = False
else:
   with open(input_file_candidates_int, 'r') as file:
      lines = file.readlines()
      if len(lines) > 1:
          cand_int = True
      else:
          cand_int = False 

if os.path.exists(input_file_candidates):
    with open(input_file_candidates, mode='r') as cand:
        csv_reader_candd = csv.reader(cand)
        header = next(csv_reader_candd)  # Skip the header row
        data = list(csv_reader_candd)
# Save the array to a binary file using numpy
    np.save('candidates_data.npy', data)

count = 0 

if cand_phase1:
    with open(input_file_candidates, mode='r') as cand:
        csv_reader_cand = csv.reader(cand)
        header = next(csv_reader_cand)  
        with open(output_file_candidates, mode='a+', buffering=1, encoding='utf-8', errors=None, newline=None, closefd=True, opener=os.open) as out_file_candidates, open(file_find_out_temp, mode='w') as find_out_temp, open(file_find_out, mode='w') as out_file_find_out, open(output_file_slopes, mode='w') as out_slopes:
            out_slopes.write('ra,dec,aox,airx,airo,aro,arx,aw1w2,arg,nu_peak,x_en_slope,gamma_sp_index,host_galaxy_sed,host_galaxy_image,blue_bump,Xray_extended\n')
            out_file_candidates.write('Number,R.A.,Dec.,Type,Name\n')
            for row in csv_reader_cand:
               count += 1 
               ra = float(row[1])
               dec = float(row[2])
               cat = row[3]
               if cat == 'type-9':
                   cat = 'NoXrayBlazarCandidate'
               obj_id = row[4]
               x_slope = get_x_slope(sed_file,ra,dec)
               extended = False 
               found = False
               nupeak = -99.
               if os.path.exists(input_file_erass1):
                   with open(input_file_erass1, mode='r') as erass1:
                       csv_reader_erass1 = csv.reader(erass1)
                       header = next(csv_reader_erass1)  # Skip the header row
                       for row_erass1 in csv_reader_erass1:
                           ra_erass1 = float(row_erass1[1])
                           dec_erass1 = float(row_erass1[2])
                           erass1_MLFlux1 = float(row_erass1[6])
                           extension = float(row_erass1[29])
                           dist = sky_distance(ra, dec, ra_erass1, dec_erass1)
                           distarcsec = dist * 3600
                           max_extension = 8
                           if erass1_MLFlux1 > 5.e-11:
                              max_extension = 20
                           if extension > max_extension and distarcsec < max(10,extension/2):
                               extended = True
               if os.path.exists(input_file_efeds):
                   with open(input_file_efeds, mode='r') as efeds:
                       csv_reader_efeds = csv.reader(efeds)
                       header = next(csv_reader_efeds)  # Skip the header row
                       for row_efeds in csv_reader_efeds:
                           ra_efeds = float(row_efeds[1])
                           dec_efeds = float(row_efeds[2])
                           extension = float(row_efeds[31])
                           dist = sky_distance(ra, dec, ra_efeds, dec_efeds)
                           distarcsec = dist * 3600
                           #if distarcsec < 7 and extension > 8:
                           if extension > 20 and distarcsec < max(10,extension/2):
                               extended = True
               if os.path.exists(input_file_4xmm):
                   with open(input_file_4xmm, mode='r') as xmm:
                       csv_reader_4xmm = csv.reader(xmm)
                       header = next(csv_reader_4xmm)  # Skip the header row
                       for row_4xmm in csv_reader_4xmm:
                           ra_4xmm = float(row_4xmm[1])
                           dec_4xmm = float(row_4xmm[2])
                           extension = float(row_4xmm[16])
                           dist = sky_distance(ra, dec, ra_4xmm, dec_4xmm)
                           distarcsec = dist * 3600
                           #if distarcsec < 7 and extension > 8:
                           if extension > 20 and distarcsec < max(10,extension/2):
                               #print("extension ",extension)
                               extended = True
               slope_4lac_4fgl = -99
               fermi4fgl_name = '--'
               if os.path.exists(input_file_4fgl):
                   with open(input_file_4fgl, mode='r') as fgl:
                       csv_reader_4fgl = csv.reader(fgl)
                       header = next(csv_reader_4fgl)  # Skip the header row
                       for row_4fgl in csv_reader_4fgl:
                           if gammaray:
                              ra_4fgl=float(row_4fgl[1])
                              dec_4fgl=float(row_4fgl[2])
                              distdeg = sky_distance(ra, dec, ra_4fgl, dec_4fgl)
                              semimajor_axis = float(row_4fgl[5])
                              if distdeg < semimajor_axis and row_4fgl[7] != '':
                                   slope_4lac_4fgl = float(row_4fgl[7])
                           if row_4fgl[10] != '':
                              ra_counterpart = float(row_4fgl[11])
                              dec_counterpart = float(row_4fgl[12])
                              dist = sky_distance(ra, dec, ra_counterpart, dec_counterpart)
                              distarcsec = dist * 3600
                              if distarcsec < 7:
                                  fermi4fgl_name = row_4fgl[12].replace(" ", "")
               if os.path.exists(input_file_4lac):
                   with open(input_file_4lac, mode='r') as lac:
                       csv_reader_4lac = csv.reader(lac)
                       header = next(csv_reader_4lac)  # Skip the header row
                       for row_4lac in csv_reader_4lac:
                           ra_4lac = float(row_4lac[1])
                           dec_4lac = float(row_4lac[2])
                           dist = sky_distance(ra, dec, ra_4lac, dec_4lac)
                           distarcsec = dist * 3600
                           if distarcsec < 7:
                               slope_4lac_4fgl = float(row_4lac[4])
                               slope_4lac_4fgl_error = float(row_4lac[5])
               min_distw = 7
               w1 = -99
               w1_err = -99
               w2 = -99
               w2_err = -99
               w3 = -99
               w3_err = -99
               w4 = -99
               w4_err = -99
               closest_source = 10000000. 
               aw1w3 = -99
               aw1w4 = -99
               aw2w4 = -99
               flux1 = -99 
               flux2 = -99
               flux3 = -99
               flux4 = -99
               if os.path.exists(input_file_wise):
                   with open(input_file_wise, mode='r') as wise:
                       csv_reader_wise = csv.reader(wise)
                       header = next(csv_reader_wise)  # Skip the header row
                       #w1 = -99
                       #w1_err = -99
                       #w2 = -99
                       #w2_err = -99
                       #w3 = -99
                       #w3_err = -99
                       #w4 = -99
                       #w4_err = -99
                       #closest_source = 10000000. 
                       #aw1w3 = -99
                       #aw1w4 = -99
                       #aw2w4 = -99
                       #flux1 = -99 
                       #flux2 = -99
                       #flux3 = -99
                       #flux4 = -99
                       for row_wise in csv_reader_wise:
                           ra_wise = float(row_wise[1])
                           dec_wise = float(row_wise[2])
                           #print("ra_wise dec_wise ",ra_wise,dec_wise)
                           if row_wise[3] != '':
                              w1 = float(row_wise[3])
                           if row_wise[4] != '':
                              w1_err = float(row_wise[4])
                           if row_wise[5] != '':
                              w2 = float(row_wise[5])
                           if row_wise[6] != '':
                              w2_err = float(row_wise[6])
                           if row_wise[7] != '':
                              w3 = float(row_wise[7])
                           if row_wise[8] != '':
                              w3_err = float(row_wise[8])
                           if row_wise[9] != '':
                              w4 = float(row_wise[9])
                           if row_wise[10] != '':
                              w4_err = float(row_wise[10])
                           variability_flags = 'nnnn'
                           if len(row_wise) > 11:
                              if row_wise[11] != '':
                                 variability_flags = row_wise[11]
                           #print("W1,W2,W3,W4 ",w1,w2,w3,w4)
                           dist = sky_distance(ra, dec, ra_wise, dec_wise)
                           distarcsec = dist * 3600
                           if w1 < 10:
                              min_distw = 10
                           if w4 < 5:
                              min_distw = 15
                           #print("distarcsec  closest_source min_distw ",distarcsec,closest_source,min_distw)
                           if distarcsec < min_distw and distarcsec < closest_source:
                               closest_source = distarcsec
                               if '7' in variability_flags or '8' in variability_flags or '9' in variability_flags:
                                  print ("------- WISE variable source ---- variability_flags:",variability_flags) 
                               if w1 > -90:
                                  fluxw1, freq_w1  = wisemag2flux(w1,'ww1')
                               if w2 > -90:
                                  fluxw2, freq_w2  = wisemag2flux(w2,'ww2')
                               if w3 > -90:
                                  fluxw3, freq_w3  = wisemag2flux(w3,'ww3')
                               if w4 > -90:
                                  fluxw4, freq_w4  = wisemag2flux(w4,'ww4')
                               if fluxw1 > 0. and fluxw3 > 1.5e-12:
                                   aw1w3 = np.log10(fluxw1/fluxw3)/np.log10(freq_w1/freq_w3)
                               else:
                                   aw1w3 = -99
                               if fluxw1 > 0. and fluxw4 > 4.e-12 and aw1w3 > -90:
                                   aw1w4 = np.log10(fluxw1/fluxw4)/np.log10(freq_w1/freq_w4)
                               if fluxw2 > 0. and fluxw4 > 4.e-12 and aw1w3 > -90:
                                   aw2w4 = np.log10(fluxw2/fluxw4)/np.log10(freq_w2/freq_w4)
                               #print("1-aw1w3 aw1w4 aw2w4 ",aw1w3,aw1w4,aw2w4)
               if cand_int:
                   with open(input_file_candidates_int, mode='r') as c_int:
                       csv_reader_cand_int = csv.reader(c_int)
                       next(csv_reader_cand_int)
                       ir = 0
                       iox = 0
                       iaro = 0
                       iro = 0
                       irx = 0
                       iw = 0
                       irg = 0
                       airxav = 0
                       aoxav = 0. 
                       aroav = 0. 
                       airoav = 0
                       arxav = 0
                       aw1w2av = 0
                       argav = 0
                       changed = False
                       for row_int in csv_reader_cand_int:
                           aox = float(row_int[0])
                           airx = float(row_int[1])
                           airo = float(row_int[2])
                           aro = float(row_int[3])
                           arx = float(row_int[4])
                           aw1w2 = float(row_int[5])
                           arg = float(row_int[6])
                           ra_int = float(row_int[7])
                           dec_int = float(row_int[8])
                           dist = sky_distance(ra, dec, ra_int, dec_int)
                           distarcsec = dist * 3600
                           if distarcsec < min_distw: 
                               #print("distarcsec min_distw ",distarcsec,min_distw)
                               found = True
                               if airx > -90:
                                   ir += 1
                                   airxav += airx
                               if aox > -90:
                                   iox += 1
                                   aoxav += aox
                               if airo > -90:
                                   iro += 1
                                   airoav += airo
                               if aro > -90:
                                   iaro += 1
                                   aroav += aro
                               if arx > -90:
                                   irx += 1
                                   arxav += arx
                               if aw1w2 > -90:
                                   iw += 1
                                   aw1w2av += aw1w2
                               if arg > -90:
                                   irg += 1
                                   argav += arg
                       blue_bump = False 
                       host_galaxy_sed = False 
                       host_galaxy_image = False 
                       if found:
                           if ir > 0:
                              airxav = airxav/ir
                           else:
                              airxav = -99
                           if iro > 0:
                              airoav = airoav/iro
                           else:
                              airoav = -99
                           if iox > 0:
                              aoxav = aoxav/iox
                           else:
                              aoxav = -99
                           if iaro > 0:
                              aroav = aroav/iaro
                           if irx > 0:
                              arxav = arxav/irx
                           else:
                              arxav = -99
                           if iw > 0:
                              aw1w2av = aw1w2av/iw
                           else:
                              aw1w2av = -99
                           if irg > 0:
                              argav = argav/irg
                           else:
                              argav = -99
                           x_slope = get_x_slope(sed_file,ra,dec)
                           newcat = 'Xray-UNCL'
                           #print("aoxav airxav airoav aoxav aroav arxav aw1w2av argav ",aoxav,airxav,airoav,aoxav,aroav,arxav,aw1w2av,argav)
                           host_galaxy_image = check_host_galaxy_sdss(ra,dec)
                           if not host_galaxy_image:
                               host_galaxy_image = check_host_galaxy_panstarrs(ra,dec)
                           if argav < 0.95 and (aw1w2av < -0.1 or abs(x_slope) < 0.9 or argav < 0.5):
                              if not host_galaxy_image:
                                 blue_bump = True
                                 #print("Possible presence of a blue bump based on r-g color",argav)
                           elif argav > 2:
                               host_galaxy_sed = True
                           if arxav > 0.38 and arxav < 0.85:
                               if host_galaxy_image: # host galaxy increases opt flux
                                  airxav_max = 1.55
                                  aoxav_max = 1.8
                               if (airoav > 0 and airoav < 0.951) or aw1w2av > 1.5:  # host galaxy increases opt flux 
                                  deltairo_irx = airxav - airoav
                                  if aw1w2av < 1.5 and deltairo_irx > 0.16 and airxav > 1:
                                     airxav_max = 1.55
                                     aoxav_max = 1.8
                                     if argav < 0.5 and (aw1w2av < .3 or abs(x_slope) < 0.9):
                                        #print("Possible presence of a blue bump 2")
                                        blue_bump = True
                                     if argav < 0.94 and (aw1w2av < -0.1 or abs(x_slope) < 0.9 or airoav < 0.9):
                                        #print("Possible presence of a blue bump 3")
                                        blue_bump = True
                                  elif aw1w2av < 2.8:
                                     airxav_max = 1.6
                                     aoxav_max = 2.0
                                     if argav < 0.5 and (aw1w2av < .5 or abs(x_slope) < 0.9):
                                        #print("Possible presence of a blue bump 4")
                                        blue_bump = True
                                     if argav < 1.0 and (aw1w2av < -0.1 or abs(x_slope) < 0.9):
                                        #print("Possible presence of a blue bump 5")
                                        blue_bump = True
                                     if aw1w2av > 1.8:
                                         host_galaxy_sed = True
                                  elif aw1w2av > 2.8:
                                     airxav_max = 1.65
                                     aoxav_max = 2.2
                                     host_galaxy_sed = True
                                  else:
                                     airxav_max = 1.5
                                     aoxav_max = 1.6
                               else:
                                  #print("Here 10")
                                  airxav_max = 1.3
                                  aoxav_max = 1.3
                               #print("airxav airxav_max aoxav aoxav_max arxav aw1w2av ",airxav,airxav_max,aoxav,aoxav_max,arxav,aw1w2av)
                               if airxav < airxav_max and aoxav < aoxav_max and arxav < 0.8 and (aw1w2av > 0 or aw1w2av < -90):
                                   #print("Here 1")
                                   if arxav < 0.8 and arxav >= 0.65:
                                      if x_slope < 0.95 and x_slope> 0.8:
                                         newcat = 'HBL-mispointed'
                                      elif x_slope < 0.8:
                                         newcat = 'LBL'
                                      else:
                                         newcat = 'HBL'
                                   elif arxav < 0.65:
                                      newcat = 'HBL'
                                   else:
                                      if x_slope > 0.95:
                                         newcat = 'HBL-mispointed?'
                                      else:
                                         newcat = 'LBL?'
                                   #if blue_bump or (aw1w2av  > -90 and aw1w2av < -0.2):
                                   if blue_bump or (aw1w2av  > -90 and aw1w2av < 0.):
                                     if arxav < 0.42 and arxav > -90:
                                        newcat = 'radio-AGN'
                                     elif arxav > 0.7:
                                        #print("Here 1")
                                        newcat = 'LBL'
                                   if ( (aw1w3 < 0 and aw1w3 > -90) or (aw1w4 < 0 and aw1w4 > -90) or (aw2w4 < 0 and aw2w4 > -90) ) and slope_4lac_4fgl < -90:
                                       #print("Here Seyf 2.5")
                                       newcat = 'Seyfert'
                               if airxav < airxav_max and aoxav < aoxav_max and arxav < 0.8 and (aw1w2av < 0 and aw1w2av > -90):
                                   #print("Here x_slope blue_bump ",x_slope,blue_bump)
                                   if (abs(x_slope) < 0.8 and aoxav > 1.0) or x_slope < -90:
                                       #print("Here 2")
                                       newcat = 'LBL'
                                   if abs(x_slope) < 0.95 and aoxav < 1.0 and ( (aw1w2av > -90 and aw1w2av < 0) or (aw1w3 < 0 and aw1w3 > -90) or (aw1w4 < 0 and aw1w4 > -90)):
                                       newcat = 'LBL'
                                   if x_slope > 0.95 and cat == 'HBL' and aw1w2av > -90 and aw1w2av < 0:
                                       if arxav < 0.7:
                                          #print("Here  Seyf 3.5")
                                          newcat = 'Seyfert'
                                       else:
                                          newcat = 'LBL'
                               elif arxav > 0.6 and arxav < 0.85 and airxav > 1 and airxav < 1.5 and aoxav > 1 and aoxav < 1.6:
                                   #print("Here 1")
                                   if blue_bump:
                                      newcat = 'LBL'
                                   else:
                                      newcat = 'IBL'
                                   if (blue_bump or aroav > 0.6) and aw1w3 < 0 and aw1w3 > -90:
                                     #print("Here 4")
                                     newcat = 'LBL'
                               elif aw1w2av < -.2 and arxav > 0.7:
                                  #print("Here 5")
                                  newcat = 'LBL'
                           elif arxav > 0.765 and (airxav > 1.0 or aoxav > 1.1 or (aw1w2av > -90 and aw1w2av < 0) or aw1w2av > 2) and abs(x_slope) < 1.1 :
                               #print("Here 6")
                               newcat = 'LBL'
                               if argav < 0.95 and (aw1w2av < -0.1 or abs(x_slope) < 0.9 or airoav < 0.9):
                                   blue_bump = True
                               if (airoav > 0 and airoav < 1.0) or aw1w2av > 2:
                                  if aw1w2av < 2 and argav < 1.0:
                                     print("Possible presence of blue bump emission")
                                     blue_bump = True
                                  else:
                                     print("Possible presence of a host galaxy or blue bump")
                           elif ((airxav > 1.05 and airoav < 0.75) or aoxav > 0.8) and arxav < 0.45 and arxav > 0: 
                               #print("Here 2")
                               newcat = 'radio-AGN'
                           if cat == 'HBL' and newcat != 'HBL' and newcat != 'LBL' and newcat != 'Seyfert':
                               if aoxav > 0.9 and arxav < 0.55 and arxav > -90: 
                                   #print("Here 3")
                                   newcat = 'radio-AGN'
                               if cat == 'HBL' and arxav < 0.42 and arxav > -90:
                                   #print("Here 4")
                                   newcat = 'radio-AGN'
                           if cat == 'HBL' or cat == 'IBL' or newcat == 'HBL' or newcat == 'IBL':
                               if aoxav > 1 and abs(x_slope) < 0.7:
                                   #print("Here 7")
                                   newcat = 'LBL'
                               if slope_4lac_4fgl > 2.4 and arxav > 0.65:
                                   #print("Here 8")
                                   newcat = 'LBL'
                           if (aoxav > 1.4 and aw1w2av > 2.0) or (aoxav > 2.0 and aoxav < 3.5) or (airoav < 0.1 and airoav != 0. and airoav > -90):
                               newcat = 'X-Star'  
                               #print("Here 1  host_galaxy_image ",host_galaxy_image)
                               if host_galaxy_image:
                                   newcat = 'Galaxy' 
                           nupeak = -99
                           if newcat == 'HBL' or newcat == 'IBL' or newcat == 'LBL':
                              if aw1w2av > -0.1 and aw1w2av < 1.6:
                                 nupeak = 3.88 * aw1w2av + 13.59
                                 if nupeak > 18:
                                     nupeak = 18.0
                              elif aw1w2av > -90 and aw1w2av < -0.1:
                                  nupeak = 1.23 * aw1w2av + 13.24
                              if host_galaxy_sed:
                                 nupeak = -99
                              #print ("W-Peak approximate estimation of nu_peak : {:.1f}".format(nupeak))
                           if obj_id == 'N.A.' and fermi4fgl_name != '--':
                               obj_id = fermi4fgl_name 
                           steepiroirx = False
                           if airoav > 1.2 or airxav > 1.:
                              steepiroirx = True
                           if newcat != 'X-Star':
                              if (newcat == 'HBL' and aw1w2av < -0.25 and aw1w2av > -90) or ((cat == 'HBL' or newcat == 'HBL') and  arxav < 0.65 and arxav > -90 and steepiroirx and abs(x_slope) < 0.9 and  slope_4lac_4fgl < -90):
                                  #print("Here 1 Seyf")
                                  if host_galaxy_image:
                                     newcat = 'Seyfert'
                                  else:
                                     newcat = 'QSO/Seyfert'
                           ir_excess = False
                           if aw1w3 < 0.3 and aw1w3 != -99:
                              ir_excess = True 
                           if aw1w4 < 0.3 and aw1w4 != -99:
                              ir_excess = True 
                           if aw2w4 < 0.3 and aw2w4 != -99:
                              ir_excess = True 
                           #print("aw1w3 aw1w4 aw2w4 newcat, cat, ir_excess slope_4lac_4fgl ",aw1w3,aw1w4,aw2w4,newcat,cat,ir_excess,slope_4lac_4fgl)
                           if (newcat == 'HBL' or cat == 'HBL') and ir_excess and airxav < 1.27 and slope_4lac_4fgl < -90 and abs(x_slope) < 1:
                               #print("Here 2 Seyf")
                               if host_galaxy_image:
                                  newcat = 'Seyfert'
                               else:
                                  newcat = 'QSO/Seyfert'
                               if aw1w4 < -0.5 and aw1w4 > -90:
                                  #print("Here 1 galaxy")
                                  newcat = 'galaxy'
                           if (cat == '3HSP' or '3HSP' in obj_id) and newcat == 'Xray-UNCL':
                                #print("Here 2")
                                newcat = 'HBL'
                           if extended:
                                if newcat == 'HBL' or newcat == 'IBL' or newcat == 'LBL' or newcat == 'radio-AGN':
                                   newcat = 'cluster?'
                                else:
                                   newcat = 'X-extended'
                           if newcat == cat:
                              print("Classification confirmed: ",cat," at position ",ra,dec)
                           else:
                              print("Classification revised to",newcat,"(from ",cat,") at position ",ra,dec)                   
                           #print("Extended 1 ",extended)
                           if newcat != 'HBL' and newcat != 'IBL' and newcat != 'LBL' and x_slope <= 0.55 and x_slope > -90:
                              newcat += '-absorbed' 
                           if newcat == 'Xray-UNCL' and blue_bump and arxav < 0.7 and x_slope < 0.7:
                              newcat = 'high-z-LBL?'
                           if ('LBL' in newcat) or ('IBL' in newcat) or ('HBL' in newcat) or ('5BZCat' in newcat): 
                              newcat +=' (blazar)'
                           if ('CRATES' in newcat):
                              newcat +=' (blazar candidate)'
                           #print("Here 1")
                           if gammaray:
                               if newcat != 'X-Star' and newcat != 'RQ-AGN' and newcat != 'Unknown' and newcat != 'type-8' and newcat != 'zw' and newcat != 'type-9' and 'zw' not in obj_id and 'mcxc' not in obj_id:
                                  out_file_candidates.write(f'{count:.0f},{ra:.5f},{dec:.5f},{newcat},{obj_id}\n')
                                  out_slopes.write(f'{ra:.5f},{dec:.5f},{aoxav:.2f},{airxav:.2f},{airoav:.2f},{aroav:.2f},{arxav:.2f},{aw1w2av:.2f},{argav:.2f},{nupeak:.1f},{x_slope:.2f},{slope_4lac_4fgl:.2f},{host_galaxy_sed},{host_galaxy_image},{blue_bump},{extended}\n')
                           else:
                               out_file_candidates.write(f'{count:.0f},{ra:.5f},{dec:.5f},{newcat},{obj_id}\n')
                               out_slopes.write(f'{ra:.5f},{dec:.5f},{aoxav:.2f},{airxav:.2f},{airoav:.2f},{aroav:.2f},{arxav:.2f},{aw1w2av:.2f},{argav:.2f},{nupeak:.1f},{x_slope:.2f},{slope_4lac_4fgl:.2f},{host_galaxy_sed},{host_galaxy_image},{blue_bump},{extended}\n')
                       else:
                           if (cat =='HBL' or cat == 'IBL') and slope_4lac_4fgl > 2.2:
                               #print("Here 7")
                               cat = 'LBL' 
                           if extended:
                                if cat == 'HBL' or cat == 'IBL' or cat == 'LBL' or cat == 'radio-AGN':
                                   cat = 'cluster?'
                                else:
                                   cat = 'X-extended'
                           #print("Here -3 cat obj_id ",cat,obj_id)
                           obj_id = 'N.A.'
                           if ('LBL' in cat) or ('IBL' in cat) or ('HBL' in cat) or ('5BZCat' in cat): 
                              cat +=' (blazar)'
                           if ('CRATES' in cat):
                              cat +=' (blazar candidate)'
                           #print("Here 2")
                           if gammaray:
                               if cat != 'X-Star' and cat != 'RQ-AGN' and cat != 'Xray-UNCL' and cat != 'Unknown' and cat != 'type-8' and cat != 'zw' and cat != 'type-9' and 'zw' not in obj_id:
                                  #print(f'{count:.0f},{ra:.5f},{dec:.5f},{cat},{obj_id}\n')
                                  out_file_candidates.write(f'{count:.0f},{ra:.5f},{dec:.5f},{cat},{obj_id}\n')
                                  out_slopes.write(f'{ra:.5f},{dec:.5f},{aoxav:.2f},{airxav:.2f},{airoav:.2f},{aroav:.2f},{arxav:.2f},{aw1w2av:.2f},{argav:.2f},{nupeak:.1f},{x_slope:.2f},{slope_4lac_4fgl:.2f},{host_galaxy_sed},{host_galaxy_image},{blue_bump},{extended}\n')
                           else:
                              out_file_candidates.write(f'{count:.0f},{ra:.5f},{dec:.5f},{cat},{obj_id}\n')
                              out_slopes.write(f'{ra:.5f},{dec:.5f},{aoxav:.2f},{airxav:.2f},{airoav:.2f},{aroav:.2f},{arxav:.2f},{aw1w2av:.2f},{argav:.2f},{nupeak:.1f},{x_slope:.2f},{slope_4lac_4fgl:.2f},{host_galaxy_sed},{host_galaxy_image},{blue_bump},{extended}\n')
               else:
                  host_galaxy_image = check_host_galaxy_sdss(ra,dec)
                  if not host_galaxy_image:
                     host_galaxy_image = check_host_galaxy_panstarrs(ra,dec)
                  if obj_id == 'type-9':
                      obj_id = 'NoXrayBlazarCandidate'
                  #print("Here no candidates in the intermediate phase host_galaxy_image ",host_galaxy_image)
                  #print("Extended 3",extended)
                  if extended:
                     if cat == 'HBL' or cat == 'IBL' or cat == 'LBL' or cat == 'radio-AGN':
                        cat = 'cluster?'
                     else:
                        cat = 'X-extended'
                  #print("Here --3 ra dec ",ra,dec)
                  if ('LBL' in cat) or ('IBL' in cat) or ('HBL' in cat) or ('5BZCat' in cat): 
                     cat +=' (blazar)'
                  if ('CRATES' in cat):
                     cat +=' (blazar candidate)'
                  #print("Here 3")
                  if gammaray:
                      if cat != 'X-Star' and cat != 'RQ-AGN' and cat != 'Xray-UNCL' and cat != 'Unknown' and cat != 'type-8' and cat != 'zw' and cat != 'type-9' and 'zw' not in obj_id:
                         out_file_candidates.write(f'{count:.0f},{ra:.5f},{dec:.5f},{cat},{obj_id}\n')
                         out_slopes.write(f'{ra:.5f},{dec:.5f},-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,False,{host_galaxy_image},False,False\n')
                  else:
                      out_file_candidates.write(f'{count:.0f},{ra:.5f},{dec:.5f},{cat},{obj_id}\n')
                      out_slopes.write(f'{ra:.5f},{dec:.5f},-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,False,{host_galaxy_image},False,False\n')
            find_out_sources  = np.load('findout_data.npy')
            out_file_candidates.seek(0)
            csv_candidates = csv.reader(out_file_candidates)
            header = next(csv_candidates)  # Skip the header row
            data = list(csv_candidates)
            np.save('candidates_d.npy', data)
            candidates  = np.load('candidates_d.npy')
            for rc in find_out_sources:
                r=rc.split()
                rafindout = float(r[0])
                decfindout = float(r[1])
                codefindout = float(r[2])
                namefindout = r[3]
                do_write_findout = True
                for row in candidates:
                    ra = float(row[1])
                    dec = float(row[2])
                    cat = row[3]
                    obj_id = row[4]
                    distf = sky_distance(ra, dec, rafindout, decfindout)
                    distf = distf * 3600.
                    if distf < 1:
                       if codefindout > 0 or cat =='X-extended':
                          c = int(codefindout/10000)
                          if cat =='X-Star':
                              #print("Here 2  host_galaxy_image ",host_galaxy_image)
                              if host_galaxy_image:
                                  cat = 'Galaxy'
                              codef= codefindout - c*10000+100000 
                              find_out_temp.write(f'{rafindout:.5f}    {decfindout:.5f}  {codef}  {namefindout}\n')  
                              do_write_findout = False
                          if cat =='radio-AGN' or cat =='Xray-UNCL' or 'Seyfert' in cat or cat =='QSO/Seyfert' or cat =='galaxy' or cat =='high-z-LBL?':
                              codef= codefindout - c*10000 + 40000
                              find_out_temp.write(f'{rafindout:.5f}    {decfindout:.5f}  {codef}  {namefindout}\n')  
                              do_write_findout = False
                          elif cat =='HBL': 
                              #print("namefindout cat distf codefindout ",namefindout,cat,distf,codefindout)
                              codef= codefindout-c*10000+10000
                              find_out_temp.write(f'{rafindout:.5f}    {decfindout:.5f}  {codef}  {namefindout}\n')  
                              #print(f'{rafindout:.5f}    {decfindout:.5f}  {codef}  {namefindout}\n')  
                              do_write_findout = False
                          elif cat =='IBL': 
                              codef= codefindout-c*10000+20000
                              find_out_temp.write(f'{rafindout:.5f}    {decfindout:.5f}  {codef}  {namefindout}\n')  
                              do_write_findout = False
                          elif cat =='LBL': 
                              codef= codefindout-c*10000+30000
                              find_out_temp.write(f'{rafindout:.5f}    {decfindout:.5f}  {codef}  {namefindout}\n')  
                              do_write_findout = False
                          elif cat =='cluster?' or cat =='X-extended': 
                              codef= -40000
                              find_out_temp.write(f'{rafindout:.5f}    {decfindout:.5f}  {codef}  {namefindout}\n')  
                              do_write_findout = False
                       else:
                          if cat =='Xray-UNCL' or cat =='radio-AGN':
                              codef= 45555
                              find_out_temp.write(f'{rafindout:.5f}    {decfindout:.5f}  {codef}  {namefindout}\n')  
                              do_write_findout = False
                       #if not do_write_findout:
                       #    print("The source code has been set to ",codef," for source ",namefindout,codefindout)
                if do_write_findout:
                    find_out_temp.write(f'{rafindout:.5f}    {decfindout:.5f}  {codefindout}  {namefindout}\n')  
                    do_write_findout = False
            #print("cand_int ",cand_int)
            if cand_int:
                  with open(input_file_optical, mode='r') as file:
                      csv_reader_opt = csv.reader(file)
                      header = next(csv_reader_opt)  # Skip the header row
                      for row in csv_reader_opt:
                         if input_file_optical == input_file_hstgsc:
                            obj_id = row[0]
                            ra = float(row[1])
                            dec = float(row[2])
                            magnitude = float(row[7])
                         elif input_file_optical == input_file_gaia:
                            obj_id = row[4]
                            ra = float(row[0])
                            dec = float(row[2])
                            magnitude = float(row[5])
                         elif input_file_optical == input_file_panstarrs:
                            obj_id = row[14]
                            ra = float(row[0])
                            dec = float(row[1])
                            magnitude = float(row[4])
                         distarcsec = 100
                         found = False
                         if magnitude < 8:
                            max_dist = 12
                         elif magnitude < 14:
                            max_dist = 8
                         else:
                            max_dist = 4
                         with open(input_file_candidates_int, mode='r') as c_int:
                            csv_reader_cand_int = csv.reader(c_int)
                            next(csv_reader_cand_int)
                            iox = 0
                            ir = 0
                            iro = 0
                            arx = 0
                            irx = 0
                            iw = 0
                            irg = 0
                            airxav = 0
                            aoxav = 0
                            airoav = 0
                            arxav = 0
                            aw1w2av = 0
                            argav = 0
                            for row_int in csv_reader_cand_int:
                               aox = float(row_int[0])
                               airx = float(row_int[1])
                               airo = float(row_int[2])
                               aro = float(row_int[3])
                               arx = float(row_int[4])
                               aw1w2 = float(row_int[5])
                               arg = float(row_int[6])
                               ra_int = float(row_int[7])
                               dec_int = float(row_int[8])
                               fx_1kev = float(row_int[9])
                               #print("ra dec ra_int dec_int ",ra,dec,ra_int,dec_int)
                               dist = sky_distance(ra, dec, ra_int, dec_int)
                               distarcsec = dist * 3600
                               type = 'OptCand' 
                               blue_bump = False 
                               host_galaxy_sed = False 
                               host_galaxy_image = False 
                               nupeak = -99
                               #print("distarcsec max_dist ",distarcsec,max_dist)
                               host_galaxy_image = check_host_galaxy_sdss(ra_int,dec_int)
                               if distarcsec < max_dist: 
                                  found = True
                                  if airx > -90:
                                      ir += 1
                                      airxav += airx
                                  if aox > -90:
                                      iox += 1
                                      aoxav += aox
                                  if airo > -90:
                                      iro += 1
                                      airoav += airo
                                  if arx > -90:
                                      irx += 1
                                      arxav += arx
                                  if aw1w2 > -90:
                                     iw += 1
                                     aw1w2av += aw1w2
                                  if arg > -90:
                                     irg += 1
                                     argav += arg
                            if ir > 0:
                               airxav = airxav/ir
                            else:
                               airxav = -99
                            if iro > 0:
                               airoav = airoav/iro
                            else:
                               airoav = -99
                            if iox > 0:
                               aoxav = aoxav/iox
                            else:
                               aoxav = -99
                            if irx > 0:
                               arxav = arxav/irx
                            else:
                               arxav = -99
                            if iw > 0:
                               aw1w2av = aw1w2av/iw
                            else:
                               aw1w2av = -99
                            if irg > 0:
                               argav = argav/irg
                            else:
                               aw1w2av = -99
                            if iox > 0:
                               code = 0
                               const = 0
                               x_slope = get_x_slope(sed_file,ra,dec)
                               #print("argav aw1w2av x_slope airoav",argav,aw1w2av,x_slope,airoav)
                               if argav < 0.95 and argav != 0 and (aw1w2av < -0.1 or abs(x_slope) < 0.9 or airoav < 0.9):
                                   blue_bump = True
                                   #print("Blue bump from r-g color",argav)
                               if aoxav > 1.2 and aoxav <= 2.0 and aw1w2av < 2.0:
                                   if arxav < 0.5:
                                      type = 'RQ-AGN' 
                                      const = 110000
                                   if airoav > 0 and airoav < 1.1:
                                       #print("Possible blue bump or strong host galaxy")
                                       if argav < 0.5 and (aw1w2av < .5 or abs(x_slope) < 0.9):
                                          blue_bump = True
                                          #print("2")
                                       if argav < 0.95 and (aw1w2av < -0.1 or abs(x_slope) < 0.9 or airoav < 0.9):
                                          blue_bump = True
                                          #print("Blue bump from r-g color",argav)
                                       elif argav > 2:
                                           host_galaxy_sed = True
                                           print("Host galaxy from r-g color",argav) 
                               elif (aoxav > 1.4 and aw1w2av > 2.0) or (aoxav > 2.0 and aoxav < 3.5):
                                   type = 'X-Star' 
                                   const = 100000
                                   if host_galaxy_sed or host_galaxy_image:
                                       type = 'Galaxy' 
                                       const = 110000
                               else:
                                   const = 70000
                               if magnitude > 19:
                                   code = const+1111
                               elif magnitude > 17:
                                   code = const+8888
                               else: 
                                   code = const+9000
                               if const == 110000:
                                   if fx_1kev > 3.e-12:
                                       code = const+8888
                                   elif fx_1kev > 3.e-13:
                                       code = const+4000
                                   elif fx_1kev > 3.e-14:
                                       code = const+2000
                                   else:
                                       code = const+1000
                               if type == 'type-9':
                                   type = 'NoXrayBlazarCandidate'
                               if flag == 'NI' and code > 100000:
                                   count += 1
                                   host_galaxy_image = check_host_galaxy_sdss(ra,dec)
                                   if not host_galaxy_image:
                                       host_galaxy_image = check_host_galaxy_panstarrs(ra,dec)
                                   do_write = check_if_incandidates(ra_int,dec_int)
                                   if do_write:
                                       do_write = False
                                       if extended:
                                          if type == 'HBL' or type == 'IBL' or type == 'LBL' or type == 'radio-AGN':
                                              type = 'cluster?'
                                          else:
                                               type = 'X-extended'
                                       #print("Extended 3 ",extended)
                                       if ('LBL' in type) or ('IBL' in type) or ('HBL' in type) or ('5BZCat' in type): 
                                          type +=' (blazar)'
                                       if ('CRATES' in type):
                                          type +=' (blazar candidate)'
                                       #print("Here 4")
                                       if gammaray:
                                           if type != 'X-Star' and type != 'RQ-AGN' and type != 'Xray-UNCL' and type != 'Unknown' and type != 'type-8' and type != 'zw' and type != 'type-9' and 'zw' not in obj_id:
                                               out_file_candidates.write(f'{count:.0f},{ra:.5f},{dec:.5f},{type},{obj_id}\n')
                                               out_file_find_out.write(f'{ra:.5f}    {dec:.5f}  {code:.0f}  {obj_id}\n')
                                               out_slopes.write(f'{ra:.5f},{dec:.5f},{aoxav:.2f},{airxav:.2f},{airoav:.2f},{aroav:.2f},{arxav:.2f},{aw1w2av:.2f},{argav:.2f},{nupeak:.1f},{x_slope:.2f},{slope_4lac_4fgl:.2f},{host_galaxy_sed},{host_galaxy_image},{blue_bump},{extended}\n')
                                       else:
                                           out_file_candidates.write(f'{count:.0f},{ra:.5f},{dec:.5f},{type},{obj_id}\n')
                                           out_file_find_out.write(f'{ra:.5f}    {dec:.5f}  {code:.0f}  {obj_id}\n')
                                           out_slopes.write(f'{ra:.5f},{dec:.5f},{aoxav:.2f},{airxav:.2f},{airoav:.2f},{aroav:.2f},{arxav:.2f},{aw1w2av:.2f},{argav:.2f},{nupeak:.1f},{x_slope:.2f},{slope_4lac_4fgl:.2f},{host_galaxy_sed},{host_galaxy_image},{blue_bump},{extended}\n')
                               elif flag != 'NI':
                                   count += 1
                                   host_galaxy_image = check_host_galaxy_sdss(ra,dec)
                                   if not host_galaxy_image:
                                       host_galaxy_image = check_host_galaxy_panstarrs(ra,dec)
                                   do_write = check_if_incandidates(ra_int,dec_int)
                                   if do_write and not gammaray:
                                       out_file_candidates.write(f'{count:.0f},{ra:.5f},{dec:.5f},HSTGSC,{obj_id}\n')
                                       out_file_find_out.write(f'{ra:.5f}    {dec:.5f}  {code:.0f}  OptCand\n')
                                       out_slopes.write(f'{ra:.5f},{dec:.5f},{aoxav:.2f},{airxav:.2f},{airoav:.2f},{aroav:.2f},{arxav:.2f},{argav:.2f},{aw1w2av:.2f},{nupeak:.1f},{x_slope:.2f},{slope_4lac_4fgl:.2f},{host_galaxy_sed},{host_galaxy_image},{blue_bump},{extended}\n')
else:
   count += 1
   with open(output_file_candidates, mode='w') as out_file_candidates, open(output_file_slopes, mode='w') as out_slopes:
      out_file_candidates.write('Number,R.A.,Dec.,Type,Name\n')
      out_slopes.write('ra,dec,aox,airx,airo,arx,aw1w2,arg,nu_peak,x_en_slope,gamma_sp_index,host_galaxy_sed,host_galaxy_image,blue_bump,Xray_extended\n')
      host_galaxy_sed = False
      blue_bump = False 
      extended = False 
      nupeak = -99
      x_slope = -99
      gamma_slope = -99
      with open(input_file_optical, mode='r') as file:
         csv_reader_opt = csv.reader(file)
         header = next(csv_reader_opt)  # Skip the header row
         with open(file_find_out, mode='w') as out_file_find_out:
            for row in csv_reader_opt:
               count += 1
               if input_file_optical == input_file_hstgsc:
                   obj_id = row[0]
                   ra = float(row[1])
                   dec = float(row[2])
                   magnitude = float(row[7])
               elif input_file_optical == input_file_gaia:
                   obj_id = row[4]
                   ra = float(row[0])
                   dec = float(row[2])
                   magnitude = float(row[5])
               elif input_file_optical == input_file_panstarrs:
                   obj_id = row[14]
                   ra = float(row[0])
                   dec = float(row[1])
                   magnitude = float(row[4])
               distarcsec = 100
               if cand_int:
                  with open(input_file_candidates_int, mode='r') as cand_int:
                     csv_reader_cand_int = csv.reader(cand_int)
                     next(csv_reader_cand_int)
                     ir = 0 
                     iaro = 0
                     iox = 0
                     iro = 0
                     iaro = 0
                     irx = 0
                     iw = 0
                     irg = 0
                     airxav = 0
                     aoxav = 0
                     aroav = 0
                     airoav = 0
                     arxav = 0
                     aw1w2av = 0
                     argav = 0
                     found = False
                     for row_int in csv_reader_cand_int:
                        aox = float(row_int[0])
                        airx = float(row_int[1])
                        airo = float(row_int[2])
                        aro = float(row_int[3]) 
                        arx = float(row_int[4])
                        aw1w2 = float(row_int[5])
                        arg = float(row_int[6])
                        ra_int = float(row_int[7])
                        dec_int = float(row_int[8])
                        fx_1kev = float(row_int[9])
                        dist = sky_distance(ra, dec, ra_int, dec_int)
                        distarcsec = dist * 3600
                        if distarcsec < 3: 
                             found = True
                             if airx > 0:
                                 ir += 1
                                 airxav += airx
                             if aox > 0:
                                 iox += 1
                                 aoxav += aox
                             if aro > 0:
                                 iaro += 1
                                 aroav += aro
                             if airo > 0:
                                 iro += 1
                                 airoav += airo
                             if arx > 0:
                                 irx += 1
                                 arxav += arx
                             if aw1w2 > 0:
                                 iw += 1
                                 aw1w2av += aw1w2
                             if arg > 0:
                                 irg += 1
                                 argav += arg
                     if ir > 0:
                         airxav = airxav/ir
                     else:
                         airxav = -99
                     if iro > 0:
                         airoav = airoav/iro
                     else:
                         airoav = -99
                     if iaro > 0:
                         aroav = aroav/iaro
                     else:
                         aroav = -99
                     if iox > 0:
                         aoxav = aoxav/iox
                     else:
                         aoxav = -99
                     if irx > 0:
                         arxav = arxav/irx
                     else:
                         arxav = -99
                     if iw > 0:
                         aw1w2av = aw1w2av/iw
                     else:
                         aw1w2av = -99
                     if irg > 0:
                         argav = argav/irg
                     else:
                         argav = -99
                     if os.path.exists(input_file_efeds):
                         with open(input_file_efeds, mode='r') as efeds:
                             csv_reader_efeds = csv.reader(efeds)
                             header = next(csv_reader_efeds)  # Skip the header row
                             for row_efeds in csv_reader_efeds:
                                 ra_efeds = float(row_efeds[1])
                                 dec_efeds = float(row_efeds[2])
                                 extension = float(row_efeds[31])
                                 dist = sky_distance(ra_int, dec_int, ra_efeds, dec_efeds)
                                 distarcsec = dist * 3600
                                 if distarcsec < 7 and extension > 20:
                                     extended = True
                     if os.path.exists(input_file_erass1):
                         with open(input_file_erass1, mode='r') as erass1:
                             csv_reader_erass1 = csv.reader(erass1)
                             header = next(csv_reader_erass1)  # Skip the header row
                             for row_erass1 in csv_reader_erass1:
                                 ra_erass1 = float(row_erass1[1])
                                 dec_erass1 = float(row_erass1[2])
                                 extension = float(row_erass1[29])
                                 dist = sky_distance(ra, dec, ra_erass1, dec_erass1)
                                 distarcsec = dist * 3600
                                 if distarcsec < 7 and extension > 20:
                                     extended = True
                         #print("eRASS1 extended ",extended)
                     if os.path.exists(input_file_4xmm):
                         with open(input_file_4xmm, mode='r') as xmm:
                             csv_reader_4xmm = csv.reader(xmm)
                             header = next(csv_reader_4xmm)  # Skip the header row
                             for row_4xmm in csv_reader_4xmm:
                                 ra_4xmm = float(row_4xmm[1])
                                 dec_4xmm = float(row_4xmm[2])
                                 extension = float(row_4xmm[16])
                                 dist = sky_distance(ra, dec, ra_4xmm, dec_4xmm)
                                 distarcsec = dist * 3600
                                 if distarcsec < 7 and extension > 20:
                                     extended = True
                     slope_4lac_4fgl = -99
                     if os.path.exists(input_file_4lac):
                         with open(input_file_4lac, mode='r') as lac:
                             csv_reader_4lac = csv.reader(lac)
                             header = next(csv_reader_4lac)  # Skip the header row
                             for row_4lac in csv_reader_4lac:
                                 ra_4lac = float(row_4lac[1])
                                 dec_4lac = float(row_4lac[2])
                                 dist = sky_distance(ra, dec, ra_4lac, dec_4lac)
                                 distarcsec = dist * 3600
                                 if distarcsec < 7:
                                     slope_4lac_4fgl = float(row_4lac[4])
                                     slope_4lac_4fgl_error = float(row_4lac[5])
               if cand_int:
                     #print("Here cand int ")
                     host_galaxy_image = check_host_galaxy_sdss(ra,dec)
                     x_slope = get_x_slope(sed_file,ra_int,dec_int)
                     if aoxav > 1.2 and aoxav <= 2.0 and aw1w2av < 2.0:
                         type = 'RQ-AGN' 
                         const = 110000
                         if airoav > 0 and airoav < 1.0:
                             print("Possible blue bump or strong host galaxy")
                     elif (aoxav > 1.4 and aw1w2av > 2.0) or (aoxav > 2.0 and aoxav < 3.5) or (airoav < 0.2 and airoav != 0. and aoxav > -90):
                         type = 'X-Star' 
                         #print("Here 4  host_galaxy_image ",host_galaxy_image)
                         const = 100000
                         if host_galaxy_image:
                            type = 'Galaxy' 
                            const = 110000
                     else:
                         const = 70000
                     if magnitude > 19:
                         code = const+1111
                     elif magnitude > 17:
                         code = const+8888
                     else: 
                         code = const+9000
                     if const == 110000:
                         if fx_1kev > 3.e-12:
                             code = const+8888
                         elif fx_1kev > 3.e-13:
                             code = const+4000
                         elif fx_1kev > 3.e-14:
                             code = const+2000
                         else:
                             code = const+1000
                     if type == 'type-9':
                         type = 'NoXrayBlazarCandidate'
                     if flag == 'NI' and code > 100000:
                         count += 1
                         host_galaxy_image = check_host_galaxy_sdss(ra,dec)
                         if not host_galaxy_image:
                             host_galaxy_image = check_host_galaxy_panstarrs(ra,dec)
                         if (type =='HBL' or type == 'IBL') and slope_4lac_4fgl > 2.2:
                             type = 'LBL' 
                         if extended:
                            if type == 'HBL' or type == 'IBL' or type == 'LBL' or type == 'radio-AGN':
                               type = 'cluster?'
                            else:
                               type = 'X-extended'
                         #print("Extended 4 ",extended)
                         if ('LBL' in type) or ('IBL' in type) or ('HBL' in type) or ('5BZCat' in type): 
                            type +=' (blazar)'
                         if ('CRATES' in type):
                             type +=' (blazar candidate)'
                         #print("Here 5")
                         if gammaray:
                             if type != 'X-Star' and type != 'RQ-AGN' and type != 'Xray-UNCL' and type != 'Unknown' and type != 'type-8' and type != 'zw' and type != 'type-9' and 'zw' not in obj_id:
                                out_file_candidates.write(f'{count:.0f},{ra:.5f},{dec:.5f},{type},{obj_id}\n')
                                out_file_find_out.write(f'{ra:.5f}    {dec:.5f}  {code:.0f}  {obj_id}\n')
                                out_slopes.write(f'{ra:.5f},{dec:.5f},{aoxav:.2f},{airxav:.2f},{airoav:.2f},{aroav:.2f},{arxav:.2f},{aw1w2av:.2f},{argav:.2f},{nupeak:.1f},{x_slope:.2f},{slope_4lac_4fgl:.2f},{host_galaxy_sed},{host_galaxy_image},{blue_bump},{extended}\n')
                         else:
                             out_file_candidates.write(f'{count:.0f},{ra:.5f},{dec:.5f},{type},{obj_id}\n')
                             out_file_find_out.write(f'{ra:.5f}    {dec:.5f}  {code:.0f}  {obj_id}\n')
                             out_slopes.write(f'{ra:.5f},{dec:.5f},{aoxav:.2f},{airxav:.2f},{airoav:.2f},{aroav:.2f},{arxav:.2f},{aw1w2av:.2f},{argav:.2f},{nupeak:.1f},{x_slope:.2f},{slope_4lac_4fgl:.2f},{host_galaxy_sed},{host_galaxy_image},{blue_bump},{extended}\n')
                     elif flag != 'NI':
                         count += 1
                         host_galaxy_image = check_host_galaxy_sdss(ra,dec)
                         if not host_galaxy_image:
                            host_galaxy_image = check_host_galaxy_panstarrs(ra,dec)
                         if not gammaray:
                            out_file_candidates.write(f'{count:.0f},{ra:.5f},{dec:.5f},HSTGSC,{obj_id}\n')
                            out_file_find_out.write(f'{ra:.5f}    {dec:.5f}  {code:.0f}  OptCand\n')
                            out_slopes.write(f'{ra:.5f},{dec:.5f},{aoxav:.2f},{airxav:.2f},{airoav:.2f},{aroav:.2f},{arxav:.2f},{aw1w2av:.2f},{argav:.2f},{nupeak:.1f},{x_slope:.2f},{slope_4lac_4fgl:.2f},{host_galaxy_sed},{host_galaxy_image},{blue_bump},{extended}\n')
               else:
                  #print("Here flag ",flag)
                  if magnitude > 19:
                     code = 71111
                  elif magnitude > 17:
                     code = 73333
                  elif magnitude > 14:
                     code = 78888
                  else: 
                     code = 79999
                  if flag != 'NI':
                     host_galaxy_image = check_host_galaxy_sdss(ra,dec)
                     if not host_galaxy_image:
                         host_galaxy_image = check_host_galaxy_panstarrs(ra,dec)
                     if not gammaray:
                         out_file_candidates.write(f'{count:.0f},{ra:.5f},{dec:.5f},HSTGSC,{obj_id}\n')
                         out_file_find_out.write(f'{ra:.5f}    {dec:.5f}  {code:.0f}  OptCand\n')
                         out_slopes.write(f'{ra:.5f},{dec:.5f},{aoxav:.2f},{airxav:.2f},{airoav:.2f},{aroav:.2f},{arxav:.2f},{aw1w2av:.2f},{argav:.2f},{nupeak:.1f},{x_slope:.2f},{slope_4lac_4fgl:.2f},{host_galaxy_sed},{host_galaxy_image},{blue_bump},{extended}\n')

