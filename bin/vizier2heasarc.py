import sys
import os
import re

def dms_to_deg(d, m, s):
    return (float(d) + float(m) / 60 + float(s) / 3600)


# Check if the file paths are provided as command-line arguments
if len(sys.argv) < 3:
    print("Please provide the input file path and output file path as command-line arguments.")
    sys.exit(1)

# Get the file paths from the command line
infile = sys.argv[1]
outfile = 'temp.txt'
catalog = sys.argv[2]

if not os.path.exists(infile):
   sys.exit(1)

# Open the input file for reading
with open(infile, 'r') as input_file:
    # Read the header line
    header = input_file.readline().strip()
    # Read the data lines
    data_lines = input_file.readlines()

# Process the data and write to the output file
with open(outfile, 'w') as output_file:
    if (catalog == 'crates'):
       # Write the header to the output file
       output_file.write('name,ctrpart_3p6_cm,flux_6_cm,ra,dec,flux_3p6_cm\n')
       for line in data_lines:
           data = line.strip().split(',')
           if len(data[2]) > 0:
              ra_parts = data[2].split()
              dec_parts = data[3].split()
              ra_deg = dms_to_deg(float(ra_parts[0]), float(ra_parts[1]), float(ra_parts[2]))
              dec_deg = dms_to_deg(abs(float(dec_parts[0])), float(dec_parts[1]), float(dec_parts[2]))
              if float(dec_parts[0]) < 0:
                 dec_deg =-dec_deg
              output_file.write('CRATES '+data[0]+',1,'+str(round(float(data[1]),1))+','+str(round(ra_deg*15.,5))+','+str(round(dec_deg,5))+','+str(round(float(data[4]),1))+'\n')
    elif (catalog == 'pulsar'):
       output_file.write('name,ra,dec\n')
       for line in data_lines:
           data = line.strip().split(',')
           ra_parts = data[1].split()
           dec_parts = data[2].split()
           ra_deg = dms_to_deg(float(ra_parts[0]), float(ra_parts[1]), float(ra_parts[2]))
           dec_deg = dms_to_deg(abs(float(dec_parts[0])), float(dec_parts[1]), float(dec_parts[2]))
           if float(dec_parts[0]) < 0:
              dec_deg =-dec_deg
           output_file.write('PSR '+data[0]+','+str(round(ra_deg*15.,5))+','+str(round(dec_deg,5))+'\n')
    elif (catalog == 'gb6'):
       output_file.write('name,ra,ra_error,dec,dec_error,flux_6_cm,flux_6_cm_error,border_flag,extent_flag,weak_flag\n')
       for line in data_lines:
           data = line.strip().split(',')
           ra_parts = data[1].split()
           dec_parts = data[3].split()
           ra_deg = dms_to_deg(float(ra_parts[0]), float(ra_parts[1]), float(ra_parts[2]))
           dec_deg = dms_to_deg(abs(float(dec_parts[0])), float(dec_parts[1]), float(dec_parts[2]))
           if float(dec_parts[0]) < 0:
              dec_deg =-dec_deg
           output_file.write('GB6 '+data[0]+','+str(round(ra_deg*15.,5))+','+data[2]+','+str(round(dec_deg,5))+','+data[4]+','+data[5]+','+data[6]+','+data[7]+','+data[8]+','+data[9]+'\n')
    elif (catalog == 'pcc'):
       output_file.write('name,ra,dec,flux,flux_error,det_flux,det_flux_error\n')
       for line in data_lines:
           data = line.strip().split(',')
           output_file.write(data[0]+','+data[1]+','+data[2]+','+data[5]+','+data[6]+','+data[3]+','+data[4]+'\n')
    elif (catalog == 'tgss150'):
       print('TGSS')
       output_file.write('name,ra,ra_error,dec,dec_error,flux_150_mhz,flux_150_mhz_error,rms_150_mhz\n')
       for line in data_lines:
           data = line.strip().split(',')
           output_file.write('TGSSADR '+data[0]+','+data[1]+','+data[2]+','+data[3]+','+data[4]+','+data[5]+','+data[6]+',3\n')
    elif (catalog == 'atpmn'):
       output_file.write('name,ra,dec,ra_accuracy,dec_accuracy,flux_3p5_cm,flux_3p5_cm_error,flux_6_cm,flux_6_cm_error,fit_flag\n')
       for line in data_lines:
          data = line.strip().split(',')
          output_file.write(data[0]+', '+data[1]+', '+data[2]+', .0001, .0001, '+data[5]+', '+data[6]+', '+data[3]+', '+data[4]+', '+data[7]+'\n')
    elif (catalog == 'nvss'):
       output_file.write('name,ra,ra_error,dec,dec_error,flux_20_cm,flux_20_cm_error,major_axis,major_axis_error,minor_axis,minor_axis_error,position_angle\n')
       for line in data_lines:
          data = line.strip().split(',')
          output_file.write(data[0]+','+data[1]+','+data[2]+','+data[3]+','+data[4]+','+data[5]+','+data[6]+','+data[7]+','+data[8]+','+data[9]+','+data[10]+','+data[11]+'\n')
    elif (catalog == 'panstarrs'):
       output_file.write('objID,RAJ2000,DEJ2000,e_RAJ2000,e_DEJ2000,gmag,e_gmag,rmag,e_rmag,imag,e_imag,zmag,e_zmag,ymag,e_ymag\n')
       for line in data_lines:
          data = line.strip().split(',')
          output_file.write(data[1]+','+data[2]+','+data[0]+','+data[3]+','+data[4]+','+data[5]+','+data[6]+','+data[7]+','+data[8]+','+data[9]+','+data[10]+','+data[11]+','+data[12]+','+data[13]+','+data[14]+'\n')
    elif (catalog == '3fhl'):
       output_file.write('name,ra,dec,error_radius,flux,flux_error,flux_10_20_gev,flux_10_20_gev_neg_err,flux_10_20_gev_pos_err,flux_20_50_gev,flux_20_50_gev_neg_err,flux_20_50_gev_pos_err,flux_50_150_gev,flux_50_150_gev_neg_err,flux_50_150_gev_pos_err,flux_150_500_gev,flux_150_500_gev_neg_err,flux_150_500_gev_pos_err,flux_0p5_2_tev,flux_0p5_2_tev_neg_err,flux_0p5_2_tev_pos_err,powerlaw_index,powerlaw_index_error,redshift\n')
       for line in data_lines:
          output_file.write('3FHL '+line)
#    elif (catalog == '4lac'):
#       output_file.write('name,ra,dec\n')
#       for line in data_lines:
#          output_file.write('4LAC-'+line)

input_file.close()
output_file.close()
os.rename(outfile,infile)
