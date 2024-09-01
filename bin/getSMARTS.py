import csv
import gzip
from astropy.coordinates import SkyCoord
import astropy.units as u
import sys

filepath = '/Users/paologiommi/app/VOU_BlazarsHybrid/data/'
if len(sys.argv) < 4:
    print("Usage: python getswiftrows.py ra(deg) dec(deg) radius(arcsec)")
    sys.exit(1)

target_ra = float(sys.argv[1])
target_dec = float(sys.argv[2])
max_distance = float(sys.argv[3])

input_file = filepath+'smarts.csv.gz'

#kev = 2.418e17
#constants = [0.5*kev, 1.0*kev, 1.5*kev , 3.0*kev, 4.5*kev]

output_file = 'matching_rows.csv'

zero = 0.0

# Pre-compute target coordinates
target_coord = SkyCoord(ra=target_ra * u.degree, dec=target_dec * u.degree)

# Store lines to write in a list
lines_to_write = []

# If the input file is gzipped, open it with gzip
if input_file.endswith('.gz'):
    with gzip.open(input_file, 'rt') as input_csv:
        reader = csv.DictReader(input_csv)
        
        header = reader.fieldnames
        for row in reader:
            try:
                ra = float(row['RA'])
                dec = float(row['DEC'])

                if abs(dec - target_dec) * 3600. < max_distance:
                    #print(row)
                    row_coord = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
                    separation = target_coord.separation(row_coord)
                    plus = float(row['FLUX']) + float(row['FLUX_ERROR'])
                    minus = float(row['FLUX']) - float(row['FLUX_ERROR'])

                    if separation < max_distance * u.arcsec:
                        columns = [float(row[column]) for column in header]
                        columns_to_write_values = [float(row[column]) for column in header[2:]]
                        line = " {:.3E}   {:.3E}   {:.3E}   {:.3E}  {:.4f}  {:.4f}   Det  SMARTS           Bonning et al. 2012, ApJ 756, 13 \n".format(columns[2],columns[4], plus , minus, columns[-2], columns[-1])
                        lines_to_write.append(line)
            except ValueError:
                pass

# Write lines to file
n = len(lines_to_write)
if n > 0:
  print(n,"SMARTS data points found")
else:
  print("No SMARTS data found")
with open(output_file, 'w') as output_txt:
    output_txt.writelines(lines_to_write)

