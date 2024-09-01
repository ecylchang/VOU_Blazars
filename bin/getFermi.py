import csv
import gzip
from astropy.coordinates import SkyCoord
import astropy.units as u
import sys

if len(sys.argv) < 3:
    print("Usage: python getFermi.py ra dec radius")
    sys.exit(1)

# Define the target sky position (RA and DEC) in degrees
target_ra = float(sys.argv[1])
target_dec = float(sys.argv[2])
max_distance = float(sys.argv[3])

# Define the input CSV file (possibly gzipped)
input_file = '/Users/paologiommi/app/VOU_BlazarsHybrid/data/Fermi_full_data.csv.gz'

# Define the output file to store matching rows
output_file = 'Fermi_matching_rows.csv'

# Open the output file in write mode
with open(output_file, 'w') as output_txt:
    writer = csv.writer(output_txt)

    # If the input file is gzipped, open it with gzip
    if input_file.endswith('.gz'):
        with gzip.open(input_file, 'rt') as input_csv:
            reader = csv.DictReader(input_csv)
            
            # Write the header row to the output
            header = reader.fieldnames
#            writer.writerow(header)
            
            for row in reader:
                try:
                    # Extract RA and DEC from the CSV rows (assuming column names 'RA' and 'DEC')
                    ra = float(row["Ra"])
                    dec = float(row["Dec"])

                    # Filter rows based on DEC difference
                    if abs(dec - target_dec) * 3600. < max_distance:
                        # Create SkyCoord objects for the target and current row
                        target_coord = SkyCoord(ra=target_ra * u.degree, dec=target_dec * u.degree)
                        row_coord = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)

                        # Calculate the angular separation between the target and current row
                        separation = target_coord.separation(row_coord)

                        # Check if the angular separation is less than the maximum allowed distance
                        if separation < max_distance * u.arcsec:
                            plus = float(row["flux"]) + float(row["err_flux"])
                            minus = float(row["flux"]) - float(row["err_flux"]) 
 
                            if row["flag"] == 'UL':
                               line = " {:.3E}   {:.3E}   {:.3E}   {:.3E}  {:.4f}  {:.4f}   UL   Fermi            Sahakyan et al. 2023 \n".format(float(row["freq."]),float(row["flux"]),plus,minus,float(row["MJD_start"]),float(row["MJD_end"]))
                            else:
                               line = " {:.3E}   {:.3E}   {:.3E}   {:.3E}  {:.4f}  {:.4f}   Det  Fermi            Sahakyan et al. 2023 \n".format(float(row["freq."]),float(row["flux"]),plus,minus,float(row["MJD_start"]),float(row["MJD_end"]))
                            output_txt.write(line)
                except ValueError:
                    # Handle cases where conversion to float fails
                    pass

