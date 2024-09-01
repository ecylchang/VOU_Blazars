import csv
from astropy.coordinates import SkyCoord
import astropy.units as u
import tempfile
import shutil

# Define the input and output file paths
input_file_path = 'tmp/output1.csv'
# Create a temporary file for the output
temp_output_file_path = tempfile.mktemp(suffix='.csv')

# Open the input CSV file and create the output CSV file
with open(input_file_path, 'r') as input_file, open(temp_output_file_path, 'w', newline='') as output_file:
    # Create CSV reader and writer objects
    csv_reader = csv.reader(input_file)
    csv_writer = csv.writer(output_file)

    # Read and write the first two lines of the header
    header_line1 = next(csv_reader)
    header_line2 = next(csv_reader)
    csv_writer.writerow(header_line1)
    csv_writer.writerow(header_line2)

    # Write the header to the output file
    #header = next(csv_reader)
    #csv_writer.writerow(header)

    # Initialize a flag to determine if the current row should be written to the output
    write_row = True

    # Initialize lists to store relevant rows
    LargeErrorCircles_rows = []
    SmallErrorCircles_rows = []

    # Iterate through the rows in the input file
    for row in csv_reader:
        source_type = row[1]
        ra = float(row[2])
        dec = float(row[3])

        # Check the source type and store in the corresponding list
        if source_type.lower().strip() in ['ipc', 'ipcsl', 'wgacat', 'rass']:
            LargeErrorCircles_rows.append((ra, dec, row))
        elif source_type.lower().strip() in ['2sxps', 'xmmsl2', 'chandracsc2' , 'xmmsl2', '4xmmdr13', 'erass1', 'bmw', 'swxcs']:
            SmallErrorCircles_rows.append((ra, dec, row))
            csv_writer.writerow(row)
        else:
            # Write other rows directly to the output file
            csv_writer.writerow(row)

    # Function to calculate the angular separation between two sky positions
    def calculate_angular_separation(coord1, coord2):
        return coord1.separation(coord2).arcminute

    # Iterate through LargeErrorCircles_rows and filter based on angular separation
    for LargeErrorCircles_row in LargeErrorCircles_rows:
        LargeErrorCircles_coord = SkyCoord(ra=LargeErrorCircles_row[0], dec=LargeErrorCircles_row[1], unit=(u.deg, u.deg))
        keep_row = True

        for SmallErrorCircles_row in SmallErrorCircles_rows:
            SmallErrorCircles_coord = SkyCoord(ra=SmallErrorCircles_row[0], dec=SmallErrorCircles_row[1], unit=(u.deg, u.deg))
            separation = calculate_angular_separation(LargeErrorCircles_coord, SmallErrorCircles_coord)
            if separation < 1.0:
                keep_row = False
                break

        # Write the row to the output file if it meets the criteria
        #if keep_row or (source_type.lower() not in ['ipcsl', 'wgacat', 'rass']):
        if keep_row:
            csv_writer.writerow(LargeErrorCircles_row[2])

# Rename the temporary output file to have the same name as the input file
shutil.move(temp_output_file_path, input_file_path)

