import csv
import sys
from astropy.coordinates import SkyCoord
import astropy.units as u
import os

def calculate_sky_distance(coord1, coord2):
    return coord1.separation(coord2).to(u.arcsecond).value

def remove_duplicates(data, distance_threshold=5.0):
    result_data = []
    
    for i in range(len(data)):
        ra_value = data[i]['R.A.']
        dec_value = data[i]['Dec.']
        name = data[i]['Name']
        type = data[i]['Type']

        if ra_value.strip() and dec_value.strip():  # Check if the values are non-empty
            try:
                ra_numeric = float(ra_value)
                dec_numeric = float(dec_value)
            except ValueError:
                print(f"Skipping row {i+2} - Non-numeric values in 'R.A.' or 'Dec.'")
                continue

            source1 = SkyCoord(ra=ra_numeric * u.deg, dec=dec_numeric * u.deg)
            duplicate = False
            for j in range(i + 1, len(data)):
                ra_value = data[j]['R.A.']
                dec_value = data[j]['Dec.']

                if ra_value.strip() and dec_value.strip():  # Check if the values are non-empty
                    try:
                        ra_numeric = float(ra_value)
                        dec_numeric = float(dec_value)
                    except ValueError:
                        print(f"Skipping row {j+2} - Non-numeric values in 'R.A.' or 'Dec.'")
                        continue

                    source2 = SkyCoord(ra=ra_numeric * u.deg, dec=dec_numeric * u.deg)
                    distance = calculate_sky_distance(source1, source2)
                    #print("distance distance_threshold ",distance,distance_threshold)
                    if distance < distance_threshold:
                        duplicate = True
                        break
            
            if not duplicate:
                result_data.append(data[i])
    
    return result_data

def renumber_column1(data):
    for i, row in enumerate(data):
        row['Number'] = i + 1
    return data

def main(input_file, output_file=None, distance_threshold=5.0):
    if not output_file:
        output_file = os.path.splitext(input_file)[0] + '_clean.csv'

    with open(input_file, mode='r') as f:
        reader = csv.DictReader(f)
        data = list(reader)
    data_with_HBL_or_LBL = [row for row in data if 'HBL' in row['Type'] or 'LBL' in row['Type'] or 'IBL' in row['Type'] or 'radio-AGN' in row['Type']]
    data_without_HBL_or_LBL = [row for row in data if 'HBL' not in row['Type'] and 'LBL' not in row['Type'] and 'IBL' not in row['Type'] and 'radio-AGN' not in row['Type']]
    rdata = data_without_HBL_or_LBL + data_with_HBL_or_LBL
    cleaned_data = remove_duplicates(rdata, distance_threshold)
    renumbered_data = renumber_column1(cleaned_data)

    with open(output_file, mode='w', newline='') as f:
        fieldnames = ['Number', 'R.A.', 'Dec.', 'Type', 'Name']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(renumbered_data)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        fov = float(sys.argv[2])
    else:
        print("Please provide an input file.")
        sys.exit(1)

    distance_threshold = min(5,max(fov*120/30,5)) #set the threshold to 1/30 of the FOV in units of arcsecs
    
    main(input_file, None, distance_threshold)

