import csv
import os

def hms_to_degrees(ra_hms, dec_dms):
    # Convert right ascension from hours, minutes, seconds to degrees
    ra_parts = ra_hms.split()
    ra_deg = float(ra_parts[0]) * 15 + float(ra_parts[1])/ 4 + float(ra_parts[2])/ 240

    # Convert declination from degrees, arcminutes, arcseconds to degrees
    dec_parts = dec_dms.split()
    dec_deg = abs(float(dec_parts[0])) + float(dec_parts[1])/ 60 + float(dec_parts[2])/ 3600

    # Check for negative declination
    if dec_dms.startswith("-"):
        dec_deg = -dec_deg

    return ra_deg, dec_deg

def is_blank(s):
    return len(s.strip()) == 0

#gammaray = True
gammaray = False 
path='tmp/'
input_files = [
    path+'cvcat_out.csv',
    path+'xrbcat_out.csv',
    path+'snrgreen_out.csv',
    path+'mwmc_out.csv',
    path+'mwsc_out.csv',
    path+'mcxc.1.csv',
    path+'zw.1.csv',
    path+'mcxc.1.csv',
    path+'abell.1.csv',
    path+'whl.1.csv',
    path+'swxcs.1.csv',
    path+'psz2.1.csv',
    path+'mquas.1.csv',
    path+'pulsar.1.csv',
    path+'f2psr.1.csv'
]

output_file_find_out = path+'find_out_temp.txt'
output_file_candidates = path+'phase1_candidates.csv'

def count_lines(file_path):
    if not os.path.exists(file_path):
        print(f'The file {file_path} does not exist.')
        return 0
    with open(file_path, 'r') as file:
        line_count = sum(1 for line in file)
    return line_count

def process_file(input_file, code):
    with open(output_file_find_out, mode='a') as out_file_find_out:
        candidates = count_lines(output_file_candidates)
        with open(output_file_candidates, mode='a') as orig_cand:
            with open(input_file, mode='r') as cand:
                cat = input_file.split('/')
                ccat = cat[1].split('.')
                csv_reader_cand = csv.reader(cand)
                next(csv_reader_cand)  # Skip header
                for row in csv_reader_cand:
                    if (input_file == path+'mcxc.1.csv') or (input_file == path+'zw.1.csv'):
                        if (not is_blank(row[1])) and (not is_blank(row[2])):
                            ra, dec = hms_to_degrees(row[1], row[2])
                            name = row[0]
                    elif input_file == path+'mquas.1.csv':
                        ra = float(row[0])
                        dec = float(row[1])
                        name = row[2]
                    else:
                        name = row[0]
                        ra = float(row[1])
                        dec = float(row[2])
                    if candidates == 0:
                        orig_cand.write('Number,R.A.,Dec.,Type,Name\n')
                        candidates = 1 
                    if (not is_blank(row[1])) and (not is_blank(row[2])):
                        out_file_find_out.write(f'{ra:.5f}   {dec:.5f}  {code:.0f}  {ccat[0]}\n')
                        if input_file != path+'mquas.1.csv':
                            if gammaray:
                                if input_file == path+'pulsar.1.csv' or input_file == path+'f2psr.1.csv':
                                   orig_cand.write(f'{candidates},{ra:.5f},{dec:.5f},{ccat[0]},{name}\n')
                                   candidates +=1
                            else:
                                orig_cand.write(f'{candidates},{ra:.5f},{dec:.5f},{ccat[0]},{name}\n')
                                candidates +=1

def main():
    for input_file, code in zip(input_files, [120000, 130000, 140000, 150000, 160000, -40000, -40000, -40000, -40000, -40000, -40000, -40000, -7000, -8888, -8888]):
        if os.path.exists(input_file):
            process_file(input_file, code)

if __name__ == "__main__":
    if not os.path.exists(output_file_find_out):
        open(output_file_find_out, 'w').close()  # Create an empty file if it doesn't exist
    main()

