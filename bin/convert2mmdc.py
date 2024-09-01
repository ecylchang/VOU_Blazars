import csv
import os
import sys

if len(sys.argv) > 1: 
    input_file = sys.argv[1]
else:       
    print("Please provide an input file.")
    sys.exit(1)     


def filter_and_write(input_file, output_file):
    with open(input_file, mode='r') as infile, open(output_file, mode='w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)

        # Write header to the output file
        header = next(reader)
        #writer.writerow(header[:3])
        outfile.write("frequency,flux,flux_err\n")

        # Iterate through rows and write selected columns
        for row in reader:
            if row[5] != 'UL':
                writer.writerow(row[:3])

if __name__ == "__main__":
    output_file = os.path.splitext(input_file)[0] + "_4mmdc.csv"
    filter_and_write(input_file, output_file)

