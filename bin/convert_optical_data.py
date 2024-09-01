import csv

# Define input and output file names
input_filename = 'optical_data.csv'
output_filename = 'out-4Sed.txt'

# Function to format numbers to scientific notation with uppercase 'E'
def format_scientific(number, precision=3):
    return f"{float(number):.{precision}E}"

# Read the input file
with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
    # Skip the first line (header)
    next(infile)
    
    reader = csv.reader(infile)
    for row in reader:
        # Convert the necessary columns to floats
        frequency = float(row[1])
        flux = float(row[2])
        flux_error = float(row[3])
        
        # Calculate the new columns
        flux_plus_error = flux + flux_error
        flux_minus_error = flux - flux_error
        detect = 'Det'
        space = '   '
        
        # Format each column according to the specified format
        formatted_row = [
            format_scientific(row[1]),           # Frequency
            format_scientific(flux),             # Flux
            format_scientific(flux_plus_error),  # Flux + Flux Error
            format_scientific(flux_minus_error), # Flux - Flux Error
        ]
        #    f"{float(row[4]):.4f}",              # Start Time
        #    f"{float(row[5]):.4f}",              # End Time
        formatted_row2 = [
            detect,                              # Detection flag
            row[7],                              # Source
            space,                               #
            row[8]                               # Reference
        ]
        # Join the formatted columns with appropriate spacing and write to the output file
        print(" " + "   ".join(formatted_row) + "  " + f"{float(row[4]):.4f}" + "  " + f"{float(row[5]):.4f}" + "   " + "  ".join(formatted_row2))
        outfile.write(" " + "   ".join(formatted_row) + "  " + f"{float(row[4]):.4f}" + "  " + f"{float(row[5]):.4f}" + "   " + "  ".join(formatted_row2) + "\n")

#print(f"Formatted data has been written to {output_filename}")

