import sys

def filter_lines(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            # Split the line into columns
            columns = line.split()

            # Check if the line has exactly 10 columns
            if len(columns) == 10 and (float(columns[9]) > 50 and float(columns[9]) < 65):
                #print(" "+line[:-4].strip()+" -"+columns[9])
                print(line.strip())

if __name__ == "__main__":
    # Check if the input file name is provided as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: python script.py input_file")
        sys.exit(1)

    # Get the input file name from the command line
    input_file = sys.argv[1]

    # Call the function with the input file name
    filter_lines(input_file)

