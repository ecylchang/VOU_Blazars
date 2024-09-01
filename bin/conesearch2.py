import argparse
import requests
import warnings
import os
import re
from io import StringIO
from bs4 import BeautifulSoup
from timeout_decorator import timeout
import warnings
import signal

class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException("Timeout occurred.")


def read_data_from_cats_ini(filecats, catalog):
    data = {}
    current_name = None

    with open(filecats, 'r') as file:
        for line in file:
            line = line.strip()
            # Check if the line contains the name in square brackets
            if line.startswith('[') and line.endswith(']'):
                current_name = line[1:-1]  # Extract the name without square brackets
            elif current_name == catalog:
                if line.startswith('url'):
                    parts = line.split('url:')
                    if len(parts) == 2:
                        data['url'] = parts[1].strip()
                elif line.startswith('columns'):
                    parts = line.split('mns:')
                    if len(parts) == 2:
                        data['columns'] = parts[1].strip()
    return data

def strip_tags(html):
    soup = BeautifulSoup(html, "lxml")
    tabledata = soup.find("tabledata")
    if not tabledata:
        print("No data found")
        exit(10) 

    table_rows = tabledata.find_all("tr")
    rows_data = []
    for row in table_rows:
        cells = row.find_all("td")
        row_data = [cell.get_text() for cell in cells]
        rows_data.append(",".join(row_data))
    return "\n".join(rows_data)

def query_VO(RA, DEC, SR, url, verbosity,print_field_names):
    
    params = {
        "RA": RA,
        "DEC": DEC,
        "SR": SR,
        "VERB": verbosity
    }
    
    response = requests.get(url, params=params)
    
    if response.status_code == 200:
        html_content = response.text
        if "System overloaded" in response.text:
            print("Remote System overloaded")
            exit (2)
        #print(html_content)
        if 'name="Error"' in response.text:
            print("Remote system error")
            exit (3)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            soup = BeautifulSoup(html_content, 'html.parser')  
        fields = soup.find_all('field')

        names_array = [field['name'] for field in fields]
        if print_field_names == 'y':
           for index, name in enumerate(names_array):
              print(f"Index: {index}, Field name: {name}")
    else:
        print("Error retrieving catalog information:", response.status_code)
        return None

    # Find the index of the starting and ending tags and extract the relevant portion
    start_tag_index = html_content.find("<TABLEDATA>")
    end_tag_index = html_content.find("</TABLEDATA>")

    if start_tag_index == -1 or end_tag_index == -1:
        print("No data found for this position.")
        exit(10)

    tabledata_content = html_content[start_tag_index:end_tag_index + len("</TABLEDATA>")]
    if "<TD>" in tabledata_content.upper():
       return names_array,tabledata_content
    else:
       print("No data found for this position.")
       exit(10)
#    pass

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Query VizieR catalog information.")
    parser.add_argument("--ra", type=str, help="RA (Right Ascension) value in decimal degrees.")
    parser.add_argument("--dec", type=str, help="DEC (Declination) value in decimal degrees.")
    parser.add_argument("--radius", type=str, help="Search radius (SR) ")
    parser.add_argument("--runit", type=str, help="units for radius (degrees, arcmin, arcsec d/f=arcmin)",default='arcmin')
    parser.add_argument("--print_field_names", type=str, help="Print all field names (d/f=no)", default='no')
    parser.add_argument("--db", type=str, help="cats1/2.ini file (d/f=none)", default='none')
    parser.add_argument("--catalog", type=str, help="catalog name  (d/f=none)", default='none')
    parser.add_argument("--output", type=str, help="output file  (d/f=none)", default='none')
    parser.add_argument("--timeout_time", type=int, help="timeout in seconds (d/f=15)", default='15')
    parser.add_argument("--verbosity", type=str, help="verbosity level of the output (d/f=3)", default='3')

    args = parser.parse_args()
    verbosity = args.verbosity
    timeout_seconds = args.timeout_time
    #print ("timeout_time ",timeout_seconds)
    if args.runit == 'arcmin':
       radius = float(args.radius)/60.
    elif args.runit == 'degrees':
       radius = float(args.radius)
    elif args.runit == 'arcsec':
       radius = float(args.radius)/3600.
    else:
       print ('Unknown radius units. Allowed values are degrees, arcmin or arcsec ')
       exit(1)

    if (args.db != 'none') & (args.catalog != 'none'):
       catalog = args.catalog
       result = read_data_from_cats_ini(args.db,catalog)
       if result:
          url = result['url']
          columns = result['columns']
       else:
          print ('Catalog not found')
          exit(1)
    try: 
       signal.signal(signal.SIGALRM, timeout_handler)
       signal.alarm(timeout_seconds)

       names_array,tabledata_content = query_VO(args.ra, args.dec, radius, url, verbosity, args.print_field_names)
    except TimeoutException:
       print("Timeout occurred.")
       exit(20)
    except requests.RequestException as e:
       print("Error retrieving catalog information:", e)
       exit(1)
    except Exception as e:
       print("An error occurred:", e)
       exit(1)

    input_names = columns.split(',')
    matching_indices = [names_array.index(name) if name in names_array else None for name in input_names]

    if tabledata_content:
       fields_separated_by_commas = strip_tags(tabledata_content)
       lines = fields_separated_by_commas.splitlines()
       columns_values = [line.split(',') for line in lines]

    if args.print_field_names == 'y':
        for j in range(0,len(columns_values)):
           if len(names_array) != len(columns_values[j]):
              print("Error: The arrays must have the same length.")
           else:
              for index, value in enumerate(columns_values[j]):
                 name = names_array[index]
                 print(f"Source nr. {j}, Field nr. {index} {name}: {value}")

    with open(args.output, "w") as f_out:
       f_out.write(columns+'\n')
       for j in range(0,len(columns_values)):
          res = ''
          for i, name in enumerate(input_names):
             if matching_indices[i] is not None:
                if res: 
                   res += ','
                res += columns_values[j][matching_indices[i]]
             else:
                print(f"Field name '{name}' not found in the names_array")
          f_out.write(res+'\n')
