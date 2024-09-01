import requests
import argparse
import re


parser = argparse.ArgumentParser(description='NH downloader')
parser.add_argument('--ra', type=float, help='R.A.', required=True)
parser.add_argument('--dec', type=float, help='Declination', required=True)

args = parser.parse_args()
ra = args.ra
dec = args.dec
entry =str(ra)+','+str(dec)

url = "https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl"

params = {
    'Entry': entry,  # Replace with your object name or coordinates
    'NR': 'GRB/SIMBAD+Sesame/NED',
    'CoordSys': 'Equatorial',
    'equinox': '2000',
    'radius': '0.1',
    'usemap': '0'
}

# Make a GET request with the specified parameters
response = requests.get(url, params=params)
#print("Response ",response)
# Check the response status code
if response.status_code == 200:
    pattern = r"Weighted average nH \(cm\*\*-2\)\s+(\d+\.\d+E[+-]\d+)"
    match = re.search(pattern, response.text)
    
    if match:
        nh_value = match.group(1)
        print(f"{nh_value}")
#        print(f"Weighted average nH (cm**-2): {nh_value}")
    else:
        print("nH value not found in the response.")
else:
    print(f"Request failed with status code: {response.status_code}")

