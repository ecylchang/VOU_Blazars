import numpy as np
import pandas as pd
import argparse
import os
import csv
import math
from astropy.table import Table
from astropy.time import Time

def find_reference(input_file_path, target_string):
    found_line = None

    with open(input_file_path, 'r') as file:
        for line in file:
            if target_string in line:
                found_line = line.strip()  # Remove leading/trailing whitespace
                break  # Stop searching after the first occurrence
        words = found_line.split()
        if len(words) >= 2:
           result_string = ' '.join(words[1:])
           return result_string
        else:
           print("Input string has only one word.")

def specfinv3(catalog,reference,input_file,output_file,refs_file):

    if input_file != '--input_file':
        data = pd.read_csv(input_file, delimiter=",", 
               names=["name", "ra", "dec", "nu", "flux_d", "e_flux_d"], skiprows=1, header=None)
        data = data[data["flux_d"]/data["e_flux_d"] > 2]  # Adjust thresholds as needed
        data["Hz"] = data["nu"] * 1.e6
        data["nufnu"] = data["Hz"] * data["flux_d"] * 1e-26
        data["e_nufnu"] = data["Hz"] * data["e_flux_d"] * 1e-26
        data["upper"] = data["nufnu"]+data["e_nufnu"] 
        data["lower"] = data["nufnu"]-data["e_nufnu"] 
        data["mjd"] = 55000.0000
        data["det"] = 'Det' 
        data["cat"] = catalog 
        data["ref"] = reference 
        numpy_array = data.values
        selected_columns = numpy_array[:, [6, 7, 9, 10, 11, 11, 12, 13, 14]]
        # Define the format for each selected column
        formats = ['%.3E','%.3E','%.3E','%.3E','%.4f','%.4f','%s','%s','%s']
        custom_delimiters = ['   ', '   ', '   ', '  ', '  ', '   ', '  ', '         ']

        formatted_rows = []
        for row in selected_columns:
            formatted_row = ' '
            for i, value in enumerate(row):
                formatted_row += formats[i] % value
                # Add custom delimiter for all columns except the last one
                if i < len(row) - 1:
                    formatted_row += custom_delimiters[i]  # Add a space as a delimiter
            formatted_rows.append(formatted_row)

        with open(output_file, 'w') as f:
           if len(formatted_rows) > 0:
              f.writelines('\n'.join(formatted_rows))
              f.writelines('\n')
        print(catalog+" data added to SED ")

def sptsz(catalog,reference,input_file,output_file,refs_file):

    if input_file != '--input_file':
        data = pd.read_csv(input_file, delimiter=",",
              names=["RAJ2000","DEJ2000","S95-best","B_S95-best","b_S95-best","S150-best","B_S150-best","b_S150-best","S220-best","B_S220-best","b_S220-best"],skiprows=1, header=None)

        selected_columns = ["S95-best","B_S95-best","b_S95-best","S150-best","B_S150-best","b_S150-best","S220-best","B_S220-best","b_S220-best"]
              
        num_subarrays = 3
        subarrays = []
        # Split the selected columns into subarrays
        for i in range(num_subarrays):
            start_index = i * 3
            end_index = i * 3 + 3
            subarray = data.loc[:, selected_columns[start_index:end_index]].copy()
            subarray["mjd"] = 55000.0000
            subarray["det"] = 'Det'
            subarray["cat"] = catalog
            subarray["ref"] = reference

            if i == 0:
                subarray["Hz"] = 95.e9
            elif i == 1:
                subarray["Hz"] = 150.e9
            elif i == 2:
                subarray["Hz"] = 220.e9

            subarray = subarray.copy()
            subarray["nufnu"] = subarray["Hz"] * subarray[selected_columns[i * 3]] * 1.e-26
            subarray["upper"] =  subarray["Hz"] * (subarray[selected_columns[i * 3 + 1]])* 1.e-26 
            subarray["lower"] =  subarray["Hz"] * (subarray[selected_columns[i * 3 + 2]])* 1.e-26 
            #print('subarray ',subarray)

            subarrays.append(subarray)

        formatted_rows = []
        custom_delimiters = [['   ', '   ', '   ', '  ', '  ', '   ', '  ', '         ']] * num_subarrays
        for i in range(num_subarrays):
           nu_array = subarrays[i].values
           #selected_columns = nu_array[:, [6, 7, 9, 10, 2, 2, 3, 4, 5]]
           selected_columns = nu_array[:, [7, 8, 9, 10, 3, 3, 4, 5, 6]]
           formats = ['%.3E', '%.3E', '%.3E', '%.3E', '%.4f', '%.4f', '%s', '%s', '%s']
           custom_delimiter = custom_delimiters[i]  # Get the delimiter list for the current 'i'
           for row in selected_columns:
               formatted_row = ' '
               for j, value in enumerate(row):
                   # Check if the value is numeric before applying formatting
                   if isinstance(value, (int, float)):
                       formatted_row += formats[j] % value
                   else:
                       formatted_row += str(value)  # Convert non-numeric values to strings
                   # Add custom delimiter for all columns except the last one
                   if j < len(row) - 1:
                       formatted_row += custom_delimiter[j]
               if "NAN" not in formatted_row:
                  formatted_rows.append(formatted_row)

        with open(output_file, 'w') as f:
           if len(formatted_rows) > 0:
              f.writelines('\n'.join(formatted_rows))
              f.writelines('\n')

        print(catalog+" data added to SED ")

def pacopccs(catalog,reference,input_file,output_file,refs_file):

    if input_file != '--input_file':
        data = pd.read_csv(input_file, delimiter=",",
              names=["RAJ2000","DEJ2000","S5500","e_S5500","S9000","e_S9000","S18000","e_S18000","S24000","e_S24000","S33000","e_S33000",
                     "S39000","e_S39000"],skiprows=1, header=None)

        selected_columns = ["S5500","e_S5500","S9000","e_S9000","S18000","e_S18000","S24000","e_S24000","S33000","e_S33000","S39000","e_S39000"]
              
        num_subarrays = 6
        subarrays = []
        # Split the selected columns into subarrays
        for i in range(num_subarrays):
            start_index = i * 2
            end_index = i * 2 + 2
            subarray = data.loc[:, selected_columns[start_index:end_index]].copy()
            subarray["mjd"] = 55000.0000
            subarray["det"] = 'Det'
            subarray["cat"] = catalog
            subarray["ref"] = reference

            if i == 0:
                subarray["Hz"] = 5.52e9
            elif i == 1:
                subarray["Hz"] = 9.0e9
            elif i == 2:
                subarray["Hz"] = 18.e9
            elif i == 3:
                subarray["Hz"] = 24.e9
            elif i == 4:
                subarray["Hz"] = 33.e9
            elif i == 5:
                subarray["Hz"] = 39.e9

            subarray = subarray.copy()
            subarray["nufnu"] = subarray["Hz"] * subarray[selected_columns[i * 2]] * 1.e-26
            subarray["e_nufnu"] = subarray["Hz"] * subarray[selected_columns[i * 2 + 1]] * 1.e-26
            subarray["upper"] = subarray["nufnu"] + subarray["e_nufnu"]
            subarray["lower"] = subarray["nufnu"] - subarray["e_nufnu"]

            subarrays.append(subarray)

        formatted_rows = []
        custom_delimiters = [['   ', '   ', '   ', '  ', '  ', '   ', '  ', '         ']] * num_subarrays
        for i in range(num_subarrays):
           nu_array = subarrays[i].values
           selected_columns = nu_array[:, [6, 7, 9, 10, 2, 2, 3, 4, 5]]
           formats = ['%.3E', '%.3E', '%.3E', '%.3E', '%.4f', '%.4f', '%s', '%s', '%s']
           custom_delimiter = custom_delimiters[i]  # Get the delimiter list for the current 'i'
           for row in selected_columns:
               formatted_row = ' '
               for j, value in enumerate(row):
                   # Check if the value is numeric before applying formatting
                   if isinstance(value, (int, float)):
                       formatted_row += formats[j] % value
                   else:
                       formatted_row += str(value)  # Convert non-numeric values to strings
                   # Add custom delimiter for all columns except the last one
                   if j < len(row) - 1:
                       formatted_row += custom_delimiter[j]
               if "NAN" not in formatted_row:
                  formatted_rows.append(formatted_row)

        with open(output_file, 'w') as f:
           if len(formatted_rows) > 0:
              f.writelines('\n'.join(formatted_rows))
              f.writelines('\n')

        print(catalog+" data added to SED ")

def paco(catalog,reference,input_file,output_file,refs_file):

    if input_file != '--input_file':
        data = pd.read_csv(input_file, delimiter=",",
              names=["RAJ2000","DEJ2000","Obs.date","S4732","e_S4732","S5244","e_S5244","S5756","e_S5756",
             "S6268","e_S6268","S8232","e_S8232","S8744","e_S8744","S9256","e_S9256","S9768","e_S9768","S17232",
             "e_S17232","S17744","e_S17744","S18256","e_S18256","S18768","e_S18768","S23232","e_S23232","S23744","e_S23744",
             "S24256","e_S24256","S24768","e_S24768","S32232","e_S32232","S32744","e_S32744","S33256","e_S33256",
             "S33768","e_S33768","S38232","e_S38232","S38744","e_S38744","S39256","e_S39256","S39768","e_S39768"],skiprows=1, header=None)

        selected_columns = ["S4732","e_S4732","S5244","e_S5244","S5756","e_S5756",
             "S6268","e_S6268","S8232","e_S8232","S8744","e_S8744","S9256","e_S9256","S9768","e_S9768","S17232",
             "e_S17232","S17744","e_S17744","S18256","e_S18256","S18768","e_S18768","S23232","e_S23232","S23744","e_S23744",
             "S24256","e_S24256","S24768","e_S24768","S32232","e_S32232","S32744","e_S32744","S33256","e_S33256",
             "S33768","e_S33768","S38232","e_S38232","S38744","e_S38744","S39256","e_S39256","S39768","e_S39768"]
              
        num_subarrays = 24
        subarrays = []

        # Split the selected columns into subarrays
        for i in range(num_subarrays):
            start_index = i * 2
            end_index = i * 2 + 2
            subarray = data.loc[:, selected_columns[start_index:end_index]].copy()
            subarray["mjd"] = 55000.0000
            subarray["det"] = 'Det'
            subarray["cat"] = catalog
            subarray["ref"] = reference

            if i == 0:
                subarray["Hz"] = 47.32e9
            elif i == 1:
                subarray["Hz"] = 52.44e9
            elif i == 2:
                subarray["Hz"] = 57.56e9
            elif i == 3:
                subarray["Hz"] = 62.68e9
            elif i == 4:
                subarray["Hz"] = 82.32e9
            elif i == 5:
                subarray["Hz"] = 87.44e9
            elif i == 6:
                subarray["Hz"] = 92.56e9
            elif i == 7:
                subarray["Hz"] = 97.68e9
            elif i == 8:
                subarray["Hz"] = 172.32e9
            elif i == 9:
                subarray["Hz"] = 177.44e9
            elif i == 10:
                subarray["Hz"] = 182.56e9
            elif i == 11:
                subarray["Hz"] = 187.68e9
            elif i == 12:
                subarray["Hz"] = 232.32e9
            elif i == 13:
                subarray["Hz"] = 327.44e9
            elif i == 14:
                subarray["Hz"] = 242.56e9
            elif i == 15:
                subarray["Hz"] = 247.68e9
            elif i == 16:
                subarray["Hz"] = 322.32e9
            elif i == 17:
                subarray["Hz"] = 327.44e9
            elif i == 18:
                subarray["Hz"] = 332.56e9
            elif i == 19:
                subarray["Hz"] = 333.68e9
            elif i == 20:
                subarray["Hz"] = 382.32e9
            elif i == 21:
                subarray["Hz"] = 387.44e9
            elif i == 22:
                subarray["Hz"] = 392.56e9
            elif i == 23:
                subarray["Hz"] = 397.68e9


            subarray = subarray.copy()
            subarray["nufnu"] = subarray["Hz"] * subarray[selected_columns[i * 2]] * 1.e-26
            subarray["e_nufnu"] = subarray["Hz"] * subarray[selected_columns[i * 2 + 1]] * 1.e-26
            subarray["upper"] = subarray["nufnu"] + subarray["e_nufnu"]
            subarray["lower"] = subarray["nufnu"] - subarray["e_nufnu"]

            subarrays.append(subarray)

        formatted_rows = []
        custom_delimiters = [['   ', '   ', '   ', '  ', '  ', '   ', '  ', '         ']] * num_subarrays
        for i in range(num_subarrays):
           nu_array = subarrays[i].values
           selected_columns = nu_array[:, [6, 7, 9, 10, 2, 2, 3, 4, 5]]
           formats = ['%.3E', '%.3E', '%.3E', '%.3E', '%.4f', '%.4f', '%s', '%s', '%s']
           custom_delimiter = custom_delimiters[i]  # Get the delimiter list for the current 'i'
           for row in selected_columns:
               formatted_row = ' '
               for j, value in enumerate(row):
                   # Check if the value is numeric before applying formatting
                   if isinstance(value, (int, float)):
                       formatted_row += formats[j] % value
                   else:
                       formatted_row += str(value)  # Convert non-numeric values to strings
                   # Add custom delimiter for all columns except the last one
                   if j < len(row) - 1:
                       formatted_row += custom_delimiter[j]
               if "NAN" not in formatted_row:
                  formatted_rows.append(formatted_row)
        with open(output_file, 'w') as f:
           if len(formatted_rows) > 0:
              f.writelines('\n'.join(formatted_rows))
              f.writelines('\n')

        print(catalog+" data added to SED ")

def eprs(catalog,reference,input_file,output_file,refs_file):

    if input_file != '--input_file':
        data = pd.read_csv(input_file, delimiter=",",
              names=["ra","dec","S200","e_S200","S76","e_S76","S84","e_S84","S92","e_S92","S99","e_S99",
              "S107","e_S107","S115","e_S115","S122","e_S122","S130","e_S130","S143","e_S143","S151", "e_S151",
              "S158","e_S158","S166","e_S166","S174","e_S174","S181","e_S181","S189","e_S189","S197","e_S197",
              "S204","e_S204","S212","e_S212","S220","e_S220","S227","e_S227"],skiprows=1, header=None)

        selected_columns = ["S200","e_S200","S76","e_S76","S84","e_S84","S92","e_S92","S99","e_S99",
              "S107","e_S107","S115","e_S115","S122","e_S122","S130","e_S130","S143","e_S143","S151", "e_S151",
              "S158","e_S158","S166","e_S166","S174","e_S174","S181","e_S181","S189","e_S189","S197","e_S197",
              "S204","e_S204","S212","e_S212","S220","e_S220","S227","e_S227"]
              
        num_subarrays = 21
        subarrays = []

        # Split the selected columns into subarrays
        for i in range(num_subarrays):
            start_index = i * 2
            end_index = i * 2 + 2
            subarray = data.loc[:, selected_columns[start_index:end_index]].copy()
            subarray["mjd"] = 55000.0000
            subarray["det"] = 'Det'
            subarray["cat"] = catalog
            subarray["ref"] = reference

            if i == 0:
                subarray["Hz"] = 200.0e6
            elif i == 1:
                subarray["Hz"] = 76.0e6
            elif i == 2:
                subarray["Hz"] = 84.0e6
            elif i == 3:
                subarray["Hz"] = 92.e6
            elif i == 4:
                subarray["Hz"] = 99.e6
            elif i == 5:
                subarray["Hz"] = 107.e6
            elif i == 6:
                subarray["Hz"] = 115.e6
            elif i == 7:
                subarray["Hz"] = 122.e6
            elif i == 8:
                subarray["Hz"] = 130.e6
            elif i == 9:
                subarray["Hz"] = 143.e6
            elif i == 10:
                subarray["Hz"] = 151.e6
            elif i == 11:
                subarray["Hz"] = 158.e6
            elif i == 12:
                subarray["Hz"] = 166.e6
            elif i == 13:
                subarray["Hz"] = 174.e6
            elif i == 14:
                subarray["Hz"] = 181.e6
            elif i == 15:
                subarray["Hz"] = 189.e6
            elif i == 16:
                subarray["Hz"] = 197.e6
            elif i == 17:
                subarray["Hz"] = 204.e6
            elif i == 18:
                subarray["Hz"] = 212.e6
            elif i == 19:
                subarray["Hz"] = 220.e6
            elif i == 20:
                subarray["Hz"] = 227.e6


            subarray = subarray.copy()
            subarray["nufnu"] = subarray["Hz"] * subarray[selected_columns[i * 2]] * 1.e-23
            subarray["e_nufnu"] = subarray["Hz"] * subarray[selected_columns[i * 2 + 1]] * 1.e-23
            subarray["upper"] = subarray["nufnu"] + subarray["e_nufnu"]
            subarray["lower"] = subarray["nufnu"] - subarray["e_nufnu"]

            subarrays.append(subarray)

        formatted_rows = []
        custom_delimiters = [['   ', '   ', '   ', '  ', '  ', '   ', '  ', '         ']] * num_subarrays
        for i in range(num_subarrays):
           nu_array = subarrays[i].values
           selected_columns = nu_array[:, [6, 7, 9, 10, 2, 2, 3, 4, 5]]
           formats = ['%.3E', '%.3E', '%.3E', '%.3E', '%.4f', '%.4f', '%s', '%s', '%s']
           custom_delimiter = custom_delimiters[i]  # Get the delimiter list for the current 'i'
           for row in selected_columns:
               formatted_row = ' '
               for j, value in enumerate(row):
                   # Check if the value is numeric before applying formatting
                   if isinstance(value, (int, float)):
                       formatted_row += formats[j] % value
                   else:
                       formatted_row += str(value)  # Convert non-numeric values to strings
                   # Add custom delimiter for all columns except the last one
                   if j < len(row) - 1:
                       formatted_row += custom_delimiter[j]
               formatted_rows.append(formatted_row)

        with open(output_file, 'w') as f:
           if len(formatted_rows) > 0:
              f.writelines('\n'.join(formatted_rows))
              f.writelines('\n')

        print(catalog+" data added to SED ")

def gleamv2(catalog,reference,input_file,output_file,refs_file):

    if input_file != '--input_file':
        data = pd.read_csv(input_file, delimiter=",",
               names=["ra","dec","Fint084","e_Fint084","Fint092","e_Fint092","Fint099","e_Fint099","Fint107",
               "e_Fint107","Fint115","e_Fint115","Fint122","e_Fint122","Fint130","e_Fint130","Fint143","e_Fint143",
               "Fint151","e_Fint151","Fint158","e_Fint158","Fint166","e_Fint166","Fint174","e_Fint174","Fint181",
               "e_Fint181","Fint189","e_Fint189","Fint197","e_Fint197","Fint204","e_Fint204","Fint212","e_Fint212",
               "Fint220","e_Fint220","Fint227","e_Fint227"],skiprows=1, header=None)

        selected_columns = ["Fint084","e_Fint084","Fint092","e_Fint092","Fint099","e_Fint099","Fint107",
               "e_Fint107","Fint115","e_Fint115","Fint122","e_Fint122","Fint130","e_Fint130","Fint143","e_Fint143",
               "Fint151","e_Fint151","Fint158","e_Fint158","Fint166","e_Fint166","Fint174","e_Fint174","Fint181",
               "e_Fint181","Fint189","e_Fint189","Fint197","e_Fint197","Fint204","e_Fint204","Fint212","e_Fint212",
               "Fint220","e_Fint220","Fint227","e_Fint227"]
        num_subarrays = 19
        subarrays = []

        # Split the selected columns into subarrays
        for i in range(num_subarrays):
            start_index = i * 2
            end_index = i * 2 + 2
            subarray = data.loc[:, selected_columns[start_index:end_index]].copy()
            subarray["mjd"] = 55000.0000
            subarray["det"] = 'Det'
            subarray["cat"] = catalog
            subarray["ref"] = reference

            if i == 0:
                subarray["Hz"] = 84.0e6
            elif i == 1:
                subarray["Hz"] = 92.0e6
            elif i == 2:
                subarray["Hz"] = 99.0e6
            elif i == 3:
                subarray["Hz"] = 107.e6
            elif i == 4:
                subarray["Hz"] = 115.e6
            elif i == 5:
                subarray["Hz"] = 122.e6
            elif i == 6:
                subarray["Hz"] = 130.e6
            elif i == 7:
                subarray["Hz"] = 143.e6
            elif i == 8:
                subarray["Hz"] = 151.e6
            elif i == 9:
                subarray["Hz"] = 158.e6
            elif i == 10:
                subarray["Hz"] = 166.e6
            elif i == 11:
                subarray["Hz"] = 174.e6
            elif i == 12:
                subarray["Hz"] = 181.e6
            elif i == 13:
                subarray["Hz"] = 189.e6
            elif i == 14:
                subarray["Hz"] = 197.e6
            elif i == 15:
                subarray["Hz"] = 204.e6
            elif i == 16:
                subarray["Hz"] = 212.e6
            elif i == 17:
                subarray["Hz"] = 220.e6
            elif i == 18:
                subarray["Hz"] = 227.e6


            subarray = subarray.copy()
            subarray["nufnu"] = subarray["Hz"] * subarray[selected_columns[i * 2]] * 1.e-23
            subarray["e_nufnu"] = subarray["Hz"] * subarray[selected_columns[i * 2 + 1]] * 1.e-23
            subarray["upper"] = subarray["nufnu"] + subarray["e_nufnu"]
            subarray["lower"] = subarray["nufnu"] - subarray["e_nufnu"]

            subarrays.append(subarray)

        formatted_rows = []
        custom_delimiters = [['   ', '   ', '   ', '  ', '  ', '   ', '  ', '         ']] * num_subarrays
        for i in range(num_subarrays):
           nu_array = subarrays[i].values
           selected_columns = nu_array[:, [6, 7, 9, 10, 2, 2, 3, 4, 5]]
           formats = ['%.3E', '%.3E', '%.3E', '%.3E', '%.4f', '%.4f', '%s', '%s', '%s']
           custom_delimiter = custom_delimiters[i]  # Get the delimiter list for the current 'i'
           for row in selected_columns:
               formatted_row = ' '
               for j, value in enumerate(row):
                   # Check if the value is numeric before applying formatting
                   if isinstance(value, (int, float)):
                       formatted_row += formats[j] % value
                   else:
                       formatted_row += str(value)  # Convert non-numeric values to strings
                   # Add custom delimiter for all columns except the last one
                   if j < len(row) - 1:
                       formatted_row += custom_delimiter[j]
               test = formatted_row.split()
               #if float(test[1]) > 0.:
               if float(test[1]) > 0. and 'NAN' not in formatted_row:
                  formatted_rows.append(formatted_row)

        with open(output_file, 'w') as f:
           if len(formatted_rows) > 0:
              f.writelines('\n'.join(formatted_rows))
              f.writelines('\n')

        print(catalog+" data added to SED ")


def ratan600(catalog, reference, input_file, output_file, refs_file):

    if input_file != '--input_file':
        data = pd.read_csv(input_file, delimiter=",", 
             names=["ra", "dec", "F21.7", "e_F21.7", "F11.2", "e_F11.2", "F7.7", "e_F7.7", "F4.8", "e_F4.8", "F2.3", "e_F2.3"], 
             skiprows=1, header=None)

        selected_columns = ["F21.7", "e_F21.7", "F11.2", "e_F11.2", "F7.7", "e_F7.7", "F4.8", "e_F4.8", "F2.3", "e_F2.3"]
        num_subarrays = 5
        subarrays = []

        # Split the selected columns into subarrays
        for i in range(num_subarrays):
            start_index = i * 2
            end_index = i * 2 + 2
            subarray = data.loc[:, selected_columns[start_index:end_index]].copy()
            subarray["mjd"] = 55000.0000
            subarray["det"] = 'Det'
            subarray["cat"] = catalog
            subarray["ref"] = reference

            if i == 0:
                subarray["Hz"] = 21.7e9
            elif i == 1:
                subarray["Hz"] = 11.2e9
            elif i == 2:
                subarray["Hz"] = 7.7e9
            elif i == 3:
                subarray["Hz"] = 4.8e9
            elif i == 4:
                subarray["Hz"] = 2.3e9
            subarray = subarray.copy()
            subarray["nufnu"] = subarray["Hz"] * subarray[selected_columns[i * 2]] * 1.e-23
            subarray["e_nufnu"] = subarray["Hz"] * subarray[selected_columns[i * 2 + 1]] * 1.e-23
            subarray["upper"] = subarray["nufnu"] + subarray["e_nufnu"]
            subarray["lower"] = subarray["nufnu"] - subarray["e_nufnu"]

            subarrays.append(subarray)

        formatted_rows = []
        custom_delimiters = [['   ', '   ', '   ', '  ', '  ', '   ', '  ', '         ']] * num_subarrays
        for i in range(num_subarrays):
           nu_array = subarrays[i].values
           selected_columns = nu_array[:, [6, 7, 9, 10, 2, 2, 3, 4, 5]]
           formats = ['%.3E', '%.3E', '%.3E', '%.3E', '%.4f', '%.4f', '%s', '%s', '%s']
           custom_delimiter = custom_delimiters[i]  # Get the delimiter list for the current 'i'
           for row in selected_columns:
               formatted_row = ' '
               ok = True
               for j, value in enumerate(row):
                   # Check if the value is numeric before applying formatting
                   if isinstance(value, (int, float)):
                       formatted_row += formats[j] % value
                       if value == 'nan' or math.isnan(value):
                         ok = False
                   else:
                       formatted_row += str(value)  # Convert non-numeric values to strings
                   # Add custom delimiter for all columns except the last one
                   if j < len(row) - 1:
                       formatted_row += custom_delimiter[j]
               if ok:
                  formatted_rows.append(formatted_row)

        with open(output_file, 'w') as f:
           if len(formatted_rows) > 0:
              f.writelines('\n'.join(formatted_rows))
              f.writelines('\n')

        print(catalog + " data added to SED")

def mmmonitoring4(catalog,reference,input_file,output_file,refs_file):

    if input_file != '--input_file':
        data = pd.read_csv(input_file, delimiter=",",names=["ra", "dec", "nu", "flux_d", "e_flux_d", "date"], skiprows=1, header=None)
        data = data[data["flux_d"]/data["e_flux_d"] > 2]  # Adjust thresholds as needed
        data["Hz"] = data["nu"] * 1.e9
        data["nufnu"] = data["Hz"] * data["flux_d"] * 1e-23
        data["e_nufnu"] = data["Hz"] * data["e_flux_d"] * 1e-23
        data["upper"] = data["nufnu"]+data["e_nufnu"] 
        data["lower"] = data["nufnu"]-data["e_nufnu"] 
#        data["date"] = pd.to_datetime(float(data["date"]), format='%Y.%f')
#        data["mjd"] = Time(data["date"]).mjd
        data["mjd"] = 55000 
        data["det"] = 'Det' 
        data["cat"] = catalog 
        data["ref"] = reference 
        numpy_array = data.values
        
        selected_columns = numpy_array[:, [6, 7, 9, 10, 11, 11, 12, 13, 14]]
        # Define the format for each selected column
        formats = ['%.3E','%.3E','%.3E','%.3E','%.4f','%.4f','%s','%s','%s']
        custom_delimiters = ['   ', '   ', '   ', '  ', '  ', '   ', '  ', '         ']

        formatted_rows = []
        for row in selected_columns:
            formatted_row = ' '
            for i, value in enumerate(row):
                formatted_row += formats[i] % value
                # Add custom delimiter for all columns except the last one
                if i < len(row) - 1:
                    formatted_row += custom_delimiters[i]  # Add a space as a delimiter
            formatted_rows.append(formatted_row)

        with open(output_file, 'w') as f:
           if len(formatted_rows) > 0:
              f.writelines('\n'.join(formatted_rows))
              f.writelines('\n')
        print(catalog+" data added to SED ")



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Vizier SED data downloader')
    parser.add_argument('--input_file', type=str, help='File name', required=True)
    parser.add_argument('--output_file', type=str, help='Output file name', default='4Sed.csv')
    parser.add_argument('--refs_file', type=str, help='File including references', default='/Users/paologiommi/app/VOU_BlazarsHybrid/bin/catrefs.txt')
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    refs_file = args.refs_file
    filename = os.path.basename(input_file)
    if not os.path.exists(input_file):
       print("Error: ",input_file+" file not found.") 
       exit()

    values = filename.strip().split(".")
    catalog=values[0].upper()
    reference=find_reference(refs_file,catalog)
    if reference is None:
       print(f"Target string '{catalog}' not found in the input file.")
    if catalog == 'SPECFINDV3':
       specfinv3(catalog,reference,input_file,output_file,refs_file)
    elif catalog == 'EPRS':
       eprs(catalog,reference,input_file,output_file,refs_file)
    elif catalog == 'RATAN600':
       ratan600(catalog,reference,input_file,output_file,refs_file)
    elif catalog == 'MMMONITORING4':
       mmmonitoring4(catalog,reference,input_file,output_file,refs_file)
    elif catalog == 'GLEAMV2':
       gleamv2(catalog,reference,input_file,output_file,refs_file)
    elif catalog == 'PACO':
       paco(catalog,reference,input_file,output_file,refs_file)
    elif catalog == 'PACOPCCS':
       pacopccs(catalog,reference,input_file,output_file,refs_file)
    elif catalog == 'SPTSZ':
       sptsz(catalog,reference,input_file,output_file,refs_file)
    else:
       exit()
