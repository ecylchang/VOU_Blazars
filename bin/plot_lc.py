import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from astropy.time import Time, TimeDelta


if len(sys.argv) < 4:
    print("Usage: python plot_lc.py input_file energy band delta-t(days) for rebinning\n")
    print("e.g.: python plot_lc.py Sed.txt 1kev 10 \n")
    print("Allowed energy bands: 3.4micron 4.6micron imag rmag vmag gmag 1kev 100mev\n")
    sys.exit(1)

input_file = sys.argv[1]
en_band = sys.argv[2]
delta_mjd = float(sys.argv[3])
if delta_mjd <= 0:
   delta_mjd = 0.1

if en_band == '1kev':
   requested_frequency = 2.418e17
   label =r'1KeV $\nu$F$_{\nu}$ flux'
   colors = '#ff00ff'
elif en_band == '5kev':
   requested_frequency =1.088E+18
   label =r'4.5KeV $\nu$F$_{\nu}$ flux'
   colors = '#ff00ff'
elif en_band == '100mev':
   requested_frequency = 2.418e24 
   label =r'100 MeV $\nu$F$_{\nu}$ flux'
   colors = '#00cccc'
elif en_band == '4.6micron':
   requested_frequency = 6.517E+13 
   label =r'4.6$\mu$  $\nu$F$_{\nu}$ flux'
   colors = '#ff0066'
elif en_band == '3.4micron':
   requested_frequency = 8.817E+13 
   label =r'3.4 $\mu$ $\nu$F$_{\nu}$ flux'
   colors = '#ff0066'
elif en_band == 'gmag':
   requested_frequency = 6.237E+14 
   label =r' 6.2E+14 Hz (gmag), $\nu$F$_{\nu}$ flux'
   colors = 'orangered'
elif en_band == 'vmag':
   requested_frequency = 5.45E+14 
   label =r' 5.4E+14 Hz (vmag), $\nu$F$_{\nu}$ flux'
   colors = '#cc3300'
elif en_band == 'rmag':
   requested_frequency = 4.85E+14 
   label =r' 4.8E+14 Hz (rmag), $\nu$F$_{\nu}$ flux'
   colors = '#559922'
elif en_band == 'imag':
   requested_frequency = 3.99E+14 
   label =r' 3.9E+14 Hz (imag), $\nu$F$_{\nu}$ flux'
   colors = '#3333ff'
else:
   print('Energy band not supported')
   print("Allowed energy bands: 3.4micron 4.6micron 1kev\n")
   exit()

start_processing = False
# Initialize empty lists to store data
x_values = []
y_values = []
y_errors = []

# Read the data from the file
with open(input_file, 'r') as file:
    for line in file:
        # Check if the line starts with "--------"
        if line.startswith("--------"):
            start_processing = True
            continue  # Skip this line and start processing the next
        if start_processing:
           # Split each line into columns
           columns = line.strip().split()
        
           # Check if the first column is equal to e.g. 2.418E+17
           freq_ok = 1000.
           if float(columns[0]) > 0:
              freq_ok = abs(1.-requested_frequency/float(columns[0]))
           if (freq_ok < 0.01) & (float(columns[4]) != 55000.0000) & (columns[6] == 'Det'):
               error=(float(columns[2])-float(columns[3]))/2
               if error > 0:
                  if float(columns[1])/error > 1.2:
                     x_values.append(float(columns[4]))  # Column 5
                     y_values.append(float(columns[1]))  # Column 2
                     y_errors.append(error)  # Error calculation
# Convert lists to NumPy arrays for plotting
x_values = np.array(x_values)
y_values = np.array(y_values)
y_errors = np.array(y_errors)
if len(x_values) == 0:
  print('Plot_lc.py: There are no data points in this energy band')
  exit(1)

sorting_indices = np.argsort(x_values)
x_values = x_values[sorting_indices]
y_values = y_values[sorting_indices]
y_errors = y_errors[sorting_indices]

# Initialize lists to store rebinned data
rebinned_x_values = []
rebinned_y_values = []
rebinned_y_errors = []
min_max_arr = []

# Calculate the weighted mean for each bin
current_bin_start = min(x_values)
current_bin_end = current_bin_start + delta_mjd
window_sum = 0.0
weight_sum = 0.0
#sum = 0.0
for x_values, y_values, y_errors in zip(x_values, y_values, y_errors):
    if x_values >= current_bin_start and x_values <= current_bin_end:
        window_sum += y_values / (y_errors ** 2)
        weight_sum += 1 / (y_errors ** 2)
    else:
        if weight_sum > 0:
            rebinned_x_values.append(current_bin_start + delta_mjd / 2)  # Use the midpoint of the bin
            rebinned_y_values.append(window_sum / weight_sum)
            rebinned_y_errors.append(1 / np.sqrt(weight_sum))
        # Move to the next bin
        current_bin_start += delta_mjd
        current_bin_end = current_bin_start + delta_mjd
        while x_values <= current_bin_start:
           current_bin_start += delta_mjd
           current_bin_end = current_bin_start + delta_mjd
        while x_values > current_bin_end:
           current_bin_start += delta_mjd
           current_bin_end = current_bin_start + delta_mjd
        window_sum = 0.0
        weight_sum = 0.0
        window_sum += y_values / (y_errors ** 2)
        weight_sum += 1 / (y_errors ** 2)
if weight_sum > 0:
    rebinned_x_values.append(current_bin_start + delta_mjd / 2)  # Use the midpoint of the bin
    rebinned_y_values.append(window_sum / weight_sum)
    rebinned_y_errors.append(1 / np.sqrt(weight_sum))
 
# Convert rebinned lists to NumPy arrays for plotting
rebinned_x_values = np.array(rebinned_x_values)
rebinned_y_values = np.array(rebinned_y_values)
rebinned_y_errors = np.array(rebinned_y_errors)
mjd_min = 54466 #2008.00
mjd_max = 60676 #2025.00
min_max_arr.append(mjd_min)
min_max_arr.append(mjd_max)
mM = np.array(min_max_arr)
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

# Assuming the necessary data arrays are already defined:
# rebinned_x_values, rebinned_y_values, rebinned_y_errors, colors, label, mjd_min, mjd_max, mM

if len(rebinned_x_values) > 0:
# Convert MJD to year.fraction
   mjd_to_year_fraction = Time(rebinned_x_values, format='mjd').decimalyear
   mjd_to_year_mM = Time(mM, format='mjd').decimalyear

# Create the plot for rebinned data
   fig, ax1 = plt.subplots(figsize=(12, 4))
   ax1.set_xlim(mjd_min, mjd_max)

# Plot the data on the primary (bottom) x-axis
   ax1.errorbar(rebinned_x_values, rebinned_y_values, yerr=rebinned_y_errors, fmt='o', markersize=2, color=colors, label='')
   ax1.set_xlabel('Time [MJD, days]', fontweight='bold')
   ax1.set_ylabel(label, fontweight='bold')

# Create a secondary (top) x-axis for year.fraction
   ax2 = ax1.twiny()
   ax2.set_xlabel('Time [years]', fontweight='bold')

# Synchronize the secondary x-axis with the primary x-axis
   def mjd_to_year(mjd):
       return Time(mjd, format='mjd').decimalyear

   def year_to_mjd(year):
       return Time(year, format='decimalyear').mjd

   ax2.set_xlim(mjd_to_year(mjd_min), mjd_to_year(mjd_max))

   year_min = mjd_to_year(mjd_min)
   year_max = mjd_to_year(mjd_max)
   deltax = year_max - year_min
   start_year = int(np.floor(year_min))
   end_year = int(np.ceil(year_max))

   if deltax > 0:
       num_intervals = 10
       num_ticks = num_intervals + 1
       #step = (year_max-year_min)/num_intervals
       step = (year_max-year_min)/num_ticks
       ax2.set_xlim(year_min, year_max)
       tick_positions = np.arange(start_year, end_year +1, step)
       tick_positions_mjd = [year_to_mjd(year) for year in tick_positions]
       ax2.set_xticks(tick_positions)
       if (end_year - start_year) >= num_intervals:
           ax2.set_xticklabels([f'{int(year)}' for year in tick_positions])
       else:
           ax2.set_xticklabels([f'{year:.1f}' for year in tick_positions])
   else:
       ax2.set_xlim(year_min - 0.1, year_max + 0.1)

   ax2.set_xticklabels([f'{year:.1f}' for year in tick_positions])
# Set y-axis scale if needed
   y_limits = ax1.get_ylim()
   if abs(y_limits[1] / y_limits[0]) > 10:
       plt.yscale('log')
   
   plt.tight_layout()
   plt.savefig('lc_plot.png', dpi=300)  # Save the plot to a file
   
