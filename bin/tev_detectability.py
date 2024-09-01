import sys
import numpy as np
import math
from ebltable.tau_from_model import OptDepth

if len(sys.argv) > 3:
    nu_peak  = float(sys.argv[1])
    nufnu_peak  = float(sys.argv[2])
    z  = float(sys.argv[3])
else:
    print("Usage: python tev_detectability.py log(nu_peak) log(nufnu_peak) redshift ")
    print("e.g. > python tev_detectability.py 16.0 -11.18 .22") 
    exit()

ref_energy_value = 0.2  # Energy in TeV

# Initialize the EBL absorption model
tau = OptDepth.readmodel(model='dominguez')

# Compute the optical depth
ref_z = 0.25
ref_tau_value = tau.opt_depth(ref_z,ref_energy_value)
ref_p_gamma_gamma = np.exp(-ref_tau_value)

tau_value = tau.opt_depth(z,ref_energy_value)
# Convert optical depth to absorption probability
p_gamma_gamma = np.exp(-tau_value)

penalty_coefficient = min(1.0,p_gamma_gamma/ref_p_gamma_gamma)
#penalty_coefficient decreases the probability of detection of high z sources
delta_gamma_slope = 0.

delta_gamma_slope = 2.+0.11*max(15.0,nu_peak)-3.64

#print("delta_gamma_slope ",delta_gamma_slope)
pivot_energy = 0.005 # 5 GEV in inuts of TeV
boost_coefficient = (0.1/pivot_energy)**delta_gamma_slope
#boost_coefficient increases the probability of detection of flatter gamma-ray spectral slope sources (estimated from the
#correlation with nu_peak)
#print("penatly_coefficient : ",penalty_coefficient)
#print("boost_coefficient   : ",boost_coefficient)

# Print the absorption probability

def check_tev_detectability (test_value):

   # Define the data as a string
   data = """-10.75 100.
   -11.00 71.
   -11.25 33.
   -11.50 26.
   -11.75 7.
   -12.00 2.
   -12.25 0.
   -15.00 0.
   """
   test = float(test_value)
   # Initialize an empty list to store the data
   data_array = []

   # Split the data string into lines
   lines = data.split('\n')

   # Iterate through the lines, split each line into values, and store them in the array
   for line in lines:
       if line:
           values = line.split()
           if len(values) == 2:
               # Convert the values to float and int, respectively, and store as a tuple
               float_value = float(values[0])
               #float_value *= penalty_coefficient
               int_value = int(values[1].rstrip('.'))
               data_array.append((float_value, int_value))
   #print("test before",test)
   test += math.log10(penalty_coefficient*boost_coefficient)
   #print("test after",test)
   for entry in data_array:
       if test > entry[0]:
           tt = entry[1]
           #if entry[1] > 95.:
           if tt > 95.:
               print(f"Probability of detection with current IACTs           : > 95%")
           else:
               #print(f"Probability of detection with current IACTs           : {entry[1]:.0f}%")
               print(f"Probability of detection with current IACTs           : {tt:.0f}%")
           break

   test += 0.5
   for entry in data_array:
       tt = entry[1]
       if test > entry[0]:
           #if entry[1] > 95.:
           if tt > 95.:
               print(f"Probability of detection in CTAO extragalactic survey : > 95%")
           else:
               print(f"Probability of detection in CTAO extragalactic survey : {tt:.0f}%")
           break

   test += 0.5
   for entry in data_array:
       tt = entry[1]
       if test > entry[0]:
           #if entry[1] > 95.:
           if tt > 95.:
               print(f"Probability of detection in a long CTAO exposure      : > 95%")
           else:
               #print(f"Probability of detection in a deep CTAO exposure      : {entry[1]:.0f}%")
               print(f"Probability of detection in a deep CTAO exposure      : {tt:.0f}%")
           break


if (nu_peak >= 13.5):
   check_tev_detectability(nufnu_peak)
else:
   print("nu_peak is too low")



