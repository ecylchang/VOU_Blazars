import csv
import math
import sys
import os

# Function to calculate the sky distance between two points
def sky_distance(ra1, dec1, ra2, dec2):
    # Convert coordinates to radians
    ra1 = math.radians(ra1)
    dec1 = math.radians(dec1)
    ra2 = math.radians(ra2)
    dec2 = math.radians(dec2)

    # Haversine formula for sky distance
    d_ra = ra2 - ra1
    d_dec = dec2 - dec1
    a = math.sin(d_dec/2)**2 + math.cos(dec1) * math.cos(dec2) * math.sin(d_ra/2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a)) * 57.32
    #57.32 converts radians to deg
    return c

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


path ='tmp/'
input_file_5bzcat = path+'5bzcat.1.csv'
input_file_3hsp = path+'3hsp.1.csv'
input_file_4lac = path+'4lacdr3.1.csv'
input_file_mquas = path+'mquas.1.csv'
input_file_sao = path+'sao.c.csv'
input_file_hd = path+'hd.c.csv'
input_file_sdss = path+'sdss_out.txt'
input_file_6df = path+'6df.1.csv'
input_file_panstarrs = path+' panstarrs.i.csv'

matched = False
input_file_candidates  = path+'phase1_candidates_clean.csv'
output_file_candidates = path+'phase1_candidates_final.csv'

if not os.path.exists(input_file_candidates):
    print("Input file ",input_file_candidates," does not exist") 
    exit

with open(output_file_candidates, mode='w') as out_cand,open(input_file_candidates, mode='r') as cand:
    out_cand.write("Number,R.A.,Dec.,Type,Name,CatalogName,redshift\n")
    candidates = csv.reader(cand)
    header = next(candidates)  # Skip the header row
    for rrow in candidates:
        written = False
        sdss = False 
        count = float(rrow[0])
        #print ("Count ",count)
        ra_cand = float(rrow[1])
        dec_cand = float(rrow[2])
        type_cand = rrow[3]
        name_cand = rrow[4]
        z_cat = -99.0
        #print("ra_cand dec_cand z_cat ",ra_cand,dec_cand,z_cat)
        if os.path.exists(input_file_mquas):
           with open(input_file_mquas, mode='r') as cat:
              catalog_sources = csv.reader(cat)
              header = next(catalog_sources)  # Skip the header row
              for rrow in catalog_sources:
                  ra_cat=float(rrow[0])
                  dec_cat=float(rrow[1])
                  name_cat= rrow[2] 
                  dist = sky_distance(ra_cand, dec_cand, ra_cat, dec_cat)
                  dist = dist* 3600 
                  if dist < 5.0:
                      if rrow[4] != '':
                         z_cat= float(rrow[4])
                      break  
        #print("mquas z_cat ra_cand, dec_cand ",z_cat,ra_cand,dec_cand)
        if os.path.exists(input_file_sdss):
           with open(input_file_sdss, mode='r') as cat:
              catalog_sdss = csv.reader(cat)
              header = next(catalog_sdss)  # Skip the header row
              for rrow in catalog_sdss:
                  #z_cat = -99 
                  ra_cat=float(rrow[1])
                  dec_cat=float(rrow[2])
                  zsp= rrow[14] 
                  zphot= rrow[15] 
                  if zsp !='':
                     z=float(zsp)
                  elif zphot !='':
                     z=float(zphot)
                  else:
                     z = -99.
                  dist = sky_distance(ra_cand, dec_cand, ra_cat, dec_cat)
                  dist = dist* 3600 
                  #print("dist z ",dist,z)
                  if dist < 4.0 and z > 0.:
                      #print("zsp ",zsp," zphot ",zphot)
                      z_cat = z 
                      sdss = True
                      break  
        #print("sdss z_cat ",z_cat)
        if os.path.exists(input_file_6df):
           with open(input_file_6df, mode='r') as cat:
              catalog_6df = csv.reader(cat)
              header = next(catalog_6df)  # Skip the header row
              for rrow in catalog_6df:
                  #z_cat = -99 
                  ra_cat,dec_cat = hms_to_degrees(rrow[1],rrow[2])
                  z= rrow[3] 
                  z_quality= rrow[4] 
                  dist = sky_distance(ra_cand, dec_cand, ra_cat, dec_cat)
                  dist = dist* 3600 
                  if dist < 5.0 and z_quality == 4:
                     z_cat = z 
        #print("6df z_cat ",z_cat)
        if os.path.exists(input_file_3hsp):
           with open(input_file_3hsp, mode='r') as cat:
              #z_cat = -99 
              catalog_sources = csv.reader(cat)
              header = next(catalog_sources)  # Skip the header row
              for rrow in catalog_sources:
                  ra_cat=float(rrow[1])
                  dec_cat=float(rrow[2])
                  name_cat= rrow[0] 
                  dist = sky_distance(ra_cand, dec_cand, ra_cat, dec_cat)
                  dist = dist* 3600 
                  #print("dist 3hsp ",dist)
                  if dist < 5.0:
                      if rrow[3] == '':
                         z_cat = -99
                      else:
                         z_cat= float(rrow[3])
                      out_cand.write(f'{count:.0f},{ra_cand:.5f},{dec_cand:.5f},{type_cand},{name_cand},{name_cat},{z_cat:.3f}\n')
                      written = True
        #print("3HSP z_cat ",z_cat)
        if os.path.exists(input_file_5bzcat):
           if not written:
              with open(input_file_5bzcat, mode='r') as cat:
                 catalog_sources = csv.reader(cat)
                 header = next(catalog_sources)  # Skip the header row
                 for rrow in catalog_sources:
                     ra_cat,dec_cat = hms_to_degrees(rrow[1],rrow[2])
                     name_cat= rrow[0] 
                     dist = sky_distance(ra_cand, dec_cand, ra_cat, dec_cat)
                     dist = dist* 3600 
                     if dist < 5.0:
                         z_cat= float(rrow[3])
                         if z_cat > -1:
                            out_cand.write(f'{count:.0f},{ra_cand:.5f},{dec_cand:.5f},{type_cand},{name_cand},{name_cat},{z_cat:.3f}\n')
                            written = True
        #print("bzcat z_cat ",z_cat)
        if os.path.exists(input_file_4lac):
           if not written:
              with open(input_file_4lac, mode='r') as cat:
                 catalog_sources = csv.reader(cat)
                 header = next(catalog_sources)  # Skip the header row
                 for row in catalog_sources:
                     name_cat= row[0] 
                     ra_cat= float(row[1])
                     dec_cat= float(row[2])
                     if row[3] != '':
                        z = float(row[3])
                     else: 
                        z = -99
                     dist = sky_distance(ra_cand, dec_cand, ra_cat, dec_cat)
                     dist = dist* 3600 
                     if dist < 5.0 and z > -1:
                         z_cat = z
                         out_cand.write(f'{count:.0f},{ra_cand:.5f},{dec_cand:.5f},{type_cand},{name_cand},{name_cat},{z_cat:.3f}\n')
                         written = True
        #print("4lac z_cat ",z_cat)
        if os.path.exists(input_file_sao):
           #z_cat = -99 
           with open(input_file_sao, mode='r') as cat:
              catalog_sources = csv.reader(cat)
              header = next(catalog_sources)  # Skip the header row
              for rrow in catalog_sources:
                  ra_cat,dec_cat = hms_to_degrees(rrow[1],rrow[2])
                  name_cat= 'SAO'+rrow[0] 
                  dist = sky_distance(ra_cand, dec_cand, ra_cat, dec_cat)
                  dist = dist* 3600 
                  if dist < 8.0:
                      out_cand.write(f'{count:.0f},{ra_cand:.5f},{dec_cand:.5f},{type_cand},{name_cand},{name_cat},-99\n')
                      written = True
        #print("sao z_cat ",z_cat)
        if os.path.exists(input_file_hd):
           if not written:
              with open(input_file_hd, mode='r') as cat:
                 catalog_sources = csv.reader(cat)
                 header = next(catalog_sources)  # Skip the header row
                 for rrow in catalog_sources:
                     ra_cat,dec_cat = hms_to_degrees(rrow[1],rrow[2])
                     name_cat= 'HD'+rrow[0] 
                     dist = sky_distance(ra_cand, dec_cand, ra_cat, dec_cat)
                     dist = dist* 3600 
                     if dist < 8.0:
                         out_cand.write(f'{count:.0f},{ra_cand:.5f},{dec_cand:.5f},{type_cand},{name_cand},{name_cat},-99\n')
                         written = True
        #print("hd z_cat ",z_cat)

        if not written:
            name_cat ='---'
            if 'mcxc' or 'zw' in name_cand:
               name_cat = name_cand
            if z_cat >= 0:
                out_cand.write(f'{count:.0f},{ra_cand:.5f},{dec_cand:.5f},{type_cand},{name_cand},{name_cat},{z_cat:.3f}\n')
            else:
                out_cand.write(f'{count:.0f},{ra_cand:.5f},{dec_cand:.5f},{type_cand},{name_cand},{name_cat},-99\n')
