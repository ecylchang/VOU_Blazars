import os 
import sys

if len(sys.argv) > 0:
    filename = sys.argv[1]

path='tmp/'

def remove_duplicates(file_path):
    unique_lines = set()
    
    if filename == path+"find_out_temp.txt":
       with open(filename, 'r') as file:
          lines = file.readlines()
       with open(filename, 'w') as file:
          for line in lines:
             columns = [col for col in line.split() if col]
             file.write(f"{float(columns[0]):.5f} {float(columns[1]):.5f} {columns[2]} {columns[3]}\n")
         
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
    with open(file_path, 'w') as file:
        for line in lines:
            if filename == path+"find_out_temp.txt":
               columns = line.split()
               #columns = [col for col in line.split() if col]
               #print(columns)
               key = (columns[0], columns[1], columns[2])
            else:
               columns = line.split(',')
               key = (columns[1], columns[2])
            
            if key not in unique_lines:
                file.write(line)
                unique_lines.add(key)

if os.path.exists(filename):
    remove_duplicates(filename)
i=-1
if filename == path+"phase1_candidates.csv":
  if os.path.exists(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    with open(filename, 'w') as file:
        for line in lines:
           i +=1
           if i == 0:
               file.write(line)
           else: 
               columns = line.split(',')
               file.write(f"{i},{float(columns[1]):.5f},{float(columns[2]):.5f},{columns[3]},{columns[4]}")
