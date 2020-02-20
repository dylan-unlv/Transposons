#this script reads a narrowPeak file and changes the range of each datum to be n bp on each side of the peak
#Usage: python peak_shift.py <input_file> n; where n is the number of bp on each side of the peak
from sys import argv

infile = argv[1]
try:
    n = int(argv[2])
except (IndexError, ValueError) as e:
    if e==IndexError:
        exit('ERROR: Specify a range value n')
    else:
        exit('ERROR: Range value must be int')

with open(infile,'r') as f:
    lines = f.readlines()

for i in range(len(lines)):
    line = lines[i].split('\t') 
    peak = int(line[1]) + int(line[-1].strip())
    line[1] = str(peak - n)
    line[2] = str(peak + n)
    lines[i] = '\t'.join(line)
    
with open(infile,'w') as f:
    f.writelines(lines)


