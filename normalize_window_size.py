#this script reads a bed file and changes the range of each datum to be n bp on each side of the middle
#Usage: python normalize_window_size.py <input_file> <optional: output_file> n; where n is the number of bp on each side of the middle
from sys import argv

if len(argv)<3:
    exit('Usage: normalizeWindow <input_file> <optional: output_file> n; where n is the number of bp on each side of the middle')

infile = argv[1]
try:
    n = int(argv[-1])
except (IndexError, ValueError) as e:
    if e==IndexError:
        exit('ERROR: Specify a range value n')
    else:
        exit('ERROR: Range value must be int')

with open(infile,'r') as f:
    lines = f.readlines()

for i in range(len(lines)):
    line = lines[i].split('\t') 
    mid = (int(line[1]) + int(line[2])) / 2
    line[1] = str(mid - n)
    line[2] = str(mid + n)
    if '\n' not in line:
        lines[i] = '\t'.join(line)
    else:
        lines[i] = '\t'.join(line)
if len(argv)==4:
    outfile = open(argv[2],'w')
    outfile.writelines(lines)
else:
    with open(infile,'w') as f:
        f.writelines(lines)


