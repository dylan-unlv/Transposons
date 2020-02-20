#restricts the number of columns in a bedGraph, writes new bedGraph
#USAGE: python cut_cols <input bedGraph> <output bedGraph> -n <number of columns in outfile>
from sys import argv
if '-n' in argv:
    try:
        numColumns = int(argv[argv.index('-n')+1])
    except: ValueError:
        exit("Number of columns must be int")
else:
    exit("Specify number of columns")

out = open(argv[2],'w')
with open(argv[1],'r') as f:
    for line in f:
        lsplit = line.strip().split("\t")
        out.write("\t".join([lsplit[0]]+lsplit[3:])+'\n')
