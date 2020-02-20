#this script takes the files '*_TE_converged.txt' and produces a bedgraph based on a family of TEs
#USAGE: python subset_by_TE.py <idr file> <TE to search for> <output bedgraph file>
from re import split
from sys import argv
infile = argv[1]
TE = argv[2]
outfile = argv[3]
out = open(outfile,'w')
for line in open(infile,'r'):
    line = line.strip().split('\t')
    te = line[-3]
    if te==TE:
        out.write('\t'.join(line)+'\n')


