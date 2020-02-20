#this script takes a transposon family name and an rmsk file and outputs a bedgraph of the transposon family specified
#USAGE: python gather_TE_family_coords.py TE <rmsk file> <output bed>
import re
from sys import argv

TE = argv[1]
rmsk = argv[2]
outfile = argv[3]

out = open(outfile,'w')
for line in open(rmsk,'r'):
    indices = [0,3,4,-1] 
    chrm, start, stop, value = [line.strip().split('\t')[x] for x in indices]
    te = value[value.find('family_id "')+len('family_id "'):value.find('";',value.find('family_id "'),len(value))] #get the value following family_id
    if te==TE:
        out.write('\t'.join([chrm, start, stop])+'\n')

