#this script takes an input data matrix and outputs a bedgraph with an optional value, whose column is denoted -n
#USAGE: python data_matrix_2_bedgraph.py <input bedGraph> <output bedGraph> -n <index of desired value>
    #-n: include the nth column from the data matrix as a value at the end of the output bedGraph 
from sys import argv
dm = open(argv[1], 'r')
header = dm.readline()
if '-n' in argv:
    value_index = int(argv[argv.index('-n')+1])
else:
    value_index = 0

outfile = open(argv[2],'w')

for line in dm:
    name = line.split()[0]
    chrm = name.split(':')[0]
    start = name.split(':')[1].split('-')[0]
    end = name.split('-')[1]
    if '-n' in argv:
        value = line.split()[motif_index].strip()
        outfile.write('\t'.join([chrm, start, end, value])+'\n')
    else:
        outfile.write('\t'.join([chrm, start, end])+'\n')
        
