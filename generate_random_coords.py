#this script generates random genome coordinates, takes a command line argument of output file, number of coordinates, and range each line should cover
#Usage: python generate_random_coords.py <output file> n_data_points range
from sys import argv
from random import choice
outfile = argv[1]
try:
    n = int(argv[2])
    rang = int(argv[3])
except e:
    if e == IndexError:
        exit("ERROR: Need number of data to generate and a range for each datum")
    elif e == ValueError:
        exit("ERROR: Number of data and range must be int")
    else:
        exit(e)

chromosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
chrange = range(10000,30000000)
out = open(outfile,'w')

for i in range(n):
    chrm = choice(chromosomes)
    start = choice(chrange)
    stop = start + rang
    out.write('\t'.join([chrm, str(start), str(stop),'a','b'])+'\n')


