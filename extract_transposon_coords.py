#this script extracts the transposon coordinates from the file 'rmsk_TE.sorted.gtf' by matching the name supplied after the -t arg
#USAGE: python extract_transposon_coords <output bedGraph name> -t <Target Transposon name, eg L1>
from sys import argv
if '-t' in argv:
    target = argv[argv.index('-t')+1]
else:
    exit('Supply Transposon name with -t argument')

g = open(argv[1],'w')
with open('rmsk_TE.sorted.gtf','r') as f:
    for line in f.readlines():
        transposon = line.split()[9].split('\"')[1]
        chromosome = line.split()[0]
        start = line.split()[3]
        end = line.split()[4]

        if transposon==target:
            g.write('\t'.join([chromosome,start,end])+'\n')


