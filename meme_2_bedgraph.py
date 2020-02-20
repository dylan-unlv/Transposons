#this script opens fimo_out/fimo.tsv and outputs a bedGraph containing the range of the motif
#USAGE: python meme_2_bedgraph.py outfile
from re import split
f = open('fimo_out/fimo.tsv','r')
g = open(argv[1],'w')

for line in f:
    info = split(r'[:-]',line.split('\t')[2])
    g.write('\t'.join(info)+'\n')   
