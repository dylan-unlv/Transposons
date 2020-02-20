#USAGE: python jsonToFasta <input json file> <output fasta file name>
from sys import argv

infile = argv[1]
outfile = argv[2]

json = open(infile,'r')
fdict = eval(json.readline())
out = open(outfile,'w')
for key, value in fdict.items():
    out.write('>'+key+'\n'+value+'\n')

