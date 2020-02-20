#this script takes a bedfile and checks for a motif (optionally with n bp up and downstream), then creates a new bedfile with the last col being 0 or 1 depending on whether the motif was found
#USAGE: python check_motif_presence.py <input bedGraph> <output bedGraph> -w <optional window size> -hg <hg19 or hg38>
        #the -w argument adds a window to the start and end of each range in the bedGraph before checking for motif
        #the -hg argument specifies a genome version, either hg19 or hg38
import pysam
import subprocess
from sys import argv

bedfile = argv[1]
outfile = argv[2]
if '-w' in argv:
    try:
        window = int(argv[argv.index('-w')+1])
    except ValueError:
        exit('-w argument must be an int')
else:
    window = 0
if '-hg' in argv:
    if argv[argv.index('-hg')+1]=='hg38':
        hg = 'hg38.fasta'
    elif argv[argv.index('-hg')+1]=='hg19':
        hg = 'hg19.fasta'
    else:
        exit('-hg must be "hg38" or "hg19"')

bed = open(bedfile,'r')
out = open(outfile,'w')

for line in bed:
    fimo = open('fimo_out/fimo.tsv','w')
    fimo.close()
    chrm, start, stop = line.split()[:3]
    start = str(int(start)+window)
    stop = str(int(stop)+window)
    temp = open('temp.fa','w')
    temp.write(pysam.faidx('hg19.fa',chrm+':'+start+'-'+stop))
    subprocess.call(['fimo','GATA1/MA0140.2_truncated.meme','temp.fa'])
    fimo = open('fimo_out/fimo.tsv','r')
    if len(fimo.readlines())>4:
        out.write(line[:-1]+'\t1\n')
    else:
        out.write(line[:-1]+'\t0\n')

