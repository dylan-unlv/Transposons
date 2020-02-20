#this script takes an input bedGraph and a genome, then produces the fastas for the ranges within, not exclusive to TE
#USAGE bedGraphToFasta <input_bedGraph> <reference_genome> <output_file.fasta>
from sys import argv
import pysam

if len(argv)<4:
    exit('USAGE: bedGraphToFasta <input_bedGraph> <reference_genome> <output_file.fasta>')

bedfile = argv[1]
genome = argv[2]
outfile = argv[3]


bed = open(bedfile,'r')
out = open(outfile,'w')

n_invalid = 0
n_valid = 0
valid_chromosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8',
'chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
'chr20','chr21','chr22','chrX','chrY'] #list of chromosomes that are valid, in order to exclude chrUn and chr1_KI270892v1_alt coords, which mess up stuff

for line in bed:
    chrm, start, stop = line.split()[:3]
    t_coords = chrm+':'+start+'-'+stop
    if chrm in valid_chromosomes:
        out.write(pysam.faidx(genome, t_coords))
        n_valid+=1
    else:
        n_invalid+=1

print(str(n_valid) + ' sequences in fasta file')
print(str(n_invalid) + ' invalid chromosome names (e.g. chrUn)')
