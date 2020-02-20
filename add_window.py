#this script adds a window to a bedfile on each line both before and after the given range
#USAGE: python add_window.py <input file> <output file>
from sys import argv

window = 100 #size of window to add / subtract from start / end

f = open(argv[1],'r')
g = open(argv[2],'w')     #fn.split('.')[0]+'_expanded.bedGraph','w')

for line in f:
    seq, start, end= line.split('\t')
    start = str(int(start) - 100)
    end = str(int(end) + 100)
    g.write('\t'.join([seq,start,end])+'\n')
