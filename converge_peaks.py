#converges the peaks of a bedgraph to make a larger range for the coordinates
#USAGE: python converge_peaks.py <input bedGraph> <output bedGraph> -p <index of value to be peaked> -s <index of value to be summed>
    #-p makes the converged ouput row contain the highest data point in the pth column of the input bedGraph (narrowPeak files for instance)
    #-s makes the converged output row contain the sum of the data points in the sth column of the input bedGraph (chip nexus data for instance)
from sys import argv

bed = argv[1]
out = argv[2]
if '-p' in argv:
    ps = '-p'
    value_col = int(argv[argv.index('-p')+1])
elif '-s' in argv:
    ps = '-s'
    value_col = int(argv[argv.index('-s'+1)])

cluster = []
last_chrm = 'chr1'
last_start = 0
last_stop = 0
outfile = open(out,'w')
for line in  open(bed,'r'):
    linesplit = line.strip().split('\t')
    if ps=='-p' or ps=='-s':
        chrm, start, stop, value = [linesplit[i] for i in [0,1,2,value_col]
        if chrm != last_chrm:
            last_chrm = chrm
            cluster = []
            last_stop = 0
            last_start = 0
        if cluster == []:
            last_start = int(start)
        if int(start) <= last_stop + 1:
            if cluster == []:
                cluster.append(last_value)
            cluster.append(float(value))
        else:
            if cluster != []:
                if ps =='-p':
                    final_value = str(max(cluster))
                elif ps == '-s':
                    final_value = str(sum(cluster))
                outfile.write(chrm+'\t'+str(last_start)+'\t'+str(last_stop)+'\t'+final_value+'\n')
                cluster = []
            else:
                if last_stop != 0:
                    outfile.write(line)

        last_stop = int(stop)
        last_value = float(value)
    
    else:
        chrm, start, stop = linesplit[0:3]
        value = 0 #now just used for the clustering of data points, value of value does not matter and is not recorded
        if chrm != last_chrm:
            last_chrm = chrm
            cluster = []
            last_stop = 0
            last_start = 0
        if cluster == []:
            last_start = int(start)
        if int(start) <= last_stop + 1:
            if cluster == []:
                cluster.append(last_value)
            cluster.append(float(value))
        else:
            if cluster != []:
                outfile.write(chrm+'\t'+str(last_start)+'\t'+str(last_stop)+'\n')
                cluster = []
            else:
                if last_stop != 0:
                    outfile.write(line)

        last_stop = int(stop)


