#this script looks for 'family_id', maps the name and the coords of the span of TE to a frequency and a summed chip_nexus value
#USAGE: python extract_TE_frequency.py <input bedfile> <output bedfile> <output freq file> -n <col number of chip_nexus value
#input bedfile must have last column from  rmsk*.gtf  with 'family_id' giving transposon name, nth col a value from chip_nexus
 
from sys import argv
import re
#input file, output file from cmd line
bed = argv[1]
out = argv[2]
freq = argv[3]
if '-n' in argv:
    ncol = int(argv[argv.index('-n')+1])

#store data with key 'chrm:start-stop' value [TE name,sum(chip_nexus)]
out_data = {}
freq_data = {}

fn = open(bed,'r')

#changes made to accomodate mouse rmsk file
#gathering the data for data-fication
for line in fn:
    chrm, start, stop = line.split('\t')[:3]
    if '-n' in argv:
        nexus = float(line.strip().split('\t')[ncol])
    species, family, genus = line.strip().split('\t')[-3:]#value = line.split('\t')[-1]     #the descriptive line from .gtf file
    TE = genus+':'+species
    #TE = value[value.find('family_id "')+len('family_id "'):value.find('";',value.find('family_id "'),len(value))] #get the value following family_id
    chrm_id = chrm+':'+start+'-'+stop
    if chrm_id not in out_data:
        if '-n' in argv:
            out_data[chrm_id] = [TE,nexus]
        else:
            out_data[chrm_id] = TE
        if TE not in freq_data:
            freq_data[TE] = 1
        else:
            freq_data[TE] += 1
    else:
        if '-n' in argv:
            out_data[chrm_id][1] += nexus
        else:
            pass

outfile = open(out,'w')
freqfile = open(freq,'w')

#writing first output, bedGraph
for key in out_data:
    if '-n' in argv:
        outfile.write('\t'.join(re.split(r'[:-]\n',key) + [str(x) for x in out_data[key]])+'\n')
    else:
        outfile.write('\t'.join(re.split(r'[:-]\n',key) + [out_data[key]]) +'\n')#+ [str(out_data[key][1])])+'\n')

#sort by frequency
freq_list = []
for key in freq_data:
    freq_list.append((key, freq_data[key]))
freq_list.sort(key=lambda x: x[1],reverse=True)

#writing second output, frequency file
for i in freq_list:
    freqfile.write(i[0]+'\t'+str(i[1])+'\n')
