'''
This script is a pipeline to gather all 5 variables with the input being a bedfile with coordinates to harvest data from
If you change variables/cell lines, you'll have to edit this script with the new file paths and columns, lines where this can be done
are preceeded by '#!'

Recommended to run sortBed on input file beforehand for the converge algorithm

command line args include:
    -c | use if input file is not yet converged
    -m | use if input file denotes motif locations
    -ch| use if input file is characterization matrix
    -tf| use if input file is transcription factor chip seq data (like GATA1)
    -dm| use if input file is dna methylation 
    -hm| use if input file is histone methylation (H3K4me3)
    -ha| use if input file is histone acetylation (H3K27ac)
    -na| use if input file is bedfile of places to harvest data from
    -cl| use to denote which column of the input file contains relevant data (starts at 0)
Usage:
    python gather_data_pipeline.py <input_bed_file> -cl <column number of value in input> [OPTIONS]
'''
from sys import argv
import subprocess
import pysam
outfile = open('15k_peaks_data_matrix.txt','w')
def converge_bed(bedfile, col): #converges bedfile
    if type(col) == int: 
        converged = []
        cluster = []
        last_chrm = 'chr1'
        last_start = 0
        last_stop = 0

        for line in  open(bedfile,'r'):
            chrm, start, stop, value = line.strip().split()[:3] + [line.strip().split()[col]]
            if chrm != last_chrm:
                last_chrm = chrm
                cluster = []
                last_stop = 0
                last_start = 0
            if cluster == []:
                last_start = int(start)
            if int(start) <= last_stop + 5:
                if cluster == []:
                    cluster.append(last_value)
                cluster.append(float(value))
            else:
                if cluster != []:
                    converged.append(chrm+'\t'+str(last_start)+'\t'+str(last_stop)+'\t'+str(max(cluster))+'\n')
                    cluster = []
                else:
                    if last_stop != 0:
                        if int(stop)-int(start)>50:
                            converged.append(line)
                        else:
                            converged.append(chrm+'\t'+str(int(start)-25)+'\t'+str(int(stop)+25)+'\t'+value+'\n')
            last_stop = int(stop)
            last_value = float(value)

    else: 
        converged = []
        last_chrm = 'chr1'
        last_start = 0
        last_stop = 0

        for line in  open(bedfile,'r'):
            chrm, start, stop= line.strip().split()[:3]
            value = 0 #used for keeping track, 0 not recorded
            if chrm != last_chrm:
                last_chrm = chrm
                cluster = []
                last_stop = 0
                last_start = 0
            if cluster == []:
                last_start = int(start)
            if int(start) <= last_stop + 5:
                if cluster == []:
                    cluster.append(last_value)
                cluster.append(value) 
            else:
                if cluster != []:
                    converged.append(chrm+'\t'+str(last_start)+'\t'+str(last_stop)+'\n')
                    cluster = []
                else:
                    if last_stop != 0:
                        if int(stop)-int(start)>50:
                            converged.append(line)
                        else:
                            converged.append(chrm+'\t'+str(int(start)-25)+'\t'+str(int(stop)+25)+'\n')
            last_stop = int(stop)
            last_value = float(value)

    #write converged file, also create dictionary used for data matrix creation, add names for each data point (converged peaks)
    bedfile_c = open('converged_'+bedfile,'w')
    for entry in converged:
        bedfile_c.write(entry)
    return converged

def create_data_dict(bedfile, col, arg):
    data = {}
    translate = {'-m':'motif', '-ch':'characterization','-tf':'GATA1','-dm':'dna_methyl','-hm':'H3K4me3','-ha':'H3K27ac','-d1':'DNAse_1','-na':''}
    data['header'] = [translate[arg]]
    for line in open(bedfile,'r'):
        if arg !='-na':
            chrm, start, stop, value = line.strip().split('\t')[:3] + [line.strip().split('\t')[col]]
            data[chrm+':'+start+'-'+stop] = {translate[arg]:value}
        else:
            chrm, start, stop= line.strip().split('\t')[0:3]
            data[chrm+':'+start+'-'+stop] = {}
    return data
 
#make sure input type is specified
args = set(['-m', '-ch','-tf','-dm','-hm','-ha','-d1','-na']) #possible input types
if len(args.intersection(set(argv)))>1:
    exit('Only one input type may be specified')
elif len(args.intersection(set(argv)))==0:
    exit('Must specify input type using -m,-na,-ch,-tf,-dm,-hm,-ha, or -d1')

print('reading input bedfile...')
#read in bedfile
bedfile = argv[1]
#find col number
try: #if there is an int for col number
    col = int(argv[argv.index('-cl')+1]) #get the command line arg following -cl to denote col number
except ValueError or IndexError: #no int for col number, or no col number supplied
    if '-na' not in argv:
        exit('Must supply column number for value of input data (eg -cl 3)')
    else:
        col = -1
#converge if necessary
if '-c' in argv: #meaning the input file is not already converged
    print('converging bedfile...')
    converged = converge_bed(bedfile, col)
    bedfile = 'converged_'+bedfile   
else:
    converged = []
    for line in open(bedfile,'r'):
        converged.append(line)

#create data structure
print('creating data structure...')
data = create_data_dict(bedfile, col, args.intersection(set(argv)).pop())
#! for each of the variables, run bedtools intersect -a -b -wb, store the output in a variable
if '-dm' not in argv:
    print('gathering dna methylation data...')
    dnaMethyl = 'GATA1/param_data/K562_MeDIP.bed'
    i_dnam = subprocess.check_output(['bedtools','intersect','-wo','-a',bedfile,'-b',dnaMethyl]).split('\n')
    data['header'].append('dna_methyl')
    dnam_col = -2
    i_dnam.remove('')
    for entry in i_dnam:
        chrm, start, stop, value = entry.split('\t')[:3] + [entry.split('\t')[dnam_col]]
        if 'dna_methyl' not in data[chrm+':'+start+'-'+stop]: #if we haven't seen this data point 
            data[chrm+':'+start+'-'+stop]['dna_methyl'] = value
        else: #if we have seen this data point, replace with the larger of the two
            data[chrm+':'+start+'-'+stop]['dna_methyl'] = str(max(float(value),float(data[chrm+':'+start+'-'+stop]['dna_methyl'])))
    for key in data:
        if 'dna_methyl' not in data[key]:
            data[key]['dna_methyl'] = '0' #maybe NA? I'm still not convinced we can treat this as 0

if '-hm' not in argv:
    print('gathering H3K4me3 data...')
    H3K4me3 = 'GATA1/param_data/K562_H3K4me3.bedGraph'
    i_H3K4 = subprocess.check_output(['bedtools','intersect','-wo','-a',bedfile,'-b',H3K4me3]).split('\n')
    data['header'].append('H3K4me3')
    H3K4_col = -2
    i_H3K4.remove('')
    for entry in i_H3K4:
        chrm, start, stop, value = entry.split('\t')[:3]+ [entry.split('\t')[H3K4_col]]
        if 'H3K3me3' not in data[chrm+':'+start+'-'+stop]: #if we haven't seen this data point 
            data[chrm+':'+start+'-'+stop]['H3K4me3'] = value
        else: #if we have seen this data point, replace with the larger of the two
            data[chrm+':'+start+'-'+stop]['H3K4me3'] = str(max(float(value), float(data[chrm+':'+start+'-'+stop]['H3K4me3'])))
    for key in data:
        if 'H3K4me3' not in data[key]:
            data[key]['H3K4me3'] = '0' #maybe NA? I'm still not convinced we can treat this as 0

if '-ha' not in argv:
    print('gathering H3K27ac data...')
    H3K27ac = 'GATA1/param_data/K562_H3K27ac.bedGraph'
    i_H3K27 = subprocess.check_output(['bedtools','intersect','-wo','-a',bedfile,'-b',H3K27ac]).split('\n')
    data['header'].append('H3K27ac')
    H3K27_col = -2
    i_H3K27.remove('')
    for entry in i_H3K27:
        chrm, start, stop, value = entry.split('\t')[:3] + [entry.split('\t')[H3K27_col]]
        if 'H3K27ac' not in data[chrm+':'+start+'-'+stop]: #if we haven't seen this data point 
            data[chrm+':'+start+'-'+stop]['H3K27ac'] = value
        else: #if we have seen this data point, replace with the larger of the two
            data[chrm+':'+start+'-'+stop]['H3K27ac'] = str(max(float(value),float(data[chrm+':'+start+'-'+stop]['H3K27ac'])))
    for key in data:
        if 'H3K27ac' not in data[key]:
            data[key]['H3K27ac'] = '0' #maybe NA? I'm still not convinced we can treat this as 0

if '-ch' not in argv:
    print('gathering characterization matrix data...')
    charMatrix = 'ENCODE/ENCODE.5groupCRE.bed'
    i_char = subprocess.check_output(['bedtools','intersect','-wo','-a',bedfile,'-b',charMatrix]).split('\n')
    data['header'].append('characterization')
    char_col = -2
    i_char.remove('')
    for entry in i_char:
        chrm, start, stop, value = entry.split('\t')[:3] + [entry.split('\t')[char_col]]
        data[chrm+':'+start+'-'+stop]['characterization'] = value #no idea how to choose here, they should be the same I'd imagine
    for key in data:
        if 'characterization' not in data[key]:
            data[key]['characterization'] = 'None' #maybe NA? I'm still not convinced we can treat this as 0

#some lines changed to process narrowpeak gata1
if '-tf' not in argv:
    print('gathering transcription factor data...')
    tf = 'GATA1_sorted.narrowPeak'#'GATA1/ENCFF616ZFJ.bedGraph'
    i_tf = subprocess.check_output(['bedtools','intersect','-wo','-a',bedfile,'-b',tf]).split('\n')
    data['header'].append('GATA1')
    tf_col = -2
    i_tf.remove('')
    for entry in i_tf:
        chrm, start, stop, value = entry.split('\t')[:3] + [entry.split('\t')[tf_col]]
        if 'GATA1' not in data[chrm+':'+start+'-'+stop]: #if we haven't seen this data point 
            data[chrm+':'+start+'-'+stop]['GATA1'] = '1'#binary value used for narrowPeak
        else: #if we have seen this data point, replace with the larger of the two
            data[chrm+':'+start+'-'+stop]['GATA1'] = '1'#str(max(float(value),float(data[chrm+':'+start+'-'+stop]['GATA1'])))
    for key in data:
        if 'GATA1' not in data[key]:
            data[key]['GATA1'] = '0' #maybe NA? I'm still not convinced we can treat this as 0


if '-d1' not in argv:
    print('gathering dnase 1 data...')
    dnase = 'K562_DNAse_1.bedGraph'
    i_dnase = subprocess.check_output(['bedtools','intersect','-wo','-a',bedfile,'-b',dnase]).split('\n')
    data['header'].append('DNAse_1')
    dnase_col = -1
    i_dnase.remove('')
    for entry in i_dnase:
        chrm, start, stop, value = entry.split('\t')[:3] + [entry.split('\t')[dnase_col]]
        if 'DNAse_1' not in data[chrm+':'+start+'-'+stop]: #if we haven't seen this data point 
            data[chrm+':'+start+'-'+stop]['DNAse_1'] = value
        else: #if we have seen this data point, replace with the larger of the two
            data[chrm+':'+start+'-'+stop]['DNAse_1'] = str(max(float(value),float(data[chrm+':'+start+'-'+stop]['DNAse_1'])))
    for key in data:
        if 'DNAse_1' not in data[key]:
            data[key]['DNAse_1'] = '0' #maybe NA? I'm still not convinced we can treat this as 0

#! also declare the meme filename for later
meme = 'GATA1/MA0140.2_truncated.meme'

if '-m' not in argv:
    print('generating fastas for motifs...')
    #for the motifs, create a temp fasta file which we will later delete
    temp = open('temp.fa','w')
    for entry in converged:
        chrm, start, stop = entry.split('\t')[:3]
        temp.write(pysam.faidx('hg19.fa',chrm+':'+start+'-'+stop))

    print('finding motifs...')
    #then, run fimo and read in the results from fimo_out/fimo.tsv 
    subprocess.call(['fimo '+meme+ ' temp.fa'], shell=True)
    #3rd col of fimo.tsv ** should ** be keys in data, store in data
    fimo = open('fimo_out/fimo.tsv','r')
    fimo.readline() #read header out
    data['header'].append('motif')
    for line in fimo: #add motif presence to data
        if line[0]!='\n': #there are some lines at the end of this file we don't need to read
            data[line.split('\t')[2]]['motif'] = line.split('\t')[7]
        else:
            break
    for key in data: #add motif absence to data
        if 'motif' not in data[key]:
            data[key]['motif'] = '1'

else: #meaning input was motif occurences
    pass

#create data matrix from these results
#write header
outfile.write('name'+'\t'.join(data['header'])+'\n')
#write body of data
if '' in data['header']:
    data['header'].remove('')
print(data['header'])
print('writing output file...')
for key in data:
    if key != 'header':
        line = [key]
        for vble in data['header']:
            line.append(data[key][vble])
        outfile.write('\t'.join(line)+'\n')
print('Done!')
#clean up? do we want to delete/move temp files?


