from os import listdir
import subprocess
gata1 = 'GATA1/GSM1003608_hg19_wgEncodeSydhTfbsK562Gata1bIggmusPk.narrowPeak'
data = []
for filename in listdir('TF_CHiPSeq'):
    if 'bigBed' not in filename:
        output = subprocess.check_output(['bedtools','intersect','-a',gata1,'-b','TF_CHiPSeq/'+filename,'-wa']).split('\n')
        output.remove('')
        output = set(output)
        data.append((filename,len(output)))
        if len(output)>7494:
            for i in output:
                print i
            exit()

data.sort(key=lambda x:x[1], reverse=True)

output = open('TF_Overlap_GATA1.txt','w')
#7494 is the num of peaks in the GATA1 file
for i in data:
    output.write("{} has a {}% overlap with GATA1 narrowpeak file\n".format(i[0],i[1]/74.94))
    
