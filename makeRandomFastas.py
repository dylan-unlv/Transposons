#not my script -dylan
import sys
import pysam
import random
import os

znfsfile = sys.argv[1]
gtffile = sys.argv[2]
sigresultsfile = sys.argv[3]
outputpath = sys.argv[4]
numberofruns = sys.argv[5]

positionsdict = {}
resultsdict = {}
znfsdict = {}
znfslist = []

def ReadZNFslist():
  #Get the list of ZNFs available in the ChIPSeq data
  with open(znfsfile, 'r') as znfs:
    print('Reading ZNFs list...')
    for line in znfs:
      znf = line.rstrip('\n')
      znfslist.append(znf)

def GetTranscriptPositions():
  #Get the positions of each of the LINE instances to be used later to grab sequences from the FASTA file
  with open(gtffile, 'r') as gtf:
    print('Grabbing transcript positions...')
    for line in gtf:
      linesp = line.rstrip('\n').split('\t')
      chrom = linesp[0].lstrip('chr')
      leftpos = linesp[3]
      rightpos = linesp[4]
      position = chrom + ':' + leftpos + '-' + rightpos
      info = linesp[8]
      #name of duplicate
      transcript = info.split(';')[1].split('"')[1]
      #only get duplicates from the LINE family
      family = info.split(';')[3].split('"')[1]
      if True: #family == 'LINE':
        #make one dictionary where the key is the instance name and the value is the position
        positionsdict[transcript] = position
        #make another dictionary where the key and value are switched
        znfsdict[position] = transcript

def GrabAllCombinations():
  #get all possible combinations of zinc finger genes and LINE instances to randomly pick from
  with open(sigresultsfile, 'r') as results:
    print('Grabbing possible combinations...')
    header = results.readline()
    for line in results:
      linesp = line.rstrip('\n').split('\t')
      transcript = linesp[0].split(':')[0]
      znfgene = linesp[1].split('|')[0]
      #first make sure that the ZNF is also in ChIPSeq data
      if znfgene in znfslist:
        if znfgene not in resultsdict:
         resultsdict[znfgene] = []
        #create a list of LINE instances that go with each ZNF gene
        resultsdict[znfgene].append(transcript)

def PickRandomCombinations(x):
  #randomly select a specified number of zinc finger gene - LINE instance pairs
  pickedcombos = set()
  combosdict = {}
  for i in range(x):
    #pick a random zinc finger gene
    znfgene = random.choice(resultsdict.keys())
    #pick a random transcript from that randomly picked gene
    transcript = random.choice(resultsdict[znfgene])
    #make sure the combination hasn't been picked already
    while (znfgene, transcript) in pickedcombos:
      transcript = random.choice(resultsdict[znfgene])
      znfgene = random.choice(resultsdict.keys())
    pickedcombos.add((znfgene, transcript))
    if znfgene not in combosdict:
      combosdict[znfgene] = []
    combosdict[znfgene].append(transcript)
  return(combosdict)

def WriteOutput(combosdict, i):
  #make new folder for each run
  if not os.path.exists(outputpath + str(i)):
    os.makedirs(outputpath + str(i))
  for znf, L1s in combosdict.iteritems():
    #make new file for each zinc finger gene
    with open(outputpath + str(i) + '/' + znf + '_LINEs.fasta', 'w+') as output:
      commandarg = ['/DataDrives/dd4/GRCH37_Reference_Fasta/human_g1k_v37.fasta']
      for transcript in L1s:
        commandarg.append(positionsdict[transcript])
      #get fasta sequences from human genome for each LINE instance
      sequences = pysam.faidx(*commandarg)
      for j in range(len(sequences)):
        #write instance name instead of position
        if sequences[j].startswith('>'):
          output.write('>' + znfsdict[sequences[j].rstrip('\n').lstrip('>')] + '\n')
        else:
          output.write(sequences[j])

ReadZNFslist()
GetTranscriptPositions()
GrabAllCombinations()
for i in range(1, int(numberofruns) + 1):
  print('Run' + str(i))
  combosdict = PickRandomCombinations(10782) #10782
  print('Done picking random combinations, outputting...')
  WriteOutput(combosdict, i)
