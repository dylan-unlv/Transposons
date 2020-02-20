import pandas as pd
#matrix = pd.read_csv('ERVK_RLTR13D6_dist_matrix_trimmed.txt',sep=' ')
''' #this bit was used to fill in the missing upper diag values with 0s
f = open('ERVK_RLTR13D6_dist_matrix_trimmed.txt','r')
lines = f.readlines()
for i in range(len(lines)):
    line = lines[i].strip().split(' ')
    if '' in line:
        line = [x for x in line if x != '']
    lines[i] = '\t'.join(line + ['0']*(len(lines[i:]) -1))+'\n'
f.close()
f = open('ERVK_RLTR13D6_dist_matrix_pandas.txt','w')
f.writelines(lines)
'''
matrix = pd.read_csv('ERVK_RLTR13D6_dist_matrix_pandas.txt',header=None,sep='\t')
qcounts = []
for index, row in matrix.iterrows():
    try:
        qcounts.append((index,matrix[index].value_counts()['?']))
    except:
        pass
qcounts.sort(key=lambda x: x[1],reverse=True)








