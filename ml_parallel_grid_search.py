import joblib
import numpy as np
from multiprocessing import Pool
import pandas
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn import svm
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


#import data with pandas
print('importing dataset...')
data = pandas.read_csv('130k_points_w_dnase_data_matrix.txt',sep='\t')

#code characterization into 0 and 1 for 3 variables promoter, enhancer, CTCF (None is 0,0,0)
print('transforming data to numerical...')
promoter = []
enhancer = []
CTCF = []
for index, row in data.iterrows():
    if row['characterization']=='promoter':
        promoter.append(1)
    else:
        promoter.append(0)
    if row['characterization']=='enhancer':
        enhancer.append(1)
    else:
        enhancer.append(0)
    if row['characterization']=='CTCF-only':
        CTCF.append(1)
    else: 
        CTCF.append(0)

data['promoter']=promoter
data['enhancer']=enhancer
data['CTCF']=CTCF

#code T/F for motif presence
motif_bin = []
for index, row in data.iterrows():
    if row['motif']=='True':
        motif_bin.append(1)
    else:
        motif_bin.append(0)

data['motif_bin']=motif_bin


#separate X, Y
print('preprocessing...')
X = StandardScaler().fit_transform(data[['motif_bin','promoter','enhancer','CTCF','DNAse_1','H3K4me3','H3K27ac','dna_methyl']])
Y = data['GATA1']

#set parameter search space, perform grid search
print('performing grid search...')
cspace = np.logspace(-2,7,10)
gspace = np.logspace(-7,2,10)
param_grid = dict(gamma=gspace, C=cspace)
#cv = StratifiedShuffleSplit(n_splits=5, test_size=0.1, random_state=42)
grid = GridSearchCV(svm.SVR(kernel='rbf'), param_grid=param_grid)
grid.fit(X, Y)

print("The best parameters are %s with a score of %0.2f"
      % (grid.best_params_, grid.best_score_))

#heatmap of results
print('creating heatmap...')
plt.figure(figsize=(8, 6))
plt.subplots_adjust(left=.2, right=0.95, bottom=0.15, top=0.95)
plt.imshow(scores, interpolation='nearest', cmap=plt.cm.hot,
           norm=MidpointNormalize(vmin=0.2, midpoint=0.75))
plt.xlabel('gamma')
plt.ylabel('C')
plt.colorbar()
plt.xticks(np.arange(len(gspace)), gspace, rotation=45)
plt.yticks(np.arange(len(cspace)), cspace)
plt.title('Validation accuracy')
plt.savefig('gridsearch_heatmap_rbf.png')
 
