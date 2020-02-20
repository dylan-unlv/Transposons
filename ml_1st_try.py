import joblib
import pandas
from sklearn.model_selection import train_test_split
from sklearn import svm
from sklearn import linear_model as lm
from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import cross_val_score

#import data with pandas
data = pandas.read_csv('130k_peaks_w_dnase_data_matrix.txt',sep='\t')

#code characterization into 0 and 1 for 3 variables promoter, enhancer, CTCF (None is 0,0,0)
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
X = data[['motif_bin','promoter','enhancer','CTCF','DNAse_1','H3K4me3','H3K27ac','dna_methyl']]
Y = data['GATA1']

X_train, X_test, Y_train, Y_test = train_test_split(X,Y,test_size= 0.3,random_state=5)

'''
clf0= GaussianNB().fit(X_train, Y_train)
y_pred=clf0.predict(X_test)
print("Number of mislabeled points out of a total %d points : %d"
       % (X_test.shape[0], (Y_test != y_pred).sum()))

clf1 = lm.LogisticRegression(C=50,penalty='elasticnet',l1_ratio=0.5,solver='saga',n_jobs=20,max_iter=10000)
scores = cross_val_score(clf1, X, Y, cv=10)
print(scores)

clf2 = lm.LogisticRegression(C=50,penalty='l2',solver='saga',n_jobs=20,max_iter=10000).fit(X_train, Y_train)
print(['l2',clf2.score(X_test,Y_test),clf2.get_params()])

clf3 = lm.LogisticRegression(C=100,penalty='elasticnet',solver='saga',n_jobs=20,l1_ratio=0.5,max_iter=10000).fit(X_train, Y_train)
print(['elasticnet',clf3.score(X_test,Y_test),clf3.get_params()])

clf4 = lm.LogisticRegression(C=5,penalty='none',solver='saga',n_jobs=20,max_iter=10000).fit(X_train, Y_train)
print(['none',clf4.score(X_test,Y_test),clf4.get_params()])
'''

score = []
mls = ['l1','l2','elasticnet','none']
for i in mls:
    print(i)
    clf = lm.LogisticRegression(penalty=i,solver='saga',n_jobs=20,l1_ratio=0.5,max_iter=10000).fit(X_train, Y_train)
    score.append(clf.score(X_test,Y_test))
    #joblib.dump(clf, i+'_GATA1_ml.joblib')
    clf = ''
for i in range(len(mls)):
    print([mls[i],score[i]])    



 
