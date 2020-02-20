import joblib
import pandas
import numpy as np 
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn import svm
from sklearn import linear_model as lm
from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import cross_val_score
from sklearn.metrics import brier_score_loss
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import classification_report

#import data with pandas
data = pandas.read_csv('15k_peaks_data_matrix.txt',sep='\t')

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
print(X)
X_train, X_test, Y_train, Y_test = train_test_split(X,Y,test_size= 0.3,random_state=5)
min_max_scaler = preprocessing.MinMaxScaler()
X_train_scaled = min_max_scaler.fit_transform(X_train)
X_test_scaled = min_max_scaler.transform(X_test)
print("train: %s 1's out of %s" % (sum(Y_train), len(Y_train)))
print("test: %s 1's out of %s" % (sum(Y_test), len(Y_test)))

'''
clf0= GaussianNB().fit(X_train, Y_train)
y_pred=clf0.predict(X_test)
print("Number of mislabeled points out of a total %d points : %d"
       % (X_test.shape[0], (Y_test != y_pred).sum()))

clf1 = lm.LogisticRegression(C=50,penalty='elasticnet',l1_ratio=0.5,solver='saga',n_jobs=20,max_iter=10000)
scores = cross_val_score(clf1, X, Y, cv=10)
print(scores)

clf2 = lm.LogisticRegression(C=50,penalty='l2',solver='saga',n_jobs=20,max_iter=10000).fit(X_train_scaled, Y_train)
print(['l2',clf2.score(X_test,Y_test),clf2.get_params()])

clf3 = lm.LogisticRegression(C=100,penalty='elasticnet',solver='saga',n_jobs=20,l1_ratio=0.5,max_iter=10000).fit(X_train_scaled, Y_train)
print(['elasticnet',clf3.score(X_test,Y_test),clf3.get_params()])

clf4 = lm.LogisticRegression(C=5,penalty='none',solver='saga',n_jobs=20,max_iter=10000).fit(X_train_scaled, Y_train)
print(['none',clf4.score(X_test,Y_test),clf4.get_params()])
'''

score = []
logmls = ['l1','l2','elasticnet','none']
svmls = ['linear','rbf','poly','sigmoid']
for C in (1, 0.1, 0.01, 0.0000001):
    for i in logmls:
        print(i)
        clf = lm.LogisticRegression(C=C, penalty=i,solver='saga',n_jobs=20,l1_ratio=0.5,max_iter=10000).fit(X_train_scaled, Y_train)

        coef_LR = clf.coef_.ravel()
        sparsity_LR = np.mean(coef_LR == 0) * 100
        print("C=%.2f" % C)
        print("{:<40} {:.2f}%".format("Sparsity:", sparsity_LR))

        y_pred = clf.predict(X_test_scaled)
        if hasattr(clf, "predict_proba"):
            prob_pos = clf.predict_proba(X_test_scaled)[:, 1]
        else:  # use decision function
            prob_pos = clf.decision_function(X_test_scaled)
            prob_pos = \
                (prob_pos - prob_pos.min()) / (prob_pos.max() - prob_pos.min())

        clf_score = brier_score_loss(Y_test, prob_pos, pos_label=1)
        print("%s:" % i)
        print("\tBrier: %1.3f" % (clf_score))
        print("\tPrecision: %1.3f" % precision_score(Y_test, y_pred))
        print("\tRecall: %1.3f" % recall_score(Y_test, y_pred))
        print("\tF1: %1.3f\n" % f1_score(Y_test, y_pred))
        print(classification_report(Y_test, y_pred))
        print('weights in trained model:\t' + str(clf.coef_))
        score.append(clf.score(X_test_scaled,Y_test))
        #joblib.dump(clf, i+'_GATA1_ml.joblib')
        clf = ''
    for i in range(len(logmls)):
        print([logmls[i],score[i]])    
    for i in svmls:
        print(i)
        clf = svm.SVC(C=C,kernel=i).fit(X_train_scaled, Y_train)

        print("C=%.2f" % C)
        print("{:<40} {:.2f}%".format("Sparsity:", sparsity_LR))

        y_pred = clf.predict(X_test_scaled)
        prob_pos = clf.decision_function(X_test_scaled)
        prob_pos = \
            (prob_pos - prob_pos.min()) / (prob_pos.max() - prob_pos.min())

        clf_score = brier_score_loss(Y_test, prob_pos, pos_label=1)
        print("%s:" % i)
        print("\tBrier: %1.3f" % (clf_score))
        print("\tPrecision: %1.3f" % precision_score(Y_test, y_pred))
        print("\tRecall: %1.3f" % recall_score(Y_test, y_pred))
        print("\tF1: %1.3f\n" % f1_score(Y_test, y_pred))
        print(classification_report(Y_test, y_pred))
        if i=='linear':
            print('weights in trained model:\t' + str(clf.coef_))
        score.append(clf.score(X_test_scaled,Y_test))
        #joblib.dump(clf, i+'_GATA1_ml.joblib')
        clf = ''
    
    


 
