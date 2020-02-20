#this program was mostly copied from travis's classification 3.0.py script, thanks travis!
import warnings
import joblib
import pandas as pd
import numpy as np 
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn import svm
from sklearn import linear_model as lm
from sklearn.model_selection import cross_val_score
from sklearn.metrics import brier_score_loss
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.metrics import classification_report, confusion_matrix

#ignore repeat User Warnings (often from sklearn about extra params being ignored)
#warnings.filterwarnings(action='once', category=UserWarning)

#import data with pandas
data = pd.read_csv('15k_peaks_data_matrix.txt',sep='\t')

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


#separate X, Y & scale X
min_max_scaler = preprocessing.MinMaxScaler()
X = data[['motif_bin','promoter','enhancer','CTCF','DNAse_1','H3K4me3','H3K27ac','dna_methyl']]
Y = data['GATA1']
#create Cross Validation Splitter
inner_k = 4
outer_k = 6
inner_cv_skf = StratifiedKFold(n_splits=inner_k, shuffle=True,  random_state=1)
outer_cv_skf = StratifiedKFold(n_splits=outer_k, shuffle=True,  random_state=1)

#outer CV to estimate generalization error
non_nested_score_svm = np.zeros(outer_k)
non_nested_score_log = np.zeros(outer_k)
nested_score_svm = np.zeros(outer_k)
nested_score_log = np.zeros(outer_k)
outer_k_index = 0
iteration_k_forSVM = 1
iteration_k_forLOG = 1

#construct pipelines for each ml algorithm
pipe_log = Pipeline([('clf', lm.LogisticRegression(solver='saga',max_iter=10000,random_state=1))])
pipe_svm = Pipeline([('clf', svm.SVC(max_iter=10000,random_state=1))])
Cs = [1e-8, 1e-6, 1e-4, 1e-2, 1, 1e2]
l1s = [0,0.25,0.5,0.75,1]
logmls = ['l1','l2','elasticnet','none']
svmls = ['linear','rbf','poly','sigmoid']
grid_params_log = [{'clf__C':Cs, 'clf__penalty':logmls,'clf__l1_ratio':l1s}]
grid_params_svm = [{'clf__C':Cs, 'clf__kernel':svmls}]


#first loop splits data up randomly into test/train
for outer_train_index, test_index in outer_cv_skf.split(X=X, y=Y):
    X_train = X.iloc[outer_train_index]
    Y_train = Y.iloc[outer_train_index]
    X_test = X.iloc[test_index]
    Y_test = Y.iloc[test_index]

    #must use DataFrame wrapper for scaling because sklearn's scaler returns numpy array with no colnames
    X_train = pd.DataFrame(min_max_scaler.fit_transform(X_train.to_numpy()), columns=X_train.columns, index=X_train.index)
    X_test = pd.DataFrame(min_max_scaler.transform(X_test.to_numpy()), columns=X_test.columns, index=X_test.index)

    #construct grids to search over in nested for loop
    gs_log = GridSearchCV(estimator=pipe_log, param_grid=grid_params_log, scoring='f1_macro',
            cv=inner_cv_skf.split(X_train,Y_train), n_jobs=20,return_train_score=True,refit=False,iid=False)
    gs_svm = GridSearchCV(estimator=pipe_svm, param_grid=grid_params_svm, scoring='f1_macro',
            cv=inner_cv_skf.split(X_train,Y_train), n_jobs=20, return_train_score=True,refit=False,iid=False)
    grids = [gs_log, gs_svm]
    grid_dict = {0:'Logistic Regression',1:'SVM'}
    best_gs_log = ''
    best_gs_svm = ''

    print("Choosing best classifier")
    for idx, gs in enumerate(grids):
        print("Classifier: {}".format(grid_dict[idx]))
        gs.fit(X_train, Y_train)
        print("Best params found on development set: {}".format(gs.best_params_))
        print("Best validation F1: %1.3f" % gs.best_score_)

        meantrainingscore = gs.cv_results_['mean_train_score'][gs.best_index_]
        stdtrainingscore = gs.cv_results_['std_train_score'][gs.best_index_]
        meanvalidationscore = gs.cv_results_['mean_test_score'][gs.best_index_]
        stdvalidationscore = gs.cv_results_['std_test_score'][gs.best_index_]
        inner_k_index=0
        #generate cross validation prediction with 10 fold stratified sampling
        for train_index, validation_index in inner_cv_skf.split(X_train,Y_train):
            gs.estimator.set_params(**gs.best_params_).fit(X_train.iloc[train_index], Y_train.iloc[train_index])
            ypred = gs.estimator.predict(X_train.iloc[validation_index])
            training_report_cv = classification_report(Y_train.iloc[validation_index], ypred)
            f1_cv_macro = f1_score(Y_train.iloc[validation_index], ypred, average='macro')
            f1_cv_weighted = f1_score(Y_train.iloc[validation_index], ypred, average='weighted')
            if grid_dict[idx] == 'SVM':
                print("SVM VALIDATION REPORT ON CV iteration #" +str(iteration_k_forSVM) + "\n " + str(training_report_cv) + "\n")
                print("SVM VALIDATION F1 score (macro) ON CV iteration #" +str(iteration_k_forSVM) + ":\t" + str(f1_cv_macro) + "\n")
                print("SVM VALIDATION F1 score (weighted) ON CV iteration #" +str(iteration_k_forSVM) + ":\t" + str(f1_cv_weighted) + "\n\n")
                iteration_k_forSVM += 1
            elif grid_dict[idx] == 'Logistic Regression':
                print("LOG REGRESSION VALIDATION REPORT ON CV iteration #" +str(iteration_k_forLOG) + "\n " + str(training_report_cv) + "\n")
                print("LOG REGRESSION VALIDATION F1 score (macro) ON CV iteration #" +str(iteration_k_forLOG) + ":\t" + str(f1_cv_macro) + "\n")
                print("LOG REGRESSION VALIDATION F1 score (weighted) ON CV iteration #" +str(iteration_k_forLOG) + ":\t" + str(f1_cv_weighted) + "\n\n")
                iteration_k_forLOG += 1
 
            inner_k_index += 1
        all_cv_results = pd.DataFrame(gs.cv_results_)
        # print(all_cv_results)
        if grid_dict[idx] == 'SVM':
            print('    [DEBUG] saving cv_results_ as csv for SVM')
            all_cv_results.to_csv('SVM_CV_Results.tsv', sep = "\t", index=True)
        elif grid_dict[idx] == 'Logistic Regression':
            print('    [DEBUG] saving cv_results_ csv for LOG')
            all_cv_results.to_csv('LR_CV_Results.tsv', sep = "\t", index=True)
        # print('  [All CV RESULTS training F1] --- END ---')
 
        if grid_dict[idx] == 'SVM':
            print("SVM mean training score   ON development DATA iteration #" +str(outer_k_index+1) + ":\t " + str(meantrainingscore) + "\n")
            print("SVM std training score    ON development DATA iteration #" +str(outer_k_index+1) + ":\t " + str(stdtrainingscore) + "\n")
            print("SVM mean validation score ON development DATA iteration #" +str(outer_k_index+1) + ":\t " + str(meanvalidationscore) + "\n")
            print("SVM std validation score  ON development DATA iteration #" +str(outer_k_index+1) + ":\t " + str(stdvalidationscore) + "\n")

        elif grid_dict[idx] == 'Logistic Regression':
            print("Log Regression mean training score    ON development DATA iteration #" +str(outer_k_index+1) + ":\t " + str(meantrainingscore) + "\n")
            print("Log Regression training score         ON development DATA iteration #" +str(outer_k_index+1) + ":\t " + str(stdtrainingscore) + "\n")
            print("Log Regression mean validation score  ON development DATA iteration #" +str(outer_k_index+1) + ":\t " + str(meanvalidationscore) + "\n")
            print("Log Regression std validation score   ON development DATA iteration #" +str(outer_k_index+1) + ":\t " + str(stdvalidationscore) + "\n")

 
        # Nested CV with parameter optimization
 
        # Predict on test data with best params
        gs.estimator.set_params(**gs.best_params_).fit(X_train, Y_train)
        Y_pred = gs.estimator.predict(X_test)
        # Test data f1 of model with best params
        f1_macro = f1_score(Y_test, Y_pred, average='macro')
        f1_weighted = f1_score(Y_test, Y_pred, average='weighted')  #nested score
 
        print('Evaluation set F1 score (macro) for best params: %.3f ' % f1_macro)
        print('Evaluation set F1 score (weighted) for best params: %.3f ' % f1_weighted)
 
        # Print and Save test result on evaluation data
        if grid_dict[idx] == 'SVM':
            svm_test_report_evaluation = classification_report(Y_test, Y_pred)              #https://joshlawman.com/metrics-classification-report-breakdown-precision-recall-f1/
            svm_test_con_matrix_evaluation = confusion_matrix(Y_test, Y_pred)                    #https://stats.stackexchange.com/questions/95209/how-can-i-interpret-sklearn-confusion-matrix
            print("SVM CLASSIFICATION REPORT ON EVALUATION DATA iteration #" +str(outer_k_index+1) + "\n " + str(svm_test_report_evaluation) + "\n" )
            print("SVM TEST ACCURACY (f1-score macro) ON EVALUATION DATA iteration #" +str(outer_k_index+1) + ": " + str(f1_macro) + "\n")
            print("SVM TEST ACCURACY (f1-score weighted) ON EVALUATION DATA iteration #" +str(outer_k_index+1) + ": " + str(f1_weighted) + "\n\n")
            print(str(svm_test_con_matrix_evaluation) + '\n\n')
            best_gs_svm = gs.estimator
            nested_score_svm[outer_k_index] = f1_macro #should be macro or weighted?
            non_nested_score_svm[outer_k_index] = meanvalidationscore
        elif grid_dict[idx] == 'Logistic Regression':
            log_test_report_evaluation = classification_report(Y_test, Y_pred)              #https://joshlawman.com/metrics-classification-report-breakdown-precision-recall-f1/
            log_test_con_matrix_evaluation = confusion_matrix(Y_test, Y_pred)                    #https://stats.stackexchange.com/questions/95209/how-can-i-interpret-sklearn-confusion-matrix
            print("LOG REGRESSION CLASSIFICATION REPORT ON EVALUATION DATA iteration #" +str(outer_k_index+1) + "\n " + str(log_test_report_evaluation) + "\n" )
            print("LOG REGRESSION TEST ACCURACY (f1-score macro) ON EVALUATION DATA iteration #" +str(outer_k_index+1) + ": " + str(f1_macro) + "\n")
            print("LOG REGRESSION TEST ACCURACY (f1-score weighted) ON EVALUATION DATA iteration #" +str(outer_k_index+1) + ": " + str(f1_weighted) + "\n\n")
            print(str(log_test_con_matrix_evaluation) + '\n\n')
            best_gs_log = gs.estimator
            nested_score_log[outer_k_index] = f1_macro
            non_nested_score_log[outer_k_index] = meanvalidationscore
    outer_k_index += 1
joblib.dump(best_gs_svm, 'gata1_svm.pkl', compress=1)
joblib.dump(best_gs_log, 'gata1_log.pkl', compress=1)

#====================================
'''


X_train, X_test, Y_train, Y_test = train_test_split(X,Y,test_size= 0.3,random_state=5)
min_max_scaler = preprocessing.MinMaxScaler()
X_train_scaled = min_max_scaler.fit_transform(X_train)
X_test_scaled = min_max_scaler.transform(X_test)
print("train: %s 1's out of %s" % (sum(Y_train), len(Y_train)))
print("test: %s 1's out of %s" % (sum(Y_test), len(Y_test)))

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
    
    
'''

 
