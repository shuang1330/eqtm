import numpy as np
import pandas as pd
import pickle
import os
from lib.read_data import dataset,Datasets
from math import copysign

from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline

# feature extractors
from sklearn.decomposition import PCA
# classifiers
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier,AdaBoostClassifier
# from sklearn.linear_model import ElasticNet
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
# finetuning
from sklearn.model_selection import GridSearchCV
# validation
from sklearn import metrics
from sklearn.metrics import confusion_matrix

# import matplotlib.pyplot as plt
# TODO: support matplotlib drawing
def read_data_set(data_table,test_size=0.25):
    '''
    convert a pandas dataframe data table into Datasets(dataset,dataset)
    '''
    train, test = train_test_split(data_table,test_size=0.25)
    train_x = train[[col for col in train.columns
                     if col not in ['zscore','direction','cpgName']]]
    features = train_x.columns
    train_x = np.array(train_x)
    test_x = np.array(test[[col for col in train.columns
                      if col not in ['zscore','direction','cpgName']]])
    train_y = np.array(train['direction'],dtype=np.int8)
    test_y = np.array(test['direction'],dtype=np.int8)

    return Datasets(train=dataset(train_x,train_y),
                    test=dataset(test_x,test_y))

def draw_roc_curve(fpr,tpr,score):
    '''
    draw roc curve
    '''
    plt.figure()
    lw = 2
    plt.plot(fpr, tpr, color='darkorange',
            lw=lw, label='ROC curve (area = %0.2f)' % score)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC curve')
    plt.legend(loc="lower right")
    plt.show()

def run_display_output(classifier,test,DRAW=False):
    '''
    get confusion matrix and auc score for test dataset
    (optional) draw roc curve
    '''
    pred = classifier.predict(test.values)
    tn, fp, fn, tp = confusion_matrix(test.labels,pred).ravel()#confusion matrix
    print(tn,fp,fn,tp)
    sensitivity = tp/(fn+tp)
    specificity = tn/(fp+tn)
    prods = classifier.predict_proba(test.values)[:,1]
    fpr, tpr, _ = metrics.roc_curve(test.labels, prods)
    score = metrics.auc(fpr,tpr) #auc score
    if DRAW:
        draw_roc_curve(fpr,tpr,score)

    return sensitivity, specificity, score

def read_gavin(gavin_res, labels):
    '''
    compare gavin results with labels for a certain subset of data
    '''
    gavin_res = gavin_res.replace('Pathogenic',1)
    gavin_res = gavin_res.replace('Benign',0)
    tn_g, fp_g, fn_g, tp_g = \
    confusion_matrix(labels, gavin_res.astype(np.int8)).ravel()
    sensitivity_g = tp_g/(fn_g+tp_g)
    specificity_g = tn_g/(fp_g+tn_g)
    return sensitivity_g, specificity_g

def display_res_gavin_and_best_model(param_grid,pipeline,mvid,filename=None):
    '''
    use model defined by pipeline to fit mvid Dataset
    gridsearchCV determine the parameters given in param_grid
    (optional) save the model in path given in filename
    '''
    classifier = GridSearchCV(estimator=pipeline,
                              param_grid=param_grid)

    print('Start training...')
    classifier.fit(mvid.train.values,mvid.train.labels)
    print('Model Description:\n',classifier.best_estimator_)
    if filename:
        pickle.dump(classifier,open(filename,'wb'))
        print('Saved model to path:',filename)
    sensitivity,specificity,score = run_display_output(classifier,mvid.test)
    print('>>> best model results: sensitivity: {:.{prec}}\tspecificity: {:.{prec}f}\tauc:{}'.\
    format(sensitivity,specificity,score,prec=3))
    return classifier

def display_res_gavin_and_elasticnet(param_grid,pipeline,mvid,filename=None):
    '''
    use model defined by pipeline to fit mvid Dataset
    gridsearchCV determine the parameters given in param_grid
    (optional) save the model in path given in filename
    '''
    classifier = GridSearchCV(estimator=pipeline,
                              param_grid=param_grid)

    print('Start training...')
    classifier.fit(mvid.train.values,mvid.train.labels)
    print('Model Description:\n',classifier.best_estimator_)
    if filename:
        pickle.dump(classifier,open(filename,'wb'))
        print('Saved model to path:',filename)
    sensitivity,specificity,score = run_display_output_elasticnet(classifier,mvid.test)
    print('>>> best model results: sensitivity: {:.{prec}}\tspecificity: {:.{prec}f}\tauc:{}'.\
    format(sensitivity,specificity,score,prec=3))
    return classifier

if __name__=='__main__':

    def load_bonderWestraData(path):
        data = pd.read_csv(bonder_path,sep=',')
        # print(data.head())
        def binarize(row):
            if row > 0:
                return 1
            else:
                return 0
        data['direction'] = data['zscore'].apply(binarize)
        dataset = read_data_set(data)
        return dataset

    bonder_path = 'data/bonder_withzscore.csv'
    westra_allFeat_path = 'data/westra_all_with_zscore.csv.csv'
    westra_bonderFeat_path = 'data/westra_bonderfeat_with_zscore.csv'

    bonder = load_bonderWestraData(bonder_path)
    westra_allFeat = load_bonderWestraData(westra_allFeat_path)
    westra_bonderFeat = load_bonderWestraData(westra_bonderFeat_path)

    print('Bonder dataset loaded.',bonder.train.values.shape)
    print('Westra with all features dataset loaded.',bonder.train.values.shape)
    print('Westra with bonder features dataset loaded.',
           bonder.train.values.shape)

    for eqtm_data in [bonder,westra_allFeat,westra_bonderFeat]:
            print('Dataset:', eqtm_data.train.values.shape)
        # ================LogisticRegression====================================
            n_components = np.arange(2,eqtm_data.train.num_features,10)
            class_weight = ['balanced',{1:4,0:1},{1:2,0:1}]
            param_grid_logr = [{'logr__penalty':['l1','l2'],
                                'logr__C':[1,2,3,4,5],
                                'logr__class_weight':class_weight}]
            # pipeline
            pipeline_logr = Pipeline(steps=[('logr',LogisticRegression())])
            classifier_logr = display_res_gavin_and_best_model(param_grid_logr,
                                             pipeline_logr,
                                             eqtm_data)
        # ====================randomForest======================================
            pipeline_ranfor = Pipeline(steps=[('ranfor',
                                       RandomForestClassifier())])
            n_estimators = [10,50,100]
            param_grid_ranfor = [{'ranfor__n_estimators':n_estimators,
                                  'ranfor__class_weight':class_weight}]
            classifier_ranfor = display_res_gavin_and_best_model(
                                             param_grid_ranfor,
                                             pipeline_ranfor,
                                             eqtm_data)
        # ====================k-nearest neighbor================================
            pipeline_kneighbor = Pipeline(steps=[('kneighbor',
                                          KNeighborsClassifier())])
            # print(pipeline_kneighbor.get_params().keys())
            n_neighbors = range(2,10)
            weights = ['uniform','distance']
            # algorithms = ['auto']
            # leaf_size = [20,30,40,50]
            # p = [1,2,3]
            param_grid_kneighbor = [{'kneighbor__n_neighbors':n_neighbors,
                                  'kneighbor__weights':weights}]#,
                                  # 'kneighbor__algorithms':algorithms}]#,
                                  # 'kneighbor__leaf_size':leaf_size,
                                  # 'kneighbor__p':p}]
            classifier_kneighbor = display_res_gavin_and_best_model(
                                             param_grid_kneighbor,
                                             pipeline_kneighbor,
                                             eqtm_data)
        # ====================gaussian_process==================================
            pipeline_gaussian = Pipeline(steps=[('gaussian',
                                         GaussianProcessClassifier())])
            param_grid_gaussian = [{'gaussian__max_iter_predict':[50]}]
            classifier_gaussian = display_res_gavin_and_best_model(
                                             param_grid_gaussian,
                                             pipeline_gaussian,
                                             eqtm_data)

        # ====================SVM===============================================
            pipeline_svc = Pipeline(steps=[('svc',SVC())])
            C = [0.5,1,2,3]
            kernel = ['linear','poly','rbf','sigmoid','precomputed']
            class_weight = ['balanced',{1:4,0:1},{1:2,0:1}]
            param_grid_svc = [{'svc__kernel':kernel,
                                  'svc__class_weight':class_weight}]
            classifier_svc = display_res_gavin_and_best_model(param_grid_svc,
                                             pipeline_svc,
                                             eqtm_data)
        # ====================Adaboost=================================================
            pipeline_ada = Pipeline(steps=[('ada',SVC())])
            base_estimator = [DecisionTreeClassifier,LogisticRegression]
            param_grid_ada = [{'ada__base_estimator':base_estimator}]
            classifier_ada = display_res_gavin_and_best_model(param_grid_ada,
                                             pipeline_ada,
                                             eqtm_data)
