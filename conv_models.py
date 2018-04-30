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
from sklearn.ensemble import RandomForestClassifier
# classifiers
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import ElasticNet
from sklearn.svm import SVC
# finetuning
from sklearn.model_selection import GridSearchCV
# validation
from sklearn import metrics
from sklearn.metrics import confusion_matrix

# import matplotlib.pyplot as plt
# TODO ElasticNet


def read_data_set(data_table,test_size=0.25):
    '''
    convert a pandas dataframe data table into Datasets(dataset,dataset)
    '''
    train, test = train_test_split(data_table,test_size=0.25)
    train_x = train[[col for col in train.columns
                     if col not in ['zscore','direction']]]
    features = train_x.columns
    train_x = np.array(train_x)
    test_x = np.array(test[[col for col in train.columns
                      if col not in ['zscore','direction']]])
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

def run_display_output_elasticnet(classifier,test,DRAW=False):
    '''
    get confusion matrix and auc score for test dataset
    (optional) draw roc curve
    '''
    pred = classifier.predict(test.values)
    pass
    # tn, fp, fn, tp = confusion_matrix(test.labels,pred).ravel()#confusion matrix
    # print(tn,fp,fn,tp)
    # sensitivity = tp/(fn+tp)
    # specificity = tn/(fp+tn)
    # prods = classifier.predict_proba(test.values)[:,1]
    # fpr, tpr, _ = metrics.roc_curve(test.labels, prods)
    # score = metrics.auc(fpr,tpr) #auc score
    # if DRAW:
    #     draw_roc_curve(fpr,tpr,score)
    #
    # return sensitivity, specificity, score

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

    features = ['X._H2A.Z_imputed_gappedPeaks','X._H2AK5ac_imputed_gappedPeaks',\
    'X._H2AK9ac_imputed_gappedPeaks','X._H2BK120ac_imputed_gappedPeaks',\
    'X._H2BK12ac_imputed_gappedPeaks','X._H2BK5ac_imputed_gappedPeaks',\
    'X._H3K14ac_imputed_gappedPeaks','X._H3K18ac_imputed_gappedPeaks',\
    'X._H3K23ac_imputed_gappedPeaks','X._H3K23me2_imputed_gappedPeaks',\
    'X._H3K27ac_imputed_gappedPeaks','X._H3K4ac_imputed_gappedPeaks',\
    'X._H3K4me1_imputed_gappedPeaks','X._H3K4me2_imputed_gappedPeaks',\
    'X._H3K4me3_imputed_gappedPeaks','X._H3K79me1_imputed_gappedPeaks',\
    'X._H3K79me2_imputed_gappedPeaks','X._H3K9ac_imputed_gappedPeaks',\
    'X._H3T11ph_imputed_gappedPeaks','X._H4K12ac_imputed_gappedPeaks',\
    'X._H4K5ac_imputed_gappedPeaks','X._H4K8ac_imputed_gappedPeaks',\
    'X._H4K91ac_imputed_gappedPeaks','X._H3K36me3_imputed_gappedPeaks',\
    'X._H3K9me3_imputed_gappedPeaks','X._H4K20me1_imputed_gappedPeaks',\
    'X._H3K27me3_imputed_gappedPeaks','X._H3K56ac_imputed_gappedPeaks',\
    'X._H2BK15ac_imputed_gappedPeaks','X._H2BK20ac_imputed_gappedPeaks',\
    'X._H3K9me1_imputed_gappedPeaks','X._H3K27me3_gapedPeaks',\
    'X._H3K36me3_gapedPeaks','X._H3K9me3_gapedPeaks',\
    'X._H3K4me1_gapedPeaks','X._H3K4me3_gapedPeaks','X._H3K27ac_gapedPeaks']
    # features.append('zscore') # output

    def read_useful_features(path,features,SAVE=False):
        all_data = pd.read_csv(path,sep='\t')
        def binarize(row):
            if row > 0:
                return 1
            else:
                return 0
        all_data['direction'] = all_data['zscore'].apply(binarize)
        features.append('direction')
        if SAVE:
            all_data[features].to_csv('mj_data/only_ratio_as_features.csv')
            print('saved data to mj_data/only_ratio_as_features.csv')
        eqtm_data = read_data_set(all_data[features])
        return eqtm_data
    eqtm_data = read_useful_features('mj_data/Anno_Value_Direction.csv',
                                     features,
                                     SAVE=True)
    # mvid,train_gavin, test_gavin = read_data_set(data,BENCHMARK=True)
    # print(data.head())
    # raise NotImplementedError # check the dataset loaded
    print('Dataset loaded.',eqtm_data.train.labels.shape)
    print(eqtm_data.train.values.shape)

    # # simple_model = LogisticRegression()
    # simple_model = RandomForestClassifier()
    # simple_model.fit(eqtm_data.train.values,eqtm_data.train.labels)
    # pred = simple_model.predict(eqtm_data.test.values)
    # tn, fp, fn, tp = confusion_matrix(eqtm_data.test.labels,pred).ravel()
    # print(tn, fp, fn, tp)


# ================model selection==========================================
    # PCA + LogisticRegression
    # Parameters
    n_components = np.arange(2,eqtm_data.train.num_features,10)
    class_weight = ['balanced',{1:4,0:1},{1:2,0:1}]
    param_grid_logr = [{'pca__n_components':n_components,
                   'logr__penalty':['l1','l2'],
                   'logr__C':[1,2,3,4,5],
                   'logr__class_weight':class_weight}]
    # pipeline
    pipeline_logr = Pipeline(steps=[('pca',PCA()),
                               ('logr',LogisticRegression())])
    # save model
    # filename = os.path.join('model','pca_logr.sav')
    # display results
    classifier_logr = display_res_gavin_and_best_model(param_grid_logr,
                                     pipeline_logr,
                                     eqtm_data)

    # # PCA + ElasticNet
    # pipeline_eln = Pipeline(steps=[('pca',PCA()),
    #                                ('eln',ElasticNet())])
    # alpha = [0.2,0.4,0.6,1]
    # l1_ratio = [0.2,0.4,0.6,0.8]
    # normalize=[True]
    # param_grid_eln = [{'pca__n_components':n_components,
    #                    'eln__alpha':alpha,
    #                    'eln__l1_ratio':l1_ratio,
    #                    'eln__normalize':normalize}]
    # display_res_gavin_and_elasticnet(param_grid_eln,
    #                                  pipeline_eln,
    #                                  mvid)

    # # PCA + RandomForest
    pipeline_ranfor = Pipeline(steps=[('pca',PCA()),
                                      ('ranfor',RandomForestClassifier())])
    n_estimators = [10,50,100]
    param_grid_ranfor = [{'pca__n_components':n_components,
                          'ranfor__n_estimators':n_estimators,
                          'ranfor__class_weight':class_weight}]
    # filename = os.path.join('model','pca_ranfor.sav')
    classifier_ranfor = display_res_gavin_and_best_model(param_grid_ranfor,
                                     pipeline_ranfor,
                                     eqtm_data)

    # # display gavin results
    # sensitivity_g,specificity_g = read_gavin(myo5b_gavin,myo5b.train.labels)
    # print('>>> gavin model results: sensitivity: {:.{prec}}\tspecificity: {:.{prec}f}'.\
    # format(sensitivity_g,specificity_g,prec=3))
# ======================================================================
