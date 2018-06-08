import __init__path
import numpy as np
import pandas as pd
import pickle
import os
from lib.read.read_data import load_data
from lib.read.dataset_class import dataset
from lib.model.modelBuild import convModelObject
from math import copysign

from sklearn.pipeline import Pipeline
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier

from sklearn import metrics
from sklearn.metrics import confusion_matrix


def examineModel(model_name,trainPath,testPath,savefilename=None,
                 keep=[],exclude=[],display=True):
    res = {'sensitivity':[],'specificity':[],'auc':[]}
    res2 = {'sensitivity':[],'specificity':[],'auc':[]}
    for i in range(4):
        data = load_data(trainPath,keep=keep,exclude=exclude)
        features_to_use = [col for col in data.train.values.columns if col not in keep]
        train_data = dataset(data.train.values[features_to_use],
                             data.train.labels.astype('int8'))
        valid_data = dataset(data.test.values[features_to_use],
                             data.test.labels.astype('int8'))

        model = choose_model(model_name,len(features_to_use),train.shape[0])
        print(model.pipeline, model.param_grid)
        model.savefilename = savefilename

        sensitivity,specificity,auc = \
        model.train_and_display_output(train_data=train_data,
                                       valid_data=valid_data,
                                       returnRes=True,
                                       display=display)
        res['sensitivity'].append(sensitivity)
        res['specificity'].append(specificity)
        res['auc'].append(auc)
        if trainPath != testPath:
            test_data = load_data(testPath,
                                  keep=keep,
                                  exclude=exclude,
                                  test_size=1).test
            sensitivity2,specificity2,auc2 = \
            model.calculate_output(model.classifier,test_data)
            print('>>> test on another dataset: sensitivity: '+
            '{:.{prec}}\tspecificity: {:.{prec}f}\tauc:{}'.\
                  format(sensitivity2,specificity2,auc2,prec=3))
            res2['sensitivity'].append(sensitivity2)
            res2['specificity'].append(specificity2)
            res2['auc'].append(auc2)
    print('Sensitivity:',
    np.mean(np.array(res['sensitivity'])),
    np.std(np.array(res['sensitivity'])))
    print('Specificity:',
    np.mean(np.array(res['specificity'])),
    np.std(np.array(res['specificity'])))
    print('AUC:',
    np.mean(np.array(res['auc'])),
    np.std(np.array(res['auc'])))
    if trainPath != testPath:
        print('Sensitivity:',
        np.mean(np.array(res2['sensitivity'])),
        np.std(np.array(res2['sensitivity'])))
        print('Specificity:',
        np.mean(np.array(res2['specificity'])),
        np.std(np.array(res2['specificity'])))
        print('AUC:',
        np.mean(np.array(res2['auc'])),
        np.std(np.array(res2['auc'])))

def compare_Tss_meanVar_exclude(filepath,testpath):
    exclude1 = ['cpgName','cpgName_split','TSS_Distance','methyMean','methyVar']
    exclude2 = ['cpgName','cpgName_split','methyMean','methyVar']
    exclude3 = ['cpgName','cpgName_split','TSS_Distance','methyVar']
    exclude4 = ['cpgName','cpgName_split','TSS_Distance','methyMean']
    exclude5 = ['cpgName','cpgName_split','methyMean']
    exclude6 = ['cpgName','cpgName_split','methyVar']
    exclude7 = ['cpgName','cpgName_split','TSS_Distance']
    exclude8 = ['cpgName','cpgName_split']

    # only check dataset with TSS and dataset with TSS & meanVariance
    for exclude in [exclude2,exclude8]:
        print('Not using features:',exclude)
        modelPerformance(filepath,testpath,exclude=exclude)

def choose_model(model_name,num_features,num_samples):
    if model_name=='ranfor':
        pipeline = Pipeline(steps=[('pca',PCA()),('ranfor',
                                           RandomForestClassifier())])
        n_estimators = [10,50,100]
        class_weight = ['balanced',{1:4,0:1},{1:2,0:1},
                        {1:1,0:2},{1:1,0:4}]
        n_components = np.arange(2,num_features,int(num_features/5))
        param_grid = [{'pca__n_components':n_components,
                       'ranfor__n_estimators':n_estimators,
                       'ranfor__class_weight':class_weight}]
    elif model_name=='knn':
        pipeline = Pipeline(steps=[('kneighbor',
                                    KNeighborsClassifier())])
        # train the model!
        n_neighbors = np.arange(5,
                      int(num_samples/2)),
                      int((num_samples/2-5)/5))
        weights = ['uniform','distance']
        param_grid = [{'kneighbor__n_neighbors':n_neighbors,
                       'kneighbor__weights':weights}]
    model = convModelObject(name='model_name',
                      pipeline=pipeline,
                      param_grid=param_grid)
    return model

if __name__=='__main__':

    PROJECT_DIR = '/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm'
    data_folder = os.path.join(PROJECT_DIR,'data','eqtmZscores','final')

    model_name = 'ranfor'
    exclude = ['SNPName','SNPChr','PValue','SNPChrPos','ProbeName','ProbeChr',
               'ProbeCenterChrPos','CisTrans', 'SNPType', 'AlleleAssessed',
               'DatasetsWhereSNPProbePairIsAvailableAndPassesQC',
               'DatasetsZScores', 'DatasetsNrSamples',
               'IncludedDatasetsMeanProbeExpression',
               'IncludedDatasetsProbeExpressionVariance', 'HGNCName',
               'IncludedDatasetsCorrelationCoefficient', 'Meta-Beta (SE)',
               'Beta (SE)','FoldChange', 'FDR','checkChr','SNPName_ProbeName']
    keep = ['OverallZScore']
    for filename in os.listdir(data_folder):
        train_path = os.path.join(data_folder,filename)
        test_path = train_path
        examineModel(model_name,train_path,test_path,
                     keep=keep,exclude=exclude,display=True)
