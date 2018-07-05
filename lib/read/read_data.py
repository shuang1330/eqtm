import numpy as np
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from lib.read.dataset_class import Datasets, dataset
import pandas as pd
import os


def check_null(datatable):
    """
    check and print columns with null values
    :param datatable: pandas dataframe
    :return: print the columns with null values
    """
    print('Here are columns with NaN values:')
    for col in datatable.columns:
        null = datatable[col].isnull().values.sum()
        if null > 0:
            print(null, col, datatable[col].dtypes)


def readDataWithRawScore(data_table,
                         label_name='direction',
                         keep=[],
                         test_size=0.25,
                         normalization=False):
    """
    convert a pandas dataframe data table into Datasets(dataset,dataset)
    :param data_table: pandas dataframe
    :param label_name: string, the name in dataframe for the label, default 'direction'
    :param keep: list of string, features will not be normalized, default None
    :param test_size: float/int, default 0.25
    :param normalization: true/false, default True
    :return: Datasets(train=dataset(x:pandas dataFrame,y:pandas dataFrame),
                test=dataset(x:pandas dataFrame,y:pandas dataFrame))
    """
    train, test = train_test_split(data_table,test_size=test_size)
    copy_train = train.copy().reset_index()
    copy_test = test.copy().reset_index()
    keep.append(label_name)
    if normalization:
        minMaxScaler = preprocessing.MinMaxScaler()
        numericalCols = [col for col in
                         train.select_dtypes(exclude=[np.object]).columns
                         if col not in keep]
        train_feature_copy = copy_train.loc[:,numericalCols]
        copy_train.loc[:,numericalCols] = \
        minMaxScaler.fit_transform(train_feature_copy)
        if copy_test.shape[0]>0:
            test_feature_copy = copy_test.loc[:,numericalCols]
            copy_test.loc[:,numericalCols] = \
            minMaxScaler.fit_transform(test_feature_copy)
        print('Data Normalized.')
    train_x = copy_train[[col for col in train.columns if col not in [label_name]]]
    test_x = copy_test[[col for col in test.columns if col not in [label_name]]]
    train_y = copy_train[label_name]
    test_y = copy_test[label_name]
    return Datasets(train=dataset(train_x,train_y),
                    test=dataset(test_x,test_y))


def load_data(path, sep=',', label_name='direction',
              exclude=[], keep=[], test_size=0.25,
              normalization=False):
    '''
    read from path and convert the csv file into dataset objects
    INPUT:
        path: string, path to data, a csv file
        sep: string, seperator, default ','
        exclude: list, features will be ignored, default []
        keep: list, features will be untouched in dataset(train)
        test_size: float/int, the train_test_split ratio, default 0.25
    OUTPUT:
        dataset(train(values,labels), test(values, labels)): self-defined object
        values, labels: pandas dataframe
    '''
    data = pd.read_csv(path, sep=sep, index_col=0)

    def binarize(row):
        if row > 0:
            return 1
        else:
            return 0
    # add a new column named direction, which binarize the zscore column
    # positive zscore gets 1 while negative gets 0
    data[label_name] = data['OverallZScore'].apply(binarize)
    print('Raw data loaded with shape:',data.shape)
    features = [feat for feat in data.columns if feat not in exclude]
    dataset = readDataWithRawScore(data[features],
                                   keep=keep,
                                   label_name=label_name,
                                   test_size=test_size,
                                   normalization=normalization)
    print('Check the null values:')
    check_null(dataset.train.values)
    print('Check the loaded dataset: \n train with shape',
          dataset.train.values.shape)
    print('test with shape:', dataset.test.values.shape)
    return dataset


def read_eqtm_names(eqtm_path):
    exclude = ['SNPChr', 'PValue', 'SNPChrPos', 'ProbeChr',
                   'ProbeCenterChrPos', 'CisTrans', 'SNPType', 'AlleleAssessed',
                   'DatasetsWhereSNPProbePairIsAvailableAndPassesQC',
                   'DatasetsZScores', 'DatasetsNrSamples',
                   'IncludedDatasetsMeanProbeExpression',
                   'IncludedDatasetsProbeExpressionVariance', 'HGNCName',
                   'IncludedDatasetsCorrelationCoefficient', 'Meta-Beta (SE)',
                   'Beta (SE)', 'FoldChange', 'FDR', 'checkChr',
                   'TssSite', 'chr', 'SNPName_ProbeName']
    keep = ['SNPName', 'ProbeName', 'OverallZScore']
    data = load_data(eqtm_path, keep=keep, exclude=exclude)
    return [data.train.values,
                data.train.labels,
                data.test.values,
                data.test.labels]


if __name__ == '__main__':
    project_rootdir = '/home/shuang/projects/development_eqtm'
    test_filepath = os.path.join(project_rootdir, 'data', 'eqtmZscores',
                             'withExpressionTSSMethyCpgOverlapGene',
                             '2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_withGeneOverlap.txt')
    train_data, train_labels, test_data, test_labels = read_eqtm_names(test_filepath)
    print(train_data.columns)
