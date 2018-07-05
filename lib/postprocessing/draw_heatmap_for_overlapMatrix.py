import os
from lib.read.read_for_cnn import read_individual_image
import numpy as np

from lib.read.read_data import load_data

import matplotlib.pyplot as plt

if __name__=='__main__':
    project_rootdir = '/home/shuang/projects/development_eqtm'
    data_dirpath = os.path.join(project_rootdir, 'data')
    cpgSites_dirpath = os.path.join(data_dirpath, 'cpgSites', 'seperate_cpgFiles')
    geneSites_dirpath = os.path.join(data_dirpath, 'eqtmZscores', 'testCNN', 'geneSites')
    eqtm_path = os.path.join(data_dirpath, 'eqtmZscores',
                             'withExpressionTSSMethyCpgOverlapGene',
                             '2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_withGeneOverlap.txt')

    def read_cnn_names(eqtm_path):
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
        return [data.train.values[['SNPName', 'ProbeName']].values,
                data.train.labels,
                data.test.values[['SNPName', 'ProbeName']].values,
                data.test.labels]

    [train_values, train_labels, _, _] = read_cnn_names(eqtm_path)
    positive_sum = np.zeros([127, 32, 1])
    negative_sum = np.zeros([127, 32, 1])
    positive_num = 0
    negative_num = 0
    for item, label in zip(train_values, train_labels):
        cpgSite_filepath = os.path.join(cpgSites_dirpath, item[0]+'.txt')
        # geneSite_filepath = os.path.join(geneSites_dirpath, item[1]+'.txt')
        # print(item)
        if label == 0:
            positive_num += 1
            negative_sum += read_individual_image(cpgSite_filepath)
        else:
            negative_num += 1
            positive_sum += read_individual_image(cpgSite_filepath)

    positive_avg = positive_sum/positive_num
    negative_avg = negative_sum/negative_num

    plt.imshow(positive_avg)
    plt.show()

    plt.imshow(negative_avg)
    plt.show()

    plt.imshow(positive_avg-negative_avg)
    plt.show()


    test_filepath = os.path.join(project_rootdir,
                                 'data/cpgSites/seperate_cpgFiles',
                                 'cg05661533.txt')
    test_image = open(test_filepath, 'r').readlines()
    print(len(test_image))
    print(len(test_image[0].split('\t')))
    test_array = [[int(num) for num in item.split('\t')] for item in test_image]
    print(np.array(test_array))
    print(np.array(test_array))
    plt.imshow(test_array)
    plt.show()