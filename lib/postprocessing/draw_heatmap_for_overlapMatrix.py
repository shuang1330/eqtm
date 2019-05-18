import os
from lib.read.read_for_cnn import read_individual_image
import numpy as np
from lib.read.read_data import load_data
import matplotlib.pyplot as plt
from matplotlib import cm
from seaborn import clustermap
from scipy.cluster.hierarchy import linkage
from lib.preprocess.selectCellType import (read_allCelltypeFile,
                                           read_bloodCelltypeFile,
                                           read_histoneList)
import pandas as pd
from matplotlib.lines import Line2D


def read_feature_importances(feature_importances_filepath):
    filecontent = open(feature_importances_filepath, 'r').readlines()
    features = [item.strip() for item in filecontent[0].split('\t')]
    importances = [float(item.strip()) for item in filecontent[1].split('\t')]
    return features, importances


def sort_feature_by_importances(feature_importances_filepath,
                                feature_importances_type,
                                feature_type_name):
    if feature_importances_type == "importances":
        features, importances = read_feature_importances(feature_importances_filepath)
    else:
        raise IOError("Needs type, options being either importances or spearcor")
    if feature_type_name == 'gene':
        histones_overlap = []
        histones_importances = []
        for ind,featurename in enumerate(features):
            if featurename.endswith('gene'):
                if featurename not in ['expressionMean', 'expressionVar',
                                       'TssDistance', 'methyMean', 'methyVar']:
                    histones_overlap.append(featurename[:-4])
                    histones_importances.append(importances[ind])
        histone_sorted = [histones_overlap[ind] for ind in np.argsort(np.array(histones_importances))[::-1]]
    elif feature_type_name == 'cpg':
        histones_overlap = []
        histones_importances = []
        for ind,featurename in enumerate(features):
            if not featurename.endswith('gene'):
                if featurename not in ['expressionMean', 'expressionVar',
                                       'TssDistance', 'methyMean', 'methyVar']:
                    histones_overlap.append(featurename)
                    histones_importances.append(importances[ind])
        histone_sorted = [histones_overlap[ind] for ind in np.argsort(np.array(histones_importances))[::-1]]
    return histone_sorted


if __name__=='__main__':
    project_rootdir = '/home/shuang/projects/development_eqtm'
    data_dirpath = os.path.join(project_rootdir, 'data')
    cpgSites_dirpath = os.path.join(data_dirpath, 'cpgSites', 'seperate_cpgFiles')
    geneSites_dirpath = os.path.join(data_dirpath, 'eqtmZscores', 'testCNN', 'geneSites')
    eqtm_path = os.path.join(data_dirpath, 'eqtmZscores',
                             'withExpressionTSSMethyCpgOverlapGene',
                             '2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_withGeneOverlap.txt')
    celltype_dir = os.path.join(project_rootdir, 'data', "features", 'celltype')
    blood_celltype_filepath = os.path.join(celltype_dir, 'roadmap-celltypes-blood-MJ.txt')
    all_celltype_filepath = os.path.join(celltype_dir, 'roadmap-celltypes.txt')
    all_celltype = read_allCelltypeFile(all_celltype_filepath)
    feature_importances_filepath = os.path.join(project_rootdir,
                                                'data',
                                                'features',
                                                'feature_importances',
                                                'featureImportances_allCellType_gt.txt')
    sorted_histone = sort_feature_by_importances(feature_importances_filepath,
                                                 'importances',
                                                 'cpg')
    sorted_celltypes = []
    all_celltype = read_allCelltypeFile(all_celltype_filepath)
    blood_celltype = read_bloodCelltypeFile(blood_celltype_filepath)
    for celltype in read_bloodCelltypeFile(blood_celltype_filepath):
        sorted_celltypes.append(celltype)
    for celltype in [item for item in all_celltype.keys() if item not in blood_celltype]:
        sorted_celltypes.append(celltype)
    print(sorted_celltypes)
    print(sorted_histone)
    # raise NotImplementedError

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
    positive_sum = np.zeros([127, 32])
    negative_sum = np.zeros([127, 32])
    positive_num = 0
    negative_num = 0
    positive_save_dirpath = os.path.join(data_dirpath, "cpgSites", "cpgImgs", "positive")
    negative_save_dirpath = os.path.join(data_dirpath, "cpgSites", "cpgImgs", "negative")
    for item, label in zip(train_values, train_labels):
        print(item[0])
        cpgSite_filepath = os.path.join(cpgSites_dirpath, item[0]+'.txt')
        # geneSite_filepath = os.path.join(geneSites_dirpath, item[1]+'.txt')
        # print(item)
        if label == 0:
            negative_num += 1
            img = read_individual_image(cpgSite_filepath).reshape([127, 32])
            # plt.imshow(img)
            # plt.savefig(os.path.join(negative_save_dirpath, item[0]+'.png'))
            negative_sum += img

        else:
            positive_num += 1
            img = read_individual_image(cpgSite_filepath).reshape([127,32])
            # plt.imshow(img)
            # plt.savefig(os.path.join(positive_save_dirpath, item[0]+'.png'))
            positive_sum += img

    positive_avg = (positive_sum-positive_sum.min())/(positive_sum.max()-positive_sum.min())
    negative_avg = (negative_sum-negative_sum.min())/(negative_sum.max()-negative_sum.min())

    # hierarchical clustering
    all_avg = (positive_sum + negative_sum) / (negative_num + positive_num)
    row_linkage = linkage(all_avg, 'ward')
    col_linkage = linkage(all_avg.T, 'ward')
    all_avg_pandas = pd.DataFrame(data=all_avg, index=sorted_celltypes, columns=sorted_histone)
    cell_type_color = pd.DataFrame(data='azure', index=sorted_celltypes, columns=["is_Blood_Cell"])
    for celltype in blood_celltype:
        cell_type_color[cell_type_color.index==celltype] = 'c'
    g = clustermap(all_avg_pandas,
                   col_linkage=col_linkage,
                   row_linkage=row_linkage,
                   row_colors=cell_type_color)
    cell_type_legends = [Line2D([0], [0], color='azure', lw=4),
                         Line2D([0], [0], color='c', lw=4)]
    # plt.legend(cell_type_legends, ['blood cells', 'non blood cells'])
    plt.savefig("../../output_fig/all_avg_clustermap.png")
    plt.show()
    positive_avg_pandas = pd.DataFrame(data=positive_avg, index=sorted_celltypes, columns=sorted_histone)
    clustermap(positive_avg_pandas, col_linkage=col_linkage, row_linkage=row_linkage, row_colors=cell_type_color)
    plt.savefig("../../output_fig/positive_avg_clustermap.png")
    plt.show()
    negative_avg_pandas = pd.DataFrame(data=negative_avg, index=sorted_celltypes, columns=sorted_histone)
    clustermap(negative_avg_pandas, col_linkage=col_linkage, row_linkage=row_linkage, row_colors=cell_type_color)
    plt.savefig("../../output_fig/negative_avg_clustermap.png")
    plt.show()
    raise NotImplementedError

    # customize the colorbar
    # Make plot with vertical (default) colorbar
    fig = plt.figure()
    ax = fig.add_subplot(121)

    cax = ax.imshow(positive_avg, interpolation='nearest', cmap=cm.coolwarm)
    ax.set_title('Positive correlated eQTMs')

    # # Add colorbar, make sure to specify tick locations to match desired ticklabels
    # cbar = fig.colorbar(cax, ticks=[-1, 0, 1], orientation='horizontal')
    # cbar.ax.set_xticklabels(['Low', 'Medium', 'High'])  # horizontal colorbar

    # Make plot with horizontal colorbar

    ax2 = fig.add_subplot(122)
    cax2 = ax2.imshow(negative_avg, interpolation='nearest', cmap=cm.coolwarm)
    ax2.set_title('Negative correlated eQTMs')

    fig.subplots_adjust(right=0.8)
    # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.75])
    fig.colorbar(cax2, ticks=[0, 0.5, 1])


    # cbar = fig.colorbar(cax, )
    # cbar.ax.set_yticklabels(['0', '0.5', '1'])  # vertically oriented colorbar
    plt.savefig("/home/shuang/projects/development_eqtm/output_fig/heatmaps.png")

    plt.show()

    # plt.imshow(positive_avg)
    # plt.show()
    #
    # plt.imshow(negative_avg)
    # plt.show()
    #
    # plt.imshow(positive_avg-negative_avg)
    # plt.show()
    #
    #
    # test_filepath = os.path.join(project_rootdir,
    #                              'data/cpgSites/seperate_cpgFiles',
    #                              'cg05661533.txt')
    # test_image = open(test_filepath, 'r').readlines()
    # print(len(test_image))
    # print(len(test_image[0].split('\t')))
    # test_array = [[int(num) for num in item.split('\t')] for item in test_image]
    # print(np.array(test_array))
    # print(np.array(test_array))
    # plt.imshow(test_array)
    # plt.show()
