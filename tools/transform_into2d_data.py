import __init__path
import os
import numpy as np
import pandas as pd

from lib.preprocess.selectCellType import read_allCelltypeFile, read_bloodCelltypeFile, read_histoneList


def read_feature_importances(feature_importances_filepath):
    filecontent = open(feature_importances_filepath,'r').readlines()
    features = [item.strip() for item in filecontent[0].split('\t')]
    importances = [float(item.strip()) for item in filecontent[1].split('\t')]
    return features, importances


def sort_feature_by_importances(feature_importances_filepath,feature_type_name):
    features, importances = read_feature_importances(feature_importances_filepath)
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


def rearrange_cpgs(overlapMatrix_filepath,
                   type_name,
                    blood_celltype_filepath,
                    all_celltype_filepath,
                    histonelist_filepath,
                    feature_importances_filepath,
                    save_dirpath):
    blood_celltype = read_bloodCelltypeFile(blood_celltype_filepath)
    all_celltype = read_allCelltypeFile(all_celltype_filepath)
    histones = read_histoneList(histonelist_filepath)
    sorted_histone = sort_feature_by_importances(feature_importances_filepath,type_name)
    sorted_features = []
    for celltype in blood_celltype:
        for histone in sorted_histone:
            feature_name = '{}-{}'.format(celltype, histone)
            sorted_features.append(feature_name)
    for celltype in [item for item in all_celltype.keys() if item not in blood_celltype]:
        for histone in sorted_histone:
            feature_name = '{}-{}'.format(celltype, histone)
            sorted_features.append(feature_name)

    overlapMatrix = pd.read_csv(overlapMatrix_filepath, sep=',', index_col=0)
    print('Loaded overlapMatrix.', flush=True)
    overlapMatrix['E005-H4K12ac'] = 0
    overlapMatrix['E008-H3T11ph'] = 0
    overlapMatrix['E017-H2AK9ac'] = 0
    num_histones = len(histones)
    num_all_celltypes = len(all_celltype.keys())

    content = overlapMatrix[sorted_features].values
    for ind, cpg_name in enumerate(overlapMatrix.index):
        print(cpg_name, flush=True)
        cpg_save_filepath = os.path.join(save_dirpath, cpg_name + '.txt')
        cpg_file = open(cpg_save_filepath, 'w')
        for i in range(0, content.shape[1], num_histones):
            line_to_write = content[ind, i:i + num_histones]
            cpg_file.write('\t'.join([str(item) for item in line_to_write]))
            cpg_file.write('\n')
        cpg_file.close()


if __name__=='__main__':
    # locally
    # project_rootdir = '/home/shuang/projects/development_eqtm'
    # calculon
    project_rootdir = '/groups/umcg-gcc/tmp04/umcg-sli/development_eqtm'
    feature_importances_filepath = os.path.join(project_rootdir,
                                                'data',
                                                'features',
                                                'feature_importances',
                                                'featureImportances_allCellType_gt.txt')
    celltype_dir = os.path.join(project_rootdir, 'data', "features", 'celltype')
    blood_celltype_filepath = os.path.join(celltype_dir, 'roadmap-celltypes-blood-MJ.txt')
    all_celltype_filepath = os.path.join(celltype_dir, 'roadmap-celltypes.txt')
    histonelist_filepath = os.path.join(project_rootdir, 'data', 'features', 'histonetype',
                                        '2017-12-09-eQTLsFDR-et0.0-flipped_histoneList.txt')

    # cpgSites
    cpg_overlapMatrix_filepath = os.path.join(project_rootdir,
                                              'data',
                                              'eqtmZscores',
                                              'complete_overlapMatrix',
                                              '2017-12-09-eQTLsFDR-gt0.0-flipped_overlapMatrixcomplete.txt')
    cpg_save_dirpath = os.path.join(project_rootdir,
                                    'data', 'cpgSites',
                                    'seperate_cpgFiles')

    gene_overlapMatrix_filepath = os.path.join(project_rootdir, 'data',
                                               'features',
                                               'geneOverlap',
                                               'gene_StartEndSite_overlapMatrix_complete.txt')
    gene_save_dirpath = os.path.join(project_rootdir,
                                     'data', 'features',
                                     'geneOverlap', 'seperate_geneFiles')

    rearrange_cpgs(cpg_overlapMatrix_filepath,
                   'cpg',
                   blood_celltype_filepath,
                   all_celltype_filepath,
                   histonelist_filepath,
                   feature_importances_filepath,
                   cpg_save_dirpath)

    rearrange_cpgs(gene_overlapMatrix_filepath,
                   'gene',
                   blood_celltype_filepath,
                   all_celltype_filepath,
                   histonelist_filepath,
                   feature_importances_filepath,
                   gene_save_dirpath)