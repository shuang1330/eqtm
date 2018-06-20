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
                if featurename not in ['expressionMean','expressionVar',
                                       'TssDistance','methyMean','methyVar']:
                    histones_overlap.append(featurename)
                    histones_importances.append(importances[ind])
        histone_sorted = [histones_overlap[ind] for ind in np.argsort(np.array(histones_importances))[::-1]]
    elif feature_type_name == 'cpg':
        histones_overlap = []
        histones_importances = []
        for ind,featurename in enumerate(features):
            if not featurename.endswith('gene'):
                if featurename not in ['expressionMean','expressionVar',
                                       'TssDistance','methyMean','methyVar']:
                    histones_overlap.append(featurename)
                    histones_importances.append(importances[ind])
        histone_sorted = [histones_overlap[ind] for ind in np.argsort(np.array(histones_importances))[::-1]]
    return histone_sorted


def transform_data(rearranged_matrix_filepath,
                   overlapMatrix_filepath,
                   blood_celltype_filepath,
                   all_celltype_filepath,
                   feature_importances_filepath,
                   feature_type_name,
                   histonelist_filepath
                   ):
    # rearranged_matrix_filepath = os.path.join(transformed_input_dirpth,'rearranged_matrix.txt')
    overlapMatrix = pd.read_csv(overlapMatrix_filepath, index_col=0)
    transposed_overlapMatrix = pd.DataFrame.transpose(overlapMatrix)
    print(transposed_overlapMatrix.head())
    rearranged_matrix = open(rearranged_matrix_filepath,'w')
    rearranged_matrix.write('\t'.join(transposed_overlapMatrix.columns))
    rearranged_matrix.write('\n')
    blood_celltypes = read_bloodCelltypeFile(blood_celltype_filepath)
    all_celltypes_info = read_allCelltypeFile(all_celltype_filepath)
    histone_sorted_by_importances = sort_feature_by_importances(feature_importances_filepath,feature_type_name)
    for celltype in blood_celltypes:
        for histone in histone_sorted_by_importances:
            index = '-'.join([celltype, histone])
            rearranged_matrix.write(index)
            rearranged_matrix.write('\t')
            try:
                rearranged_matrix.write('\t'.join([str(item) for item in transposed_overlapMatrix.loc[index].values]))
            except:
                rearranged_matrix.write('\t'.join(['0' for item in np.arange(transposed_overlapMatrix.shape[1])]))
            rearranged_matrix.write('\n')
    other_celltype = [col for col in all_celltypes_info.keys() if col not in blood_celltypes]
    for celltype in other_celltype:
        for histone in histone_sorted_by_importances:
            index = '-'.join([celltype, histone])
            rearranged_matrix.write(index)
            rearranged_matrix.write('\t')
            try:
                rearranged_matrix.write('\t'.join([str(item) for item in transposed_overlapMatrix.loc[index].values]))
            except:
                rearranged_matrix.write('\t'.join(['0' for item in np.arange(transposed_overlapMatrix.shape[1])]))
            rearranged_matrix.write('\n')
    rearranged_matrix.close()
    print('Wrote the re-arranged complete matrix to path:',rearranged_matrix_filepath)


def write_inidividualFiles(histonelist_filepath,
                           all_celltype_filepath,
                           rearranged_matrix_filepath,
                           transformed_input_dirpth):

    # write inividual files
    histone_list = read_histoneList(histonelist_filepath)
    num_histone = len(histone_list)
    all_celltypes_info = read_allCelltypeFile(all_celltype_filepath)
    num_allcelltypes = len(all_celltypes_info.keys())
    final_gt_rearranged_matrix = pd.read_csv(rearranged_matrix_filepath, sep='\t')
    print('Read matrix file.')
    for col in final_gt_rearranged_matrix.columns:
        file_to_write = open(os.path.join(transformed_input_dirpth,col+'.txt'),'w')
        for i in range(0,final_gt_rearranged_matrix.shape[0],num_allcelltypes):
            content_to_write = [str(row) for row in final_gt_rearranged_matrix[col][i:(i+num_allcelltypes)].values]
            file_to_write.write('\t'.join(content_to_write))
            file_to_write.write('\n')
        file_to_write.close()
    print('Finished writing all individual files to path:',transformed_input_dirpth)


if __name__=='__main__':
    project_rootdir = '/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm'
    feature_importances_filepath = os.path.join(project_rootdir,
                                                'data',
                                                'features',
                                                'feature_importances',
                                                'featureImportances_allCellType_gt.txt')
    celltype_dir = os.path.join(project_rootdir, 'data', 'features', 'celltype')
    blood_celltype_filepath = os.path.join(celltype_dir, 'roadmap-celltypes-blood-MJ.txt')
    all_celltype_filepath = os.path.join(celltype_dir, 'roadmap-celltypes.txt')
    histonelist_filepath = os.path.join(project_rootdir, 'data', 'features', 'histonetype',
                                        '2017-12-09-eQTLsFDR-et0.0-flipped_histoneList.txt')

    # # cpgSites
    # cpg_rearranged_matrix_filepath = os.path.join(project_rootdir, 'data',
    #                                               'eqtmZscores',
    #                                               'testCNN',
    #                                               'cpgSites',
    #                                               'rearranged_matrix.txt')
    # cpg_overlapMatrix_filepath = os.path.join(project_rootdir, 'data', 'eqtmZscores',
    #                                           'complete_overlapMatrix',
    #                                           '2017-12-09-eQTLsFDR-gt0.0-flipped_overlapMatrixcomplete.txt')
    # cpg_transformed_input_dirpth = os.path.join(project_rootdir, 'data', 'eqtmZscores', 'testCNN', 'cpgSites')
    # transform_data(cpg_rearranged_matrix_filepath,
    #                cpg_overlapMatrix_filepath,
    #                blood_celltype_filepath,
    #                all_celltype_filepath,
    #                feature_importances_filepath,
    #                'cpg',
    #                histonelist_filepath
    #                )
    # write_inidividualFiles(histonelist_filepath,
    #                        all_celltype_filepath,
    #                        cpg_rearranged_matrix_filepath,
    #                        cpg_transformed_input_dirpth)

    # geneSites
    gene_rearranged_matrix_filepath = os.path.join(project_rootdir, 'data', 'eqtmZscores',
                                                   'testCNN', 'geneSites',
                                                   'gene_rearranged_overlapMatrix.txt')
    gene_overlapMatrix_filepath = os.path.join(project_rootdir, 'data', 'features',
                                               'geneOverlap',
                                               'gene_StartEndSite_overlapMatrix.txt')
    transformed_input_dirpth = os.path.join(project_rootdir, 'data', 'eqtmZscores', 'testCNN')
    gene_overlapMatrix_save_dirpath = os.path.join(transformed_input_dirpth, 'geneSites')
    gene_transformed_input_dirpth = os.path.join(project_rootdir, 'data',
                                                 'eqtmZscores', 'testCNN', 'geneSites')
    transform_data(gene_rearranged_matrix_filepath,
                   gene_overlapMatrix_filepath,
                   blood_celltype_filepath,
                   all_celltype_filepath,
                   feature_importances_filepath,
                   'gene',
                   histonelist_filepath
                   )
    write_inidividualFiles(histonelist_filepath,
                           all_celltype_filepath,
                           gene_rearranged_matrix_filepath,
                           gene_transformed_input_dirpth)
