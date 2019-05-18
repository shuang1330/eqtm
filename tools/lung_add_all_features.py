import __init__path
import os
import numpy as np
import pandas as pd
from lib.preprocess.add_emptyFeature_nonOverlappingCpg_toOverlapMatrix import (
    add_emptyFeature_nonOverlappingCpg)
import argparse


def featureList(feature_folder):
    feature_list = []
    for name in os.listdir(feature_folder):
        if name.endswith('gz'):
            feature_name = '.'.join(name.split('.')[:-5])
            feature_list.append(feature_name)
    print('featurelist:',len(feature_list))
    return feature_list


def histoneList(feature_list):
    histone_list = []
    for col in feature_list:
        histone = col.split('-')[1]
        if histone not in histone_list:
            histone_list.append(histone)
    print('Histone list:',len(histone_list))
    return histone_list


def cellList(feature_list):
    cell_list = []
    for col in feature_list:
        cell = col.split('-')[0]
        if cell not in cell_list:
            cell_list.append(cell)
    print('Cell list:',len(cell_list))
    return cell_list


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dirname", dest="input_dirname", type=str)
    parser.add_argument("--cpg_filename", dest="cpg_filename", type=str)
    parser.add_argument("--original_eqtm_filepath", dest="original_eqtm_filepath", type=str)
    parser.add_argument("--name", dest="name", type=str)
    args = parser.parse_args()

    # input_dirname = 'cpg_largerthan0.05'
    # cpg_filename = 'eqtm_FDR_larger_than_0.05_bedtoolsFormat.txt'
    # input_dirname = 'cpg'
    # cpg_filename = 'UCEC_bedtoolsFormat.txt'
    # original_eqtm_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/tcga/cis_results/TCGA-UCEC-mval/eQTLSNPsFDR0.05.txt"

    cpg_filename = args.cpg_filename
    input_dirname = args.input_dirname
    original_eqtm_filepath = args.original_eqtm_filepath
    name = args.name

    # all the directories
    PROJECT_DIR = '/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm'
    feature_folder = '/groups/umcg-bios/tmp03/projects/2018-methylation' \
                     '/input/features/Roadmap/consolidatedImputedGappedPeak'

    data_folder = os.path.join(PROJECT_DIR, 'data', 'temp',
                               'intersect_features', input_dirname)
    cpg_path = os.path.join(PROJECT_DIR,
                            'data',
                            'eqtmZscores',
                            'allCpgs',
                            cpg_filename)
    save_path = os.path.join(PROJECT_DIR,
                             'data',
                             'eqtmZscores',
                             'allCpgs',
                             name+'_overlapMatrix.txt')
    complete_overlapMatrix_savepath = os.path.join(PROJECT_DIR,
                                                   "data",
                                                   "eqtmZscores",
                                                   "allCpgs",
                                                   name + '_overlapMatrix_complete.txt')
    ratio_table_savepath = os.path.join(PROJECT_DIR,
                                        'data',
                                        'eqtmZscores',
                                        'allCpgs',
                                        name+'_overlapRatio.txt')

    # create the result overlap_matrix
    feature_list = featureList(feature_folder)
    histone_list = histoneList(feature_list)
    cell_type_list = cellList(feature_list)
    print("Cell type number: ", len(cell_type_list))
    raise NotImplementedError

    cpg_colnames = ["chr", "startSite", "endSite", "SNP", 'allele', 'minor_allele', 'traits']
    cpg = pd.read_csv(cpg_path,
                      sep='\t',
                      names=cpg_colnames)
    print(cpg.head())
    overlapTable = pd.DataFrame(data=0,
                                index=cpg['SNP'].unique(),
                                columns=feature_list,
                                dtype=np.int8)

    print('With shape: ', overlapTable.shape)
    names = ['chr_feat', 'str_feat', 'end_feat', 'type',
             'chr_cpg', 'str_cpg', 'end_cpg', 'cpg_name', 'allele', 'minor_allele', 'traits', 'pubids']

    # start processing
    print('Process features in folder: ', data_folder)
    for filename in os.listdir(data_folder):
        print("Processing feature:", filename, flush=True)
        filepath = os.path.join(data_folder, filename)
        if os.stat(filepath).st_size == 0:
            print('Empty file:', filename, flush=True)
            continue
        else:
            bedtoolsRes = pd.read_csv(filepath, sep='\t', names=names)
            # print(bedtoolsRes.head())
            bedtoolsRes['featureName'] = '.'.\
                join(filename.split('.')[:-2])
            row_index = [overlapTable.index.get_loc(i)
                         for i in bedtoolsRes['cpg_name'].unique()]
            feature_name = '.'.join(filename.split('.')[:-2])
            col_index = overlapTable.columns.get_loc(feature_name)
            overlapTable.values[row_index, col_index] = 1

    overlapTable.to_csv(save_path)
    print('The overlapTable looks like: \n', overlapTable.head(), flush=True)
    print('With shape: ', overlapTable.shape, flush=True)
    print('OverlapMatrix Finished.', flush=True)

    # add empty features
    overlapMatrix_filepath = save_path


    empty_featureList_filepath = os.path.join(PROJECT_DIR, "data",
                                              "features",
                                              "emptyFeatureFiles",
                                              "emptyFeatureList.txt")
    add_emptyFeature_nonOverlappingCpg(empty_featureList_filepath,
                                       overlapMatrix_filepath,
                                       original_eqtm_filepath,
                                       complete_overlapMatrix_savepath,
                                       cpgcolname="SNPName")

    # make the overlapRatio table
    overlapTable = pd.read_csv(complete_overlapMatrix_savepath, index_col=0)
    print(overlapTable.head(), flush=True)
    print(overlapTable.columns, flush=True)
    ratio_table = pd.DataFrame(data=0,
                               index=overlapTable.index,
                               columns=histone_list,
                               dtype=np.int8)
    for histone in histone_list:
        print('Processing histone:', histone, flush=True)
        for feature in overlapTable.columns:
            if histone in feature:
                ratio_table[histone] += overlapTable[feature]/len(cell_type_list)

    print(ratio_table.values[ratio_table.values > 0], flush=True)
    ratio_table.to_csv(ratio_table_savepath, sep='\t')
    print("Ratio table saved in path: \n", ratio_table_savepath, flush=True)
