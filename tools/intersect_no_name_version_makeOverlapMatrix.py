import __init__path
import os
import numpy as np
import pandas as pd
from lib.read.read_roadmap_features import roadmap


def create_overlapMatrix(cpg_path,
                         feature_list,
                         data_folder,
                         save_path):
    cpg = pd.read_csv(cpg_path, sep=',')
    print(cpg.head())
    print('Cpg file with shape:', cpg.shape)
    cpg['chr_pos'] = (['chr{}_{}'.format(i[0], i[1]) for i in
                       zip(cpg["chr"].map(str), (cpg["startSite"] - 100).map(str))])
    cpg_dic = cpg[['chr_pos', 'geneName']].set_index('chr_pos').T.to_dict()

    overlapTable = pd.DataFrame(data=0,
                                index=cpg['geneName'].unique(),
                                columns=feature_list,
                                dtype=np.int8)
    print('Process features in folder: ', data_folder)
    names = ['chr_feat', 'str_feat', 'end_feat', 'type', 'chr_cpg', 'str_cpg', 'end_cpg']
    for filename in os.listdir(data_folder):
        print("Processing feature:", filename, flush=True)
        filepath = os.path.join(data_folder, filename)
        if os.stat(filepath).st_size == 0:
            print('Empty file:', filename, flush=True)
            continue
        else:
            bedtoolsRes = pd.read_csv(filepath, sep='\t', names=names)
            bedtoolsRes['featureName'] = '.'.join(filename.split('.')[:-2])
            bedtoolsRes['chr_pos'] = ['_'.join(i) for i in
                                      zip(bedtoolsRes['chr_cpg'],
                                          bedtoolsRes['str_cpg'].map(str))]
            bedtoolsRes['cpgName'] = [cpg_dic[i]['geneName'] for i in bedtoolsRes['chr_pos']]
            row_index = [overlapTable.index.get_loc(i) for i in bedtoolsRes['cpgName'].unique()]
            feature_name = '.'.join(filename.split('.')[:-2])
            col_index = overlapTable.columns.get_loc(feature_name)
            overlapTable.values[row_index, col_index] = 1
    overlapTable.to_csv(save_path)
    print('The overlapTable looks like: \n', overlapTable.head(), flush=True)
    print('With shape: ', overlapTable.shape, flush=True)
    return overlapTable


def create_overlapRatio(overlapTable,
                        histone_list,
                        cell_list,
                        ratio_table_savepath):
    ratio_table = pd.DataFrame(data=0,
                               index=overlapTable.index,
                               columns=histone_list,
                               dtype=np.int8)
    for histone in histone_list:
        print('Processing histone:', histone, flush=True)
        for feature in overlapTable.columns:
            if histone in feature:
                ratio_table[histone] += overlapTable[feature] / len(cell_list)

    ratio_table.to_csv(ratio_table_savepath, sep='\t')
    print("Ratio table saved in path: \n", ratio_table_savepath, flush=True)
    return ratio_table


def create_promoter_overlapMatrix(cpg_path,
                                  feature_list,
                                  data_folder,
                                  save_path):
    cpg_names = ["chr", "startSite", "endSite", "geneName"]
    cpg = pd.read_csv(cpg_path, sep='\t', names=cpg_names)
    print('Gene file with shape:', cpg.shape)
    print(cpg.head())
    cpg['chr_pos'] = (['{}_{}'.format(i[0], i[1]) for i in
                       zip(cpg["chr"].map(str), (cpg["startSite"]).map(str))])
    cpg_dic = cpg[['chr_pos', 'geneName']].set_index('chr_pos').T.to_dict()
    for key in cpg_dic:
        if "chr1_" in key:
            print(key, cpg_dic[key])

    overlapTable = pd.DataFrame(data=0,
                                index=cpg['geneName'].unique(),
                                columns=feature_list,
                                dtype=np.int8)
    print('Process features in folder: ', data_folder)
    names = ['chr_feat', 'str_feat', 'end_feat', 'type', 'chr_cpg', 'str_cpg', 'end_cpg']
    for filename in os.listdir(data_folder):
        print("Processing feature:", filename, flush=True)
        filepath = os.path.join(data_folder, filename)
        if os.stat(filepath).st_size == 0:
            print('Empty file:', filename, flush=True)
            continue
        else:
            bedtoolsRes = pd.read_csv(filepath, sep='\t', names=names)
            bedtoolsRes['featureName'] = '.'.join(filename.split('.')[:-2])
            bedtoolsRes['chr_pos'] = ['_'.join(i) for i in
                                      zip(bedtoolsRes['chr_cpg'],
                                          bedtoolsRes['str_cpg'].map(str))]
            bedtoolsRes['cpgName'] = [cpg_dic[i]['geneName'] for i in bedtoolsRes['chr_pos']]
            row_index = [overlapTable.index.get_loc(i) for i in bedtoolsRes['cpgName'].unique()]
            feature_name = '.'.join(filename.split('.')[:-2])
            col_index = overlapTable.columns.get_loc(feature_name)
            overlapTable.values[row_index, col_index] = 1
    overlapTable.to_csv(save_path)
    print('The overlapTable looks like: \n', overlapTable.head(), flush=True)
    print('With shape: ', overlapTable.shape, flush=True)
    return overlapTable


if __name__=='__main__':
    project_rootdir = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
    cpg_path = os.path.join(project_rootdir, "data", "features",
                            "TSSDistance",
                            "promoter_startEndSite_bedtoolsFormat.txt")
    roadmap_features = roadmap()
    feature_list = roadmap_features.readFeatureList(
        feature_list_filepath=os.path.join(project_rootdir, "data", "temp", "meta",
                                           "feature_list.txt"))
    data_folder = os.path.join(project_rootdir, "data",
                               "temp", "intersect_features",
                               "promoter")
    save_path = os.path.join(project_rootdir, "data", "features", "geneOverlap",
                             "promoter_overlapMatrix_complete.txt")
    # overlapTable = create_promoter_overlapMatrix(cpg_path,
    #                                              feature_list,
    #                                              data_folder,
    #                                              save_path)

    overlapTable = pd.read_csv(save_path, index_col=0)
    histone_list = roadmap_features.readHistoneList(os.path.join(project_rootdir, "data",
                                                                 "temp", "meta",
                                                                 "histone_list.txt"))
    cell_list = roadmap_features.readCellList(os.path.join(project_rootdir,
                                                           "data", "temp", "meta",
                                                           "cell_list.txt"))
    ratio_table_savepath = os.path.join(project_rootdir, "data",
                                        "features", "geneOverlap",
                                        "promoter_overlapRatio.txt")
    _ = create_overlapRatio(overlapTable,
                            histone_list,
                            cell_list,
                            ratio_table_savepath)

