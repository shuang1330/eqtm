import __init__path
import os
import numpy as np
import pandas as pd
import argparse
from time import time


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

    ratio_table_savepath = os.path.join(PROJECT_DIR,
                                        'data',
                                        'eqtmZscores',
                                        'allCpgs',
                                        name+'_overlapRatio.txt')
    ratio_table_complete_savepath = os.path.join(PROJECT_DIR,
                                        'data',
                                        'eqtmZscores',
                                        'allCpgs',
                                        name + '_overlapRatio_withoriginalInfo.txt')

    # create the result overlap_matrix
    feature_list = featureList(feature_folder)
    histone_list = histoneList(feature_list)
    cell_type_list = cellList(feature_list)
    print("Cell type number: ", len(cell_type_list))
    # raise NotImplementedError

    original_eqtm = pd.read_csv(original_eqtm_filepath, sep='\t')
    original_eqtm['chromosome'] = ['chr{}'.format(chrom) for chrom in original_eqtm['chr']]
    original_eqtm["index"] = ["_".join([str(ele) for ele in item]) for item in
                              original_eqtm[["chromosome", "pos", 'rsID']].values]

    original_eqtm["index_nr"] = [i for i in range(original_eqtm.shape[0])]
    print(max(original_eqtm['index_nr']))
    print("Original data: ")
    print(original_eqtm.head())
    startime = time()
    index_nr = 0
    original_eqtm_index_dict = {}
    for index in original_eqtm['index'].unique():
        original_eqtm_index_dict[index] = index_nr
        index_nr += 1
    # original_eqtm_index_dict = original_eqtm[["index", "index_nr"]].set_index("index").T.to_dict()
    print("Converting to dictionary time:", time() - startime)
    print(len(original_eqtm_index_dict.keys()))
    print(original_eqtm.shape)

    names = ['chr_feat', 'str_feat', 'end_feat', 'type',
             'chr_cpg', 'str_cpg', 'end_cpg', 'rs']

    # start processing
    print('Process features in folder: ', data_folder)

    ratio_data = np.zeros([len(original_eqtm_index_dict.keys()), len(histone_list)])
    histone_dict = {}
    for i in range(len(histone_list)):
        histone_dict[histone_list[i]] = i
    file_index = 0
    for filename in os.listdir(data_folder):
        file_index += 1
        if file_index%200 ==0:
            print("Processing feature:", filename, file_index, flush=True)
        filepath = os.path.join(data_folder, filename)
        if os.stat(filepath).st_size == 0:
            print('Empty file:', filename, flush=True)
            continue
        else:
            bedtoolsRes = pd.read_csv(filepath, sep='\t', names=names)
            bedtoolsRes["original_pos"] = bedtoolsRes["str_cpg"] + 25
            bedtoolsRes["index_name"] = ["_".join([str(ele) for ele in item]) for item
                                         in bedtoolsRes[['chr_cpg', 'original_pos', 'rs']].values]
            bedtoolsRes['featureName'] = '.'.\
                join(filename.split('.')[:-2])
            row_index = [original_eqtm_index_dict[item] for item in bedtoolsRes["index_name"].values]
            feature_name = '.'.join(filename.split('.')[:-2])
            histone_name = feature_name.split("-")[1]
            histone_index = histone_dict[histone_name]
            ratio_data[row_index, histone_index] += 1.0/len(cell_type_list)

    ratio_dataframe = pd.DataFrame(data=ratio_data, index=original_eqtm['index'].unique(), columns=histone_list)
    print(ratio_dataframe.head())
    ratio_dataframe.to_csv(ratio_table_savepath, sep="\t")

    # merged_data = pd.concat([ratio_dataframe, original_eqtm], axis=1)
    # print(merged_data.head())
    # merged_data.to_csv(ratio_table_complete_savepath, sep="\t", index=False)


