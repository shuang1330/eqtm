import os
import numpy as np
import pandas as pd
import gzip
import shlex, subprocess
import time

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

if __name__=='__main__':
    PROJECT_DIR = '/groups/umcg-gcc/tmp03/umcg-sli/eqtm_project'
    data_folder = os.path.join(PROJECT_DIR,'data','temp',
                               'overlap','random20k_gt0.5')
    cpg_path = os.path.join(PROJECT_DIR,'data','eqtmZscores','random20k_gt0.5.txt')
    feature_folder = '/groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap/consolidatedImputedGappedPeak'
    save_path = os.path.join(PROJECT_DIR,'data','output','random20k_gt0.5','random20k_gt0.5_overlapMatrix.txt')
    ratio_table_savepath = os.path.join(PROJECT_DIR,'data','output','random20k_gt0.5','random20k_gt0.5_overlapRatio.txt')
    # create the result overlap_matrix
    feature_list = featureList(feature_folder)
    histone_list = histoneList(feature_list)
    cell_type_list = cellList(feature_list)
    cpg = pd.read_csv(cpg_path,sep='\t')
    overlapTable = pd.DataFrame(data=0,
                                index=cpg['SNPName'],
                                columns=feature_list,
                                dtype=np.int8)
    names = ['chr_feat','str_feat','end_feat', 'type',
             'chr_cpg','str_cpg','end_cpg','featureName','cpgName']
    dur1 = []
    dur2 = []
    for filename in os.listdir(data_folder)[:10]:
        print("Processing feature:",filename)
        filepath = os.path.join(data_folder,filename)
        try:
            start_time = time.time()
            bedtoolsRes = pd.read_csv(filepath,sep=',',names=names)
            time1 = time.time()
            dur1.append(time1-start_time)
            overlapTable.loc[bedtoolsRes['cpgName'].unique(),
                             bedtoolsRes['featureName']]=1
            dur2.append(time.time()-time1)
        except:
            print('Empty file:',filename)
            continue
    print(dur1)
    print(dur2)
    raise NotImplementedError
    overlapTable.to_csv(save_path)
    print(overlapTable.head())
    print(overlapTable.shape)
    print('Finished.')

    # overlapTable = pd.read_csv(save_path)
    print(overlapTable.head())
    print(overlapTable.columns)
    ratio_table = pd.DataFrame(data=0,
                               index=overlapTable.index,
                               columns = histone_list,
                               dtype=np.int8)
    ratio_table['SNPName'] = overlapTable['SNPName']
    for histone in histone_list:
        print('Processing histone:',histone)
        for feature in overlapTable.columns:
            if histone in feature:
                # print(ratio_table.index.is_unique,overlapTable.index.is_unique)
                # raise NotImplementedError
                ratio_table[histone] += overlapTable[feature]/len(cell_type_list)

    ratio_table.to_csv(ratio_table_savepath,sep='\t')
    print("Ratio table saved in path: \n",ratio_table_savepath)
