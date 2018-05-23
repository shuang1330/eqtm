import os
import numpy as np
import pandas as pd
import gzip
import shlex, subprocess

def featureList(feature_folder):
    feature_list = []
    for name in os.listdir(feature_folder):
        if name.endswith('gz'):
            feature_name = '.'.join(name.split('.')[:-5])
            feature_list.append(feature_name)
    print('featurelist:',feature_list)
    return feature_list

def histoneList(feature_list):
    histone_list = []
    for col in feature_list:
        histone = col.split('-')[1]
        if histone not in histone_list:
            histone_list.append(histone)
    print('Histone list:',histone_list)
    return histone_list

def cellList(feature_list):
    cell_list = []
    for col in feature_list:
        cell = col.split('-')[0]
        if cell not in cell_list:
            cell_list.append(cell)
    print('Cell list:',cell_list)
    return cell_list

if __name__=='__main__':
    PROJECT_DIR = '/groups/umcg-gcc/tmp03/umcg-sli/eqtm_project'
    data_folder = os.path.join(PROJECT_DIR,'data','temp',
                               'features','random20k_gt0.1_sm0.2')
    cpg_path = os.path.join(PROJECT_DIR,'data','eqtmZscores','random20k_gt0.1_sm0.2.txt')
    feature_folder = '/groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap/consolidatedImputedGappedPeak'
    save_path = os.path.join(PROJECT_DIR,'data','output','random20k_gt0.1_sm0.2','random20k_gt0.1_sm0.2_overlapMatrix.txt')
    # create the result overlap_matrix
    feature_list = featureList(feature_folder)
    histone_list = histoneList(feature_list)
    cell_type_list = cellList(feature_list)
    cpg = pd.read_csv(cpg_path,sep='\t')
    overlapTable = pd.DataFrame(data=0,
                                index=cpg['SNPName'],
                                columns=feature_list,
                                dtype=np.int8)
    names = ['chr_feat','str_feat','end_feat', 'type','chr_cpg','str_cpg','end_cpg']

    def mapStartEndPosition(row):
        return cpg['SNPName'][(cpg['SNPChr']==int(row[4][3:]))&\
        (cpg['SNPChrPos']==row[5]+25)].values[0]

    print('Start Processing.')
    # print(os.listdir(data_folder))
    for filename in os.listdir(data_folder):
        print("Processing feature:",filename)
        filepath = os.path.join(data_folder,filename)
        try:
            bedtoolsRes = pd.read_csv(filepath,sep='\t',names=names)
            # print(bedtoolsRes.head())
            bedtoolsRes['featureName'] = '.'.join(filename.split('.')[:-2])
            bedtoolsRes['cpgName'] = bedtoolsRes.apply(mapStartEndPosition,axis=1)
            # print(bedtoolsRes.head())
            # raise NotImplementedError
            overlapTable.loc[bedtoolsRes['cpgName'],bedtoolsRes['featureName']]=1
        except:
            print('Empty File:',filename)
            continue
    overlapTable.to_csv(save_path)
    print(overlapTable.head())
    print(overlapTable.shape)
    print('Finished.')
