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
                               'features','2017-12-09-eQTLsFDR-gt0.0-flipped')
    cpg_path = os.path.join(PROJECT_DIR,'data','eqtmZscores','2017-12-09-eQTLsFDR-gt0.0-flipped.txt')
    feature_folder = '/groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap/consolidatedImputedGappedPeak'
    save_path = os.path.join(PROJECT_DIR,'data','output','2017-12-09-eQTLsFDR-gt0.0-flipped','2017-12-09-eQTLsFDR-gt0.0-flipped_overlapMatrix.txt')
    ratio_table_savepath = os.path.join(PROJECT_DIR,'data','output','2017-12-09-eQTLsFDR-gt0.0-flipped','2017-12-09-eQTLsFDR-gt0.0-flipped_overlapRatio.txt')
    # create the result overlap_matrix
    feature_list = featureList(feature_folder)
    histone_list = histoneList(feature_list)
    cell_type_list = cellList(feature_list)
    cpg = pd.read_csv(cpg_path,sep='\t')
    print('Cpg file with shape:', cpg.shape)

    def mapStartEndPosition(row):
        return cpg['SNPName'][(cpg['SNPChr']==int(row[4][3:]))&\
        (cpg['SNPChrPos']==row[5]+25)].values[0]

    overlapTable = pd.DataFrame(data=0,
                                index=cpg['SNPName'].unique(),
                                columns=feature_list,
                                dtype=np.int8)
    print('With shape: ',overlapTable.shape)
    names = ['chr_feat','str_feat','end_feat', 'type','chr_cpg','str_cpg','end_cpg']

    for filename in os.listdir(data_folder):
        print("Processing feature:",filename)
        filepath = os.path.join(data_folder,filename)
        if os.stat(filepath).st_size == 0:
            print('Empty file:',filename)
            continue
        else:
            bedtoolsRes = pd.read_csv(filepath,sep='\t',names=names)
            bedtoolsRes['featureName'] = '.'.join(filename.split('.')[:-2])
            bedtoolsRes['cpgName'] = bedtoolsRes.apply(mapStartEndPosition,axis=1)
            overlapTable.loc[bedtoolsRes['cpgName'],bedtoolsRes['featureName']]=1


    overlapTable.to_csv(save_path)
    print('The overlapTable looks like: \n',overlapTable.head())
    print('With shape: ',overlapTable.shape)
    print('Finished.')

    # overlapTable = pd.read_csv(save_path)
    print(overlapTable.head())
    print(overlapTable.columns)
    ratio_table = pd.DataFrame(data=0,
                               index=overlapTable.index,
                               columns = histone_list,
                               dtype=np.int8)
    for histone in histone_list:
        print('Processing histone:',histone)
        for feature in overlapTable.columns:
            if histone in feature:
                # print(ratio_table.index.is_unique,overlapTable.index.is_unique)
                # raise NotImplementedError
                ratio_table[histone] += overlapTable[feature]/len(cell_type_list)

    ratio_table.to_csv(ratio_table_savepath,sep='\t')
    print("Ratio table saved in path: \n",ratio_table_savepath)
