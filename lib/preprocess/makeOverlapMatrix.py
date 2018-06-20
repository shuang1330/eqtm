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

    cpg_name = 'GenomicCoordinates' # probably the only thing needs changing

    # all the directories
    PROJECT_DIR = '/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm'
    feature_folder = '/groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap/consolidatedImputedGappedPeak'

    data_folder = os.path.join(PROJECT_DIR,'data','temp',
                               'features',cpg_name)
    cpg_path = os.path.join(PROJECT_DIR,'data','eqtmZscores','ORIGIN',
                            cpg_name+'.txt')
    save_path = os.path.join(PROJECT_DIR,'data','output',
                             cpg_name,
                             cpg_name+'_overlapMatrix.txt')
    ratio_table_savepath = os.path.join(PROJECT_DIR,'data','output',
                             cpg_name,
                             cpg_name+'_overlapRatio.txt')
    # create the result overlap_matrix
    feature_list = featureList(feature_folder)
    histone_list = histoneList(feature_list)
    cell_type_list = cellList(feature_list)

    cpg = pd.read_csv(cpg_path,sep='\t')
    print('Cpg file with shape:', cpg.shape)
    cpg['chr_pos'] = (['chr{}_{}'.format(i[0],i[1]) for i in
                      zip(cpg["SNPChr"].map(str),(cpg["SNPChrPos"]-25).map(str))])
    cpg_dic = cpg[['chr_pos', 'SNPName']].set_index('chr_pos').T.to_dict()

    def mapChr_pos(row):
        print(row)
        return cpg_dic[row]

    overlapTable = pd.DataFrame(data=0,
                                index=cpg['SNPName'].unique(),
                                columns=feature_list,
                                dtype=np.int8)
    print('With shape: ',overlapTable.shape)
    names = ['chr_feat','str_feat','end_feat', 'type','chr_cpg','str_cpg','end_cpg']

    # start processing
    print('Process features in folder: ',data_folder)
    record_path = os.path.join(PROJECT_DIR,'sbatch','record',cpg_name+'.txt')
    record_file = open(record_path, 'w')

    for filename in os.listdir(data_folder):
        print("Processing feature:",filename)
        record_file.write("Processing feature: {}\n".format(filename))
        record_file.flush()
        filepath = os.path.join(data_folder,filename)
        if os.stat(filepath).st_size == 0:
            print('Empty file:',filename)
            continue
        else:
            bedtoolsRes = pd.read_csv(filepath,sep='\t',names=names)
            bedtoolsRes['featureName'] = '.'.join(filename.split('.')[:-2])
            bedtoolsRes['chr_pos'] = ['_'.join(i) for i in
                                      zip(bedtoolsRes['chr_cpg'],
                                          bedtoolsRes['str_cpg'].map(str))]
            bedtoolsRes['cpgName'] = [cpg_dic[i]['SNPName'] for i in bedtoolsRes['chr_pos']]
            row_index = [overlapTable.index.get_loc(i) for i in bedtoolsRes['cpgName'].unique()]
            feature_name = '.'.join(filename.split('.')[:-2])
            col_index = overlapTable.columns.get_loc(feature_name)
            overlapTable.values[row_index,col_index] = 1
    record_file.close()
    overlapTable.to_csv(save_path)
    print('The overlapTable looks like: \n',overlapTable.head())
    print('With shape: ',overlapTable.shape)
    print('Finished.')

    # make the overlapRatio table
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
                ratio_table[histone] += overlapTable[feature]/len(cell_type_list)

    print(ratio_table.values[ratio_table.values>0])
    ratio_table.to_csv(ratio_table_savepath,sep='\t')
    print("Ratio table saved in path: \n",ratio_table_savepath)
