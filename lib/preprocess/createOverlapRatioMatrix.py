import os
import numpy as np
import pandas as pd
import gzip
import shlex, subprocess

def isEmpty(path):
    return os.stat(path).st_size==0

def readFeatureList(feature_folder,feature_list_filename):
    feature_list = []
    # feature_list_filename = 'feature_list'
    feature_list_path = os.path.join(feature_folder,
                                     feature_list_filename+'.txt')
    if os.path.exists(feature_list_filename):
        with open(feature_list_path,'r') as f:
            for line in f.readlines():
                feature_list.append(line.strip())
            f.close()
        if 'feature_list' in feature_list:
            print('There should not be a item called feature list.')
    else:
        # create and save a feature list
        feature_list_savepath = os.path.join(feature_folder,
                                             feature_list_filename+'.txt')
        feature_f = open(feature_list_savepath,'w')
        for name in os.listdir(feature_folder):
            if name.endswith('gz'):
                feature_name = '.'.join(name.split('.')[:-5])
                feature_list.append(feature_name)
                feature_f.write('%s\n'%feature_name) # write to a feature list file
        feature_f.close()
        print('Saved feature list to path: \n%s'%feature_list_savepath)

    return feature_list

def readHistoneList(histone_folder,histone_list_filename,feature_list):
    histone_list = []
    # histone_list_filename = 'histone_list'
    histone_list_path = os.path.join(histone_folder,'histone_list.txt')
    if os.path.exists(histone_list_path):
        with open(histone_list_path,'r') as f:
            for line in f.readlines():
                histone_list.append(line.strip())
            f.close()
        print('Read histone list from path: \n',histone_list_path)
        print(histone_list)
    else:
        histone_list_savepath = os.path.join(histone_folder,
                                             histone_list_filename+'.txt')
        for col in feature_list:
            histone = col.split('-')[1]
            if histone not in histone_list:
                histone_list.append(histone)
        with open(histone_list_savepath,'w+') as f:
            for item in histone_list:
                f.write('%s\n'%item)
            f.close()
        print('Saved histone list to path:\n',histone_list_savepath)

    return histone_list

def readCellList(cell_folder,cell_list_filename,feature_list):
    cell_list = []
    # cell_list_filename = 'cell_list'
    cell_list_path = os.path.join(cell_folder,'cell_list.txt')
    if os.path.exists(cell_list_filename):
        with open(cell_list_path,'r') as f:
            for line in f.readlines():
                cell_list.append(line.strip())
            f.close()
        print('Read cell list from path: \n',cell_list_path)
        print(cell_list)
    else:
        # create and save a cell type list
        cell_list_savepath = os.path.join(cell_folder,
                                          cell_list_filename+'.txt')
        for col in feature_list:
            cell = col.split('-')[0]
            if cell not in cell_list:
                cell_list.append(cell)
        with open(cell_list_savepath,'w+') as f:
            for item in cell_list:
                f.write('%s\n'%item)
            f.close()
        print('Saved cell list to path:\n',cell_list_savepath)

    return cell_list

def readCpgFileInbedToolsFormat(temp_cpgFolder,cpg_folder,cpg_filename):
    # cpg_filename = 'bonder-eQTMsFDR0.0-CpGLevel-split'
    cpg_path = os.path.join(cpg_folder,
                            cpg_filename+'.txt')
    cpg = pd.read_csv(cpg_path,sep='\t')
    cpg['startSite'] = cpg['SNPChrPos']-25
    cpg['endSite'] = cpg['SNPChrPos']+25
    cpg['chr'] = cpg['SNPChr'].apply(lambda x:'chr{}'.format(str(x)))

    # create the cpg file with bedtools format
    cpg_bedtoolFormat_name = cpg_filename+'_bedtoolsFormat'
    cpg_bedtoolFormat_path = os.path.join(temp_cpgFolder,
    cpg_bedtoolFormat_name+'.tsv')
    cpg = cpg.sort_values(by=['SNPChr','startSite'])
    cpg[['chr','startSite','endSite']].to_csv(cpg_bedtoolFormat_path,
                                              header=False,index=False,
                                              sep='\t')
    print('Saved the cpg file in bedtools format in path: \n%s'%cpg_bedtoolFormat_path)
    # check the saved cpg file
    def check_cpg():
        check_cpg = pd.read_csv(cpg_bedtoolFormat_path,sep='\t',header=None)
        print('Hey, please check the saved file:\n',check_cpg.head())
    check_cpg()
    return cpg[['chr','startSite','endSite']].copy()

def callBedtools(PROJECT_DIR,execute_bedtoolFilename,
                 feature_list,cpg_filename):
    # call bedtools for intersect files
    # execute_bedtoolFilename = 'findOverlap.sh'
    cpg_bedtoolFormat_name = cpg_filename+'_bedtoolsFormat'
    findOverlap_path = os.path.join(PROJECT_DIR,execute_bedtoolFilename+'.sh')
    for feature in feature_list:
        subprocess.check_call([findOverlap_path,
                               '%s'%feature,
                               '%s'%cpg_bedtoolFormat_name])
        print('Built the overlap matrix for feature: %s'%feature)

def buildRatioMatrix(cpg_folder,cpg_filename,
                     temp_output,temp_output_withcpgName,
                     histone_list,feature_list,cell_list,
                     ratio_folder):
    # build the overlap ratio matrix
    cpg_path = os.path.join(cpg_folder,
                            cpg_filename+'.txt')
    cpg = pd.read_csv(cpg_path,sep='\t')
    def mapStartEndPosition(row):
        return cpg['SNPName'][(cpg['SNPChr']==int(row[4][3:]))&\
        (cpg['SNPChrPos']==row[5]+25)].values[0]

    res = pd.DataFrame(data=0,
                       index=cpg['SNPName'],
                       columns=feature_list,
                       dtype=np.int8)

    for intersect_raw in os.listdir(temp_output):
        intersect_rawPath = os.path.join(temp_output,
                                         intersect_raw)
        if not isEmpty(intersect_rawPath):
            print('Processing feature:',intersect_raw)
            test_intersectFile = pd.read_csv(intersect_rawPath,
                                             sep='\t',
                                             header=None)
            ############################################
            # print('.'.join(intersect_raw.split('.')[:-2]))
            test_intersectFile['featureName'] = '.'.join(intersect_raw.split('.')[:-2])
            ############################################
            test_intersectFile['cpgName'] = \
            test_intersectFile.apply(mapStartEndPosition,axis=1)
            res.loc[test_intersectFile['cpgName'],
                    test_intersectFile['featureName']]=1
            save_path = os.path.join(temp_output_withcpgName,
                                     intersect_rawPath)
            test_intersectFile.to_csv(save_path,
                                      index=False,
                                      header=False)
        else:
            print('Empty File: ',intersect_raw)

    # create the overlapping ratio
    ratio_table = pd.DataFrame(data=0,
                               index=cpg['SNPName'],
                               columns = histone_list,
                               dtype=np.int8)
    for histone in histone_list:
        for feature in res.columns:
            if histone in feature:
                ratio_table[histone] += res[feature]/len(cell_list)

    ratio_table_name = cpg_filename+'_ratioTable'
    ratio_table_savepath = os.path.join(ratio_folder,ratio_table_name+'.csv')
    ratio_table.to_csv(ratio_table_savepath,
                       sep='\t')
    print("Ratio table saved in path: \n",ratio_table_savepath)

def overlapMatrixPipeline(PROJECT_DIR,execute_bedtoolFilename,
             feature_folder,feature_list_filename,
             histone_folder,histone_list_filename,
             temp_cpgFolder,cpg_folder,cpg_filename,
             temp_output,temp_output_withcpgName,ratio_folder,
             cell_list_filename):
    feature_list = readFeatureList(feature_folder,feature_list_filename)
    print(feature_list)
    histone_list = readHistoneList(feature_folder,histone_list_filename,feature_list)
    cell_list = readCellList(feature_folder,cell_list_filename,feature_list)
    _ = readCpgFileInbedToolsFormat(temp_cpgFolder,cpg_folder,cpg_filename)
    callBedtools(PROJECT_DIR,execute_bedtoolFilename,
                 feature_list,cpg_filename)
    buildRatioMatrix(cpg_folder,cpg_filename,
                     temp_output,temp_output_withcpgName,
                     histone_list,feature_list,cell_list,
                     ratio_folder)
    print("Finished.")
