import os
import numpy as np
import pandas as pd
import gzip
import shlex, subprocess

def isEmpty(path):
    return os.stat(path).st_size==0

def readFeatureList(dirs):
    feature_list = []
    feature_list_path = os.path.join(dirs.temp_meta,
                                     dirs.feature_list_filename+'.txt')
    if os.path.exists(dirs.feature_list_filename):
        with open(feature_list_path,'r') as f:
            for line in f.readlines():
                feature_list.append(line.strip())
            f.close()
        if 'feature_list' in feature_list:
            print('There should not be a item called feature list.')
    else:
        feature_f = open(feature_list_path,'w')
        for name in os.listdir(dirs.feature_folder):
            if name.endswith('gz'):
                feature_name = '.'.join(name.split('.')[:-5])
                feature_list.append(feature_name)
                feature_f.write('%s\n'%feature_name)
        feature_f.close()
        print('Saved feature list to path: \n%s'%feature_list_path)

    return feature_list

def readHistoneList(temp_meta,histone_list_filename,feature_list):
    histone_list = []
    histone_list_path = os.path.join(temp_meta,histone_list_filename+'.txt')
    if os.path.exists(histone_list_path):
        with open(histone_list_path,'r') as f:
            for line in f.readlines():
                histone_list.append(line.strip())
            f.close()
        print('Read histone list from path: \n',histone_list_path)
        print(histone_list)
    else:
        for col in feature_list:
            histone = col.split('-')[1]
            if histone not in histone_list:
                histone_list.append(histone)
        with open(histone_list_path,'w+') as f:
            for item in histone_list:
                f.write('%s\n'%item)
            f.close()
        print('Saved histone list to path:\n',histone_list_path)

    return histone_list

def readCellList(temp_meta,cell_list_filename,feature_list):
    cell_list = []
    cell_list_filename = 'cell_list'
    cell_list_path = os.path.join(temp_meta,cell_list_filename+'.txt')
    if os.path.exists(cell_list_filename):
        with open(cell_list_path,'r') as f:
            for line in f.readlines():
                cell_list.append(line.strip())
            f.close()
        print('Read cell types from path: \n',cell_list_path)
        print(cell_list)
    else:
        for col in feature_list:
            cell = col.split('-')[0]
            if cell not in cell_list:
                cell_list.append(cell)
        with open(cell_list_path,'w+') as f:
            for item in cell_list:
                f.write('%s\n'%item)
            f.close()
        print('Saved cell list to path:\n',cell_list_path)

    return cell_list

def readCpgFileInbedToolsFormat(dirs):
    # cpg_filename = 'bonder-eQTMsFDR0.0-CpGLevel-split'
    cpg_path = os.path.join(dirs.cpg_folder,
                            dirs.cpg_filename+'.txt')
    cpg = pd.read_csv(cpg_path,sep='\t')
    cpg['startSite'] = cpg['SNPChrPos']-25
    cpg['endSite'] = cpg['SNPChrPos']+25
    cpg['chr'] = cpg['SNPChr'].apply(lambda x:'chr{}'.format(str(x)))

    # create the cpg file with bedtools format
    cpg_bedtoolFormat_name = dirs.cpg_filename+'_bedtoolsFormat'
    cpg_bedtoolFormat_path = os.path.join(dirs.temp_cpgFolder,
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

def callBedtools(dirs,execute_bedtoolFilename,feature_list):
    # call bedtools for intersect files
    cpg_bedtoolFormat_name = dirs.cpg_filename+'_bedtoolsFormat'
    findOverlap_path = os.path.join(dirs.project_rootdir,
                                    execute_bedtoolFilename+'.sh')
    for feature in feature_list:
        subprocess.check_call([findOverlap_path,
                               '%s'%dirs.project_rootdir,
                               '%s'%dirs.feature_folder,
                               '%s'%dirs.cpg_folder,
                               '%s'%dirs.temp_cpgFolder,
                               '%s'%dirs.temp_featureFolder,
                               '%s'%dirs.temp_output,
                               '%s'%feature,
                               '%s'%cpg_bedtoolFormat_name])
        print('Built the overlap matrix for feature: %s'%feature)

def buildRatioMatrix(dirs,histone_list,feature_list,cell_list):
    # build the overlap ratio matrix
    cpg_path = os.path.join(dirs.cpg_folder,
                            dirs.cpg_filename+'.txt')
    cpg = pd.read_csv(cpg_path,sep='\t')
    def mapStartEndPosition(row):
        return cpg['SNPName'][(cpg['SNPChr']==int(row[4][3:]))&\
        (cpg['SNPChrPos']==row[5]+25)].values[0]
    res = pd.DataFrame(data=0,
                       index=cpg['SNPName'],
                       columns=feature_list,
                       dtype=np.int8)
    for intersect_raw in os.listdir(dirs.temp_featureFolder):
        intersect_rawPath = os.path.join(dirs.temp_featureFolder,
                                         intersect_raw)
        if not isEmpty(intersect_rawPath):
            print('Processing feature:',intersect_raw)
            test_intersectFile = pd.read_csv(intersect_rawPath,
                                             sep='\t',
                                             header=None)
            test_intersectFile['featureName'] = '.'.join(intersect_raw.split('.')[:-2])
            test_intersectFile['cpgName'] = \
            test_intersectFile.apply(mapStartEndPosition,axis=1)
            res.loc[test_intersectFile['cpgName'],
                    test_intersectFile['featureName']]=1
            save_path = os.path.join(dirs.temp_overlap,
                                     intersect_raw)
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

    ratio_table_name = dirs.cpg_filename+'_ratioTable'
    ratio_table_savepath = os.path.join(dirs.output_ratio,ratio_table_name+'.csv')
    ratio_table.to_csv(ratio_table_savepath,
                       sep='\t')
    print("Ratio table saved in path: \n",ratio_table_savepath)

def overlapMatrixPipeline(dirs,execute_bedtoolFilename):
    feature_list = readFeatureList(dirs)
    histone_list = readHistoneList(dirs.temp_meta,dirs.histone_list_filename,
                                   feature_list)
    cell_list = readCellList(dirs.temp_meta,dirs.cell_list_filename,
                             feature_list)
    _ = readCpgFileInbedToolsFormat(dirs)
    callBedtools(dirs,execute_bedtoolFilename,feature_list)
    buildRatioMatrix(dirs,histone_list,feature_list,cell_list)
    print("Finished.")
