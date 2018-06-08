#TODO looking into the /dev/shm

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
    if os.path.exists(feature_list_path):
        with open(feature_list_path,'r') as f:
            for line in f.readlines():
                feature_list.append(line.strip())
            f.close()
        print('Read feature list from path: \n',feature_list_path)
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
    # cell_list_path = os.path.join(temp_meta,'cell_list.txt')
    if os.path.exists(cell_list_path):
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
    
def readCpgFileInBedToolsFormat(dirs):
    cpg_filepath = os.path.join(dirs.cpgDir,dirs.cpgname+'.txt')
    dirs.cpg_bedtoolFormat_filepath = os.path.join(dirs.cpgDir,
                                      dirs.cpgname+'_bedtoolsFormat.txt')
    cpg_bedtoolFormat = open(dirs.cpg_bedtoolFormat_filepath,'w')
    with open(cpg_filepath,'r') as cpgFile:
        cpgInfo = [row.strip().split('\t') for row in cpgFile.readlines()[1:]]
        for cpg in cpgInfo:
            cpg_name,cpg_chr,start,end = cpg[0],'chr'+str(cpg[1]),int(cpg[2])-25,int(cpg[2])+25
            cpg_bedtoolFormat.write('{}\t{}\t{}\t{}\n'.format(cpg_chr,start,end,cpg_name))
    cpg_bedtoolFormat.close()
    return cpg_bedtoolFormat

def callBedtools(dirs,execute_bedtoolFilename,feature_list):
    # call bedtools for intersect files
    findOverlap_path = os.path.join(dirs.project_rootdir,'tools',
                                    execute_bedtoolFilename+'.sh')
    print(dirs.cpg_bedtoolFormat_filepath)
    for feature in feature_list:
        feature_filename = feature+dirs.feature_filetype
        subprocess.check_call([findOverlap_path,
                               '%s'%dirs.project_rootdir,
                               '%s'%dirs.feature_folder,
                               '%s'%dirs.temp_featureFolder,
                               '%s'%dirs.cpgname,
                               '%s'%dirs.cpg_bedtoolFormat_filepath,
                               '%s'%feature,
                               '%s'%feature_filename,
                               '%s'%dirs.feature_type])
        print('Built the overlap matrix for feature: %s'%feature)

def overlapMatrixPipeline(dirs,execute_bedtoolFilename):
    feature_list = readFeatureList(dirs)
    histone_list = readHistoneList(dirs.temp_meta,dirs.histone_list_filename,
                                   feature_list)
    cell_list = readCellList(dirs.temp_meta,dirs.cell_list_filename,
                             feature_list)
    _ = readCpgFileInBedToolsFormat(dirs)
    callBedtools(dirs,execute_bedtoolFilename,feature_list)
    print("Bedtools results saved in feature folder.")
