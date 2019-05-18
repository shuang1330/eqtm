import os
import pandas as pd
import numpy as np


def read_bloodCelltypeFile(blood_celltype_filepath):
    with open(blood_celltype_filepath,'r') as f:
        return list(filter(None,[item.strip() for item in f.readlines()]))


def read_allCelltypeFile(all_celltype_filepath):
    dic = {} # cell type number as key, description and other stuff as values
    with open(all_celltype_filepath,'r',encoding="ISO-8859-1") as f:
        content = f.readlines()
        for line in content:
            line_list = line.strip().split('\t')
            dic[line_list[0]] = '_'.join(line_list[1:])
    return dic


def read_histoneList(histonelist_filepath):
    with open(histonelist_filepath,'r') as f:
        return [row.strip() for row in f.readlines()]


def calculate_overlapRatio(overlapMatrix,
                           histonelist,
                           num_celltype,
                           overlapRatio_savepath):
    overlapRatio = pd.DataFrame(
    data=0,
    index=overlapMatrix.index,
    columns=histonelist
    )
    for col in overlapMatrix.columns:
        histone_name = col.split('-')[1]
        overlapRatio[histone_name] += overlapMatrix[col]/num_celltype
    overlapRatio.to_csv(overlapRatio_savepath)
    print(overlapRatio.head())


def select_bloodNonBloodCelltypes(blood_celltype_filepath,
                                  all_celltype_filepath,
                                  histonelist_filepath,
                                  overlapMatrix_name,
                                  overlapMatrix_filepath,
                                  overlapMatrix_save_dirpath,
                                  overlapRatio_name,
                                  overlapRatio_save_dirpath
                                  ):

    blood_celltype = read_bloodCelltypeFile(blood_celltype_filepath)
    all_celltype_dic = read_allCelltypeFile(all_celltype_filepath)

    # read overlapMatrix
    overlapMatrix = pd.read_csv(overlapMatrix_filepath,index_col=0)

    # blood_celltype
    blood_features = []
    for featurename in overlapMatrix.columns:
        for celltype in blood_celltype:
            if celltype in featurename:
                blood_features.append(featurename)
    overlapMatrix_bloodCelltype = overlapMatrix[blood_features].copy()
    overlapMatrix_nonBlood = overlapMatrix[[col for col in overlapMatrix.columns if col not in blood_features]].copy()
    overlapMatrix_blood_savepath = os.path.join(overlapMatrix_save_dirpath,
                                overlapMatrix_name + '_overlapMatrixblood.txt')
    overlapMatrix_nonBlood_savepath = os.path.join(overlapMatrix_save_dirpath,
                                overlapMatrix_name + '_overlapMatrixnonBlood.txt')

    overlapMatrix_bloodCelltype.to_csv(overlapMatrix_blood_savepath)
    overlapMatrix_nonBlood.to_csv(overlapMatrix_nonBlood_savepath)
    print('Blood celltype Matrix saved to {}\nNon-BloodCelltype saved to {}'.format(
    overlapMatrix_blood_savepath, overlapMatrix_nonBlood_savepath
    ))

    # calculate the overlapratio
    histonelist = read_histoneList(histonelist_filepath)
    num_allCellType = len(all_celltype_dic.keys())
    num_bloodCellType = len(blood_celltype)
    num_nonBlood = num_allCellType-num_bloodCellType
    ratio_allcelltype_savepath = os.path.join(
    overlapRatio_save_dirpath,
    overlapRatio_name+'allCelltype.txt'
    )
    ratio_blood_savepath = os.path.join(overlapRatio_save_dirpath,
                               overlapRatio_name+'blood.txt')
    ratio_nonBlood_savepath = os.path.join(overlapRatio_save_dirpath,
                                  overlapRatio_name +'nonBlood.txt')
    calculate_overlapRatio(overlapMatrix,
                           histonelist,
                           num_allCellType,
                           ratio_allcelltype_savepath)
    print('OverlapRatio for all celltypes saved in:',ratio_allcelltype_savepath)
    calculate_overlapRatio(overlapMatrix_bloodCelltype,
                           histonelist,
                           num_bloodCellType,
                           ratio_blood_savepath)
    print('OverlapRatio for blood celltypes saved in :',ratio_blood_savepath)
    calculate_overlapRatio(overlapMatrix_nonBlood,
                           histonelist,
                           num_nonBlood,
                           ratio_nonBlood_savepath)
    print('OverlapRatio for non-Blood celltypes saved in :',ratio_nonBlood_savepath)
    return None


if __name__ == '__main__':

    celltype_dir = '/home/shuang/projects/development_eqtm/data/features/celltype'
    blood_celltype_filepath = os.path.join(celltype_dir,'roadmap-celltypes-blood-MJ.txt')
    all_celltype_filepath = os.path.join(celltype_dir,'roadmap-celltypes.txt')

    histonelist_filepath = '/home/shuang/projects/development_eqtm/data/features/histonetype/2017-12-09-eQTLsFDR-et0.0-flipped_histoneList.txt'

    # calculate for SNP ratio
    overlapMatrix_dirpath = '/home/shuang/projects/development_eqtm/data/eqtmZscores/complete_overlapMatrix'
    overlapMatrix_dic2 = {
    'GenomicCoordinates_overlapMatrixcomplete.txt':
    'GenomicCoordinates'
    }

    # calculate for cpg ratio
    overlapMatrix_dic = {
    '2017-12-09-eQTLsFDR-et0.0-flipped_overlapMatrixcomplete.txt':
    '2017-12-09-eQTLsFDR-et0.0-flipped',
    '2017-12-09-eQTLsFDR-gt0.0-flipped_overlapMatrixcomplete.txt':
    '2017-12-09-eQTLsFDR-gt0.0-flipped',
    'random20k_gt0.5_overlapMatrixcomplete.txt':
    'random20k_gt0.5'
    }
    overlapRatio_save_dirpath = overlapMatrix_dirpath

    for filename in overlapMatrix_dic2.keys():
        overlapMatrix_name = overlapMatrix_dic2[filename]
        overlapMatrix_filepath = os.path.join(overlapMatrix_dirpath,
                                 filename)
        overlapMatrix_save_dirpath = overlapMatrix_dirpath
        overlapRatio_name = overlapMatrix_name+'_overlapRatiocomplete'

        select_bloodNonBloodCelltypes(blood_celltype_filepath,
                                      all_celltype_filepath,
                                      histonelist_filepath,
                                      overlapMatrix_name,
                                      overlapMatrix_filepath,
                                      overlapMatrix_save_dirpath,
                                      overlapRatio_name,
                                      overlapRatio_save_dirpath
                                      )
