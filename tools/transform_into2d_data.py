import __init__path
import os
import numpy as np
import pandas as pd
from lib.preprocess.selectCellType import (read_allCelltypeFile,
                                           read_bloodCelltypeFile,
                                           read_histoneList)
from lib.read.read_for_cnn import read_individual_image
import matplotlib.pyplot as plt


def read_feature_importances(feature_importances_filepath):
    filecontent = open(feature_importances_filepath, 'r').readlines()
    print(filecontent[:2])
    features = [item.strip() for item in filecontent[0].split('\t')]
    importances = [float(item.strip()) for item in filecontent[1].split('\t')]
    return features, importances


def read_spearcor(spear_cor_filepath):
    with open(spear_cor_filepath, "r") as f:
        all_content = f.readlines()
        features = [line.split("\t")[0] for line in all_content]
        importances = [float(line.split("\t")[1]) for line in all_content]
        f.close()
    return features, importances


def sort_feature_by_importances(feature_importances_filepath,
                                feature_importances_type,
                                feature_type_name):
    if feature_importances_type == "importances":
        features, importances = read_feature_importances(feature_importances_filepath)
    elif feature_importances_type == "spearcor":
        features, importances = read_spearcor(feature_importances_filepath)
    else:
        raise IOError("Needs type, options being either importances or spearcor")
    if feature_type_name == 'gene':
        histones_overlap = []
        histones_importances = []
        for ind,featurename in enumerate(features):
            if featurename.endswith('gene'):
                if featurename not in ['expressionMean', 'expressionVar',
                                       'TssDistance', 'methyMean', 'methyVar']:
                    histones_overlap.append(featurename[:-4])
                    histones_importances.append(importances[ind])
        histone_sorted = [histones_overlap[ind] for ind in np.argsort(np.array(histones_importances))[::-1]]
    elif feature_type_name == 'cpg':
        histones_overlap = []
        histones_importances = []
        for ind,featurename in enumerate(features):
            if not featurename.endswith('gene'):
                if featurename not in ['expressionMean', 'expressionVar',
                                       'TssDistance', 'methyMean', 'methyVar']:
                    histones_overlap.append(featurename)
                    histones_importances.append(importances[ind])
        histone_sorted = [histones_overlap[ind] for ind in np.argsort(np.array(histones_importances))[::-1]]
    return histone_sorted


def rearrange_cpgs(overlapMatrix_filepath,
                   type_name,
                   blood_celltype_filepath,
                   all_celltype_filepath,
                   histonelist_filepath,
                   feature_importances_filepath,
                   feature_importances_type,
                   save_dirpath):
    blood_celltype = read_bloodCelltypeFile(blood_celltype_filepath)
    all_celltype = read_allCelltypeFile(all_celltype_filepath)
    histones = read_histoneList(histonelist_filepath)
    sorted_histone = sort_feature_by_importances(feature_importances_filepath,
                                                 feature_importances_type,
                                                 type_name)
    sorted_features = []
    for celltype in blood_celltype:
        for histone in sorted_histone:
            feature_name = '{}-{}'.format(celltype, histone)
            sorted_features.append(feature_name)
    for celltype in [item for item in all_celltype.keys() if item not in blood_celltype]:
        for histone in sorted_histone:
            feature_name = '{}-{}'.format(celltype, histone)
            sorted_features.append(feature_name)

    # # save the sorted cell types to filepath,
    # with open("/home/shuang/projects/development_eqtm/data/clustermap_meta/sorted_celltypes.txt") as f:
    #     for celltype in blood_celltype:
    #         f.write(str(celltype))
    #         f.write("\n")
    #     for celltype in [item for item in all_celltype.keys() if item not in blood_celltype]:
    #         f.write(str(celltype))
    #         f.write("\n")
    #     f.close()
    # with open("/home/shuang/projects/development_eqtm/data/clustermap_meta/sorted_histones.txt") as f:
    #     for histone in sorted_histone:
    #         f.write(str(histone))
    #         f.write("\n")
    #     f.close()
    # raise NotImplementedError

    overlapMatrix = pd.read_csv(overlapMatrix_filepath, sep=',', index_col=0)
    print('Loaded overlapMatrix.', flush=True)
    overlapMatrix['E005-H4K12ac'] = 0
    overlapMatrix['E008-H3T11ph'] = 0
    overlapMatrix['E017-H2AK9ac'] = 0
    num_histones = len(histones)
    num_all_celltypes = len(all_celltype.keys())

    content = overlapMatrix[sorted_features].values
    for ind, cpg_name in enumerate(overlapMatrix.index):
        print(cpg_name, flush=True)
        cpg_save_filepath = os.path.join(save_dirpath, cpg_name + '.txt')
        cpg_file = open(cpg_save_filepath, 'w')
        for i in range(0, content.shape[1], num_histones):
            line_to_write = content[ind, i:i + num_histones]
            cpg_file.write('\t'.join([str(item) for item in line_to_write]))
            cpg_file.write('\n')
        cpg_file.close()


if __name__=='__main__':
    # locally
    project_rootdir = '/home/shuang/projects/development_eqtm'
    # calculon
    # project_rootdir = '/groups/umcg-gcc/tmp04/umcg-sli/development_eqtm'
    feature_importances_filepath = os.path.join(project_rootdir,
                                                'data',
                                                'features',
                                                'feature_importances',
                                                'spearcor.txt')
    celltype_dir = os.path.join(project_rootdir, 'data', "features", 'celltype')
    blood_celltype_filepath = os.path.join(celltype_dir, 'roadmap-celltypes-blood-MJ.txt')
    all_celltype_filepath = os.path.join(celltype_dir, 'roadmap-celltypes.txt')
    histonelist_filepath = os.path.join(project_rootdir, 'data', 'features', 'histonetype',
                                        '2017-12-09-eQTLsFDR-et0.0-flipped_histoneList.txt')

    # cpgSites
    cpg_overlapMatrix_filepath = os.path.join(project_rootdir,
                                              'data',
                                              'eqtmZscores',
                                              'complete_overlapMatrix',
                                              '2017-12-09-eQTLsFDR-et0.0-flipped_overlapMatrixcomplete.txt')
    cpg_save_dirpath = os.path.join(project_rootdir,
                                    'data',
                                    'cpgSites',
                                    'et')

    gene_overlapMatrix_filepath = os.path.join(project_rootdir,
                                               'data',
                                               'features',
                                               'geneOverlap',
                                               'gene_StartEndSite_overlapMatrix_complete.txt')
    gene_save_dirpath = os.path.join(project_rootdir,
                                     'data',
                                     'features',
                                     'geneOverlap',
                                     'seperate_geneFiles')

    # rearrange_cpgs(cpg_overlapMatrix_filepath,
    #                'cpg',
    #                blood_celltype_filepath,
    #                all_celltype_filepath,
    #                histonelist_filepath,
    #                feature_importances_filepath,
    #                "spearcor",
    #                cpg_save_dirpath)

    # rearrange_cpgs(gene_overlapMatrix_filepath,
    #                'gene',
    #                blood_celltype_filepath,
    #                all_celltype_filepath,
    #                histonelist_filepath,
    #                feature_importances_filepath,
    #                gene_save_dirpath)

    # save them in img format
    eqtm_filepath = os.path.join(project_rootdir,
                                 "data",
                                 "eqtmZscores",
                                 "ORIGIN",
                                 "2017-12-09-eQTLsFDR-et0.0-flipped.txt")
    eqtm = pd.read_csv(eqtm_filepath, sep="\t")
    positive_cpgnames = eqtm["SNPName"][eqtm["OverallZScore"] > 0]
    negative_cpgnames = eqtm["SNPName"][eqtm["OverallZScore"] < 0]

    def save_in_imgFormat(matrix_dirpath,
                          cpgnames,
                          save_dirpath):
        for cpgname in cpgnames:
            cpg_filepath = os.path.join(matrix_dirpath,
                                           cpgname+".txt")
            content = read_individual_image(cpg_filepath)
            img = content.reshape([content.shape[0], content.shape[1]])
            plt.imshow(img)
            plt.savefig(os.path.join(save_dirpath,"{}.png".format(cpgname)))

    matrix_dirpath = os.path.join(project_rootdir,
                                  "data",
                                  "cpgSites",
                                  "et")
    positive_save_dirpath = os.path.join(project_rootdir,
                                         "data",
                                         "cpgSites",
                                         "et_pos_img")
    negative_save_dirpath = os.path.join(project_rootdir,
                                         "data",
                                         "cpgSites",
                                         "et_neg_img")
    save_in_imgFormat(matrix_dirpath, positive_cpgnames, positive_save_dirpath)
    print("Processed all positive samples.")
    save_in_imgFormat(matrix_dirpath, negative_cpgnames, negative_save_dirpath)
