import __init__path
import numpy as np
import pandas as pd
import os
from seaborn import clustermap
from scipy.cluster.hierarchy import linkage
from lib.preprocess.selectCellType import read_allCelltypeFile
from time import time
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def read_histone_list(histone_filepath):
    with open(histone_filepath, "r") as f:
        histone_list = [item.strip() for item in f.readlines()]
        f.close()
    return histone_list


def calculate_heatmap_matrix(dataset_path, histone_list, cell_type_list):
    # timeit
    startTime = time()
    dataset= pd.read_csv(dataset_path, sep=",", index_col=0)
    print("Dataset loaded. Time: ", time()-startTime)
    heatmap_matrix = pd.DataFrame(data=0, index=cell_type_list, columns=histone_list)
    for histone_ind in tqdm(range(len(histone_list))):
        histone = histone_list[histone_ind]
        for cell_type in cell_type_list:
            if cell_type:
                feature = "{}-{}".format(cell_type, histone)
                heatmap_matrix.loc[cell_type, histone] = dataset[feature].mean()
    return heatmap_matrix


if __name__=='__main__':
    project_rootdir = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
    # project_rootdir = "/home/shuang/projects/development_eqtm"
    data_dir = os.path.join(project_rootdir, "data/eqtmZscores/allCpgs")
    heatmaps_dir = os.path.join(project_rootdir, "data/heatmaps")
    origin_name = "eqtm_FDR_smaller_than_0.05_bedtoolsFormat_overlapMatrix.txt"
    new_blood_name = "cpg_fibroblast_overlapMatrix_complete.txt"
    tcga_name = "cpg_TCGA_overlapMatrix_complete.txt"
    cell_type_list = read_allCelltypeFile(os.path.join(project_rootdir,
                                                       "data/features/celltype/roadmap-celltypes-blood-MJ.txt")).keys()
    print(cell_type_list)
    histone_list = read_histone_list(os.path.join(project_rootdir,
                                                  "data/features/histonetype/2017-12-09-eQTLsFDR-et0.0-flipped_histoneList.txt"))

    # for dataset_name in [origin_name, new_blood_name, tcga_name]:
    #     dataset_path = os.path.join(data_dir, dataset_name)
    #     heatmap = calculate_heatmap_matrix(dataset_path, histone_list, cell_type_list)
    #     heatmap.to_csv(os.path.join(heatmaps_dir, "{}.txt".format(dataset_name[:-4])))

    origin_heatmap = pd.read_csv(os.path.join(heatmaps_dir, origin_name[:-4]+".txt"), index_col=0)
    new_blood_heatmap = pd.read_csv(os.path.join(heatmaps_dir, new_blood_name[:-4] + ".txt"), index_col=0)
    tcga_heatmap = pd.read_csv(os.path.join(heatmaps_dir, tcga_name[:-4] + ".txt"), index_col=0)

    row_linkage = linkage(origin_heatmap.values, 'ward')
    col_linkage = linkage(origin_heatmap.values.T, 'ward')
    # cell_type_color = pd.DataFrame(data='azure', index=origin_heatmap.index, columns=["is_Blood_Cell"])
    # for celltype in blood_celltype:
    #     cell_type_color[cell_type_color.index==celltype] = 'c'
    g = clustermap(origin_heatmap,
                   col_linkage=col_linkage,
                   row_linkage=row_linkage)
    # cell_type_legends = [Line2D([0], [0], color='azure', lw=4),
    #                      Line2D([0], [0], color='c', lw=4)]
    # plt.legend(cell_type_legends, ['blood cells', 'non blood cells'])
    plt.savefig(os.path.join(project_rootdir, "output_fig/heatmaps/origin.png"))
    g2 = clustermap(new_blood_heatmap,
                   col_linkage=col_linkage,
                   row_linkage=row_linkage)
    plt.savefig(os.path.join(project_rootdir, "output_fig/heatmaps/new_blood.png"))
    g3 = clustermap(tcga_heatmap,
                   col_linkage=col_linkage,
                   row_linkage=row_linkage)
    plt.savefig(os.path.join(project_rootdir, "output_fig/heatmaps/tcga.png"))