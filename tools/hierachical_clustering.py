"""
clustering for features
"""
import __init__path
import os
import numpy as np
# from scipy.cluster.hierarchy import dendrogram, linkage
# from  scipy.cluster.hierarchy import cophenet
# from scipy.spatial.distance import pdist
from lib.read.read_data import load_data
from seaborn import clustermap
# import matplotlib as mpl
# mpl.use('Agg')
# import matplotlib.pyplot as plt


if __name__ == "__main__":
    project_rootdir = "/home/shuang/projects/development_eqtm"
    # cluster
    # project_rootdir = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
    eqtm_dirpath = os.path.join(project_rootdir,
                                "data",
                                "eqtmZscores",
                                "withExpressionTSSMethyCpgOverlapGenePromoter")
    eqtm_names = os.listdir(eqtm_dirpath)
    eqtm_filepath = os.path.join(eqtm_dirpath, "")
    exclude = ['SNPName', 'SNPChr', 'PValue', 'SNPChrPos', 'ProbeName', 'ProbeChr',
               'ProbeCenterChrPos', 'CisTrans', 'SNPType', 'AlleleAssessed',
               'DatasetsWhereSNPProbePairIsAvailableAndPassesQC',
               'DatasetsZScores', 'DatasetsNrSamples',
               'IncludedDatasetsMeanProbeExpression',
               'IncludedDatasetsProbeExpressionVariance', 'HGNCName',
               'IncludedDatasetsCorrelationCoefficient', 'Meta-Beta (SE)',
               'Beta (SE)', 'FoldChange', 'FDR', 'checkChr', 'SNPName_ProbeName',
               'TssSite', 'chr']
    keep = ['OverallZScore']
    input_data = load_data(eqtm_filepath, exclude=exclude, keep=keep, normalization=False)

