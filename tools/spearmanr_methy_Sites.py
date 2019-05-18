import __init__path
import os
from lib.read.read_data import load_data
import numpy as np
from scipy.stats import spearmanr

if __name__ == "__main__":
    project_rootdir = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
    # project_rootdir = "/home/shuang/projects/development_eqtm"
    cpg_overlapRatio_filepath = os.path.join(project_rootdir,
                                             "data",
                                             "eqtmZscores",
                                             "withExpressionTSSMethyCpgOverlapGene",
                                             "2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_withGeneOverlap.txt")
    exclude = ['SNPChr', 'PValue', 'SNPChrPos', 'ProbeChr',
               'ProbeCenterChrPos', 'CisTrans', 'SNPType', 'AlleleAssessed',
               'DatasetsWhereSNPProbePairIsAvailableAndPassesQC',
               'DatasetsZScores', 'DatasetsNrSamples',
               'IncludedDatasetsMeanProbeExpression',
               'IncludedDatasetsProbeExpressionVariance', 'HGNCName',
               'IncludedDatasetsCorrelationCoefficient', 'Meta-Beta (SE)',
               'Beta (SE)', 'FoldChange', 'FDR', 'checkChr',
               'TssSite', 'chr', 'SNPName_ProbeName']
    keep = ['SNPName', 'ProbeName', 'OverallZScore']
    data = load_data(cpg_overlapRatio_filepath,
                     keep=keep,
                     exclude=exclude,
                     test_size=0)

    features = [col for col in data.train.values.columns if col not in
                ["SNPName", "ProbeName", "OverallZScore"]]

    positive_samples = data.train.values[features][data.train.values["OverallZScore"] > 0]
    negative_samples = data.train.values[features][data.train.values["OverallZScore"] < 0]

    # # calculate the spearmanr between positive and negative samples
    # pos_neg_spearmanr = np.zeros([positive_samples.shape[0], negative_samples.shape[0]])
    # for i1, row1 in enumerate(positive_samples.itertuples()):
    #     for i2, row2 in enumerate(negative_samples.itertuples()):
    #         pos_neg_spearmanr[i1][i2] = spearmanr(row1[1], row2[1])[0]
    # np.save(os.path.join(project_rootdir,
    #                      "corr_save",
    #                      "pos_neg_spearmanr.npy"),
    #         pos_neg_spearmanr)
    # print("Saved pos_neg_spearmanr to path: ", os.path.join(project_rootdir,
    #                                                "corr_save",
    #                                                "pos_neg_spearmanr.npy"),
    #       flush=True)
    #
    # pos_spearmanr = positive_samples.T.corr(method='spearman')
    # pos_spearmanr.to_csv(os.path.join(project_rootdir,
    #                                   "corr_save",
    #                                   "pos_spearmanr.csv"))
    # print("Saved pos_corr to path: ", os.path.join(project_rootdir,
    #                                                "corr_save",
    #                                                "neg_spearmanr.csv"),
    #       flush=True)

    neg_spearmanr = negative_samples.T.corr(method='spearman')
    neg_spearmanr.to_csv(os.path.join(project_rootdir,
                                      "corr_save",
                                      "neg_spearmanr.csv"))
    print("Saved neg_corr to path: ", os.path.join(project_rootdir,
                                      "corr_save",
                                      "neg_spearmanr.csv"),
          flush=True)
