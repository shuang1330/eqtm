import __init__path
import os
import pandas as pd, seaborn as sns
from lib.read.read_data import load_data
import numpy as np
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

sns.set(font="monospace")

# project_rootdir = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/"
project_rootdir = "/home/shuang/projects/development_eqtm"
overlap_ratio_filepath = os.path.join(project_rootdir,
                                      "data/eqtmZscores",
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
data = load_data(overlap_ratio_filepath,
                 keep=keep,
                 exclude=exclude,
                 test_size=0)

names = ["SNPName", "ProbeName"]
additional_features = ["OverallZScore", "expressionMean", "expressionVar",
                       "TssDistance", "methyMean", "methyVar"]
all_features = [col for col in data.train.values.columns
                if col not in additional_features and col not in names]
histone_features = [col for col in all_features if not col.endswith("gene")]
data_features = data.train.values[histone_features]

all_zeros_cpgsites = {}
for row in data_features.iterrows():
    if len(pd.unique(row[1].values)) == 1:
        print(row[0], data.train.values.loc[row[0]].values)
        all_zeros_cpgsites[row[0]] = data.train.values[additional_features].loc[row[0]]


def detect_sign(value):
    if value < 0:
        return 0
    else:
        return 1


all_additional_features = pd.DataFrame.from_dict(all_zeros_cpgsites, orient="index")
all_additional_features["direction"] = [detect_sign(row) for row in all_additional_features["OverallZScore"].values]

for feature in [col for col in all_additional_features.columns if col not in ["SNPName", "ProbeName"]]:
    min_x = all_additional_features[feature].min()
    max_x = all_additional_features[feature].max()
    spear = spearmanr(all_additional_features[feature], all_additional_features["direction"])
    corr = spear[0]
    p = spear[1]
    bins = np.linspace(min_x, max_x, 50)
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.hist(all_additional_features[feature][all_additional_features["OverallZScore"] > 0],
            bins=bins, normed=True, stacked=True, alpha=0.2,
            label="Positive")
    ax.hist(all_additional_features[feature][all_additional_features["OverallZScore"] < 0],
            bins=bins, normed=True, stacked=True, alpha=0.2,
            label="Negative")
    ax.legend(loc="upper right")
    display_s = "{}\ncorr:{}\npvalue:{}".format(feature,
                                                round(corr, 2),
                                                round(p, 2))
    ax.set_xlabel(feature)

    plt.text(0.1, 0.9, display_s,
             horizontalalignment='left',
             verticalalignment='center',
             transform=ax.transAxes)

    plt.savefig("./development_eqtm/output_fig/all_zeros_cpgsites/{}.png".format(str(feature)))
    plt.show()

# TODO: find nearby cpgsites

positive_histone_features = data.train.values[histone_features][data.train.values.OverallZScore > 0]
negative_histone_features = data.train.values[histone_features][data.train.values.OverallZScore < 0]

corr_positive = positive_histone_features.T.corr()
corr_negative = negative_histone_features.T.corr()

sns.clustermap(corr_positive)
plt.savefig("./development_eqtm/output_fig/corr_positive.png")
plt.show()
sns.clustermap(corr_negative)
plt.savefig("./development_eqtm/output_fig/corr_negative.png")
plt.show()
