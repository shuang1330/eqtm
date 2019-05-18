import __init__path
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import pickle
from seaborn import clustermap
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist

if __name__ == "__main__":
    # project_rootdir = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
    # locally
    project_rootdir = "/home/shuang/projects/development_eqtm"
    exclude = ['SNPName', 'SNPChr', 'PValue', 'SNPChrPos', 'ProbeName', 'ProbeChr',
               'ProbeCenterChrPos', 'CisTrans', 'SNPType', 'AlleleAssessed',
               'DatasetsWhereSNPProbePairIsAvailableAndPassesQC',
               'DatasetsZScores', 'DatasetsNrSamples',
               'IncludedDatasetsMeanProbeExpression',
               'IncludedDatasetsProbeExpressionVariance', 'HGNCName',
               'IncludedDatasetsCorrelationCoefficient', 'Meta-Beta (SE)',
               'Beta (SE)', 'FoldChange', 'FDR', 'checkChr', 'SNPName_ProbeName',
               'TssSite', 'chr', 'OverallZScore']
    data_filepath = os.path.join(project_rootdir,
                                 "data",
                                 "eqtmZscores",
                                 "withExpressionTSSMethyCpgOverlapPromoter",
                                 "et_gt_all_celltypes.txt")
    data = pd.read_csv(data_filepath, sep=",")
    features = [feature for feature in data.columns if feature not in exclude]
    for feature in features:
        if data[feature].describe()["std"] == 0:
            print(feature)
    # hierarchical clustering on histone markers for positive eQTMs
    data_positive = data[data["OverallZScore"] > 0]
    data_negative = data[data["OverallZScore"] < 0]
    # hierarchical clustering
    Z = linkage(data[features].values, 'ward')
    print(Z)
    c, coph_dists = cophenet(Z, pdist(data[features].values))
    den = dendrogram(
        Z,
        truncate_mode='lastp',  # show only the last p merged clusters
        p=12,  # show only the last p merged clusters
        show_leaf_counts=False,  # otherwise numbers in brackets are counts
        leaf_rotation=90,
        leaf_font_size=12,
        show_contracted=True  # to get a distribution impression in truncated branches
    )

    raise NotImplementedError
    # seaborn draw clustermaps
    g_positive = clustermap(data_positive[features].corr())
    plt.savefig(os.path.join(project_rootdir, "output_fig", "cluster_map_positive.png"))
    g_negative = clustermap(data_negative[features].corr())
    plt.savefig(os.path.join(project_rootdir, "output_fig", "cluster_map_negative.png"))
    g = clustermap(data[features].corr())
    plt.savefig(os.path.join(project_rootdir, "output_fig", "cluster_map_histones.png"))


    raise NotImplementedError
    # draw model trained on different cell types feature importances
    # load data
    data_filepath = os.path.join(project_rootdir,
                                 "data", "eqtmZscores",
                                 "withExpressionTSSMethyCpgOverlapPromoter",
                                 "2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_PromoterOverlap.txt")
    # load model
    model_dirpath = os.path.join(project_rootdir, "model")
    model_names = {"Blood cell types": "2017-12-09-eQTLsFDR-et0.0-flipped_withExpressionTssMethy_blood_PromoterOverlap.txt.pkl",
                   "Non-blood cell types": "2017-12-09-eQTLsFDR-et0.0-flipped_withExpressionTssMethy_nonBlood_PromoterOverlap.txt.pkl",
                   "All cell types": "2017-12-09-eQTLsFDR-et0_withExpressionTssMethyOverlap_PromoterOverlap.txt.pkl"}

    # load data and find feature names
    def find_feature_names():
        input_data = pd.read_csv(data_filepath, index_col=0)
        print(input_data.head())
        return [item for item in input_data.columns if item not in exclude]

    data_feature_names = find_feature_names()

    # draw feature importances
    save_fig_filepath = os.path.join(project_rootdir, "output_fig", "feature_importances_celltype.png")
    fig = plt.gcf()
    fig.set_size_inches(18.5, 12.5)
    for label in model_names.keys():
        model = pickle.load(open(os.path.join(model_dirpath, model_names[label]), "rb"))
        importances = model.best_estimator_.steps[0][1].feature_importances_
        # print(model.best_estimator_.steps[0][1])
        # print(importances)
        # print(importances.shape)
        # print(data_feature_names)
        assert len(importances) == len(data_feature_names)
        plt.bar(data_feature_names, importances, alpha=0.2, label=label)
    plt.xticks(rotation=90)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_fig_filepath, dpi=100)

    print("Fin.")
