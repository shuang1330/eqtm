import __init__path
import os
from lib.read.read_data import load_data
import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import pandas as pd

if __name__ == "__main__":
    # project_rootdir = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
    project_rootdir = "/home/shuang/projects/development_eqtm"
    cpg_overlapRatio_filepath = os.path.join(project_rootdir,
                                             "data",
                                             "eqtmZscores",
                                             "withExpressionTSSMethyCpgOverlapGene",
                                             "2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_withGeneOverlap.txt")
    et = "2017-12-09-eQTLsFDR-et0_withExpressionTssMethyOverlap_PromoterOverlap.txt"
    et_blood = "2017-12-09-eQTLsFDR-et0.0-flipped_withExpressionTssMethy_blood_withGeneOverlap_PromoterOverlap.txt"
    et_nonblood = "2017-12-09-eQTLsFDR-et0.0-flipped_withExpressionTssMethy_nonBlood_withGeneOverlap_PromoterOverlap.txt"
    gt = "2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_withGeneOverlap_PromoterOverlap.txt"
    gt_blood = "2017-12-09-eQTLsFDR-gt0.0-flipped_withExpressionTssMethy_blood_withGeneOverlap_PromoterOverlap.txt"
    gt_nonblood = "2017-12-09-eQTLsFDR-gt0.0-flipped_withExpressionTssMethy_nonBlood_withGeneOverlap_PromoterOverlap.txt"

    all_data = {"et":et,
                "gt":gt,
                "et_blood":et_blood,
                "et_nonblood":et_nonblood,
                "gt_blood":gt_blood,
                "gt_nonblood":gt_nonblood}

    exclude = ['SNPChr', 'PValue', 'SNPChrPos', 'ProbeChr',
               'ProbeCenterChrPos', 'CisTrans', 'SNPType', 'AlleleAssessed',
               'DatasetsWhereSNPProbePairIsAvailableAndPassesQC',
               'DatasetsZScores', 'DatasetsNrSamples',
               'IncludedDatasetsMeanProbeExpression',
               'IncludedDatasetsProbeExpressionVariance', 'HGNCName',
               'IncludedDatasetsCorrelationCoefficient', 'Meta-Beta (SE)',
               'Beta (SE)', 'FoldChange', 'FDR', 'checkChr',
               'TssSite', 'chr', 'SNPName_ProbeName', 'SNPName', 'ProbeName',
               'expressionMean', 'expressionVar', 'methyMean', 'methyVar']
    keep = ['OverallZScore']

    et_data = load_data(os.path.join(project_rootdir, "data",
                                     "eqtmZscores", "withExpressionTSSMethyCpgOverlapPromoter", et),
                        exclude=exclude, keep=keep, test_size=0)
    features = [feat for feat in et_data.train.values.columns if feat not in keep]

    new_blood_eqtm_filepath = os.path.join(project_rootdir,
                                           "data/eqtmZscores/withExpressionTSSMethyCpgOverlapPromoter/eqtm_fibroblast_allFeaturesAdded.txt")
    new_blood_eqtm = pd.read_csv(new_blood_eqtm_filepath, index_col=0)
    # tcga_ucec = load_data(os.path.join(project_rootdir,
    #                                    "data/eqtmZscores/withExpressionTSSMethyCpgOverlapPromoter/TCGA-UCEC-mval_all_features_added.txt"),
    #                       exclude=exclude, keep=keep, test_size=0)


    def tell_direction(item):
        return 1 if item > 0 else 0
    new_blood_eqtm["label"] = [tell_direction(item) for item in new_blood_eqtm["rvalue"]]
    print("New Blood Positive", new_blood_eqtm[new_blood_eqtm["label"] == 1].shape)
    print("New Blood Negative", new_blood_eqtm[new_blood_eqtm["label"] == 0].shape)

    # draw histogram
    for feature in features:
        # if feature.endswith("_promoter"):
        if feature:
            spear_Corr_original = spearmanr(et_data.train.values[feature], et_data.train.labels)
            spear_Corr_new = spearmanr(new_blood_eqtm[feature], new_blood_eqtm["label"])
            # spear_Corr_tcga_ucec = spearmanr(tcga_ucec.train.values[feature], tcga_ucec.train.labels)
            f = plt.figure()
            ax = f.add_subplot(111)
            min_x = min(new_blood_eqtm[feature].min(), et_data.train.values[feature].min())
            max_x = max(new_blood_eqtm[feature].max(), et_data.train.values[feature].max())
            bins = np.linspace(min_x, max_x, 100)
            n1, _, _ = ax.hist(et_data.train.values[feature],
                               bins=bins, alpha=0.2,
                               density=True, stacked=True,
                               label="Original")
            # n2, _, _ = ax.hist(new_blood_eqtm[feature],
            #                    bins=bins, alpha=0.2,
            #                    density=True, stacked=True,
            #                    label="NewBlood")
            n2, _, _ = ax.hist(new_blood_eqtm[feature],
                               bins=bins, alpha=0.2,
                               density=True, stacked=True,
                               label="TCGA_UCEC")
            plt.legend(loc="upper right")
            display_s = "{}\nOriginal:{}\nTCGA_UCEC:{}".format(feature,
                                                        round(spear_Corr_original[0], 2),
                                                        round(spear_Corr_new[0], 2))

            plt.text(0.1, 0.9, display_s,
                     horizontalalignment='left',
                     verticalalignment='center',
                     transform=ax.transAxes)

            plt.savefig("./development_eqtm/output_fig/new_blood_vs_original_et/{}.png".format(feature))
            plt.show()


    raise NotImplementedError
    print("TCGA UCEC Positive", tcga_ucec.train.values[tcga_ucec.train.labels==1].shape)
    print("TCGA UCEC Negative", tcga_ucec.train.values[tcga_ucec.train.labels==0].shape)

    # draw histogram
    for feature in features:
        # if feature.endswith("_promoter"):
        if feature:
            spear_Corr_original = spearmanr(et_data.train.values[feature], et_data.train.labels)
            # spear_Corr_new = spearmanr(new_blood_eqtm[feature], new_blood_eqtm["label"])
            spear_Corr_tcga_ucec = spearmanr(tcga_ucec.train.values[feature], tcga_ucec.train.labels)
            f = plt.figure()
            ax = f.add_subplot(111)
            min_x = min(tcga_ucec.train.values[feature].min(), et_data.train.values[feature].min())
            max_x = max(tcga_ucec.train.values[feature].max(), et_data.train.values[feature].max())
            bins = np.linspace(min_x, max_x, 100)
            n1, _, _ = ax.hist(et_data.train.values[feature],
                               bins=bins, alpha=0.2,
                               density=True, stacked=True,
                               label="Original")
            # n2, _, _ = ax.hist(new_blood_eqtm[feature],
            #                    bins=bins, alpha=0.2,
            #                    density=True, stacked=True,
            #                    label="NewBlood")
            n2, _, _ = ax.hist(tcga_ucec.train.values[feature],
                               bins=bins, alpha=0.2,
                               density=True, stacked=True,
                               label="TCGA_UCEC")
            plt.legend(loc="upper right")
            display_s = "{}\nOriginal:{}\nTCGA_UCEC:{}".format(feature,
                                                        round(spear_Corr_original[0], 2),
                                                        round(spear_Corr_tcga_ucec[0], 2))

            plt.text(0.1, 0.9, display_s,
                     horizontalalignment='left',
                     verticalalignment='center',
                     transform=ax.transAxes)

            plt.savefig("./development_eqtm/output_fig/tcga_ucec_vs_original_et/{}.png".format(feature))
            plt.show()

    raise NotImplementedError
    # pandas description
    # def read_return_description(filename):
    #     filepath = os.path.join(project_rootdir,
    #                             "data",
    #                             "eqtmZscores",
    #                             "withExpressionTSSMethyCpgOverlapGenePromoter",
    #                             filename)
    #     data = load_data(filepath, exclude=exclude + keep, keep=[], test_size=0)
    #     return data.train.values.describe()
    # # naively use describe for different datasets
    # describe_dic = {}
    # for file_prefix in all_data.keys():
    #     filename = all_data[file_prefix]
    #     describe_dic[file_prefix] = read_return_description(filename)

    # the following only examined gt dataset

    data = load_data(cpg_overlapRatio_filepath,
                     keep=keep,
                     exclude=exclude,
                     test_size=0)

    features = [col for col in data.train.values.columns if col not in
                ["SNPName", "ProbeName", "OverallZScore"]]

    positive_samples = data.train.values[features][data.train.values["OverallZScore"] > 0]
    negative_samples = data.train.values[features][data.train.values["OverallZScore"] < 0]

    # # draw scatter plot for overallzscore and each feature
    # for feature in features:
    #     spear_Corr = spearmanr(data.train.values[feature], data.train.values["OverallZScore"])
    #     f = plt.figure()
    #     ax = f.add_subplot(111)
    #     min_x = data.train.values[feature].min()
    #     max_x = data.train.values[feature].max()
    #     bins = np.linspace(min_x, max_x, num=100)
    #     ax.scatter(x=data.train.values[feature], y=data.train.values["OverallZScore"],
    #                alpha=0.5)
    #     ax.set_xlabel(xlabel=feature)
    #     ax.set_ylabel(ylabel="OverallZscore")
    #     display_s = "Spearman Corr:{}\nSpearman P:{}".format(spear_Corr[0], spear_Corr[1])
    #     ax.text(0.1, 0.9, display_s, horizontalalignment='left',
    #             verticalalignment='center', transform=ax.transAxes)
    #     plt.savefig("./development_eqtm/output_fig/unitvariate_test/scatter/{}_Zscore_scatter.png")
    #     plt.show()

    raise NotImplementedError
    # draw histogram
    for feature in features:
        if feature.endswith("_promoter"):
            spear_Corr = spearmanr(data.train.values[feature], data.train.labels)
            f = plt.figure()
            ax = f.add_subplot(111)
            min_x = data.train.values[feature].min()
            max_x = data.train.values[feature].max()
            bins = np.linspace(min_x, max_x, 100)
            n1, _, _ = ax.hist(positive_samples[feature],
                               bins=bins, alpha=0.2,
                               density=True, stacked=True,
                               label="Positive")
            n2, _, _ = ax.hist(negative_samples[feature],
                               bins=bins, alpha=0.2,
                               density=True, stacked=True,
                               label="Negative")
            plt.legend(loc="upper right")
            display_s = "{}\ncorr:{}\npvalue:{}".format(feature,
                                                        round(spear_Corr[0], 2),
                                                        round(spear_Corr[1], 2))

            plt.text(0.1, 0.9, display_s,
                     horizontalalignment='left',
                     verticalalignment='center',
                     transform=ax.transAxes)

            plt.savefig("./development_eqtm/output_fig/unitvariate_test/{}.png".format(feature))
            plt.show()

    def autolabel(rects, value):
        """
        Attach a text label above each bar displaying its height
        """
        for ind, rect in enumerate(rects):
            height = round(value[ind], 2)
            ax.text(rect.get_x() + rect.get_width() / 2., 0,
                    '%d' % int(height),
                    ha='center', va='bottom').set_fontsize(2)
        return None

    spear_corr = []
    spear_p = []
    for feature in features:
        corr = spearmanr(data.train.values[feature], data.train.labels)
        spear_corr.append(corr[0])
        spear_p.append(corr[1])
    f = plt.figure()
    ax = f.add_subplot(111)
    rect = ax.bar(features, spear_corr)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(5)
        tick.label.set_rotation("vertical")
    plt.savefig("./development_eqtm/output_fig/spearman_corr_for_all_features.png")
    autolabel(rect, spear_p)
    plt.show()

    save_filepath = os.path.join(project_rootdir,
                                 "data",
                                 "features",
                                 "feature_importances",
                                 "spearcor.txt")
    with open(save_filepath, "w") as f:
        for ind in range(len(features)):
            f.write("{}\t{}\t{}\n".format(features[ind],
                                        spear_corr[ind],
                                        spear_p[ind]))
        f.close()



