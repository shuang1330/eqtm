import os
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
import pandas as pd
from sklearn.model_selection import train_test_split
import numpy as np
import matplotlib.pyplot as plt


def measure_auc(classifier, test_data_values, test_data_labels):
    prods = classifier.predict_proba(test_data_values)[:, 1]
    fpr, tpr, _ = metrics.roc_curve(test_data_labels, prods)
    score = metrics.auc(fpr, tpr)
    return score


def tell_direction(item):
    if item > 0:
        return 1
    else:
        return 0


if __name__ == '__main__':
    # locally
    # project_rootdir = '/home/shuang/projects/development_eqtm'
    # calculon project root dir
    project_rootdir = '/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm'
    bios_part2_path = os.path.join(project_rootdir, "data",
                                   "eqtmZscores", "withExpressionTSSMethyCpgOverlapPromoter",
                                   "2017-12-09-eQTLsFDR-et0_withExpressionTssMethyOverlap_PromoterOverlap.txt")
                                   # "BIOS_resplit_Part2_flipped_all_features_added.txt")
                                   # "BIOS_Part2_all_features_added.txt")
    bios_part1_path = os.path.join(project_rootdir, "data",
                                   "eqtmZscores", "withExpressionTSSMethyCpgOverlapPromoter",
                                   "meta_analysis_tcga_flipped_all_features_added.txt")
                                   # "BIOS_resplit_Part2_flipped_all_features_added.txt")

    bios_part2 = pd.read_csv(bios_part2_path, index_col=0)
    if "flippedZscore" not in bios_part2.columns:
        bios_part2["flippedZscore"] = bios_part2["OverallZScore"]
    bios_part2["eqtm_combi"] = ["{}_{}".format(item[0], item[1]) for item in
                                zip(bios_part2["SNPName"], bios_part2["ProbeName"])]
    bios_part2["label"] = [tell_direction(item) for item in bios_part2["flippedZscore"]]
    bios_part2_eqtms_names = bios_part2["eqtm_combi"].unique()
    exclude_columns = ['SNPName', 'SNPChr', 'PValue', 'SNPChrPos', 'ProbeName', 'ProbeChr',
               'ProbeCenterChrPos', 'CisTrans', 'SNPType', 'AlleleAssessed',
               'DatasetsWhereSNPProbePairIsAvailableAndPassesQC',
               'DatasetsZScores', 'DatasetsNrSamples',
               'IncludedDatasetsMeanProbeExpression',
               'IncludedDatasetsProbeExpressionVariance', 'HGNCName',
               'IncludedDatasetsCorrelationCoefficient', 'Meta-Beta (SE)',
               'Beta (SE)', 'FoldChange',  'checkChr', 'SNPName_ProbeName', 'chr', 'TssSite', "sign",
                       "OverallZScore", "FDR",
                       "expressionMean", "expressionVar",
                       "methyMean", "methyVar", "H3K9me3",
                       "H3K9me1", "H3K27me3", "H3K36me3",
                       "H4K20me1", "flippedZscore", "label", "eqtm_combi"]
    bios_part1 = pd.read_csv(bios_part1_path, index_col=0)
    if "flippedZscore" not in bios_part1.columns:
        bios_part1["flippedZscore"] = bios_part1["OverallZScore"]
    bios_part1["label"] = [tell_direction(item) for item in bios_part1["flippedZscore"]]
    bios_part1["eqtm_combi"] = ["{}_{}".format(item[0], item[1]) for item in
                                zip(bios_part1["SNPName"], bios_part1["ProbeName"])]


    def isCommon(item):
        return True if item in bios_part2_eqtms_names else False
    bios_part1["isCommon"] = [isCommon(item) for item in bios_part1["eqtm_combi"]]

    features = [item for item in bios_part2.columns
                if item not in exclude_columns and not item.endswith("promoter")]
    print(features)
    # features = ['TssDistance', 'H4K5ac', 'H2A.Z', 'H2BK120ac', 'H3K79me2', 'H3K27ac',
    #             'H2BK20ac', 'H3K14ac', 'H3K9ac', 'H3K4ac', 'H2AK5ac',
    #             'H3K4me1', 'H3K18ac', 'H3K23ac', 'H2BK5ac', 'H3K4me3',
    #             'H2BK12ac', 'H3K23me2', 'H4K12ac', 'DNase', 'H2BK15ac',
    #             'H3K4me2', 'H3K79me1', 'H2AK9ac', 'H3T11ph', 'H4K8ac',
    #             'H4K91ac', 'H3K56ac']
    print("Length of train data features: ", len(features))
    assert "OverallZScore" not in features
    model = RandomForestClassifier(n_estimators=50)
    train_x, test_x, train_y, test_y = train_test_split(bios_part2[features], bios_part2["label"], random_state=42)
    print(train_y.head())
    model.fit(train_x, train_y)
    train_score = measure_auc(model, test_x, test_y)
    print("Train score", train_score)

    # all test data
    test_score = measure_auc(model, bios_part1[features], bios_part1["label"])
    print("Test score", test_score)

    # bios part 2 data
    test_score = measure_auc(model, bios_part1[features][bios_part1["isCommon"] == False],
                             bios_part1["label"][bios_part1["isCommon"] == False])
    print(" Take all samples seenin the training set out. \nTest score", test_score)


    # draw feature importances
    importances = model.feature_importances_
    std = np.std([tree.feature_importances_ for tree in model.estimators_],
                 axis=0)
    indices = np.argsort(importances)[::-1]
    feature_names = [features[ind] for ind in indices]
    print(feature_names)

    # Print the feature ranking
    print("Feature ranking:")

    for f in range(train_x.shape[1]):
        print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))

    # Plot the feature importances of the forest
    plt.figure()
    plt.title("Feature importances")
    plt.bar(range(train_x.shape[1]), importances[indices],
            color="r", yerr=std[indices], align="center", alpha=0.5)
    plt.xticks(range(train_x.shape[1]), feature_names, rotation='vertical')

    plt.xlim([-1, train_x.shape[1]])
    plt.tight_layout()
    plt.savefig("/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/output_fig/feature_importances.png")
