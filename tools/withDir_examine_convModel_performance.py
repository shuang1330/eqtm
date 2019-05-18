import __init__path
from lib.read.read_data import load_data
import os
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import permutation_test_score
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
from sklearn import metrics
import pandas as pd


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
    eqtm_path = os.path.join(project_rootdir, "data",
                             "eqtmZscores", "withExpressionTSSMethyCpgOverlapPromoter",
                             "2017-12-09-eQTLsFDR-et0_withExpressionTssMethyOverlap_PromoterOverlap.txt")
    tcga_path = os.path.join(project_rootdir,  "data",
                             "eqtmZscores", "withExpressionTSSMethyCpgOverlapPromoter",
                             "TCGA-BRCA-mval_all_features_added.txt")
    bios_part1_path = os.path.join(project_rootdir, "data",
                                   "eqtmZscores", "withExpressionTSSMethyCpgOverlapPromoter",
                                   "BIOS_Part1_all_features_added.txt")
    bios_part2_path = os.path.join(project_rootdir, "data",
                                   "eqtmZscores", "withExpressionTSSMethyCpgOverlapPromoter",
                                   "BIOS_Part2_all_features_added.txt")
    model_name = 'ranfor'
    exclude = ['SNPName', 'SNPChr', 'PValue', 'SNPChrPos', 'ProbeName', 'ProbeChr',
               'ProbeCenterChrPos', 'CisTrans', 'SNPType', 'AlleleAssessed',
               'DatasetsWhereSNPProbePairIsAvailableAndPassesQC',
               'DatasetsZScores', 'DatasetsNrSamples',
               'IncludedDatasetsMeanProbeExpression',
               'IncludedDatasetsProbeExpressionVariance', 'HGNCName',
               'IncludedDatasetsCorrelationCoefficient', 'Meta-Beta (SE)',
               'Beta (SE)', 'FoldChange',  'checkChr', 'SNPName_ProbeName', 'chr', 'TssSite', "sign"]
    keep = ['OverallZScore', 'FDR']
    train_path = eqtm_path
    data_names = ["TCGA-BLCA-mval", "TCGA-BRCA-mval", "TCGA-HNSC-mval",
                  "TCGA-PRAD-mval", "TCGA-THCA-mval", "TCGA-UCEC-mval",
                  "TCGA-CHOL-mval", "TCGA-KIRC-mval"]

    train_data = load_data(train_path, exclude=exclude, keep=keep)
    train_columns = train_data.train.values.columns
    features = [item for item in train_data.train.values.columns
                if item not in ["OverallZScore", "FDR",
                                "expressionMean", "expressionVar",
                                "methyMean", "methyVar", "H3K9me3",
                                "H3K9me1", "H3K27me3", "H3K36me3",
                                "H4K20me1"] and not item.endswith("promoter")]
    # features = ['TssDistance', 'H4K5ac', 'H2A.Z', 'H2BK120ac', 'H3K79me2', 'H3K27ac',
    #             'H2BK20ac', 'H3K14ac', 'H3K9ac', 'H3K4ac', 'H2AK5ac',
    #             'H3K4me1', 'H3K18ac', 'H3K23ac', 'H2BK5ac', 'H3K4me3',
    #             'H2BK12ac', 'H3K23me2', 'H4K12ac', 'DNase', 'H2BK15ac',
    #             'H3K4me2', 'H3K79me1', 'H2AK9ac', 'H3T11ph', 'H4K8ac',
    #             'H4K91ac', 'H3K56ac']
    print("Length of train data features: ", len(features))
    print(features)
    assert "OverallZScore" not in features
    model = RandomForestClassifier(n_estimators=50)
    model.fit(train_data.train.values[features], train_data.train.labels)
    train_score = measure_auc(model, train_data.test.values[features], train_data.test.labels)
    print("CV score", train_score)


    # bios part 2 data
    bios_part1 = load_data(bios_part2_path, exclude=exclude, keep=keep)



    # liver dataset
    liver_eqtms_filepath = os.path.join(project_rootdir,
                                        "data","eqtmZscores","withExpressionTSSMethyCpgOverlapPromoter",
                                        "lung_eqtms_with_all_features.txt")
    liver_eqtms = pd.read_csv(liver_eqtms_filepath, sep=",", index_col=0)
    liver_eqtms["label"] = [tell_direction(item) for item in liver_eqtms["Spearmanr"]]
    new_eqtm_score = measure_auc(model, liver_eqtms[features], liver_eqtms["label"])
    print("Liver eqtm score", new_eqtm_score)
    # 0.629736354001

    # new blood eqtm dataset
    new_blood_eqtm = pd.read_csv(
        os.path.join(project_rootdir,
                     "data","eqtmZscores","withExpressionTSSMethyCpgOverlapPromoter",
                     "eqtm_fibroblast_allFeaturesAdded.txt"),
        index_col=0
    )
    no_features_in_newData = ["METHYL_label", "GENE_ID,METHYL_chromosome",
                              "GENE_chromosome", "METHYL_location", "GENE_location",
                              "rvalue", "pvalue", "-log10(pvalue)"]

    new_blood_eqtm["direction"] = [tell_direction(item) for item in new_blood_eqtm["rvalue"].values]
    new_eqtm_score = measure_auc(model, new_blood_eqtm[features], new_blood_eqtm["direction"])
    print("New blood eqtm score", new_eqtm_score)
    # 0.651939146843


    # TCGA test dataset
    for filename in os.listdir(os.path.join(project_rootdir, "data", "eqtmZscores", "withExpressionTSSMethyCpgOverlapPromoter")):
        if filename.startswith("TCGA"):
            tcga_filepath = os.path.join(project_rootdir, "data", "eqtmZscores", "withExpressionTSSMethyCpgOverlapPromoter",
                                         filename)
            tcga_data = pd.read_csv(tcga_filepath, index_col=0)
            tcga_data_labels = [tell_direction(item) for item in tcga_data["flippedZscore"]]
            tcga_score = measure_auc(model, tcga_data[features], tcga_data_labels)
            print(filename, " Test score: ", tcga_score)
    # TCGA - CHOL - mval_all_features_added.txt Test score: nan
    # TCGA - BLCA - mval_all_features_added.txt Test score: 0.496517056943
    # TCGA - KIRC - mval_all_features_added.txt Test score: nan
    # TCGA - HNSC - mval_all_features_added.txt Test score: 0.45
    # TCGA - THCA - mval_all_features_added.txt Test score: 0.503546910755
    # TCGA - PRAD - mval_all_features_added.txt Test score: 0.474358974359
    # TCGA - UCEC - mval_all_features_added.txt Test score: 0.514509242092
    # TCGA - BRCA - mval_all_features_added.txt Test score: 0.549476888306

    raise NotImplementedError
    # permutation tests
    train_data = load_data(train_path, exclude=exclude, keep=keep)
    # raise NotImplementedError
    features = [item for item in train_data.train.values.columns if item not in ["OverallZScore"]]
    model = RandomForestClassifier()
    cv = StratifiedKFold(2)
    model.fit(train_data.train.values[features], train_data.train.labels)
    print(model.score(train_data.test.values[features], train_data.test.labels))
    score, permutation_scores, pvalue = permutation_test_score(
        model, train_data.train.values[features], train_data.train.labels,
        scoring="roc_auc", cv=cv, n_permutations=100, n_jobs=4)

    print("Classification score %s (pvalue : %s)" % (score, pvalue))

    # #############################################################################
    # View histogram of permutation scores
    plt.hist(permutation_scores, 20, label='Permutation scores',
             edgecolor='black')
    ylim = plt.ylim()
    # BUG: vlines(..., linestyle='--') fails on older versions of matplotlib
    # plt.vlines(score, ylim[0], ylim[1], linestyle='--',
    #          color='g', linewidth=3, label='Classification Score'
    #          ' (pvalue %s)' % pvalue)
    # plt.vlines(1.0 / n_classes, ylim[0], ylim[1], linestyle='--',
    #          color='k', linewidth=3, label='Luck')
    plt.plot(2 * [score], ylim, '--g', linewidth=3,
             label='Classification Score'
                   ' (pvalue %s)' % pvalue)

    plt.ylim(ylim)
    plt.legend()
    plt.xlabel('Score')
    plt.show()
