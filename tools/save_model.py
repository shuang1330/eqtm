import __init__path
import pickle
import os
from sklearn.ensemble import RandomForestClassifier
from lib.read.read_data import load_data
from sklearn.metrics import auc, roc_curve
import matplotlib.pyplot as plt
import pandas as pd

if __name__ == "__main__":
    project_rootdir = '/home/shuang/projects/development_eqtm'
    eqtm_dirpath = os.path.join(project_rootdir,
                                "data",
                                "eqtmZscores",
                                "withExpressionTSSMethyCpgOverlapGene")
    eqtm_names = ["2017-12-09-eQTLsFDR-et0.0-flipped_withExpressionTssMethy_blood_withGeneOverlap.txt",
                  "2017-12-09-eQTLsFDR-et0.0-flipped_withExpressionTssMethy_nonBlood_withGeneOverlap.txt",
                  "2017-12-09-eQTLsFDR-et0_withExpressionTssMethyOverlap_withGeneOverlap.txt",
                  "2017-12-09-eQTLsFDR-gt0.0-flipped_withExpressionTssMethy_blood_withGeneOverlap.txt",
                  "2017-12-09-eQTLsFDR-gt0.0-flipped_withExpressionTssMethy_nonBlood_withGeneOverlap.txt",
                  "2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_withGeneOverlap.txt",
                  "et_gt_allcell_withExpressionTssMethyOverlap_withGeneOverlap.txt",
                  "et_gt_blood_withExpressionTssMethy_blood_withGeneOverlap.txt",
                  "et_gt_nonBlood_withExpressionTssMethy_withGeneOverlap.txt"]
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

    train_path = os.path.join(project_rootdir, "data",
                              "eqtmZscores",
                              "withExpressionTSSMethyCpgOverlapGene",
                              "2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_withGeneOverlap.txt")
    test_path = os.path.join(project_rootdir, "data",
                              "eqtmZscores",
                             "withExpressionTSSMethyCpgOverlapGene",
                              "2017-12-09-eQTLsFDR-et0_withExpressionTssMethyOverlap_withGeneOverlap.txt")
    train_data = load_data(train_path, exclude=exclude, keep=keep, test_size=0)
    test_data = load_data(test_path, exclude=exclude, keep=keep, test_size=0)

    features = [cor for cor in train_data.train.values.columns if cor not in ["OverallZScore"]]
    model = pickle.load(open(os.path.join(project_rootdir,
                                          "model",
                                          "train_gt.pkl"), "rb"))
    print(model.feature_importances_)
    print(features)
    print(len(features), len(model.feature_importances_))
    assert len(features) == len(model.feature_importances_)

    model = RandomForestClassifier(n_jobs=8, n_estimators=300, class_weight="balanced")
    model.fit(train_data.train.values[features], train_data.train.labels)

    all_data_values = pd.concat([train_data.train.values, test_data.train.values], axis=0)
    all_data_labels = pd.concat([train_data.train.labels, test_data.train.labels], axis=0)
    print(all_data_values.shape)
    all_positive_data = all_data_values[all_data_labels > 0]
    probs = model.predict_proba(all_positive_data[features])[:, 1]
    all_negative_data = all_data_values[all_data_labels == 0]
    probs_neg = model.predict_proba(all_negative_data[features])[:, 0]
    zscores = all_positive_data["OverallZScore"]
    zscores_neg = all_negative_data["OverallZScore"]

    f = plt.figure()
    ax = f.add_subplot(111)
    ax.scatter(probs, zscores, alpha=0.2, color="b", label="Positive")
    ax.scatter(probs_neg, zscores_neg, alpha=0.2, color="r", label="Negative")
    ax.set_xlabel("Confidence Score")
    ax.set_ylabel("Zscores")
    ax.legend(loc=2)
    plt.savefig("./development_eqtm/output_fig/confidenceScore_Zscore.png")
    plt.show()

    # fpr, tpr, _ = roc_curve(test_data.train.labels, prods)
    # score = auc(fpr, tpr)
    # print(score)
    # pickle.dump(model, open(os.path.join(project_rootdir, "model", "train_gt.pkl"), "wb"))

    # FDR and confidence score
    exclude2 = ['SNPName', 'SNPChr', 'PValue', 'SNPChrPos', 'ProbeName', 'ProbeChr',
               'ProbeCenterChrPos', 'CisTrans', 'SNPType', 'AlleleAssessed',
               'DatasetsWhereSNPProbePairIsAvailableAndPassesQC',
               'DatasetsZScores', 'DatasetsNrSamples',
               'IncludedDatasetsMeanProbeExpression',
               'IncludedDatasetsProbeExpressionVariance', 'HGNCName',
               'IncludedDatasetsCorrelationCoefficient', 'Meta-Beta (SE)',
               'Beta (SE)', 'FoldChange', 'checkChr', 'SNPName_ProbeName',
               'TssSite', 'chr']
    keep2 = ['FDR', 'OverallZScore']
    train_data = load_data(train_path, exclude=exclude2, keep=keep2, test_size=0)
    test_data = load_data(test_path, exclude=exclude2, keep=keep2, test_size=0)

    features = [cor for cor in train_data.train.values.columns
                if cor not in ["OverallZScore", "FDR"]]

    all_data_values = pd.concat([train_data.train.values, test_data.train.values], axis=0)
    all_data_labels = pd.concat([train_data.train.labels, test_data.train.labels], axis=0)
    print(all_data_values.shape)
    all_positive_data = all_data_values[all_data_labels > 0]
    probs = model.predict_proba(all_positive_data[features])[:, 1]
    zscores = all_positive_data["FDR"]

    all_negative_data = all_data_values[all_data_labels == 0]
    probs_neg = model.predict_proba(all_negative_data[features])[:, 0]
    zscores_neg = all_negative_data["FDR"]

    f = plt.figure()
    ax = f.add_subplot(111)
    ax.scatter(probs, zscores, alpha=0.2, color="b", label="Positive")
    ax.scatter(probs_neg, zscores_neg, alpha=0.2, color="r", label="Negative")
    ax.set_xlabel("ConfidenceScore")
    ax.set_ylabel("FDR")
    ax.legend(loc=2)
    plt.savefig("./development_eqtm/output_fig/confidenceScore_FDR.png")
    plt.show()

    # find all tp and tn samples
    exclude3 = ['SNPChr', 'PValue', 'SNPChrPos', 'ProbeName', 'ProbeChr',
                'ProbeCenterChrPos', 'CisTrans', 'SNPType', 'AlleleAssessed',
                'DatasetsWhereSNPProbePairIsAvailableAndPassesQC',
                'DatasetsZScores', 'DatasetsNrSamples',
                'IncludedDatasetsMeanProbeExpression',
                'IncludedDatasetsProbeExpressionVariance', 'HGNCName',
                'IncludedDatasetsCorrelationCoefficient', 'Meta-Beta (SE)',
                'Beta (SE)', 'FoldChange', 'checkChr', 'SNPName_ProbeName',
                'TssSite', 'chr']
    keep3 = ['SNPName', 'FDR', 'OverallZScore']
    train_data = load_data(train_path, exclude=exclude3, keep=keep3, test_size=0)
    test_data = load_data(test_path, exclude=exclude3, keep=keep3, test_size=0)
    features = [item for item in train_data.train.values.columns
                if item not in ['SNPName', 'FDR', 'OverallZScore']]
    model = pickle.load(open(os.path.join(project_rootdir,
                                          "model",
                                          "train_gt.pkl"), "rb"))
    all_samples = pd.concat([train_data.train.values, test_data.train.values], axis=0)
    all_labels = pd.concat([train_data.train.labels, test_data.train.labels], axis=0)
    all_probs = model.predict_proba(all_samples)

