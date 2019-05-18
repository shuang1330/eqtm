import pandas as pd
import os
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.model_selection import train_test_split
import pickle


def measure_auc(classifier, test_data_values, test_data_labels):
    prods = classifier.predict_proba(test_data_values)[:,1]
    fpr, tpr, _ = metrics.roc_curve(test_data_labels, prods)
    score = metrics.auc(fpr, tpr)
    return score


if __name__ == "__main__":
    DATA_DIR = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/withExpressionTSSMethyCpgOverlapPromoter/"
    bios_sig_path = os.path.join(DATA_DIR, "eqtms_biosComplete_sig_all_features_added.txt")
    bios_insig_path = os.path.join(DATA_DIR, "eqtms_biosComplete_insig_all_features_added.txt")
    bios_between_path = os.path.join(DATA_DIR, "eqtms_biosComplete_between0.05_0.5_all_features_added.txt")

    tcga_sig_path = os.path.join(DATA_DIR, "eqtms_metaTCGA_sig_all_features_added.txt")
    tcga_insig_path = os.path.join(DATA_DIR, "eqtms_metaTCGA_insig_all_features_added.txt")
    tcga_between_path = os.path.join(DATA_DIR, "eqtms_metaTCGA_between0.05_0.5_all_features_added.txt")

    bios_sig = pd.read_csv(bios_sig_path)
    bios_insig = pd.read_csv(bios_insig_path)
    bios_between = pd.read_csv(bios_between_path)
    tcga_sig = pd.read_csv(tcga_sig_path)
    tcga_insig = pd.read_csv(tcga_insig_path).sample(923, random_state=6)
    tcga_between = pd.read_csv(tcga_between_path)
    print("Loaded all datasets.")
    print(bios_sig.head())
    print(tcga_sig.head())

    features = ['H4K5ac', 'H2A.Z', 'H2BK120ac', 'H3K79me2', 'H3K27ac',
                'H2BK20ac', 'H3K14ac', 'H3K9ac', 'H4K20me1', 'H3K4ac',
                'H2AK5ac', 'H3K9me3', 'H3K36me3', 'H3K4me1', 'H3K18ac',
                'H3K23ac', 'H2BK5ac', 'H3K4me3', 'H2BK12ac', 'H3K23me2',
                'H4K12ac', 'DNase', 'H2BK15ac', 'H3K9me1', 'H3K4me2',
                'H3K27me3', 'H3K79me1', 'H2AK9ac', 'H3T11ph', 'H4K8ac',
                'H4K91ac', 'H3K56ac', 'TssDistance','expressionMean',
                'expressionVar', 'methyMean', 'methyVar']
    bios_dataset = pd.concat([bios_sig, bios_insig])
    determineLabel = lambda x:1 if x < 0.05 else 0
    bios_dataset["label"] = [determineLabel(fdr) for fdr in bios_dataset["FDR"]]
    bios_dataset["eqtm_combi"] = ["{}_{}".format(item[0], item[1])
                                  for item in zip(bios_dataset["SNPName"], bios_dataset["ProbeName"])]
    bios_eqtm_combis = set(bios_dataset["eqtm_combi"])
    tcga_dataset = pd.concat([tcga_sig, tcga_insig])
    tcga_dataset["label"] = [determineLabel(fdr) for fdr in tcga_dataset["FDR"]]
    tcga_dataset["eqtm_combi"] = ["{}_{}".format(item[0], item[1]) for item in
                                  zip(tcga_dataset["SNPName"], tcga_dataset["ProbeName"])]
    isSeen = lambda x: True if x in bios_eqtm_combis else False
    tcga_dataset["isSeen"] = [isSeen(item) for item in tcga_dataset["eqtm_combi"]]
    print("How many eqtms are already seen in the bios dataset?\n", tcga_dataset["isSeen"].value_counts())

    train_x, test_x, train_y, test_y = train_test_split(bios_dataset[features],
                                                        bios_dataset["label"],
                                                        test_size=0.2,
                                                        random_state=4)
    valid_x = tcga_dataset[features][tcga_dataset["isSeen"] == False]
    valid_y = tcga_dataset["label"][tcga_dataset["isSeen"] == False]

    print("Start to train the model.")
    model = RandomForestClassifier(n_estimators=300)
    model.fit(train_x, train_y)
    pickle.dump(model,
                open("/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/model/stage1.pkl", "wb"))
    print("Model saved in /groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/model/stage1.pkl")
    test_auc = measure_auc(model, test_x, test_y)
    valid_auc = measure_auc(model, valid_x, valid_y)
    print(test_auc, valid_auc)

    test_preds = model.predict(test_x[features])
    tn, fp, fn, tp = metrics.confusion_matrix(test_y, test_preds).ravel()
    print("Test dataset: \n")
    print("tn\tfp\tfn\ttp\n")
    print(tn, fp, fn, tp)
    sig_bios_index = set(test_y[test_y == 1].index)
    right_prediction_index = set(test_y[test_y == test_preds].index)
    true_positive_index = list(sig_bios_index & right_prediction_index)
    truePositive_eqtms = bios_dataset.iloc[true_positive_index, :]
    print(truePositive_eqtms.shape)
    truePositive_eqtms.to_csv(os.path.join(DATA_DIR,
                                           "eqtms_biosComplete_truePositives.txt"),
                              index=False)

    valid_preds = model.predict(valid_x[features])
    tn, fp, fn, tp = metrics.confusion_matrix(valid_y, valid_preds).ravel()
    print("Validation dataset: \n")
    print("tn\tfp\tfn\ttp\n")
    print(tn, fp, fn, tp)

    sig_tcga_index = set(valid_y[valid_y == 1].index)
    right_prediction_index = set(valid_y[valid_y == valid_preds].index)
    true_positive_index = list(sig_tcga_index & right_prediction_index)
    truePositive_eqtms = tcga_dataset.iloc[true_positive_index, :]
    print(truePositive_eqtms.shape)
    truePositive_eqtms.to_csv(os.path.join(DATA_DIR,
                                           "eqtms_metaTCGA_truePositives_notSeenInBIOSTrainingSet.txt"),
                              index=False)