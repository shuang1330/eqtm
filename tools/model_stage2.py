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
    DATA_DIR = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/withExpressionTSSMethyCpgOverlapPromoter"
    bios_sig_path = os.path.join(DATA_DIR, "eqtms_biosComplete_sig_all_features_added.txt")
    tcga_sig_path = os.path.join(DATA_DIR, "eqtms_metaTCGA_truePositives_notSeenInBIOSTrainingSet.txt")
    bios_sig = pd.read_csv(bios_sig_path)
    tcga_sig = pd.read_csv(tcga_sig_path)
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

    determineLabel = lambda x: 1 if x > 0 else 0
    bios_sig["label"] = [determineLabel(zscore) for zscore in bios_sig["OverallZScore"]]
    tcga_sig["label"] = [determineLabel(zscore) for zscore in tcga_sig["flippedZscore"]]

    train_x, test_x, train_y, test_y = train_test_split(bios_sig[features], bios_sig["label"],
                                                        test_size=0.2,
                                                        random_state=5)
    valid_x, valid_y = tcga_sig[features], tcga_sig["label"]

    model = RandomForestClassifier(n_estimators=100)
    model.fit(train_x, train_y)
    pickle.dump(model,
                open("/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/model/stage2.pkl", "wb"))
    print("Model saved in /groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/model/stage2.pkl")
    test_auc = measure_auc(model, test_x, test_y)
    valid_auc = measure_auc(model, valid_x, valid_y)

    print(test_auc, valid_auc)

    test_preds = model.predict(test_x[features])
    tn, fp, fn, tp = metrics.confusion_matrix(test_y, test_preds).ravel()
    print("Test dataset: \n")
    print("tn\tfp\tfn\ttp\n")
    print(tn, fp, fn, tp)

    valid_preds = model.predict(valid_x[features])
    tn, fp, fn, tp = metrics.confusion_matrix(valid_y, valid_preds).ravel()
    print("Validation dataset: \n")
    print("tn\tfp\tfn\ttp\n")
    print(tn, fp, fn, tp)
