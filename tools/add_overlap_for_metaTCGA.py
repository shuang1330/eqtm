import os
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics


def measure_auc(classifier, test_data_values, test_data_labels):
    prods = classifier.predict_proba(test_data_values)[:, 1]
    fpr, tpr, _ = metrics.roc_curve(test_data_labels, prods)
    score = metrics.auc(fpr, tpr)
    return score

if __name__ == "__main__":
    tcga_eqtm_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/tcga/meta_analysis/output/eQTMs_meta.txt"
    tcga_overlap_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/allCpgs/cpg_metaTCGA_overlapRatio.txt"

    tcga_eqtm = pd.read_csv(tcga_eqtm_filepath, sep="\t")
    tcga_overlapRatio = pd.read_csv(tcga_overlap_filepath, sep="\t", index_col=0)
    tcga_overlapRatio_dict = tcga_overlapRatio.T.to_dict()

    # print(tcga_overlapRatio_dict)
    print(tcga_eqtm.head())

    def findOverlapRatio(x, col):
        if x in tcga_overlapRatio_dict:
            return tcga_overlapRatio_dict[x][col]
        else:
            print(x)
            return 0
    for col in tcga_overlapRatio.columns:
        tcga_eqtm[col] = [findOverlapRatio(item, col) for item in tcga_eqtm["SNPName"]]

    tss = pd.read_csv("/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/features/TSSDistance/gene_startEndSite.txt")
    tss_dict = tss[["geneName", "TssSite"]].set_index("geneName").T.to_dict()
    print(tss_dict["ENSG00000258484"])
    findDistance = lambda name, pos: abs(tss_dict[name]["TssSite"]-pos)
    tcga_eqtm["TssDistance"] = [findDistance(item[0], item[1]) for item in zip(tcga_eqtm["ProbeName"], tcga_eqtm["ProbeCenterChrPos"])]
    flipDirection = lambda zscore, allele:zscore if allele == "C" else -zscore
    tcga_eqtm["flippedZscore"] = [flipDirection(item[0], item[1]) for item in zip(tcga_eqtm["OverallZScore"], tcga_eqtm["AlleleAssessed"])]
    tellDirection = lambda x:0 if x <0 else 1
    tcga_eqtm["label"] = [tellDirection(item) for item in tcga_eqtm["flippedZscore"]]
    print(tcga_eqtm.head())
    features = ['TssDistance', 'H4K5ac', 'H2A.Z', 'H2BK120ac', 'H3K79me2', 'H3K27ac',
                'H2BK20ac', 'H3K14ac', 'H3K9ac', 'H3K4ac', 'H2AK5ac',
                'H3K4me1', 'H3K18ac', 'H3K23ac', 'H2BK5ac', 'H3K4me3',
                'H2BK12ac', 'H3K23me2', 'H4K12ac', 'DNase', 'H2BK15ac',
                'H3K4me2', 'H3K79me1', 'H2AK9ac', 'H3T11ph', 'H4K8ac',
                'H4K91ac', 'H3K56ac']


    tcga_eqtm.to_csv("/groups/umcg-gcc/tmp03/umcg-sli/temp_data/added_overlapRatio_metaTCGA.txt")

    train_data = pd.read_csv("/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/withExpressionTSSMethyCpgOverlapPromoter/2017-12-09-eQTLsFDR-et0_withExpressionTssMethyOverlap_PromoterOverlap.txt")
    train_data["label"] = [tellDirection(item) for item in train_data["OverallZScore"]]
    model = RandomForestClassifier(n_estimators=100)
    model.fit(train_data[features], train_data["label"])
    preds = model.predict(tcga_eqtm[features])
    print(preds)
    print(tcga_eqtm[tcga_eqtm["label"] == 0].shape)
    print(tcga_eqtm[tcga_eqtm["label"] == 1].shape)
    auc = measure_auc(model, tcga_eqtm[features], tcga_eqtm["label"])
    print(auc)