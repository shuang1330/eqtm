import os
import gzip
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics


def read_eqtm_as_dict(eqtm_filepath, sep="\t",zscore_col_name="OverallZScore"):
    eqtm = pd.read_csv(eqtm_filepath, sep=sep)
    eqtm["eqtm_combi"] = ["{}_{}".format(item[0], item[1]) for
                          item in zip(eqtm["SNPName"], eqtm["ProbeName"])]
    eqtm_dict = eqtm[["eqtm_combi", zscore_col_name]].set_index("eqtm_combi").T.to_dict()
    return eqtm_dict


def read_and_check(eqtm_dict, tcga_filepath, eqtm_dict_name="OverallZScore"):
    eqtms = []
    zscores = [[],[]]
    with gzip.open(tcga_filepath, "rb") as f:
        info = f.readlines()
        # col_names = info[0].decode("utf-8").split("\t")
        for ind in tqdm(range(1, len(info))):
            line = [item.strip() for item in info[ind].decode("utf-8").split("\t")]
            eqtm_combi = "{}_{}".format(line[1], line[4])
            zscore = float(line[10])
            allele = line[9]
            if eqtm_combi in eqtm_dict:
                eqtms.append(eqtm_combi)
                if allele == "T":
                    zscores[1].append(-zscore)
                else:
                    zscores[1].append(zscore)
                zscores[0].append(eqtm_dict[eqtm_combi][eqtm_dict_name])

    return eqtms, zscores


def check_and_return_common_and_unique_eqtms(eqtm_dict, tcga_filepath):
    common_eqtms = []
    unique_to_tcga_eqtms = []
    with gzip.open(tcga_filepath, "rb") as f:
        info = f.readlines()
        for ind in tqdm(range(1, len(info))):
            line = [item.strip() for item in info[ind].decode("utf-8").split("\t")]
            eqtm_combi = "{}_{}".format(line[1], line[4])
            if eqtm_combi in eqtm_dict:
                common_eqtms.append(eqtm_combi)
            else:
                unique_to_tcga_eqtms.append(eqtm_combi)
    unique_to_ori = list(set(eqtm_dict.keys())-set(common_eqtms))

    return common_eqtms, unique_to_tcga_eqtms, unique_to_ori


def compare_two_datasets(et_gt_dict, new_gz_dataset_filepath,
                         new_dataset_name, plot_savepath):
    eqtms, zscores = read_and_check(et_gt_dict, new_gz_dataset_filepath)
    print("Found in common: ", len(eqtms), " For example: ",
          eqtms[0], zscores[0][0], zscores[1][0])

    f = plt.figure()
    ax = f.add_subplot(111)
    ax.scatter(zscores[0], zscores[1])
    ax.set_xlabel("zscores in BIOS")
    ax.set_ylabel("zscores in %s"%new_dataset_name)
    ax.text(5, 4, "%d eqtms"%len(zscores[0]))
    plt.savefig(plot_savepath)
    print("saved figure to %s"%plot_savepath)
    return None


def drawing(eqtms, zscores, dataset1_name,
            dataset2_name, plot_savepath,
            ratio=0.0):
    print("Found in common: ", len(eqtms), " For example: ",
          eqtms[0], zscores[0][0], zscores[0][0])
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.scatter(zscores[0], zscores[1], alpha=0.2)
    ax.set_xlabel("zscores in %s" % dataset1_name)
    ax.set_ylabel("zscores in %s" % dataset2_name)
    ax.text(25, 4, "%d eqtms" % len(zscores[0]))
    if ratio > 0:
        ax.text(25, 3, "{} in concordance".format(ratio))
    ax.axvline(x=0)
    ax.axhline(y=0)
    plt.savefig(plot_savepath)
    print("saved figure to %s" % plot_savepath)
    return None


def drawing_for_cord_blood_and_bios(eqtms, zscores, dataset1_name,
                                    dataset2_name, plot_savepath,
                                    num_diff=0):
    print("Found in common: ", len(eqtms), " For example: ",
          eqtms[0], zscores[0][0], zscores[0][0])
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.scatter(zscores[0], zscores[1], alpha=0.2)
    ax.set_xlabel("rvalue in %s" % dataset1_name)
    ax.set_ylabel("zscores in %s" % dataset2_name)
    ax.text(0.4, 30, "%d eqtms" % len(zscores[0]))
    if num_diff > 0:
        ax.text(0.4, 26, "%d in different direction" % num_diff)
    ax.axvline(x=0)
    ax.axhline(y=0)
    plt.savefig(plot_savepath)
    print("saved figure to %s" % plot_savepath)
    return None


def drawing_for_liver_and_bios(eqtms, zscores, dataset1_name,
                               dataset2_name, plot_savepath,
                               num_diff=0):
    print("Found in common: ", len(eqtms), " For example: ",
          eqtms[0], zscores[0][0], zscores[0][0])
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.scatter(zscores[0], zscores[1], alpha=0.2)
    ax.set_xlabel("Spearmanr in %s" % dataset1_name)
    ax.set_ylabel("zscores in %s" % dataset2_name)
    ax.text(0.4, 30, "%d eqtms" % len(zscores[0]))
    if num_diff > 0:
        ax.text(0.4, 26, "%d in different direction" % num_diff)
    ax.axvline(x=0)
    ax.axhline(y=0)
    plt.savefig(plot_savepath)
    print("saved figure to %s" % plot_savepath)
    return None


def compare_two_dicts_return_common_eqtms(dic1, dic2,
                                          zscore1_col_name="OverallZScore",
                                          zscore2_col_name="OverallZScore"):
    eqtms = list(set(dic1.keys()) & set(dic2.keys()))
    print(dic1[eqtms[0]])
    zscores = [[dic1[item][zscore1_col_name] for item in eqtms],
               [dic2[item][zscore2_col_name] for item in eqtms]]
    return eqtms, zscores


def find_diff_direction_bios_split(eqtms, zscores):
    diff_eqtms = []
    for ind in range(len(eqtms)):
        zscore1 = zscores[0][ind]
        zscore2 = zscores[1][ind]
        if zscore1*zscore2 < 0:
            diff_eqtms.append(eqtms[ind])
    print("%d number of eQTMs with opposite directions."%len(diff_eqtms))
    return diff_eqtms


def examine_common_gene_involved(dataset1, dataset2, genename_col_1, genename_col_2):
    common_genes = set(dataset1[genename_col_1]) & set(dataset2[genename_col_2])
    unique_to_dataset1 = set(dataset1[genename_col_1]) - common_genes
    unique_to_dataset2 = set(dataset2[genename_col_2]) - common_genes
    print("Dataset1 %d, Dataset2 %d"%(len(dataset1[genename_col_1].unique()),
                                      len(dataset2[genename_col_2].unique())))
    print("Common genes %d"%len(common_genes))
    print("Unique to dataset1: %d"%len(unique_to_dataset1))
    print("Unique to dataset2: %d"%len(unique_to_dataset2))


def examine_all_three_datasets_common_genes(dataset1, dataset2, dataset3,
                                            genename_col1, genename_col2, genename_col3):
    common_genes = set(dataset1[genename_col1]) & set(dataset2[genename_col2]) & set(dataset3[genename_col3])
    print("Genes common to all three datasets: %d"%len(common_genes))


def compare_three_datasets_return_4_subsets(dataset1_dic, dataset2_dic, dataset3_dic,
                                            zname1):
    dat1 = dataset1_dic.keys()
    dat2 = dataset2_dic.keys()
    dat3 = dataset3_dic.keys()

    subset1 = [dataset1_dic[name][zname1]
               for name in dat1 if name not in dat2 and name not in dat3]
    subset2 = [dataset1_dic[name][zname1]
               for name in dat1 if name in dat2 and name not in dat3]
    subset3 = [dataset1_dic[name][zname1]
               for name in dat1 if name not in dat2 and name in dat3]
    subset4 = [dataset1_dic[name][zname1]
               for name in dat1 if name in dat2 and name in dat3]

    print("For example")
    print("Subset1", len(subset1), subset1[0])
    print("Subset2", len(subset2), subset2[0])
    print("Subset3", len(subset3), subset3[0])
    print("Subset4", len(subset4), subset4[0])
    return subset1, subset2, subset3, subset4


def draw_histogram_for_4_different_subsets(savepath, subset1, subset2, subset3, subset4,
                                          name1="unique", name2="commonTdata1",
                                          name3="commonTdata2", name4="allCommon"):
    f = plt.figure()
    ax = f.add_subplot(111)
    min_value = min(min(subset1), min(subset2),
                    min(subset3), min(subset4))
    max_value = max(max(subset1), max(subset2),
                    max(subset3), max(subset4))
    # min_value = min(min(subset3), min(subset4))
    # max_value = max(max(subset3), max(subset4))
    ax.set_xlim(min_value, max_value)
    ax.hist(subset1, bins=50, label=name1, alpha=0.2)
    ax.hist(subset2, bins=50, label=name2, alpha=0.2)
    ax.hist(subset3, bins=50, label=name3, alpha=0.5)
    ax.hist(subset4, bins=50, label=name4, alpha=0.5)
    plt.legend()

    plt.savefig(savepath)
    print("Plot saved in %s"%savepath)


def measure_auc(classifier, test_data_values, test_data_labels):
    prods = classifier.predict_proba(test_data_values)[:, 1]
    fpr, tpr, _ = metrics.roc_curve(test_data_labels, prods)
    score = metrics.auc(fpr, tpr)
    return score


if __name__ == "__main__":
    meta_tcga_path = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/withExpressionTSSMethyCpgOverlapPromoter/meta_analysis_tcga_flipped_all_features_added.txt"
    et_gt_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/withExpressionTSSMethyCpgOverlapPromoter/et_gt_all_celltypes.txt"

    meta_tcga = pd.read_csv(meta_tcga_path, sep=",")
    et_gt = pd.read_csv(et_gt_filepath, sep=",")
    meta_tcga["eqtm_combi"] = ["{}_{}".format(item[0], item[1]) for
                          item in zip(meta_tcga["SNPName"], meta_tcga["ProbeName"])]
    et_gt["eqtm_combi"] = ["{}_{}".format(item[0], item[1]) for
                               item in zip(et_gt["SNPName"], et_gt["ProbeName"])]
    getDirection = lambda x: 0 if x < 0 else 1
    et_gt["label"] = [getDirection(item) for item in et_gt["OverallZScore"]]
    common_eqtms = set(meta_tcga["eqtm_combi"].values) & set(et_gt["eqtm_combi"])


    isCommon = lambda x:True if x in common_eqtms else False
    et_gt["isCommon"] = [isCommon(item) for item in et_gt["eqtm_combi"]]
    meta_tcga["isCommon"] = [isCommon(item) for item in meta_tcga["eqtm_combi"]]
    meta_tcga["label"] = [getDirection(item) for item in meta_tcga["flippedZscore"]]

    print(len(common_eqtms))
    for eqtm_combi in common_eqtms:
        if et_gt["label"][et_gt["eqtm_combi"] == eqtm_combi].values!=meta_tcga["label"][meta_tcga["eqtm_combi"] == eqtm_combi].values:
            print(et_gt[et_gt["eqtm_combi"] == eqtm_combi])
            print(meta_tcga[meta_tcga["eqtm_combi"] == eqtm_combi])
    # cg05388281 ENSG00000250312
    # cg13084669 ENSG00000244165

    features = ['TssDistance', 'H4K5ac', 'H2A.Z', 'H2BK120ac', 'H3K79me2', 'H3K27ac',
                'H2BK20ac', 'H3K14ac', 'H3K9ac', 'H3K4ac', 'H2AK5ac',
                'H3K4me1', 'H3K18ac', 'H3K23ac', 'H2BK5ac', 'H3K4me3',
                'H2BK12ac', 'H3K23me2', 'H4K12ac', 'DNase', 'H2BK15ac',
                'H3K4me2', 'H3K79me1', 'H2AK9ac', 'H3T11ph', 'H4K8ac',
                'H4K91ac', 'H3K56ac']

    model = RandomForestClassifier(n_estimators=100)
    model.fit(meta_tcga[features][meta_tcga["isCommon"] == False], meta_tcga["label"][meta_tcga["isCommon"] == False])
    auc1 = measure_auc(model, meta_tcga[features][meta_tcga["isCommon"] == False], meta_tcga["label"][meta_tcga["isCommon"] == False])
    auc2 = measure_auc(model, et_gt[features][et_gt["isCommon"] == True], et_gt["label"][et_gt["isCommon"] == True])
    auc3 = measure_auc(model, et_gt[features][et_gt["isCommon"] == False], et_gt["label"][et_gt["isCommon"] == False])

    print(auc1)
    print(auc2)
    print(auc3)