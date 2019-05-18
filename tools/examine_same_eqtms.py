import os
import gzip
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt


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


if __name__ == "__main__":
    et_gt_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/ORIGIN/et_gt.txt"  # 36k eqtms
    et_gt_dict = read_eqtm_as_dict(et_gt_filepath, sep="\t")
    # et_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/withExpressionTSSMethyCpgOverlapPromoter/2017-12-09-eQTLsFDR-et0_withExpressionTssMethyOverlap_PromoterOverlap.txt"  # 36k eqtms
    # et_dict = read_eqtm_as_dict(et_filepath, sep=",")
    bios_replicates_in_tcga_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/replication_output/test/eQTLs.txt.gz"
    eqtms_bios_replicates_in_tcga, zscores_bios_replicates_in_tcga = read_and_check(et_gt_dict,
                                                                                    bios_replicates_in_tcga_filepath,
                                                                                    eqtm_dict_name="OverallZScore")
    diff_eqtms = find_diff_direction_bios_split(eqtms_bios_replicates_in_tcga, zscores_bios_replicates_in_tcga)
    bios_replicates_in_tcga_save_figpath = "/groups/umcg-gcc/tmp03/umcg-sli/replication_output/figs/bios_replicates_in_tcga.png"
    # drawing(eqtms_bios_replicates_in_tcga, zscores_bios_replicates_in_tcga, "BIOS", "TCGA",
    #         bios_replicates_in_tcga_save_figpath, ratio=(1 - round(len(diff_eqtms)/len(set(et_dict.keys())), 2)))
    meta_tcga_path = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/withExpressionTSSMethyCpgOverlapPromoter/meta_analysis_tcga_flipped_all_features_added.txt"
    meta_tcga_dic = read_eqtm_as_dict(meta_tcga_path, sep=",", zscore_col_name="flippedZscore")
    common_eqtms = set(eqtms_bios_replicates_in_tcga) & set(meta_tcga_dic.keys())
    common_zscores = [[],[]]
    for eqtm_name in eqtms_bios_replicates_in_tcga:
        if eqtm_name in common_eqtms:
            common_zscores[0].append(et_gt_dict[eqtm_name]["OverallZScore"])
            common_zscores[1].append(meta_tcga_dic[eqtm_name]["flippedZscore"])
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.scatter(common_zscores[0], common_zscores[1], alpha=0.2, color="r")
    ax.scatter(zscores_bios_replicates_in_tcga[0], zscores_bios_replicates_in_tcga[1], alpha=0.2, color="b")
    ax.set_xlabel("zscores in BIOS")
    ax.set_ylabel("zscores in TCGA")
    ax.text(20, 7.5, "%d eqtms" % len(zscores_bios_replicates_in_tcga[0]))
    ax.text(20, 5.5, "{} in concordance".format(1-round(len(diff_eqtms)/len(eqtms_bios_replicates_in_tcga), 2)))
    ax.axvline(x=0)
    ax.axhline(y=0)
    plt.savefig("/groups/umcg-gcc/tmp03/umcg-sli/replication_output/figs/bios_replicates_in_tcga_significant_in_red.png")
    raise NotImplementedError
    # check for meta TCGA
    meta_tcga_path = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/withExpressionTSSMethyCpgOverlapPromoter/meta_analysis_tcga_flipped_all_features_added.txt"
    meta_tcga_dic = read_eqtm_as_dict(meta_tcga_path, sep=",", zscore_col_name="flippedZscore")
    meta_replicates_in_bios_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/replication_output/concat_replicate_results/liver_metaTCGA_replicates_in_bios_flipped.txt"
    meta_replicates_in_bios_dic = read_eqtm_as_dict(meta_replicates_in_bios_filepath,
                                                     sep="\t",
                                                     zscore_col_name="flippedZscore")
    eqtms_meta_replicates_in_bios, zscores_meta_replicates_In_bios = \
        compare_two_dicts_return_common_eqtms(meta_tcga_dic, meta_replicates_in_bios_dic,
                                              zscore1_col_name="flippedZscore",
                                              zscore2_col_name="flippedZscore")

    # check for liver
    liver_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/withExpressionTSSMethyCpgOverlapPromoter/lung_eqtms_with_all_features.txt"
    liver = pd.read_csv(liver_filepath, sep=",")
    liver["eqtm_combi"] = ["{}_{}".format(item[0], item[1].split(".")[0]) for
                           item in zip(liver["CpgSite"], liver["GeneName"])]
    liver_dic = liver[["eqtm_combi", "Spearmanr"]].set_index("eqtm_combi").T.to_dict()
    liver_replicates_in_bios_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/replication_output/concat_replicate_results/liver_metaTCGA_replicates_in_bios_flipped.txt"
    liver_replicates_in_bios_dic = read_eqtm_as_dict(liver_replicates_in_bios_filepath,
                                                          sep="\t",
                                                          zscore_col_name="flippedZscore")
    eqtms_liver_replicates_in_bios, zscores_liver_replicates_In_bios = \
        compare_two_dicts_return_common_eqtms(liver_dic, liver_replicates_in_bios_dic,
                                              zscore1_col_name="Spearmanr",
                                              zscore2_col_name="flippedZscore")

    # check for cord_blood
    cord_blood_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/withExpressionTSSMethyCpgOverlapPromoter/eqtm_fibroblast_allFeaturesAdded.txt"

    cord_blood_replicates_in_bios_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/replication_output/concat_replicate_results/cordBlood_replicates_in_bios_flipped.txt"

    cord_blood = pd.read_csv(cord_blood_filepath, sep=",")
    cord_blood["eqtm_combi"] = ["{}_{}".format(item[0], item[1].split(".")[0]) for
                          item in zip(cord_blood["METHYL_label"], cord_blood["GENE_ID"])]
    cord_blood_dic = cord_blood[["eqtm_combi", "rvalue"]].set_index("eqtm_combi").T.to_dict()
    cord_blood_replicates_in_bios_dic = read_eqtm_as_dict(cord_blood_replicates_in_bios_filepath,
                                                          sep="\t",
                                                          zscore_col_name="flippedZscore")
    print(set(cord_blood_dic.keys())- set(cord_blood_replicates_in_bios_dic.keys()))
    eqtms_cord_blood_replicates_in_bios, zscores_cord_blood_replicates_In_bios = \
        compare_two_dicts_return_common_eqtms(cord_blood_dic, cord_blood_replicates_in_bios_dic,
                                              zscore1_col_name="rvalue",
                                              zscore2_col_name="flippedZscore")

    list_of_eqtms_zscores = [[eqtms_meta_replicates_in_bios, zscores_meta_replicates_In_bios, "metaTCGA"],
                             [eqtms_liver_replicates_in_bios, zscores_liver_replicates_In_bios, "liver"],
                             [eqtms_cord_blood_replicates_in_bios, zscores_cord_blood_replicates_In_bios, "cord blood"]]
    diff_eqtms_dic = {}
    for eqtms_zscores in list_of_eqtms_zscores:
        print(eqtms_zscores[2])
        diff_eqtms = find_diff_direction_bios_split(eqtms_zscores[0], eqtms_zscores[1])
        diff_eqtms_dic[eqtms_zscores[2]]=len(diff_eqtms)

    # metaTCGA
    print(eqtms_meta_replicates_in_bios[0], zscores_meta_replicates_In_bios[0][0])
    plot_save_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/replication_output/figs/metaTCGA_replicates_in_bios.png"
    drawing(eqtms_meta_replicates_in_bios, zscores_meta_replicates_In_bios,
                               "metaTCGA", "replicates_in_BIOS", plot_save_filepath, diff_eqtms_dic["metaTCGA"])
    # liver
    print(eqtms_liver_replicates_in_bios[0], zscores_liver_replicates_In_bios[0][0])
    plot_save_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/replication_output/figs/liver_replicates_in_bios.png"
    drawing_for_liver_and_bios(eqtms_liver_replicates_in_bios, zscores_liver_replicates_In_bios,
                                    "Liver", "replicates_in_BIOS", plot_save_filepath, diff_eqtms_dic["liver"])

    # cord blood
    print(eqtms_cord_blood_replicates_in_bios[0], zscores_cord_blood_replicates_In_bios[0][0])
    plot_save_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/replication_output/figs/cord_blood_replicates_in_bios.png"
    drawing_for_cord_blood_and_bios(eqtms_cord_blood_replicates_in_bios, zscores_cord_blood_replicates_In_bios,
            "CordBlood", "replicates_in_BIOS", plot_save_filepath, diff_eqtms_dic["cord blood"])

    raise NotImplementedError
    replicates_in_meta_tcga = "/groups/umcg-gcc/tmp03/umcg-sli/replication_output/test/eQTLs.txt.gz"

    # read et_gt as dictionary
    et_gt_dict = read_eqtm_as_dict(et_gt_filepath)
    et_replicates_in_tcga_eqtms, et_replicates_in_tcga_zscores = read_and_check(et_gt_dict, replicates_in_meta_tcga)
    save_figpath = "/groups/umcg-gcc/tmp03/umcg-sli/replication_output/figs/bios_replicate_in_tcga.png"
    drawing(et_replicates_in_tcga_eqtms, et_replicates_in_tcga_zscores, "BIOS", "meta TCGA", save_figpath)
    print("Figure saved in %s"%save_figpath)

    raise NotImplementedError
    # check significant eqtms from bios and tcga
    bios_part1_path = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/withExpressionTSSMethyCpgOverlapPromoter/BIOS_resplit_Part1_flipped_all_features_added.txt"
    bios_part2_path = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/withExpressionTSSMethyCpgOverlapPromoter/BIOS_resplit_Part2_flipped_all_features_added.txt"
    meta_tcga_path = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/allCpgs/meta_analysis_tcga_flipped.txt"


    # load as dataframe
    et_gt = pd.read_csv(et_gt_filepath, sep="\t")
    bios_part1 = pd.read_csv(bios_part1_path, sep=",", index_col=0)
    bios_part2 = pd.read_csv(bios_part2_path, sep=",", index_col=0)
    et_gt["eqtm_combi"] = ["{}_{}".format(item[0], item[1]) for
                          item in zip(et_gt["SNPName"], et_gt["ProbeName"])]
    bios_part1["eqtm_combi"] = ["{}_{}".format(item[0], item[1]) for
                          item in zip(bios_part1["SNPName"], bios_part1["ProbeName"])]
    bios_part2["eqtm_combi"] = ["{}_{}".format(item[0], item[1]) for
                          item in zip(bios_part2["SNPName"], bios_part2["ProbeName"])]

    # print(et_gt.head())
    # print(bios_part1.head())
    # print(bios_part2.head())

    # et_gt_dict = read_eqtm_as_dict(et_gt_filepath)

    # # examine meta analysis
    # meta_dict = read_eqtm_as_dict(meta_tcga_path, sep=",", zscore_col_name="flippedZscore")
    # eqtms_et_meta, zscores_et_meta = compare_two_dicts_return_common_eqtms(et_gt_dict, meta_dict,
    #                                                                        zscore1_col_name="OverallZScore",
    #                                                                        zscore2_col_name="flippedZscore")
    # save_et_bios2_figpath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/output_fig/BIOSComplete_metaTCGA.png"
    # drawing(eqtms_et_meta, zscores_et_meta, "BIOScomplete", "meta_analysis", save_et_bios2_figpath)
    # raise NotImplementedError


    # # examine bios part1 and part2
    # bios_part1_dict = read_eqtm_as_dict(bios_part1_path, sep=",",
    #                                     zscore_col_name="flippedZscore")
    # bios_part2_dict = read_eqtm_as_dict(bios_part2_path, sep=",",
    #                                     zscore_col_name="flippedZscore")
    print("Dictionaries loaded.")

    # # get four subsets from et_gt dataset
    # subset1, subset2, subset3, subset4 = \
    #     compare_three_datasets_return_4_subsets(
    #         et_gt_dict, bios_part2_dict, bios_part1_dict,
    #     "OverallZScore"
    # )
    # histogram_plot_savepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/output_fig/et_biosResplit_4subsets.png"
    # draw_histogram_for_4_different_subsets(histogram_plot_savepath,
    #                                        subset1, subset2, subset3, subset4)
    # raise NotImplementedError

    # examine common and unique cpgs in each pair of two datasets
    pairs = [[et_gt, bios_part1, "eqtm_combi", "eqtm_combi", "complete", "part1"],
             [et_gt, bios_part2, "eqtm_combi", "eqtm_combi", "complete", "part2"],
             [bios_part1, bios_part2, "eqtm_combi", "eqtm_combi", "part1", "part2"]]

    for pair in pairs:
        print(pair[4], pair[5])
        examine_common_gene_involved(pair[0], pair[1], pair[2], pair[3])
    examine_all_three_datasets_common_genes(et_gt, bios_part1, bios_part2,
                                            "eqtm_combi", "eqtm_combi", "eqtm_combi")

    raise NotImplementedError

    # examine common and unique genes in each pair of two datasets
    pairs = [[et_gt, bios_part1, "ProbeName", "ProbeName", "complete", "part1"],
             [et_gt, bios_part2, "ProbeName", "ProbeName", "complete", "part2"],
             [bios_part1, bios_part2, "ProbeName", "ProbeName", "part1", "part2"]]

    for pair in pairs:
        print(pair[4], pair[5])
        examine_common_gene_involved(pair[0], pair[1], pair[2], pair[3])
    examine_all_three_datasets_common_genes(et_gt, bios_part1, bios_part2,
                                            "ProbeName", "ProbeName", "ProbeName")


    # load bios complete and bios split part1 and part2 as dictionary
    bios_part1_dict = read_eqtm_as_dict(bios_part1_path, sep=",",
                                        zscore_col_name="flippedZscore")
    bios_part2_dict = read_eqtm_as_dict(bios_part2_path, sep=",",
                                        zscore_col_name="flippedZscore")
    et_gt_dict = read_eqtm_as_dict(et_gt_filepath)

    # examine bios complete and bios split part2
    eqtms_et_bios1, zscores_et_bios1 = \
        compare_two_dicts_return_common_eqtms(et_gt_dict,
                                              bios_part1_dict,
                                              zscore1_col_name="OverallZScore",
                                              zscore2_col_name="flippedZscore")
    save_et_bios1_figpath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/output_fig/BIOScomplete_part1.png"
    drawing(eqtms_et_bios1, zscores_et_bios1, "BIOScomplete", "part1", save_et_bios1_figpath)
    # examine bios complete and bios split part 2
    eqtms_et_bios2, zscores_et_bios2 = \
        compare_two_dicts_return_common_eqtms(et_gt_dict,
                                              bios_part2_dict,
                                              zscore1_col_name="OverallZScore",
                                              zscore2_col_name="flippedZscore")
    save_et_bios2_figpath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/output_fig/BIOScomplete_part2.png"
    drawing(eqtms_et_bios2, zscores_et_bios2, "BIOScomplete", "part2", save_et_bios2_figpath)
    # examine bios split part 1 and bios split part 2
    eqtms, zscores = compare_two_dicts_return_common_eqtms(bios_part1_dict,
                                                           bios_part2_dict,
                                                           zscore1_col_name="flippedZscore",
                                                           zscore2_col_name="flippedZscore")
    i = 0
    for eqtm in eqtms:
        if eqtm in eqtms_et_bios1 and eqtm in eqtms_et_bios2:
            i += 1
    print("Common eqtm number:", i)
    print("Eqtms with different directions", find_diff_direction_bios_split(eqtms, zscores))
    save_figpath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/output_fig/BIOS_part1_part2.png"
    drawing(eqtms, zscores, "part1", "part2", save_figpath)

    raise NotImplementedError
    # examine bios complete and tcga datasets
    et_gt_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/ORIGIN/et_gt.txt"  # 36k eqtms
    tcga_brca_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/tcga/test_traitfile/output/eQTLsFDR.txt.gz" # 300k
    new_bios_filepath = "/groups/umcg-bios/tmp03/projects/2018-methylation/output/eqtm/2018-04-10-cis-1000k/eQTLsFDR0.05.txt.gz"

    et_gt_dict = read_eqtm_as_dict(et_gt_filepath) #{"cgs*_ENSG*":Zscore}
    # get unique eqtms
    common_eqtm_combis, unique_to_new, unique_to_ori = \
        check_and_return_common_and_unique_eqtms(et_gt_dict, new_bios_filepath)
    with open("/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/extras/common_eqtms_56k_36k.txt", "w") as f:
        f.write("\n".join(common_eqtm_combis))
    with open("/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/extras/unique_to_56k.txt", "w") as f:
        f.write("\n".join(unique_to_new))
    with open("/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/extras/unique_to_36k.txt", "w") as f:
        f.write("\n".join(unique_to_ori))
    plot_save_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/output_fig/etGt_56kBIOS.png"
    compare_two_datasets(et_gt_dict, new_bios_filepath, "20180410 56kBIOS", plot_save_filepath)

    # check bios complete and tcga breast cancer dataset
    eqtms, zscores = read_and_check(et_gt_dict, tcga_brca_filepath)
    print("Found in common: ", len(eqtms), " For example: ", eqtms[0], zscores["et_gt"][0], zscores["tcga"][0])

    f = plt.figure()
    ax = f.add_subplot(111)
    ax.scatter(zscores["et_gt"], zscores["tcga"])
    ax.set_xlabel("zscores in BIOS")
    ax.set_ylabel("zscores in TCGA-BRCA")
    ax.text(5, 4, "200 eqtms")
    plt.savefig("/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/output_fig/etGt_brca.png")
    print("saved figure to /groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/output_fig/etGt_brca.png")
