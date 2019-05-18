import pandas as pd
import os


if __name__ == '__main__':
    basepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
    original_eqtm_filepath = os.path.join(basepath, "data", "eqtmZscores", "allCpgs",
                                          "cpg_combined_significant_SNPs_location_2019-02-21.txt")
    overlapRatio_filepath = os.path.join(basepath, "data", "eqtmZscores", "allCpgs",
                                         "cpg_combined_significant_SNPs_location_2019-02-21_overlapRatio_blood.txt")
    overlapRatio = pd.read_csv(overlapRatio_filepath, sep="\t", index_col=0)
    features = overlapRatio.columns
    overlapratio_dict = overlapRatio.T.to_dict()

    original_eqtm = pd.read_csv(original_eqtm_filepath, sep='\t')
    original_eqtm['chromosome'] = ['chr{}'.format(chrom) for chrom in original_eqtm['chr']]
    original_eqtm["index"] = ["_".join([str(ele) for ele in item]) for item in
                              original_eqtm[["chromosome", "pos", 'rsID']].values]

    for feature in features:
        print("Adding feature:", feature)
        findFeaturevalue = lambda x:overlapratio_dict[x][feature] if x in overlapratio_dict else 0
        original_eqtm[feature] = [findFeaturevalue(variant) for variant in original_eqtm['index']]

    original_eqtm.to_csv(os.path.join(basepath, "data", "eqtmZscores", "allCpgs",
                                      "cpg_combined_significant_SNPs_location_2019-02-21_overlapRatio_blood_withOriginalInfo.txt"),
                         sep="\t")
    print("Saved in :", os.path.join(basepath, "data", "eqtmZscores", "allCpgs",
                                      "cpg_combined_significant_SNPs_location_2019-02-21_overlapRatio_blood_withOriginalInfo.txt"))


    overlapRatio_filepath = os.path.join(basepath, "data", "eqtmZscores", "allCpgs",
                                         "cpg_combined_significant_SNPs_location_2019-02-21_overlapRatio_nonBlood.txt")
    overlapRatio = pd.read_csv(overlapRatio_filepath, sep="\t", index_col=0)
    features = overlapRatio.columns
    overlapratio_dict = overlapRatio.T.to_dict()

    original_eqtm = pd.read_csv(original_eqtm_filepath, sep='\t')
    original_eqtm['chromosome'] = ['chr{}'.format(chrom) for chrom in original_eqtm['chr']]
    original_eqtm["index"] = ["_".join([str(ele) for ele in item]) for item in
                              original_eqtm[["chromosome", "pos", 'rsID']].values]

    for feature in features:
        print("Adding feature:", feature)
        findFeaturevalue = lambda x: overlapratio_dict[x][feature] if x in overlapratio_dict else 0
        original_eqtm[feature] = [findFeaturevalue(variant) for variant in original_eqtm['index']]

    original_eqtm.to_csv(os.path.join(basepath, "data", "eqtmZscores", "allCpgs",
                                      "cpg_combined_significant_SNPs_location_2019-02-21_overlapRatio_nonBlood_withOriginalInfo.txt"),
                         sep="\t")
    print("Saved in :", os.path.join(basepath, "data", "eqtmZscores", "allCpgs",
                                     "cpg_combined_significant_SNPs_location_2019-02-21_overlapRatio_nonBlood_withOriginalInfo.txt"))