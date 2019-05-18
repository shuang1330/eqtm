import pandas as pd
import os


def add_cpg_overlapRatio(eqtm, overlap_dict):
    a_key = list(overlap_dict.keys())[0]
    features = overlap_dict[a_key].keys()
    for col in features:
        def find_overlap_vector(item):
            if item in overlap_dict:
                return overlap_dict[item][col]
            else:
                return 0
        col_name = col
        eqtm[col_name] = [find_overlap_vector(item) for item in eqtm['SNPName'].values]
    print(eqtm.shape, eqtm.head())
    return eqtm


if __name__ == "__main__":
    project_rootdir = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
    cpg_overlapRatio_filepath = os.path.join(project_rootdir, "data", "eqtmZscores",
                                             "allCpgs", "cpg_methy_anno_overlapRatio.txt")
    splits_dirpath = os.path.join(project_rootdir, "data", "splits")

    cpg_overlapRatio_dict = pd.read_csv(cpg_overlapRatio_filepath, sep="\t", index_col=0).T.to_dict()

    for split_name in os.listdir(splits_dirpath):
        print("Processing ", split_name)
        split_filepath = os.path.join(splits_dirpath, split_name)
        split_data = pd.read_csv(split_filepath, sep="\t")
        split_data_with_ratio = add_cpg_overlapRatio(split_data, cpg_overlapRatio_dict)
        split_data_with_ratio.to_csv(os.path.join(project_rootdir, "data", "splits_with_features",
                                                  split_name))