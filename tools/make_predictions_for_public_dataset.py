import os
import pickle
import pandas as pd
import numpy as np


if __name__ == "__main__":
    project_rootdir = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
    model1_path = os.path.join(project_rootdir, "model", "stage1.pkl")
    model2_path = os.path.join(project_rootdir, "model", "stage2.pkl")
    model1 = pickle.load(open(model1_path, "rb"))
    model2 = pickle.load(open(model2_path, "rb"))

    splits_with_features_dirpath = os.path.join(project_rootdir, "data",
                                                "splits_with_features")
    with open(os.path.join(project_rootdir, "extras", "model1_2_features.txt"), "r") as f:
        features = [item.strip() for item in f.readlines()]

    for split_name in os.listdir(splits_with_features_dirpath):
        print("Processing ", split_name)
        split_filepath = os.path.join(splits_with_features_dirpath, split_name)
        split = pd.read_csv(split_filepath, sep=",", index_col=0)
        print(split[features].head())
        for col in features:
            null = split[col].isnull().sum()
            if null > 0:
                print(col, null)
        valid_sample_index = split[split[features].isna().any(axis=1)==False].index
        split["isComplete"] = ~split[features].isna().any(axis=1)
        conf_score1 = model1.predict_proba(split[features][split["isComplete"] == True])
        scores = np.zeros([split.shape[0], ])
        scores[valid_sample_index] = conf_score1[:, 1]
        split["model1_score"] = scores

        conf_score2 = model2.predict_proba(split[features][split["isComplete"] == True])
        scores2 = np.zeros([split.shape[0], ])
        scores2[valid_sample_index] = conf_score2[:, 1]
        split["model2_score"] = scores2
        save_path = os.path.join(project_rootdir, "data", "splits_with_predictions", split_name)
        print(split[["model1_score", "model2_score", "isComplete"]].head())
        split.to_csv(save_path, index=False)