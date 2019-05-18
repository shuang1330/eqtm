import os
import pandas as pd
import numpy as np


if __name__ == "__main__":
    et_gt_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/ORIGIN/et_gt.txt"
    et_gt = pd.read_csv(et_gt_filepath, sep="\t")

    def isConsistent(item):
        scores = [float(score) for score in item.split(";")]
        pos_num = 0
        neg_num = 0
        for score in scores:
            if score > 0:
                pos_num += 1
            else:
                neg_num += 1
        if pos_num*neg_num == 0:
            return True
        else:
            return False

    et_gt["isConsistent"] = [isConsistent(item) for item in et_gt["DatasetsZScores"]]
    print(et_gt["isConsistent"].describe())
    print(et_gt["DatasetsZScores"][et_gt["isConsistent"]==False].values)
