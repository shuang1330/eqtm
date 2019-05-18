import pandas as pd
import os

if __name__ == '__main__':
    project_rootdir = "/home/shuang/projects/development_eqtm"
    data_dirpath = os.path.join(project_rootdir, "data", "eqtmZscores", "ORIGIN")
    et_filepath = os.path.join(data_dirpath, "2017-12-09-eQTLsFDR-et0.0-flipped.txt")
    gt_filepath = os.path.join(data_dirpath, "2017-12-09-eQTLsFDR-gt0.0-flipped.txt")

    et_file = pd.read_csv(et_filepath, sep="\t")
    gt_file = pd.read_csv(gt_filepath, sep="\t")

    print(et_file.head())
    print(et_file["FDR"][et_file["OverallZScore"]>0].count())