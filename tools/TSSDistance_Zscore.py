import os
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":
    project_rootdir = "/home/shuang/projects/development_eqtm"
    eqtm_dirpath = os.path.join(project_rootdir,
                                 "data",
                                 "eqtmZscores",
                                 "withExpressionTSS")
    zscores = []
    tss_distances = []
    for eqtm_filename in os.listdir(eqtm_dirpath):
        eqtm_filepath = os.path.join(eqtm_dirpath,
                                     eqtm_filename)
        eqtm = pd.read_csv(eqtm_filepath)
        zscores.append(eqtm["OverallZScore"])
        tss_distances.append(eqtm["TssDistance"])

    all_zscores = [item for sublist in zscores for item in sublist]
    all_tss = [item for sublist in tss_distances for item in sublist]

    plt.scatter(all_zscores, all_tss, alpha=0.5)
    plt.xlabel("Zscores")
    plt.ylabel("TSS Distance")
    plt.savefig("./development_eqtm/output_fig/zscore_tssDistance.png")
    plt.show()
