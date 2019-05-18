import gzip
import os
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd


def read_tcga_file(filepath):
    eqtm_zscores = {}
    with gzip.open(filepath, "r") as f:
        info = f.readlines()[1:]
        for ind in range(len(info)):
            line = info[ind].decode("utf-8").strip().split("\t")
            eqtm = "{}_{}".format(line[1], line[4])
            zscore = float(line[10])
            assert eqtm not in eqtm_zscores
            eqtm_zscores[eqtm] = zscore
        f.close()
    return eqtm_zscores


def find_common_eqtms_zscores(origin_eqtms, tcga_eqtm_zscores):
    zscores_origin = []
    zscores_tcga = []
    with open(origin_eqtms, "r") as origin_f:
        info = origin_f.readlines()[1:]
        for ind in tqdm(range(len(info))):
            origin_line = info[ind].strip().split("\t")
            origin_eqtm = "{}_{}".format(origin_line[1], origin_line[4])
            if origin_eqtm in tcga_eqtm_zscores:
                zscores_origin.append(float(origin_line[10]))
                zscores_tcga.append(tcga_eqtm_zscores[origin_eqtm])
            else:
                continue
        origin_f.close()
    return zscores_origin, zscores_tcga


if __name__ == "__main__":
    origin_eqtms = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/ORIGIN/full_dataset.txt"
    tcga_dir = "/groups/umcg-gcc/tmp03/umcg-sli/tcga/cis_results"
    figure_dir = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/output_fig/common_eqtms"
    et_gt = pd.read_csv("/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/ORIGIN/et_gt.txt",
                        sep="\t")
    et_gt["eqtm_combi"] = ["{}_{}".format(item[0], item[1])
                           for item in zip(et_gt["SNPName"].values, et_gt["ProbeName"].values)]
    et_gt_dict = et_gt[["eqtm_combi", "OverallZScore"]].set_index("eqtm_combi").T.to_dict()
    print("Loaded et_gt dict.", et_gt_dict["cg07772999_ENSG00000189223"])

    for tcga_dataname in os.listdir(tcga_dir):
        tcga_filepath = os.path.join(tcga_dir, tcga_dataname, "eQTLsFDR.txt.gz")
        # try:
        if os.path.exists(tcga_filepath):
            print("Processing: ", tcga_dataname)
            tcga_eqtm_zscores = read_tcga_file(tcga_filepath)
            print("Loaded tcga dataset.")
            zscores_origin, zscores_tcga = find_common_eqtms_zscores(origin_eqtms, tcga_eqtm_zscores)
            print("Found common eqtms number: ", len(zscores_tcga))
            for eqtm_combi in tcga_eqtm_zscores.keys():
                if eqtm_combi in et_gt_dict:
                    zscores_origin.append(et_gt_dict[eqtm_combi]["OverallZScore"])
                    zscores_tcga.append(tcga_eqtm_zscores[eqtm_combi])
                else:
                    continue
            f = plt.figure()
            ax = f.add_subplot(111)
            ax.scatter(zscores_origin, zscores_tcga)
            plt.savefig(os.path.join(figure_dir, "{}.png".format(tcga_dataname)))
            plt.close()
        # except:
        #     print("Not Processing dataset: ", tcga_dataname)
        #     continue