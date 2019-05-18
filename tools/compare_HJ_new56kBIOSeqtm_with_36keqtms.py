import os
import pandas as pd
import matplotlib.pyplot as plt


def compare_common_eqtms(data1_dict, data2_dict):

    data1_keys = set(data1_dict.keys())
    data2_keys = set(data2_dict.keys())
    common_eqtms_keys = list(data1_keys & data2_keys)
    print("In total there are %d common eqtms."%len(common_eqtms_keys))
    data1_eqtms_keys = list(data1_dict-common_eqtms_keys)
    data2_eqtms_keys = list(data2_dict-common_eqtms_keys)

    common_zscores_data1 = [item["OverallZScore"] for item in data1_dict[common_eqtms_keys]]
    common_zscores_data2 = [item["OverallZScore"] for item in data2_dict[common_eqtms_keys]]

    return [common_zscores_data1, common_zscores_data2], data1_eqtms_keys, data2_eqtms_keys


if __name__ == "__main__":
    data_dir = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/ORIGIN"
    new_eqtms_filepath = os.path.join(data_dir, "new_bios_eqtms0.05.txt")
    old_eqtms_filepath = os.path.join(data_dir, "et_gt.txt")

    new_eqtms = pd.read_csv(new_eqtms_filepath, sep="\t")
    old_eqtms = pd.read_csv(old_eqtms_filepath, sep="\t")

    new_eqtms["eqtm_combi"] = ["{}_{}".format(item[0], item[1]) for item in
                               zip(new_eqtms["SNPName"], new_eqtms["ProbeName"])]
    old_eqtms["eqtm_combi"] = ["{}_{}".format(item[0], item[1]) for item in
                               zip(old_eqtms["SNPName"], old_eqtms["ProbeName"])]

    new_eqtms_dict = new_eqtms[["eqtm_combi", "OverallZScore"]].set_index("eqtm_combi").T.to_dict()
    old_eqtms_dict = old_eqtms[["eqtm_combi", "OverallZScore"]].set_index("eqtm_combi").T.to_dict()

    [new_zscores, old_zscores], new_unique_eqtms, old_unique_eqtms = compare_common_eqtms(new_eqtms_dict,
                                                                                          old_eqtms_dict)

    # draw
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.scatter(new_zscores, old_zscores)
