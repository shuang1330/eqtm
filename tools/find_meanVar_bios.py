import gzip
import os
from tqdm import tqdm
import numpy as np


def read_gz_file_and_calculate_meanVar(gz_path, save_path):
    with open(save_path, "w") as write_f:
        with gzip.open(gz_path, "rb") as read_f:
            info = read_f.readlines()[1:]
            for ind in range(len(info)):
                line = info[ind].decode("utf-8").split("\t")
                cpg = line[0]
                numbers = [float(item.strip()) for item in line[1:]]
                mean = np.mean(numbers)
                var = np.var(numbers)
                write_f.write("{}\t{}\t{}\n".format(cpg, mean, var))
    print("Saved meanVar in: ", save_path)


if __name__ == "__main__":
    methylation_path = "/groups/umcg-bios/tmp03/projects/2017-meQTL-BonderEtAl/input/methylationdata/CODAM_NTR_LLS_LLD_RS_BBMRI_450K.QuantileNormalized.ProbesCentered.SamplesZTransformed.22PCAsOverSamplesRemoved.txt.gz"
    save_path = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/features/meanVar/bios_methylation_meanVar.txt"

    read_gz_file_and_calculate_meanVar(methylation_path, save_path)