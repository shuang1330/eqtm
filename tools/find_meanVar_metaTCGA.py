import gzip
import os
import numpy as np
import pandas as pd
from tqdm import tqdm


def read_gzip_file(filepath, common_individuals=None):
    meanVar_dic = {}
    with gzip.open(filepath, "r") as f:
        samples=[item.strip().split(".")[0] for item in f.readline().decode("utf-8").split("\t")]
        print(samples)
        info = f.readlines()[1:]
        for ind in tqdm(range(len(info))):
            line = info[ind]
            values = line.decode("utf-8").split("\t")
            probename = values[0]
            if common_individuals:
                common_inds = [samples.index(individual) for individual in common_individuals]
                rna_counts = [float(values[ind].strip()) for ind in common_inds]
            else:
                rna_counts = [float(value.strip()) for value in values[1:]]
            meanVar_dic[probename] = rna_counts
        return meanVar_dic


def rna():
    rna_dir = "/groups/umcg-gcc/tmp03/umcg-sli/tcga/meta_analysis/rna"
    rna_suffix = "-beta.SampleSelection.ProbesWithZeroVarianceRemoved.TMM.txt.Log2Transformed.ProbesCentered.SamplesZTransformed.5PCAsOverSamplesRemoved.txt.gz"
    common_inds_record_dir = "/groups/umcg-gcc/tmp03/umcg-sli/tcga/meta_analysis/common_individuals"
    datasets_meanVar = {}
    probenames = []
    for filename in os.listdir(rna_dir):
        if os.path.exists(os.path.join(common_inds_record_dir, filename[:-len(rna_suffix)]+".txt")):
            print(filename[:-len(rna_suffix)])
            with open(os.path.join(common_inds_record_dir, filename[:-len(rna_suffix)]+".txt"), "r") as f:
                common_inds = [item.strip() for item in f.readlines()]
            rna_path = os.path.join(rna_dir, filename)
            dataset_meanVar = read_gzip_file(rna_path, common_inds)
            datasets_meanVar[filename.split(".")[0]] = dataset_meanVar
            probenames.append(set(dataset_meanVar.keys()))

    common_probenames = set()
    for i in range(len(probenames)):
        if i == 0:
            common_probenames = probenames[i]
        else:
            common_probenames = common_probenames | probenames[i]
    probe_meanVar = {}

    for probename in common_probenames:
        rna_counts = []
        for dataset in datasets_meanVar:
            i = 0
            if probename in datasets_meanVar[dataset]:
                if i == 0:
                    rna_counts = datasets_meanVar[dataset][probename]
                elif i > 0:
                    for count in datasets_meanVar[dataset][probename]:
                        rna_counts.append(count)
        # rna_counts = [item for dataset in datasets_meanVar for item in datasets_meanVar[dataset][probename]]
        mean = np.mean(rna_counts)
        var = np.var(rna_counts)
        probe_meanVar[probename] = [mean, var]

    probe_meanVar_pd = pd.DataFrame.from_dict(probe_meanVar,
                                              orient="index",
                                              columns=["Mean", "Var"])
    probe_meanVar_pd.to_csv("/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/features/meanVar/metaTCGA_rna_meanVar.txt")
    return None


def methy():
    methy_dir = "/groups/umcg-gcc/tmp03/umcg-sli/tcga/methy_log_norm_pca"
    methy_suffix = "-mval.QuantileNormalized.ProbesCentered.SamplesZTransformed.5PCAsOverSamplesRemoved.txt.gz"
    common_inds_record_dir = "/groups/umcg-gcc/tmp03/umcg-sli/tcga/meta_analysis/common_individuals"
    datasets_meanVar = {}
    probenames = []
    for filename in os.listdir(methy_dir):
        if os.path.exists(os.path.join(common_inds_record_dir, filename[:-len(methy_suffix)]+".txt")):
            print(filename[:-len(methy_suffix)])
            if filename.endswith(".QuantileNormalized.ProbesCentered.SamplesZTransformed.5PCAsOverSamplesRemoved.txt.gz"):
                with open(os.path.join(common_inds_record_dir, filename[:-len(methy_suffix)] + ".txt"), "r") as f:
                    common_inds = [item.strip() for item in f.readlines()]
                methy_path = os.path.join(methy_dir, filename)
                dataset_meanVar = read_gzip_file(methy_path, common_inds)
                datasets_meanVar[filename.split(".")[0]] = dataset_meanVar
                probenames.append(set(dataset_meanVar.keys()))


    common_probenames = set()
    for i in range(len(probenames)):
        if i == 0:
            common_probenames = probenames[i]
        else:
            common_probenames = common_probenames | probenames[i]
    probe_meanVar = {}

    for probename in common_probenames:
        rna_counts = []
        for dataset in datasets_meanVar:
            i = 0
            if probename in datasets_meanVar[dataset]:
                if i == 0:
                    rna_counts = datasets_meanVar[dataset][probename]
                elif i > 0:
                    for count in datasets_meanVar[dataset][probename]:
                        rna_counts.append(count)
        # rna_counts = [item for dataset in datasets_meanVar for item in datasets_meanVar[dataset][probename]]
        mean = np.mean(rna_counts)
        var = np.var(rna_counts)
        probe_meanVar[probename] = [mean, var]

    probe_meanVar_pd = pd.DataFrame.from_dict(probe_meanVar,
                                              orient="index",
                                              columns=["Mean", "Var"])
    probe_meanVar_pd.to_csv("/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/features/meanVar/metaTCGA_methy_meanVar.txt")
    return None


if __name__ == "__main__":
    rna()
    methy()