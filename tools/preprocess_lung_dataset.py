import os
import pandas as pd
import numpy as np


if __name__ == "__main__":
    project_rootdir = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
    cpg_anno = pd.read_csv("/groups/umcg-gcc/tmp03/umcg-sli/tcga/anno_files/methy_anno.txt",
                           sep="\t")
    # gene_anno = pd.read_csv("/groups/umcg-gcc/tmp03/umcg-sli/tcga/anno_files/gene_anno.txt",
    #                         sep="\t")[["Symbol", "ArrayAddress"]]
    # gene_anno["Symbol_name"] = [item.strip() for item in gene_anno["Symbol"]]
    # gene_anno_dict = gene_anno[["Symbol_name","ArrayAddress"]].set_index("Symbol_name").T.to_dict()
    # print(gene_anno_dict["SORT1"])
    # # raise NotImplementedError
    # gene_Sites = pd.read_csv("/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/geneSites/all_geneSites.txt")[["geneName", "TssSite"]].set_index("geneName").T.to_dict()
    # print(gene_Sites["ENSG00000134243"])
    # lung_eqtm_filepath = os.path.join(project_rootdir, "data", "eqtmZscores",
    #                                   "ORIGIN", "lung_eqtms.csv")
    # lung_columns = ["chr", "HT12v3_ArrayAddress", "CpgSite", "GeneName_ExpressionProbe",
    #                 "GeneName_CpgSite", "Spearmanr", "PermutationP"]
    # lung_eqtms = pd.read_csv(lung_eqtm_filepath, sep=",", names=lung_columns,
    #                          dtype={"GeneName_ExpressionProbe": np.object})
    # print(lung_eqtms.head())
    # print(lung_eqtms["GeneName_ExpressionProbe"].head())
    #
    # def find_symbol(item):
    #     return str(item).split(";")[0] if item else np.nan
    # lung_eqtms["Symbol"] = [find_symbol(item)
    #                         for item in lung_eqtms["GeneName_ExpressionProbe"]]
    # print(lung_eqtms["Symbol"].head())
    # # print(lung_eqtms[lung_eqtms["GeneName_ExpressionProbe"].isnull()])
    #
    # def find_geneName(symbol):
    #     return gene_anno_dict[symbol]["ArrayAddress"] if symbol in gene_anno_dict else np.nan
    # lung_eqtms["GeneName"] = [find_geneName(symbol)
    #                           for symbol in lung_eqtms["Symbol"]]
    # print(lung_eqtms["GeneName"].head())
    #
    # def find_tssSite(name):
    #     return gene_Sites[name]["TssSite"] if name in gene_Sites else np.nan
    # lung_eqtms["TssSite"] = [find_tssSite(name) for name in lung_eqtms["GeneName"]]
    # print(lung_eqtms.head())
    # lung_eqtms.to_csv(os.path.join(project_rootdir, "data", "eqtmZscores",
    #                                "ORIGIN", "lung_eqtms_annotated.csv"),
    #                   sep="\t")
    lung_eqtms_annotated = pd.read_csv(os.path.join(project_rootdir, "data", "eqtmZscores",
                                                    "ORIGIN", "lung_eqtms_annotated.csv"),
                                       sep="\t", index_col=0)
    print(lung_eqtms_annotated.head())
    print(lung_eqtms_annotated.shape)
    lung_eqtms_annotated = lung_eqtms_annotated.dropna()
    print("After dropping nan values, dataset left with: ", lung_eqtms_annotated.shape)
    methy_anno_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/tcga/anno_files/methy_probecenter.txt"
    methy_anno = pd.read_csv(methy_anno_filepath, sep="\t", index_col=0)[["HT12v4-ArrayAddress", "ChrStart"]].set_index("HT12v4-ArrayAddress").T.to_dict()
    print(methy_anno["cg13912858"])
    lung_eqtms_annotated["probeCenter"] = [methy_anno[cpgname]["ChrStart"] for cpgname in lung_eqtms_annotated["CpgSite"]]
    lung_eqtms_annotated["chromosome"] = ["chr{}".format(item) for item in lung_eqtms_annotated["chr"]]
    lung_eqtms_annotated["startSite"] = lung_eqtms_annotated["probeCenter"] - 25
    lung_eqtms_annotated["endSite"] = lung_eqtms_annotated["probeCenter"] + 25
    print(lung_eqtms_annotated[["chromosome", "startSite", "endSite", "CpgSite"]].head())
    save_bedtoolsFormat_path = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/allCpgs/lung_cpg_bedtoolsFormat.txt"
    lung_eqtms_annotated[["chromosome", "startSite", "endSite", "CpgSite"]].to_csv(save_bedtoolsFormat_path, index=False, header=False, sep="\t")
    print("Finished.")
