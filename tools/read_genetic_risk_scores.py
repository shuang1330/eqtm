import pandas as pd
import os

if __name__ == '__main__':
    project_rootdir = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
    data_path = os.path.join(project_rootdir, "data", "eqtmZscores", "allCpgs", "cpg_genetic_risk_factors_cleaned_traits_standardized.txt")
    data = pd.read_csv(data_path, sep="\t")
    print(data.head())
    print(data.columns)
    print(data["traits"].head())

    bedtools_format_data = pd.DataFrame(data=0, columns=["chr", "start", "end", "SNP",
                                                         "allele", "minor_allele", "traits",
                                                         "pubids"],
                                        index=data.index)
    bedtools_format_data["chr"] =