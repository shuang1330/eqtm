import __init__path
import os
from lib.read.read_roadmap_features import roadmap
import numpy as np
from lib.read.read_data import load_data
from lib.read.read_allCpgFile import cpgFile


if __name__ == "__main__":
    project_rootdir = "/groups/umcg-gcc/tmp04/umcg-sli/development_eqtm"
    test_eqtm_ratio_filepath = os.path.join(project_rootdir, "data",
                                            "eqtmZscores",
                                            "withExpressionTSSMethyCpgOverlapGene",
                                            "2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_withGeneOverlap.txt")
    feature_list_filepath = os.path.join(project_rootdir,
                                         "data",
                                         "temp",
                                         "meta",
                                         "feature_list.txt")
    histone_list_filepath = os.path.join(project_rootdir,
                                         "data",
                                         "temp",
                                         "meta",
                                         "histone_list.txt")

    roadmap_features = roadmap()
    feature_list = roadmap_features.readFeatureList(feature_list_filepath)
    histone_list = roadmap_features.readHistoneList(histone_list_filepath)
    # test_eqtm_ratio = pd.read_csv(test_eqtm_ratio_filepath, index_col=0, sep="\t")

    test_eqtm_ratio = cpgFile(allCpg_filepath=test_eqtm_ratio_filepath)
    exclude, keep = cpgFile.cpgFile_DefaultExcludeAndKeep()
    data = load_data(test_eqtm_ratio.allCpg_filepath, keep=keep, exclude=exclude)
