import __init__path
import os
import argparse
from lib.read.read_roadmap_features import roadmap
from lib.read.read_tss_file import TSS_file
from lib.read.read_eqtm import eqtm
from lib.read.read_allCpgFile import cpgFile
from lib.preprocess.call_bedtools import callBedtools

if __name__ == "__main__":
    project_rootdir = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
    data_rootdir = "/groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap"
    feature_filetype = "consolidatedImputedGappedPeak"
    roadmap_feature = roadmap(data_rootdir=data_rootdir, feature_type=feature_filetype)
    temp_res_dirpath = os.path.join(project_rootdir, "data", "temp", "meta")
    featureList_filepath = os.path.join(temp_res_dirpath, "feature_list.txt")

    # input file as TSS
    tss_dirpath = os.path.join(project_rootdir, "data", "features", "TSSDistance")
    tss_raw_filepath = os.path.join(tss_dirpath,
                                    "Homo_sapiens.GRCh37.71.gtf")
    tss_file = TSS_file(tss_raw_filepath, sep="\t")
    gene_startEndSite_savepath = os.path.join(tss_dirpath,
                                              'gene_startEndSite.txt')
    _ = tss_file.find_allGene_startEndSite(gene_startEndSite_savepath)

    promoter_bedtoolsformat_filepath = os.path.join(tss_dirpath,
                                                    "promoter_startEndSite_bedtoolsFormat.txt")
    _ = tss_file.promoter_into_bedtoolsFormat(promoter_bedtoolsformat_filepath)

    # bedtools bash file
    findOverlap_path = os.path.join(project_rootdir,
                                    "tools",
                                    "findOverlap.sh")
    callBedtools(findOverlap_path,
                 project_rootdir,
                 roadmap_feature.feature_dirpath,
                 roadmap_feature.readFeatureList(featureList_filepath),
                 feature_filetype,
                 promoter_bedtoolsformat_filepath,
                 "promoter")
