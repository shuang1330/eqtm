import __init__path
import os
from lib.read.read_roadmap_features import roadmap
from lib.read.read_tss_file import TSS_file
from lib.read.read_eqtm import eqtm
from lib.read.read_allCpgFile import cpgFile
from lib.preprocess.call_bedtools import callBedtools
import argparse


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project_rootdir", dest="project_rootdir", type=str)
    parser.add_argument("--data_rootdir", dest="data_rootdir", type=str)
    parser.add_argument("--feature_filetype", dest="feature_filetype", type=str)
    parser.add_argument("--input_type", dest="input_type", type=str,
                        help="Either promoter or cpg*")
    parser.add_argument("--cpg_filename", dest="cpg_filename", type=str)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse()
    print("Called with arguments: ", args)
    project_rootdir = args.project_rootdir
    data_rootdir = args.data_rootdir
    feature_filetype = args.feature_filetype
    roadmap_feature = roadmap(data_rootdir=data_rootdir, feature_type=feature_filetype)
    temp_res_dirpath = os.path.join(project_rootdir, "data", "temp", "meta")
    featureList_filepath = os.path.join(temp_res_dirpath, "feature_list.txt")

    if args.input_type == "promoter":
        # input file as T
        tss_dirpath = os.path.join(project_rootdir, "data", "features", "TSSDistance")
        tss_raw_filepath = os.path.join(tss_dirpath, "Homo_sapiens.GRCh37.71.gtf")
        tss_file = TSS_file(tss_raw_filepath, sep="\t")
        gene_startEndSite_savepath = os.path.join(tss_dirpath, 'gene_startEndSite.txt')
        _ = tss_file.find_allGene_startEndSite(gene_startEndSite_savepath)

        input_bedtoolsformat_filepath = os.path.join(tss_dirpath, "promoter_startEndSite_bedtoolsFormat.txt")
        _ = tss_file.promoter_into_bedtoolsFormat(input_bedtoolsformat_filepath)
    else:
        cpg_filename = args.cpg_filename
        input_bedtoolsformat_filepath = os.path.join(project_rootdir, "data", "eqtmZscores", "allCpgs", cpg_filename)

    # bedtools bash file
    findOverlap_path = os.path.join(project_rootdir, "tools", "findOverlap.sh")
    callBedtools(findOverlap_path,
                 project_rootdir,
                 roadmap_feature.feature_dirpath,
                 roadmap_feature.readFeatureList(featureList_filepath),
                 feature_filetype,
                 input_bedtoolsformat_filepath,
                 args.input_type)
