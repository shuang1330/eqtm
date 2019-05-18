import __init__path
from lib.read.read_eqtm import eqtm
import os
import argparse


def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project_rootdir", dest="project_rootdir", type=str)
    parser.add_argument("--eqtm_filename", dest="eqtm_origin_filename", type=str)
    parser.add_argument("--cond_colName", dest="cond_colName", type=str)
    parser.add_argument("--cond_operator", dest="cond_operator", type=str)
    parser.add_argument("--cond_threshold", dest="cond_threshold", type=float)
    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = parse()
    print("Call with args:", args)
    project_rootdir = args.project_rootdir
    eqtm_origin_filepath = os.path.join(project_rootdir,
                                        "data",
                                        "eqtmZscores",
                                        "ORIGIN",
                                        args.eqtm_origin_filename)
    cond_colName = args.cond_colName
    cond_operator = args.cond_operator
    cond_threshold = args.cond_threshold
    cond_name = "eqtm_{}_{}_{}_bedtoolsFormat.txt".\
        format(cond_colName, cond_operator, cond_threshold)
    bedtools_savepath = os.path.join(project_rootdir,
                                     "data",
                                     "eqtmZscores",
                                     "allCpgs",
                                     cond_name)

    eqtm_file = eqtm(eqtm_origin_filepath)
    eqtm_file.read_eqtmFile(
        cond_colName=cond_colName,
        cond_threshold=cond_threshold,
        cond_operator=cond_operator
    )
    print("Read eqtmfile with condition: ",
          cond_colName, cond_operator, cond_threshold)
    print("For examination: ",
          eqtm_file.data[:2])
    eqtm_file.save_in_bedtoolsFormat(bedtools_savepath)
    print("Selected eqtm file transformed into bedtools format and saved in {}".
          format(bedtools_savepath))
