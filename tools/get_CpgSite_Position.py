import pandas as pd
import numpy as np
import argparse


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_path", dest='input_path', type=str)
    parser.add_argument("--bedtoolFormat_savepath", dest='bedtoolFormat_savepath', type=str)
    parser.add_argument("--annotation_path", dest="annotation_path",
                        default='/groups/umcg-bios/tmp03/projects/2018-methylation/input/helpfiles/Illumina450K_MQtlMappingFile_MJB.txt.gz')
    return parser.parse_args()


def add_position(input, savepath, annotation_dict):
    findPosition = lambda x:[annotation_dict[x]['ChrStart'], annotation_dict[x]['ChrEnd']] if x in annotation_dict else np.nan
    start_end_positions = [findPosition(cpgsite) for cpgsite in input['SNPName']]
    input['start_position'] = [item[0] for item in start_end_positions]
    input['end_position'] = [item[1] for item in start_end_positions]
    input['chromosome'] = ['chr{}'.format(item) for item in input['SNPChr']]
    input[['chromosome', 'start_position', 'end_position', 'SNPName']].to_csv(savepath, header=False, index=False, sep='\t')
    return input


def main():
    arguments = parser()
    annotation_file = arguments.annotation_path
    annotation = pd.read_csv(annotation_file, sep='\t', compression='gzip')
    annotation_dict = annotation.set_index('HT12v4-ArrayAddress').T.to_dict()
    input = pd.read_csv(arguments.input_path, sep='\t')
    input_annotated = add_position(input, savepath=arguments.bedtoolFormat_savepath, annotation_dict=annotation_dict)
    return input_annotated


if __name__ == '__main__':
    _ = main()
