import os
import argparse
from lib.preprocess.createOverlapRatioMatrix import overlapMatrixPipeline

def initFolders(cpg_filename,
                PROJECT_DIR='/groups/umcg-gcc/tmp03/umcg-sli/eqtm_project/',
                DATA_ROOTDIR='/groups/umcg-wijmenga/tmp03/projects/eQTMPrediction'):
    """
    Initiate the file structures for saving results
    """
    print('The project root directory is: ',PROJECT_DIR)
    feature_folder = os.path.join(DATA_ROOTDIR,
                                  'features',
                                  'Roadmap',
                                  'consolidatedNarrowPeak')
    cpg_folder = os.path.join(DATA_ROOTDIR,'data','eqtmZscores')
    ratio_folder = os.path.join(DATA_ROOTDIR,'data','overlapRatio')
    # folder for intermediate results
    temp_folder = os.path.join(PROJECT_DIR,'data','temp')
    temp_cpgFolder = os.path.join(temp_folder,'cpg')
    # temp_featureFolder = os.path.join(temp_folder,'features')
    temp_output = os.path.join(temp_folder,'output')
    output_withcpgName = os.path.join(PROJECT_DIR,'data','output',cpg_filename)
    folder_lists = [ratio_folder,
                    temp_folder,
                    temp_cpgFolder,
                    temp_output,
                    output_withcpgName]
    for folder in folder_lists:
        if not os.path.exists(folder):
            os.makedirs(folder)
    print('Init folders.')
    return PROJECT_DIR,feature_folder,cpg_folder,ratio_folder,temp_cpgFolder,\
    temp_output,output_withcpgName

def parse_args():
    """
    Parse input arguments
    """
    parser = argparse.ArgumentParser(description=\
    'Build overlapping ratio matrix for all Roadmap features.')
    parser.add_argument('--rootDir', dest='PROJECT_DIR',
                        help='the project root directory',
                        default=None, type=str)
    parser.add_argument('--featList', dest='feature_list_filename',
                        help='a txt file with feature names',
                        default=None, type=str)
    parser.add_argument('--histList', dest='histone_list_filename',
                        help='a txt file with histone names',
                        default=None, type=str)
    parser.add_argument('--cellList', dest='cell_list_filename',
                        help='a txt file with all cell type names',
                        default=None, type=str)
    parser.add_argument('--cpg', dest='cpg_filename',
                        help='a csv file with cpg sites and eqtm zscores',
                        default=None, type=str)
    parser.add_argument('--bedtools', dest='execute_bedtoolFilename',
                        help='a bash command file for determining the overlap',
                        default=None, type=str)
    args = parser.parse_args()
    return args

def startBuild(PROJECT_DIR,
               feature_folder,
               histone_folder,
               temp_output,
               output_withcpgName,
               temp_cpgFolder,
               cpg_folder,
               ratio_folder,
               feature_list_filename='feature_list',
               histone_list_filename='histone_list',
               cell_list_filename='cell_list',
               execute_bedtoolFilename='findOverlap',
               cpg_filename='bonder-eQTMsFDR0.0-CpGLevel-split'):
    # print(cpg_folder)
    overlapMatrixPipeline(PROJECT_DIR,
                          execute_bedtoolFilename,
                          feature_folder,
                          feature_list_filename,
                          histone_folder,
                          histone_list_filename,
                          temp_cpgFolder,
                          cpg_folder,
                          cpg_filename,
                          temp_output,
                          output_withcpgName,
                          ratio_folder,
                          cell_list_filename)


if __name__=='__main__':

    args = parse_args()

    print('Called with args:')
    print(args)

    PROJECT_DIR,feature_folder,cpg_folder,ratio_folder,temp_cpgFolder,\
    temp_output,output_withcpgName = \
    initFolders(args.cpg_filename)

    # print(cpg_folder)
    startBuild(PROJECT_DIR,
               feature_folder,
               feature_folder,
               temp_output,
               output_withcpgName,
               temp_cpgFolder,
               cpg_folder,
               ratio_folder)
