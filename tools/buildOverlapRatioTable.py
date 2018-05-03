import os
import argparse
from lib.preprocess.createOverlapRatioMatrix import overlapMatrixPipeline
from lib.path.localPaths import directoryLocal

def parse_args():
    """
    Parse input arguments
    """
    parser = argparse.ArgumentParser(description=
    'Build overlapping ratio matrix for all Roadmap features.')
    parser.add_argument('--rootDir', dest='PROJECT_DIR',
                        help='the project root directory',
                        default=None, type=str)
    parser.add_argument('--cpg', dest='cpg_filename',
                        help='a csv file with cpg sites and eqtm zscores',
                        default=None, type=str)
    parser.add_argument('--bedtools', dest='execute_bedtoolFilename',
                        help='a bash command file for determining the overlap',
                        default=None, type=str)
    args = parser.parse_args()
    return args

if __name__=='__main__':

    args = parse_args()

    print('Called with args:\n')
    print(args)

    dirs = directoryLocal(args.PROJECT_DIR,args.cpg_filename)
    dirs.initializeFolders()

    overlapMatrixPipeline(dirs, args.execute_bedtoolFilename)
