import pandas as pd
import os
import numpy as np
import sklearn
from math import copysign
PROJECT_DIR = '/home/shuang/projects/eqtm'

if __name__=='__main__':
    annotation_file = os.path.join(PROJECT_DIR,
                                   'mj_data',
                                   'Annotation450k_AdditionMJ_v10.txt')
    cg_level_file = os.path.join(PROJECT_DIR,
                                 'mj_data',
                                 'eQTLSNPsFDR0.05-SNPLevel-flipped.txt')

    # read zscore from cg file
    def read_zscore(cg_level_file):
        cg_level = pd.read_csv(cg_level_file,sep='\t')
        return cg_level.set_index(['SNPName'])['OverallZScore'].to_dict()
    zscore_dict = read_zscore(cg_level_file)

    # # overlap check
    # def annotation_cg_set():
    #     annotation = pd.read_csv(annotation_file,sep='\t')
    #     return annotation['HT12v4.ArrayAddress'].unique()
    # annotation_cgs = annotation_cg_set()
    # print(len(list(set(zscore_dict.keys()).intersection(set(annotation_cgs)))))

    # add direction to annotation
    annotation = pd.read_csv(annotation_file,sep='\t')
    def find_zscore(row):
        if row in zscore_dict:
            return zscore_dict[row]
        else:
            return None
    annotation['zscore'] = annotation['HT12v4.ArrayAddress'].apply(find_zscore)
    annotation = annotation[pd.notnull(annotation['zscore'])]
    print(annotation.shape)
    def sign(row):
        return copysign(1,row)
        annotation['direction'] = annotation['zscore'].apply(sign)

    # save the rows with direction and ascore
    annotation.to_csv(os.path.join(PROJECT_DIR,
                                   'mj_data',
                                   'Anno_Value_Direction.csv'),
                      sep='\t',index=False)
