import pandas as pd
import numpy as np
import os

PROJECT_ROOTDIR = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
INPUT_DIR = os.path.join(PROJECT_ROOTDIR, "data/eqtmZscores/withExpressionTSSMethyCpgOverlap")
OUTPUT_DIR = os.path.join(PROJECT_ROOTDIR, "data/eqtmZscores/withExpressionTSSMethyCpgOverlapPromoter")
feature_filepath = os.path.join(PROJECT_ROOTDIR, "data/features/geneOverlap/promoter_overlapRatio.txt")

feature = pd.read_csv(feature_filepath,sep='\t',index_col=0)


def add_overlap_to_eqtm(eqtm, overlap):
    overlap_dict = overlap.T.to_dict()
    print(eqtm.columns)
    for col in overlap.columns:
        col_name = col+'gene'
        eqtm[col_name] = [overlap_dict[row][col] for row in eqtm['ProbeName'].values]
    return eqtm


if __name__=='__main__':
    # for filename in ['2017-12-09-eQTLsFDR-et0.0-flipped_withExpressionTssMethy_blood.txt',
    #                  '2017-12-09-eQTLsFDR-et0.0-flipped_withExpressionTssMethy_nonBlood.txt',
    #                  '2017-12-09-eQTLsFDR-gt0.0-flipped_withExpressionTssMethy_blood.txt',
    #                  '2017-12-09-eQTLsFDR-gt0.0-flipped_withExpressionTssMethy_nonBlood.txt']:
    for filename in ["2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap.txt"]:
        filepath = os.path.join(INPUT_DIR, filename)
        if 'multi_class' not in filename:
            filecontent = pd.read_csv(filepath)
            eqtm = add_overlap_to_eqtm(filecontent, feature)
            eqtm.to_csv(os.path.join(OUTPUT_DIR, '.'.join(filename.split('.')[:-1])+'_withPromoterOverlap.txt'))
            print(eqtm.shape)
