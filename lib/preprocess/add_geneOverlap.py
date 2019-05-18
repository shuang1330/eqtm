import pandas as pd
import os

# PROJECT_ROOTDIR = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
PROJECT_ROOTDIR = "/home/shuang/projects/development_eqtm"
INPUT_DIR = os.path.join(PROJECT_ROOTDIR, "data/eqtmZscores/withExpressionTSSMethyCpgOverlap")
OUTPUT_DIR = os.path.join(PROJECT_ROOTDIR, "data/eqtmZscores/withExpressionTSSMethyCpgOverlapPromoter")
feature_filepath = os.path.join(PROJECT_ROOTDIR, "data/features/geneOverlap/promoter_overlapRatio.txt")

feature = pd.read_csv(feature_filepath, sep='\t', index_col=0)


def add_overlap_to_eqtm(eqtm, overlap):
    overlap_dict = overlap.T.to_dict()
    print(eqtm.columns)
    for col in overlap.columns:
        col_name = col+'_promoter'
        eqtm[col_name] = [overlap_dict[row][col] for row in eqtm['ProbeName'].values]
    return eqtm


if __name__=='__main__':
    for filename in os.listdir(INPUT_DIR):
        filepath = os.path.join(INPUT_DIR, filename)
        if 'multi_class' not in filename:
            filecontent = pd.read_csv(filepath, index_col=0)
            eqtm = add_overlap_to_eqtm(filecontent, feature)
            eqtm.to_csv(os.path.join(OUTPUT_DIR, filename[:-4]+'_PromoterOverlap.txt'))
            print(eqtm.shape)
