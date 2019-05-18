import os
import pandas as pd
import numpy as np

def add_OverlapRatio_basedOnSNPName_toEQTMFile(eqtm_filepath,
                                               overlap_filepath,
                                               eqtm_savepath):
    '''
    read overlapRatio file and add them to eqtm file
    '''
    overlap = pd.read_csv(overlap_filepath,sep=',',index_col=0)
    overlap_dict = overlap.T.to_dict()
    eqtm = pd.read_csv(eqtm_filepath,index_col=0)
    for col in overlap.columns:
        eqtm[col] = [overlap_dict[row][col] for row in eqtm['SNPName'].values]
    eqtm.to_csv(eqtm_savepath)
    print(eqtm.columns)
    print('Saved to path:',eqtm_savepath)
    return eqtm

if __name__=='__main__':
    inputdir = '/home/shuang/projects/development_eqtm/data/eqtmZscores/withExpressionTSSMethy'
    overlapDir = '/home/shuang/projects/development_eqtm/data/eqtmZscores/overlapMatrix'
    saveDir = '/home/shuang/projects/development_eqtm/data/eqtmZscores/withExpressionTSSMethyCpgOverlap'
    save_geneOverlapDir = '/home/shuang/projects/development_eqtm/data/eqtmZscores/withExpressionTSSMethyCpgOverlapGene'
    eqtm_name = ['2017-12-09-eQTLsFDR-et0_withExpressionTssMethy',
                 '2017-12-09-eQTLsFDR-gt0_withExpressionTssMethy']

    _ = add_OverlapRatio_basedOnSNPName_toEQTMFile(
    os.path.join(inputdir,'2017-12-09-eQTLsFDR-et0_withExpressionTssMethy.txt'),
    os.path.join(overlapDir,'2017-12-09-eQTLsFDR-et0.0-flipped_blood.txt'),
    os.path.join(saveDir,'2017-12-09-eQTLsFDR-et0.0-flipped_withExpressionTssMethy_blood.txt'))

    _ = add_OverlapRatio_basedOnSNPName_toEQTMFile(
    os.path.join(inputdir,'2017-12-09-eQTLsFDR-et0_withExpressionTssMethy.txt'),
    os.path.join(overlapDir,'2017-12-09-eQTLsFDR-et0.0-flipped_nonBlood.txt'),
    os.path.join(saveDir,'2017-12-09-eQTLsFDR-et0.0-flipped_withExpressionTssMethy_nonBlood.txt'))

    _ = add_OverlapRatio_basedOnSNPName_toEQTMFile(
    os.path.join(inputdir,'2017-12-09-eQTLsFDR-gt0_withExpressionTssMethy.txt'),
    os.path.join(overlapDir,'2017-12-09-eQTLsFDR-gt0.0-flipped_blood.txt'),
    os.path.join(saveDir,'2017-12-09-eQTLsFDR-gt0.0-flipped_withExpressionTssMethy_blood.txt'))

    _ = add_OverlapRatio_basedOnSNPName_toEQTMFile(
    os.path.join(inputdir,'2017-12-09-eQTLsFDR-gt0_withExpressionTssMethy.txt'),
    os.path.join(overlapDir,'2017-12-09-eQTLsFDR-gt0.0-flipped_nonBlood.txt'),
    os.path.join(saveDir,'2017-12-09-eQTLsFDR-gt0.0-flipped_withExpressionTssMethy_nonBlood.txt'))
