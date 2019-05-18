#TODO:add directories to startup files
import __init__path
from lib.preprocess.add_ExpressionTSSMethylation import (
read_expression_data,
read_methylation_data,
add_Expression_basedOnProbeName_toEQTMFile,
add_TSS_basedOnProbeName_toEQTMFile,
add_Methy_basedOnSNPName_toEQTMFile
)
import os
import argparse


def parse_args():
    '''
    Parse input arguments
    '''
    parser = argparse.ArgumentParser(description=
    'Add Expression an TSS and Methylation data to eqtm file.')
    parser.add_argument('--dataDir', dest='DATA_DIR',
                        help='the data root directory',
                        default=None, type=str)
    args = parser.parse_args()
    return args


if __name__=='__main__':

    args = parse_args()

    print('Called with args:\n')
    print(args)

    DATA_FOLDER = args.DATA_DIR
    OVERLAP_DIR = '/groups/umcg-gcc/tmp03/umcg-sli/eqtm/data/output'
    EQTM_DATADIR = os.path.join(DATA_FOLDER,'eqtmZscores')
    INPUT_FOLDER = os.path.join(EQTM_DATADIR,'ORIGIN')
    OUTPUT_EXPRESS = os.path.join(EQTM_DATADIR,'withExpressionData')
    OUTPUT_EXPRESS_TSS = os.path.join(EQTM_DATADIR,'withExpressionTSS')
    OUTPUT_EXPRESS_TSS_METHY = os.path.join(EQTM_DATADIR,
                                            'withExpressionTSSMethy')
    OUTPUT_EXPRESS_TSS_METHY_OVERLAP = os.path.join(EQTM_DATADIR,
                                            'final')

    # inputfiles
    eqtm_names = {'2017-12-09-eQTLsFDR-et0.0-flipped.txt':'2017-12-09-eQTLsFDR-et0.0-flipped',
                  '2017-12-09-eQTLsFDR-gt0.0-flipped.txt':'2017-12-09-eQTLsFDR-gt0.0-flipped',
                  'random20k_gt0.5.txt':'random20k_gt0.5'}
    expression_filepath = os.path.join(DATA_FOLDER,'features','meanVar',
                                       'expression-MeanAndVarianceRows.txt')
    expression = read_expression_data(expression_filepath)
    tss_filepath = os.path.join(DATA_FOLDER,'features','TSSDistance',
                                'Homo_sapiens.GRCh37.71.gtf')
    tss_raw = read_tss_data(tss_filepath)
    methy_filepath = os.path.join(DATA_FOLDER,'features','meanVar',
                                  'methylation-MeanAndVarianceRows.txt')
    methylation = read_expression_data(methy_filepath)

    # read expression data and add expression to eqtm files
    for filename in eqtm_names.keys():
        print("Processing file:",filename)
        eqtm_name = eqtm_names[filename]
        eqtm_filepath = os.path.join(INPUT_FOLDER,filename)
        output_express_filepath = os.path.join(OUTPUT_EXPRESS,
                                  eqtm_name+'_withExpression.txt')
        output_expressTss_filepath = os.path.join(OUTPUT_EXPRESS_TSS,
                                  eqtm_name+'_withExpressionTSS.txt')
        output_expressTssMethy_filepath = os.path.join(OUTPUT_EXPRESS_TSS_METHY,
                                  eqtm_name+'_withExpressionTSSMethy.txt')
        output_expressTssMethyOverlap_filepath = os.path.join(
                                  OUTPUT_EXPRESS_TSS_METHY_OVERLAP,
                                  eqtm_name+'_withExpressionTSSMethyOverlap.txt')
        overlapRatio_filepath = os.path.join(OVERLAP_DIR,
                                             eqtm_name,
                                             eqtm_name+'overlapRatio.txt')
        _ = add_Expression_basedOnProbeName_toEQTMFile(eqtm_filepath,
                                                    expression,
                                                    output_express_filepath)
        _ = add_TSS_basedOnProbeName_toEQTMFile(output_express_filepath,
                                                tss_raw,
                                                output_expressTss_filepath)
        _ = add_Methy_basedOnSNPName_toEQTMFile(output_expressTss_filepath,
                                                methylation,
                                                output_expressTssMethy_filepath)
        # _ = add_OverlapRatio_basedOnSNPName_toEQTMFile(
        #                                 output_expressTssMethy_filepath,
        #                                 overlapRatio_filepath,
        #                                 output_expressTssMethyOverlap_filepath)
