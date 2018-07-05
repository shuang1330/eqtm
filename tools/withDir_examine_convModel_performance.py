import __init__path
from lib.model.finetune_model import examineModel
import os
import argparse

if __name__ == '__main__':
    # locally
    # project_rootdir = '/home/shuang/projects/development_eqtm'
    # calculon project root dir
    project_rootdir = '/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm'

    eqtm_path = os.path.join(project_rootdir, "data",
                             "eqtmZscores", "withExpressionTSSMethyCpgOverlapPromoter",
                             "2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_withPromoterOverlap.txt")
    model_name = 'ranfor'
    exclude = ['SNPName', 'SNPChr', 'PValue', 'SNPChrPos', 'ProbeName', 'ProbeChr',
               'ProbeCenterChrPos', 'CisTrans', 'SNPType', 'AlleleAssessed',
               'DatasetsWhereSNPProbePairIsAvailableAndPassesQC',
               'DatasetsZScores', 'DatasetsNrSamples',
               'IncludedDatasetsMeanProbeExpression',
               'IncludedDatasetsProbeExpressionVariance', 'HGNCName',
               'IncludedDatasetsCorrelationCoefficient', 'Meta-Beta (SE)',
               'Beta (SE)', 'FoldChange', 'FDR', 'checkChr', 'SNPName_ProbeName']
    keep = ['OverallZScore']
    train_path = eqtm_path
    test_path = train_path
    examineModel(model_name, train_path, test_path,
                 keep=keep, exclude=exclude, display=True)