
def cpgFile_DefaultExcludeAndKeep():
    exclude = ['SNPName','SNPChr','PValue','SNPChrPos','ProbeName','ProbeChr',
               'ProbeCenterChrPos','CisTrans', 'SNPType', 'AlleleAssessed',
               'DatasetsWhereSNPProbePairIsAvailableAndPassesQC',
               'DatasetsZScores', 'DatasetsNrSamples',
               'IncludedDatasetsMeanProbeExpression',
               'IncludedDatasetsProbeExpressionVariance', 'HGNCName',
               'IncludedDatasetsCorrelationCoefficient', 'Meta-Beta (SE)',
               'Beta (SE)','FoldChange', 'FDR','checkChr','SNPName_ProbeName',
               'TssSite', 'chr']
    keep = ['OverallZScore']
    return exclude,keep
