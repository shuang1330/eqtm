class cpgFile(object):

    def __init__(self,
                 allCpg_filepath=""):
        self.allCpg_filepath = allCpg_filepath
        self.allCpg_bedtoolsFormat_filepath = ""

    def readCpgFileInBedToolsFormat(self, allCpg_bedtoolsFormat_filepath):
        self.allCpg_bedtoolsFormat_filepath = allCpg_bedtoolsFormat_filepath
        cpg_bedtoolFormat = open(self.allCpg_bedtoolsFormat_filepath, 'w')
        with open(self.allCpg_filepath, 'r') as cpgFile:
            cpgInfo = [row.strip().split('\t') for row in cpgFile.readlines()[1:]]
            for cpg in cpgInfo:
                cpg_name = cpg[0]
                cpg_chr, start, end = 'chr'+str(cpg[1]), int(cpg[2])-25, int(cpg[2])+25
                cpg_bedtoolFormat.write('{}\t{}\t{}\t{}\n'.format(cpg_chr, start, end, cpg_name))
        cpg_bedtoolFormat.close()
        return cpg_bedtoolFormat

    @staticmethod
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
        return exclude, keep

