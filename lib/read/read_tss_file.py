import os
import pandas as pd
import numpy as np


class TSS_file(object):

    def __init__(self, filepath, sep='\t'):
        self.tss_filepath = filepath
        self._sep=sep
        self.gene_startEndSite_filepath = ""
        self.gene_bedtoolsFormat_filepath = ""
        self.promoter_bedtoolsFormat_filepath = ""

    def read_tss_data(self):
        """
        read tss file from tss_filepath
        :param tss_filepath:
        :param sep:
        :return:
        """

        colnames = ['chr', 'regionFunction', 'regionType', 'startSite',
                    'endSite', 'score', 'strand', 'sthunknown', 'geneInfo']
        dtype = {'chr': object, 'regionFunction': object, 'regionType': object,
                 'startSite': int, 'endSite': int, 'score': object,
                 'strand': object, 'sthunknown': object, 'geneInfo': object}
        tss_raw = pd.read_csv(self.tss_filepath,
                              sep=self._sep,
                              header=None,
                              names=colnames,
                              dtype=dtype)
        print("Loaded tss raw file", flush=True)
        print(tss_raw.head(), flush=True)
        return tss_raw

    def find_allGene_startEndSite(self, save_filepath):
        tss_raw = self.read_tss_data()

        # extract gene name from geneInfo for tss file
        def findGeneName(item):
            item = [thing for thing in list(filter(None, item.strip().split(";")))][0]
            name = item.replace('"', '').replace(';', '').strip().split(' ')[1]
            return name
        tss_raw['geneName'] = tss_raw['geneInfo'].apply(findGeneName)

        # find the tss sites for each gene in the tss file
        groupbyTss = tss_raw.groupby('geneName').agg({
            'chr': lambda x: x.unique(),
            'startSite': np.min,
            'endSite': np.max,
            'strand': lambda x: x.unique()
        })

        def findTssSite(series):
            if series[3] == '-':
                return series[2]
            else:
                return series[1]
        groupbyTss['TssSite'] = groupbyTss.apply(findTssSite, axis=1)
        groupbyTss.to_csv(save_filepath)
        self.gene_startEndSite_filepath = save_filepath
        print("Gene-start-end position file saved to path %s" % save_filepath, flush=True)
        return groupbyTss

    def gene_into_bedtoolsFormat(self, save_filepath):
        self.gene_bedtoolsFormat_filepath = save_filepath
        if self.gene_startEndSite_filepath == "":
            raise IOError("Try using find_allGene_startEndSite first.")
        gene_startEndSite = [item.split(',') for item in
                             open(self.gene_startEndSite_filepath, 'r').readlines()[1:]]
        bedtoolsformat_file = open(self.gene_bedtoolsFormat_filepath, 'w')
        for item in gene_startEndSite:
            if item[4] == '-':
                startSite = int(item[3])
                endSite = int(item[2])
            elif item[4] == '+':
                startSite = int(item[2])
                endSite = int(item[3])
            else:
                raise IOError("wrong input file %d" % str(item))
            bedtoolsformat_file.write("chr{}\t{}\t{}\t{}".format(item[1],
                                                             startSite,
                                                             endSite,
                                                             item[0]))
            bedtoolsformat_file.write("\n")
        bedtoolsformat_file.close()
        print("Saved bedtools Format file in path %s" % self.gene_bedtoolsFormat_filepath)
        return bedtoolsformat_file

    def promoter_into_bedtoolsFormat(self, save_filepath):
        self.promoter_bedtoolsFormat_filepath = save_filepath
        if self.gene_startEndSite_filepath == "":
            raise IOError("Try using find_allGene_startEndSite first.")
        gene_startEndSite = [item.split(',') for item in
                             open(self.gene_startEndSite_filepath, 'r').readlines()[1:]]
        bedtoolsformat_file = open(self.promoter_bedtoolsFormat_filepath, 'w')
        for item in gene_startEndSite:
            if item[4] == '-':
                promoter_start = item[3]
                promoter_end = int(item[3])+100
            elif item[4] == '+':
                promoter_start = int(item[2])-100
                promoter_end = item[2]
            else:
                raise IOError("wrong input file %d"%str(item))
            bedtoolsformat_file.write("chr{}\t{}\t{}\t{}".format(item[1],
                                                             str(promoter_start),
                                                             str(promoter_end),
                                                             str(item[0])))
            bedtoolsformat_file.write("\n")
        bedtoolsformat_file.close()
        print("Saved bedtools Format file in path %s" % self.promoter_bedtoolsFormat_filepath)
        return bedtoolsformat_file


if __name__ == '__main__':
    project_rootdir = '/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm'
    tss_dirpath = os.path.join(project_rootdir,
                               'data',
                               'features',
                               'TSSDistance')
    tss_raw_filepath = os.path.join(tss_dirpath,
                                    'Homo_sapiens.GRCh37.71.gtf')
    tss_file = TSS_file(tss_raw_filepath)
    gene_startEndSite_savepath = os.path.join(tss_dirpath,
                                              'gene_startEndSite.txt')
    gene_startEndSite = tss_file.find_allGene_startEndSite(gene_startEndSite_savepath)
    _ = tss_file.gene_into_bedtoolsFormat(gene_startEndSite_savepath[:-4]+"_bedtoolsFormat.txt")
    promoter_bedtoolsformat_filepath = os.path.join(tss_dirpath,
                                                    'promoter_startEndSite_bedtoolsFormat.txt')
    promoter_bedtoolsFormatfile = tss_file.promoter_into_bedtoolsFormat(promoter_bedtoolsformat_filepath)
