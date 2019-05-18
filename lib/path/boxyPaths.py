import os

# boxy

class directoryBoxy:
    def __init__(self,
                 PROJECT_DIR,
                 FEATURE_ROOTDIR,
                 ):
        self.project_rootdir = PROJECT_DIR

        # feature files
        self._feature_rootdir = FEATURE_ROOTDIR
        self.feature_folder = os.path.join(self._feature_rootdir,
                                           'features',
                                           'Roadmap',
                                           'consolidatedImputedGappedPeak')
        self.feature_type = 'imputedGappedPeak'
        self.feature_filetype = '.imputed.gappedPeak.bed.gPk.gz'

        # self.eqtm_folder = os.path.join(self.project_rootdir,
        #                                 'data','eqtmZscores','ORIGIN')
        # self.eqtmname = eqtmname
        # self.feature_list_filename = self.feature_type+'_featureList'
        # self.histone_list_filename = self.feature_type+'_histoneList'
        # self.cell_list_filename = self.feature_type+'_cellList'

        # self.cpgDir = os.path.join(self.project_rootdir,'data',
                                   # 'eqtmZscores','ORIGIN')
        # self.cpgname = 'GenomicCoordinates'
        # self.geneDir = os.path.join(self.project_rootdir,'data','geneSites')
        # self.genename = 'all_geneSites'
        # self.cpg_bedtoolFormat_filepath = None

        # for saving intermediate results
        self.temp_output_dirpath = os.path.join(self.project_rootdir, 'data', 'temp')
        self.temp_meta = os.path.join(self.temp_output_dirpath, 'meta')
        self.temp_featureFolder = os.path.join(self.temp_output_dirpath, 'features',
                                               self.feature_type)

    def initializeFolders(self):

        folder_list = [self.cpgDir, self.geneDir,
                       self.temp_output_dirpath, self.temp_meta,
                       self.temp_featureFolder]
        for folder in folder_list:
            if not os.path.exists(folder):
                os.makedirs(folder)
        print('Initiated all folders on boxy.')
