import os

# boxy

class directoryBoxy:
    def __init__(self,PROJECT_DIR,DATA_ROOTDIR,cpg_filename):
        self.project_rootdir = PROJECT_DIR
        self._data_rootdir = DATA_ROOTDIR
        self.feature_folder = os.path.join(self._data_rootdir,
                                           'features',
                                           'Roadmap',
                                           'consolidatedImputedGappedPeak')
        self.cpg_folder = os.path.join(self.project_rootdir,'data','eqtmZscores')
        self.cpg_filename = cpg_filename
        self.feature_list_filename = self.cpg_filename+'_featureList'
        self.histone_list_filename = self.cpg_filename+'_histoneList'
        self.cell_list_filename = self.cpg_filename+'_cellList'

        # for saving intermediate results
        self.temp_output = os.path.join(self.project_rootdir,'data','temp')
        self.temp_meta = os.path.join(self.temp_output,'meta')
        self.temp_cpgFolder = os.path.join(self.temp_output,'cpg')
        self.temp_featureFolder = os.path.join(self.temp_output,'features',
                                               self.cpg_filename)
        self.temp_overlap = os.path.join(self.temp_output,'overlap',
                                         self.cpg_filename)
        # final results
        self.output_ratio = os.path.join(self.project_rootdir,
                                         'data','output',
                                         self.cpg_filename)

    def initializeFolders(self):
        folder_list = [self.temp_output,self.temp_meta,
                       self.temp_cpgFolder,self.temp_overlap,
                       self.temp_featureFolder,self.output_ratio]
        for folder in folder_list:
            if not os.path.exists(folder):
                os.makedirs(folder)
        print('Initiated all folders on boxy.')
