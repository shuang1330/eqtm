import os

class directoryLocal:
    def __init__(self,PROJECT_DIR,eqtm_filename):
        self.project_rootdir = PROJECT_DIR
        self.eqtm_filename = eqtm_filename
        self.feature_list_filename = self.eqtm_filename+'_featureList'
        self.histone_list_filename = self.eqtm_filename+'_histoneList'
        self.cell_list_filename = self.eqtm_filename+'_cellList'
        self.feature_folder = os.path.join(self.project_rootdir,
                                           'data','RoadmapFeatureSamples')
        self.cpg_folder = os.path.join(self.project_rootdir,
                                       'data','eqtmZscores')
        # for saving intermediate results
        self.temp_output = os.path.join(self.project_rootdir,'data','temp')
        self.temp_meta = os.path.join(self.temp_output,'meta')
        self.temp_cpgFolder = os.path.join(self.temp_output,'cpg')
        self.temp_featureFolder = os.path.join(self.temp_output,'features')
        self.temp_overlap = os.path.join(self.temp_output,'overlap',
                                         self.eqtm_filename)
        # final results
        self.output_ratio = os.path.join(self.project_rootdir,
                                         'data','output',
                                         self.eqtm_filename)

    def initializeFolders(self):
        folder_list = [self.temp_output,self.temp_meta,
                       self.temp_cpgFolder,self.temp_overlap,
                       self.temp_featureFolder,self.output_ratio]
        for folder in folder_list:
            if not os.path.exists(folder):
                os.makedirs(folder)
        print('Initiated all folders locally.')
