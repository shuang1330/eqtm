import os


class roadmap(object):

    def __init__(self,
                 data_rootdir="/groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap",
                 feature_type="consolidatedImputedGappedPeak"):
        self._data_rootdir = data_rootdir
        self.feature_type = feature_type
        self.feature_dirpath = os.path.join(data_rootdir, feature_type)
        self.feature_list_filepath = ""
        self.histone_list_filepath = ""
        self.cell_list_filepath = ""
        self.feature_list = []
        self.histone_list = []
        self.cell_list = []

    def readFeatureList(self, feature_list_filepath):
        self.feature_list_filepath = feature_list_filepath
        # feature_list_path = os.path.join(dirs.temp_meta,
        #                                  dirs.feature_list_filename+'.txt')
        if os.path.exists(self.feature_list_filepath):
            with open(self.feature_list_filepath, 'r') as f:
                for line in f.readlines():
                    self.feature_list.append(line.strip())
                f.close()
            print('Read feature list from path: \n', self.feature_list_filepath)
            if 'feature_list' in self.feature_list:
                print('There should not be a item called feature list.')
        else:
            feature_f = open(self.feature_list_filepath, 'w')
            for name in os.listdir(self.feature_dirpath):
                if name.endswith('gz'):
                    feature_name = '.'.join(name.split('.')[:-5])
                    self.feature_list.append(feature_name)
                    feature_f.write('%s\n' % feature_name)
            feature_f.close()
            print('Saved feature list to path: \n%s' % self.feature_list_filepath)

        return self.feature_list

    def readHistoneList(self, histone_list_filepath):
        self.histone_list_filepath = histone_list_filepath
        if os.path.exists(self.histone_list_filepath):
            with open(self.histone_list_filepath, 'r') as f:
                for line in f.readlines():
                    self.histone_list.append(line.strip())
                f.close()
            print('Read histone list from path: \n', self.histone_list_filepath)
            print(self.histone_list)
        else:
            if len(self.feature_list) < 1:
                raise IOError("No feature list provided.\nTry call readFeatureList first.")
            for col in self.feature_list:
                histone = col.split('-')[1]
                if histone not in self.histone_list:
                    self.histone_list.append(histone)
            with open(self.histone_list_filepath, 'w+') as f:
                for item in self.histone_list:
                    f.write('%s\n' % item)
                f.close()
            print('Saved histone list to path:\n', self.histone_list_filepath)

        return self.histone_list

    def readCellList(self, cell_list_filepath):
        self.cell_list_filepath = cell_list_filepath
        if os.path.exists(self.cell_list_filepath):
            with open(self.cell_list_filepath, 'r') as f:
                for line in f.readlines():
                    self.cell_list.append(line.strip())
                f.close()
            print('Read cell types from path: \n', self.cell_list_filepath)
            print(self.cell_list)
        else:
            if len(self.feature_list) < 1:
                raise IOError("No feature list provided.\nTry call readFeatureList first.")
            for col in self.feature_list:
                cell = col.split('-')[0]
                if cell not in self.cell_list:
                    self.cell_list.append(cell)
            with open(self.cell_list_filepath, 'w+') as f:
                for item in self.cell_list:
                    f.write('%s\n' % item)
                f.close()
            print('Saved cell list to path:\n', self.cell_list_filepath)

        return self.cell_list

