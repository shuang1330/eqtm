import os
from sklearn.ensemble import RandomForestClassifier
from lib.read.read_data import read_eqtm_names

if __name__ == "__main__":
    project_rootdir = '/home/shuang/projects/development_eqtm'
    data_dirpath = os.path.join(project_rootdir, 'data')
    cpgSites_dirpath = os.path.join(data_dirpath, 'cpgSites', 'seperate_cpgFiles')
    geneSites_dirpath = os.path.join(data_dirpath, 'eqtmZscores', 'testCNN', 'geneSites')
    eqtm_path = os.path.join(data_dirpath, 'eqtmZscores',
                             'withExpressionTSSMethyCpgOverlapGene',
                             '2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_withGeneOverlap.txt')

    # model
    model = RandomForestClassifier()