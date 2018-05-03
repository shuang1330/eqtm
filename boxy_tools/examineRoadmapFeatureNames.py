import os

ROADMAP_DIR = '/groups/umcg-wijmenga/tmp03/projects/eQTMPrediction/features/Roadmap'
folders = ['consolidatedImputedGappedPeak','consolidatedImputedGappedPeak']

name_list_lengthDic = {}

for featureFolder in folders:
    featureFolderPath= os.path.join(ROADMAP_DIR,featureFolder)
    for featureFile in os.listdir(featureFolderPath):
        if featureFile.endswith('.gz'):
            name_list = featureFile.split('.')
            if len(name_list) in name_list_lengthDic:
                name_list_lengthDic[len(name_list)] += 1
                if name_list_lengthDic[len(name_list)] < 5:
                    print(name_list)
            else:
                print(name_list)
                name_list_lengthDic[len(name_list)] = 1

print(name_list_lengthDic)
