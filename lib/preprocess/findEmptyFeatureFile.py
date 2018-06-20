import os
import gzip


def gz_size(fname):
    # https://stackoverflow.com/a/37875919
    with gzip.open(fname, 'rb') as f:
        return f.seek(0, whence=2)

def findEmptyFeature(feature_dir,empty_featureList_savepath):
    empty_featureList = open(empty_featureList_savepath,'w')
    i = 0
    for filename in os.listdir(feature_dir):
        if filename.endswith('.gz'):
            filepath = os.path.join(feature_dir,filename)
            if gz_size(filepath)==0:
                i += 1
                empty_featureList.write('{}\n'.format(filename))
                print('Empty feature:',filename)
    print('Finished checking all features. {} out of {} are empty.'.format(
    i,len(os.listdir(feature_dir))
    ))

if __name__=='__main__':
    feature_dir = '/groups/umcg-bios/tmp03/projects/2018-methylation/input/features/Roadmap/consolidatedImputedGappedPeak'
    empty_featureList_savepath = '/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/features/emptyFeatureFiles/emptyFeatureList.txt'
    findEmptyFeature(feature_dir,empty_featureList_savepath)
