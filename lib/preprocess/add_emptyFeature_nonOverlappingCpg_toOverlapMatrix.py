import pandas as pd


def read_emptyFeatureList(empty_featureList_filepath):
    with open(empty_featureList_filepath, 'r') as f:
        return [row.strip() for row in f.readlines()]


def read_overlapMatrix(overlapMatrix_filepath, sep=','):
    return pd.read_csv(overlapMatrix_filepath, sep=sep, index_col=0)


def read_completeCpgs(original_eqtm_filepath,
                      sep='\t',
                      cpgname='SNPName',
                      header='infer',
                      names=None):
    return set(pd.read_csv(original_eqtm_filepath,
                           sep=sep,
                           header=header,
                           names=names)[cpgname].unique())


def add_emptyFeature_nonOverlappingCpg(empty_featureList_filepath,
                                       overlapMatrix_filepath,
                                       original_eqtm_filepath,
                                       complete_overlapMatrix_savepath,
                                       cpgcolname='SNPName',
                                       cpgheader='infer',
                                       cpgnames=None,
                                       geneSiteSep='\t'):
    empty_featureList = read_emptyFeatureList(empty_featureList_filepath)  # list
    print('Empty features: ', len(empty_featureList))
    complete_cpgs = read_completeCpgs(original_eqtm_filepath,
                                      cpgname=cpgcolname,
                                      header=cpgheader,
                                      names=cpgnames,
                                      sep=geneSiteSep)  # set
    overlapMatrix = read_overlapMatrix(overlapMatrix_filepath) #df

    for empty_feature in empty_featureList:
        featurename = '.'.join(empty_feature.split('.')[:-5])
        overlapMatrix[featurename] = 0

    print(len(complete_cpgs), len(overlapMatrix.index))
    nonOverlappingCpgs = complete_cpgs - set(overlapMatrix.index)
    print('In total, there are {} non-overlapping cpg sites.'.format(len(nonOverlappingCpgs)))
    for cpg in nonOverlappingCpgs:
        new_CpgRow = pd.DataFrame(data=0,
                                  index=[cpg],
                                  columns=overlapMatrix.columns)
        overlapMatrix.append(new_CpgRow)
    overlapMatrix.to_csv(complete_overlapMatrix_savepath)
    print(
    'Complete overlapMatrix with empty feature and non-overlapping cpg sites saved in',
    complete_overlapMatrix_savepath
    )
    return overlapMatrix


if __name__ == '__main__':

    empty_featureList_filepath = '/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/features/emptyFeatureFiles/emptyFeatureList.txt'
    # overlapMatrix_dirpath = '/home/shuang/projects/development_eqtm/data/eqtmZscores/overlapMatrix'
    # original_eqtm_dirpath = '/home/shuang/projects/development_eqtm/data/eqtmZscores/ORIGIN'
    # complete_overlapMatrix_dirpath = '/home/shuang/projects/development_eqtm/data/eqtmZscores/complete_overlapMatrix'
    # overlapMatrix_dic2 = {
    # 'GenomicCoordinates_overlapMatrix.txt':
    # 'GenomicCoordinates.txt'
    # }
    # overlapMatrix_dic = {
    # '2017-12-09-eQTLsFDR-gt0.0-flipped_overlapMatrix.txt':
    # '2017-12-09-eQTLsFDR-gt0.0-flipped.txt',
    # '2017-12-09-eQTLsFDR-et0.0-flipped_overlapMatrix.txt':
    # '2017-12-09-eQTLsFDR-et0.0-flipped.txt',
    # 'random20k_gt0.5_overlapMatrix.txt':
    # 'random20k_gt0.5.txt'
    # }

    overlapMatrix_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/allCpgs/eqtm_FDR_larger_than_0.05_bedtoolsFormat_overlapMatrix.txt"
    original_eqtm_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/allCpgs/eqtm_FDR_larger_than_0.05_bedtoolsFormat.txt"
    complete_overlapMatrix_savepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/allCpgs/eqtm_FDR_larger_than_0.05_bedtoolsFormat_overlapMatrix_complete.txt"
    add_emptyFeature_nonOverlappingCpg(empty_featureList_filepath,
                                       overlapMatrix_filepath,
                                       original_eqtm_filepath,
                                       complete_overlapMatrix_savepath,
                                       cpgcolname='SNP',
                                       cpgheader=None,
                                       cpgnames=["chr", "startSite",
                                                 "endSite", "SNP"]
                                       )
    #
    # for overlapMatrix_filename in overlapMatrix_dic2.keys():
    #     overlapMatrix_name = overlapMatrix_filename[:-4]
    #     overlapMatrix_filepath = os.path.join(overlapMatrix_dirpath,
    #                                          overlapMatrix_filename)
    #     original_eqtm_filename = overlapMatrix_dic2[overlapMatrix_filename]
    #     original_eqtm_filepath = os.path.join(original_eqtm_dirpath,
    #                                       original_eqtm_filename)
    #     complete_overlapMatrix_savepath = os.path.join(complete_overlapMatrix_dirpath,
    #                                       overlapMatrix_name+'complete.txt')
    #     add_emptyFeature_nonOverlappingCpg(empty_featureList_filepath,
    #                                        overlapMatrix_filepath,
    #                                        original_eqtm_filepath,
    #                                        complete_overlapMatrix_savepath,
    #                                        cpgcolname='SNP')
