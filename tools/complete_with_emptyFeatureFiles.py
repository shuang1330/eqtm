import __init__path
import os
from lib.preprocess.add_emptyFeature_nonOverlappingCpg_toOverlapMatrix import add_emptyFeature_nonOverlappingCpg

if __name__=='__main__':
    # locally
    # project_rootdir = '/home/shuang/projects/development_eqtm'
    # calculon
    project_rootdir = '/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm'
    empty_featureList_filepath = os.path.join(project_rootdir,
                                              'data',
                                              'features',
                                              'emptyFeatureFiles',
                                              'emptyFeatureList.txt')
    # # cpg sites
    # cpg_overlapMatrix_filepath = os.path.join(project_rootdir,
    #                                           'data',
    #                                           'cpgSites',
    #                                           'all_cpgSites_overlapMatrix.txt')
    # cpg_original_eqtm_filepath = os.path.join(project_rootdir,
    #                                          'data',
    #                                          'cpgSites',
    #                                          'all_cpgSites.txt')
    # cpg_complete_overlapMatrix_savepath = os.path.join(project_rootdir,
    #                                                   'data',
    #                                                   'cpgSites',
    #                                                   'all_cpgSites_overlapMatrix_complete.txt')
    #
    # add_emptyFeature_nonOverlappingCpg(empty_featureList_filepath,
    #                                    cpg_overlapMatrix_filepath,
    #                                    cpg_original_eqtm_filepath,
    #                                    cpg_complete_overlapMatrix_savepath,
    #                                    cpgcolname='SNPName')

    # gene sites
    gene_overlapMatrix_filepath = os.path.join(project_rootdir,
                                              'data', 'features',
                                              'geneOverlap',
                                              'promoter_overlapMatrix.txt')
    gene_original_eqtm_filepath = os.path.join(project_rootdir,
                                              'data',
                                              'geneSites',
                                              'all_geneSites.txt')
    gene_complete_overlapMatrix_savepath = os.path.join(project_rootdir,
                                                       'data', 'features',
                                                       'geneOverlap',
                                                       'promoter_overlapMatrix_complete.txt')

    add_emptyFeature_nonOverlappingCpg(empty_featureList_filepath,
                                       gene_overlapMatrix_filepath,
                                       gene_original_eqtm_filepath,
                                       gene_complete_overlapMatrix_savepath,
                                       cpgcolname='geneName',
                                       geneSiteSep=',')