import pandas as pd
import os
# import matplotlib.pyplot as plt
import numpy as np

def addTssDistance2eQTMwithZscoreFile(eQTM_path,tss_path,save_path = None):
    '''
    add TssDistance to eqtmZscore file
    INPUT:
        eQTM_path, string, path to eqtmZscore file
        tss_path, string, path to TssSite file
        save_path, string, path to save the new eqtmZscore_withTSSDistance file
    OUTPUT:
        eQTMs, pandas dataframe
    '''
    # read tss file from tss_path
    colnames = ['chr','regionFunction','regionType','startSite',
                'endSite','score','strand','sthunknown','geneInfo']
    dtype = {'chr':object,'regionFunction':object,'regionType':object,
             'startSite':int,'endSite':int,'score':object,
             'strand':object,'sthunknown':object,'geneInfo':object}
    tss_raw = pd.read_csv(tss_path,sep='\t',header=None,names=colnames,dtype=dtype)

    # reading the eQTMs
    eQTMs = pd.read_csv(eQTM_path,sep='\t')

    # extract gene name from geneInfo for tss file
    def findGeneName(item):
        item = [thing for thing in list(filter(None,item.strip().split(";")))][0]
        name = item.replace('"','').replace(';','').strip().split(' ')[1]
        return name
    tss_raw['geneName'] = tss_raw['geneInfo'].apply(findGeneName)

    # find the tss sites for each gene in the tss file
    groupbyTss = tss_raw.groupby('geneName').agg({
        'chr':lambda x: x.unique(),
        'startSite':np.min,
        'endSite':np.max,
        'strand':lambda x: x.unique()
    })
    def findTssSite(series):
        if series[3] == '-':
            return series[2]
        else:
            return series[1]
    groupbyTss['TssSite'] = groupbyTss.apply(findTssSite,axis=1)

    # add tss sites and tss distance to the eqtm file
    def mapSite(row):
        return groupbyTss.loc[row]['TssSite']
    def calculateDis(row):
        return abs(row[0]-row[1])
    def findChr(row):
        return groupbyTss.loc[row]['chr']
    def checkChr(row):
        if str(row[0])==str(row[1]):
            return True
        else:
            return False
    eQTMs['TssSite'] = eQTMs['ProbeName'].apply(mapSite)
    eQTMs['chr'] = eQTMs['ProbeName'].apply(findChr)
    eQTMs['TssDistance'] = eQTMs[['SNPChrPos','TssSite']].apply(calculateDis,axis=1)
    eQTMs['checkChr'] = eQTMs[['chr','SNPChr']].apply(checkChr,axis=1)
    # check whether they are from the same chromosome
    assert len(eQTMs['checkChr'].unique()) == 1

    if save_path:
        # save the eQTM file
        eQTMs.to_csv(save_path,index=False)
        print('Saved eQTM file to: ',save_path)

    return eQTMs

def changeCpGColumnName_addZscoreTss_forRandom(feature_file_path,
                                               feature_sep,
                                               zscore_path,
                                               zscore_sep,
                                               save_path = None):
    '''
    add Zscore and TSSDistance to overlapRatio file
    INPUT:
        feature_file_path, string, path to overlapRatio file
        feature_sep, string of seperator, '\t' or ','
        zscore_path, string, path to eqtmZscore_withTSSDistance file
        zscore_sep, string of seperator
        save_path, string, path to save the overlapRatioZscoreTss file,
            default None
    OUTPUT:
        feature_file, pandas dataframe, the overlapRatioZscoreTss
    '''
    feature_file = pd.read_csv(feature_file_path,sep=feature_sep,index_col=0)
    feature_file = feature_file.rename(index=str,columns={'SNPName':'cpgName'})
    def read_zscore(zscore_path):
        zscore_file = pd.read_csv(zscore_path,sep=zscore_sep,index_col=0)

        zscore_file['cpgName'] = zscore_file[['SNPName']]
        zscore_dic = zscore_file.set_index('cpgName')['OverallZScore'].to_dict()
        tss_dic = zscore_file.set_index('cpgName')['TssDistance'].to_dict()
        return zscore_dic,tss_dic
    zscore_dic,tss_dic = read_zscore(zscore_path)
    def map_zscore(row):
        if row in zscore_dic:
            return zscore_dic[row]
        return None
    def map_tss(row):
        if row in tss_dic:
            return tss_dic[row]
        return None
    feature_file['zscore'] = feature_file['cpgName'].apply(map_zscore)
    feature_file['TssDistance'] = feature_file['cpgName'].apply(map_tss)

    if save_path:
        feature_file.to_csv(save_path)
        print('overlapRatioZscoreTss saved to path:',save_path)

    return feature_file

def addMeanVar(meanVar_filepath, cpg_filepath, final_res_filepath=None):

    def read_meanVar_asDic(meanVar_filepath):
        meanVar = pd.read_csv(meanVar_filepath,sep='\t')
        meanVar_dic = meanVar[['ID','Mean','Var']].set_index('ID').T.to_dict('list')
        return meanVar_dic

    dic = read_meanVar_asDic(meanVar_filepath)
    def findMean(row):
        if row in dic:
            return dic[row][0]
        else:
            print('Cpg not found. Return None.')
            return None

    def findVar(row):
        if row in dic:
            return dic[row][1]
        else:
            print('Cpg not found. Return None.')
            return None

    print("Processing the datafile %s."%cpg_filename)
    cpg = pd.read_csv(cpg_filepath,sep=',')
    print("Adding the mean value.")

    cpg['methyMean'] = cpg['cpgName'].apply(findMean)
    print("Adding the variance.")
    cpg['methyVar'] = cpg['cpgName'].apply(findVar)
    print("Done adding the mean and variance. Please check below:\n",
          cpg[['cpgName','zscore','methyMean','methyVar']].head())

    if final_res_filepath:
        cpg.to_csv(final_res_filepath,index=False)
        print('Saved overlapRatioZscoreTssMeanVar file to path %s.'%(final_res_filepath))

    return cpg

if __name__=='__main__':

    PROJECT_DIR='/groups/umcg-gcc/tmp03/umcg-sli/boxy_eqtm'
    eqtm_name = 'random20k_gt0.5'

    # gene position
    tss_raw_path = os.path.join(PROJECT_DIR,'data','features','TSSDistance',
                                'Homo_sapiens.GRCh37.71.gtf')

    # eqtmZscore file path
    eqtmZscore_folder = os.path.join(PROJECT_DIR,'data','eqtmZscores')
    eqtmZscore_path = os.path.join(eqtmZscore_folder,
                                   'fdr_gt0.05',eqtm_name+'.txt')
    eqtmZscoreTss_savepath = os.path.join(eqtmZscore_folder,
                                   'fdr_gt0.05',eqtm_name+'withZscoreTSS.txt')
    _ = addTssDistance2eQTMwithZscoreFile(eqtmZscore_path,
                                          tss_raw_path,
                                          eqtmZscoreTss_savepath)

    # add to overlapRatio the zscore and Tss
    overlapRatio_path = os.path.join(PROJECT_DIR,'data',
                                     'dataReadyForModeling',
                                     'overlapRatio',
                                     eqtm_name+'_overlapRatio.txt')
    overlapRatioZscoreTss_path = os.path.join(PROJECT_DIR,'data',
                                     'dataReadyForModeling',
                                     'overlapRatioTss',
                                     eqtm_name+'_overlapRatioTss.txt')
    _ = changeCpGColumnName_addZscoreTss_forRandom(overlapRatio_path,
                                     '\t',eqtmZscoreTss_savepath,',',
                                     save_path=overlapRatioZscoreTss_path)

    # add mean Variance
    overlapRatioTss_folder = os.path.join(PROJECT_DIR,'data',
                              'dataReadyForModeling','overlapRatioTss')
    overlapRatioTssMeanVar_folder = os.path.join(PROJECT_NEWDIR,'data',
                              'dataReadyForModeling','overlapRatioTssMeanVar')
    meanVar_filepath = os.path.join(PROJECT_DIR,'data',
                                    'features','meanVar',
                                    'methylation-MeanAndVarianceRows.txt')
    overlapRatioTss_filename = eqtm_name+'_withZscoreTss'
    overlapRatioTss_filepath = os.path.join(overlapRatioTss_folder,
                                      overlapRatioTss_filename+'.csv')
    overlapRatioTssMeanVar_filepath = os.path.join(overlapRatioTss_folder,
                                      overlapRatioTss_filename+'MeanVar.csv')
    _ =addMeanVar(meanVar_filepath,overlapRatioTssMeanVar_filepath,
                  overlapRatioTssMeanVar_filepath)
