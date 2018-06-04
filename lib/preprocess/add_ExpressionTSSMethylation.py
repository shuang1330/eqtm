import os
import pandas as pd
import numpy as np

def read_expression_data(expression_filepath,sep=','):
    '''
    read gene expression file:
    expression-MeanAndVarianceRows.txt
    '''
    expression = pd.read_csv(expression_filepath,sep=sep,index_col=0)
    expression_dict = expression_data[['ID','Mean','Var']]
    .set_index('ID')
    .T.to_dict('list')
    return expression_dict

def add_Expression_basedOnProbeName_toEQTMFile(eqtm_filepath,
                                                 expression_filepath,
                                                 eqtm_savepath):
    '''
    read expression data and add them to the eqtm file based on
    the methylation sites and gene name (ProbeName)
    INPUT:
        eqtm_filepath: str
        expression_filepath: str
        eqtm_savepath: str
    OUTPUT:
        eqtm file with gene expression data included: pandas dataFrame
        eqtm file saved in eqtm_savepath
    '''
    expression_dict = read_expression_data(expression_filepath)
    eqtm = pd.read_csv(eqtm_filepath,sep='\t',index_col=0)
    eqtm['expressionMean'] = [expression_dict[key][0] for key in eqtm['ProbeName'].values]
    eqtm['expressionVar'] = [expression_dict[key][1] for key in eqtm['ProbeName'].values]
    eqtm.to_csv(eqtm_savepath)
    print("Saved to path:",eqtm_savepath)
    return eqtm

def read_tss_data(tss_filepath):
    '''
    read tss file from tss_filepath
    '''
    colnames = ['chr','regionFunction','regionType','startSite',
                'endSite','score','strand','sthunknown','geneInfo']
    dtype = {'chr':object,'regionFunction':object,'regionType':object,
             'startSite':int,'endSite':int,'score':object,
             'strand':object,'sthunknown':object,'geneInfo':object}
    tss_raw = pd.read_csv(tss_filepath,sep='\t',header=None,
                          names=colnames,dtype=dtype)
    return tss_raw

def add_TSS_basedOnProbeName_toEQTMFile(eqtm_filepath,tss_filepath,
                                      eqtm_savepath=None):
    '''
    add TssDistance to eqtm file
    INPUT:
        eQTM_path, string, path to eqtmZscore file
        tss_filepath, string, path to TssSite file
        save_path, string, path to save the new eqtmZscore_withTSSDistance file
    OUTPUT:
        eQTMs, pandas dataframe
    '''
    # read the tss file
    tss_raw = read_tss_data(tss_filepath)
    # reading the eQTMs
    eQTMs = pd.read_csv(eqtm_filepath,sep=',')

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

    if eqtm_savepath:
        # save the eQTM file
        eQTMs.to_csv(eqtm_savepath,index=False)
        print('Saved eQTM file to: ',eqtm_savepath)

    return eQTMs

def read_methylation_data(methy_filepath):
    '''
    read methylation level files
    '''
    methy = pd.read_csv(methy_filepath,sep='\t')
    methy_dict = methy[['ID','Mean','Var']].set_index('ID').T.to_dict('list')
    return methy_dict

def add_Methy_basedOnSNPName_toEQTMFile(eqtm_filepath,
                                        methy_filepath,
                                        eqtm_savepath):
    '''
    add methylation level to each site
    '''
    methy_dict = read_methylation_data(methy_filepath)
    inputfile = pd.read_csv(eqtm_filepath,index_col=0)
    inputfile['SNPName_ProbeName'] = ['{}_{}'.format(row[0],row[1])
                                      for row in inputfile[['SNPName','ProbeName']].values]
    inputfile['methyMean'] = [methy_dict[row][0] for row in inputfile['SNPName'].values]
    inputfile['methyVar'] = [methy_dict[row][1] for row in inputfile['SNPName'].values]
    inputfile.to_csv(eqtm_savepath)
    print('Eqtm file with gene expressino data, TSS distance and methylation data saved to:',
          eqtm_savepath)
    return inputfile


if __name__=='__main__':

    DATA_FOLDER = '/home/shuang/projects/boxy_eqtm/data'
    EQTM_DATADIR = os.path.join(DATA_FOLDER,'eqtmZscores')
    INPUT_FOLDER = os.path.join(EQTM_DATADIR,'ORIGIN')
    OUTPUT_EXPRESS = os.path.join(EQTM_DATADIR,'withExpressionData')
    OUTPUT_EXPRESS_TSS = os.path.join(EQTM_DATADIR,'withExpressionTSS')
    OUTPUT_EXPRESS_TSS_METHY = os.path.join(EQTM_DATADIR,
                                            'withExpressionTSSMethy')

    # inputfiles
    eqtm_names = {'2017-12-09-eQTLsFDR-et0.0-flipped.txt':'2017-12-09-eQTLsFDR-et0.0',
                  '2017-12-09-eQTLsFDR-gt0.0-flipped.txt':'2017-12-09-eQTLsFDR-gt0.0',
                  'random20k_gt0.5_withTss.txt':'random20k_gt0.5'}

    # expression data
    expression_filepath = os.path.join(DATA_FOLDER,'features','meanVar',
                                       'expression-MeanAndVarianceRows.txt')

    # read expression data and add expression to eqtm files
    for filename in eqtm_names.keys():
        print("Processing file:",filename)
        eqtm_name = eqtm_names[filename]
        eqtm_filepath = os.path.join(INPUT_FOLDER,filename)
        output_express_filepath = os.path.join(OUTPUT_EXPRESS,
                                  eqtm_name+'_withExpression.txt')
        output_expressTss_filepath = os.path.join(OUTPUT_EXPRESS_TSS,
                                  eqtm_name+'_withExpressionTSS.txt')
        output_expressTssMethy_filepath = os.path.join(OUTPUT_EXPRESS_TSS_METHY,
                                  eqtm_name+'_withExpressionTSSMethy.txt')
        _ = add_Expression_basedOnProbeName_toEQTMFile(eqtm_filepath,
                                                    expression_filepath,
                                                    output_express_filepath)
        _ = add_TSS_basedOnProbeName_toEQTMFile(output_express_filepath,
                                                tss_filepath,
                                                output_expressTss_filepath)
        _ = add_Methy_basedOnSNPName_toEQTMFile(output_expressTss_filepath,
                                                tss_filepath,
                                                output_expressTss_filepath)
