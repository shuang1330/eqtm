import os
import pandas as pd
import numpy as np

eqtmfolder = '/groups/umcg-gcc/tmp03/umcg-sli/eqtm_project/data/eqtmZscores'
for eqtm_name in os.listdir(eqtmfolder):
    eqtm_path = os.path.join(eqtmfolder,eqtm_name)
    eqtmfile = pd.read_csv(eqtm_path,sep='\t')
    print(eqtm_name)
    print('With shape: ',eqtmfile.shape)
    # print(eqtmfile.columns)
    print('Unique eqtm names:',len(eqtmfile['SNPName'].unique()))
    groupbyEqtm = eqtmfile[['OverallZScore','SNPName']].groupby(['SNPName']).agg(['count'])
    eqtmDup = groupbyEqtm['OverallZScore'][groupbyEqtm['OverallZScore']['count']>1]
    # print(eqtmDup.index.values)
    print(eqtmfile[eqtmfile['SNPName']==eqtmDup.index.values[0]][['SNPName','SNPChrPos','OverallZScore','ProbeName','FDR']])
