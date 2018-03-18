import pandas as pd
import numpy as np
import os
PROJECT_DIR='/home/shuang/projects/eqtm'

if __name__=='__main__':
    bedtools_res_file = os.path.join(PROJECT_DIR,'bedtools_intersect_output.txt')
    bedtools_res = pd.read_csv(bedtools_res_file,sep='\t',header=None)

    bedtools_res[3] = bedtools_res[3].apply(lambda row:row.split('/')[-1].split('.')[0])
    # print(bedtools_res[3].head(5))

    def create_list(row):
        return list(row)
    grouped_df = pd.DataFrame(bedtools_res.groupby(7)[3].apply(create_list).reset_index())
    grouped_df = grouped_df.rename(index=str,columns={7:'cg',3:'feature'})
    overlap_matrix = grouped_df[['cg']].join(grouped_df['feature'].str.join('|').str.get_dummies())

    # add the value to the overlap_matrix
    cg_value_file = os.path.join(PROJECT_DIR, 'eQTMsUniqueCGs-FDR0.05.txt')
    cg_value_dataframe = pd.read_csv(cg_value_file,sep='\t')
    cg_value_dict = cg_value_dataframe.set_index(['SNPName'])['OverallZScore'].to_dict()

    def find_value(row):
        return cg_value_dict[row]
    overlap_matrix['cg_value'] = grouped_df['cg'].apply(find_value)
    print(overlap_matrix.head())
    # save the overlap_matrix
    overlap_matrix.to_csv(PROJECT_DIR+'/output/overlap_matrix.csv',index=False)
