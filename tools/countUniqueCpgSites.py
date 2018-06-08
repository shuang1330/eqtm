import os
import numpy as np
import pandas as pd

CPG_DATADIR = '/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/ORIGIN'
all_cpg_filepath = '/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/cpgSites/all_cpgSites.txt'

dic = {}
pos_dic = {}
cpgname_set = set()
probename_set = set()
all_cpg = open(all_cpg_filepath,'w')
all_cpg.write('SNPName\tChr\tPos\n')

for filename in os.listdir(CPG_DATADIR):
    if filename.endswith('.txt'):
        print('Processing file:',filename)
        with open(os.path.join(CPG_DATADIR,filename),'r') as eqtm:
            for line in eqtm.readlines():
                if line[0].isdigit():
                    if len(line.strip().split('\t'))<2 and len(line)>2:
                        infos = line.strip().split(',')
                    else:
                        infos = line.strip().split('\t')
                    cpgname = infos[1]
                    probename = infos[4]
                    cpg_chr = infos[2]
                    cpg_pos = infos[3]
                    if cpgname not in dic:
                        dic[cpgname] = [probename]
                        pos_dic[cpgname] = [cpg_chr,cpg_pos]
                    else:
                        dic[cpgname].append(probename)
for key in pos_dic:
    # print(pos_dic[key])
    all_cpg.write(key)
    all_cpg.write('\t')
    all_cpg.write('\t'.join(pos_dic[key]))
    all_cpg.write('\n')

all_cpg.close()
print("Finished, let's take a look at the results!")
print(len(dic.keys()))
gene_num = []
for key in dic:
    gene_num.append(len(dic[key]))
hist,bins = np.histogram(gene_num)
print('Hist',hist)
print('Bins',bins)
