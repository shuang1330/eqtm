# eqtm
## file description:
./data/eqtmZscores/
2017-12-09-eQTLsFDR-et0.0-flipped.txt contains the eQTMs with FDR==0.0
2017-12-09-eQTLsFDR-gt0.0-flipped.txt contains the eQTMs with 0.05>FDR>0.0
2017-12-09-eQTLsFDR-gt0.05-flipped.txt with absolutely insignificant eqtms with FDR>0.05
random20000-eQTLsFDR-gt0.05-flipped.txt with a random subset of 20000 insignificant eqtms


## Bonder file with 8 features

KNN results:
Sensitivity: 0.79 0.032
Specificity: 0.82 0.013
AUC: 0.87 0.011

Ranfor results:
Sensitivity: 0.80 0.068
Specificity: 0.79 0.135
AUC: 0.88 0.0187

PCA+ranfor results:
Sensitivity: 0.80 0.022
Specificity: 0.82 0.014
AUC: 0.88 0.010


## eQTM - et 0.0
### withought TSS distance
#### test on itself
PCA+ranfor results:
Sensitivity: 0.58 0.018
Specificity: 0.96 0.005
AUC: 0.92 0.002
#### test on gt
PCA+ranfor results:
Sensitivity: 0.23 0.009
Specificity: 0.92 0.005
AUC: 0.71 0.004

### with TSS Distance
#### test on itself
Sensitivity:
Specificity:
AUC: 0.9255365184536088 0.0025035726727016317
#### test on gt
Sensitivity: 0.231861246847401 0.01792657087569302
Specificity: 0.9192671144331843 0.011409685020798133
AUC: 0.71 0.004
### with mean and var
#### test on itself
Sensitivity:
Specificity:
AUC: 0.9220396512916833 0.007815859259167493
#### test on gt
Sensitivity: 0.2528913084878358 0.012664690196804896
Specificity: 0.9072004016868995 0.007673719631246879
AUC: 
### with tss+mean var
#### test on itself
#### test on gt

## eQTM - gt 0.0
### without TSS distance
PCA+ranfor results:
Sensitivity: 0.5119424094973878 0.015669131423281392
Specificity: 0.8603116428673548 0.0024415519870774026
AUC: 0.8078494912694445 0.004256705611510144
### with TSS Distance
### with mean and var
### with tss+mean var

# train on each other datasets and test


# removing features from blood cell types


## eQTM - gt 0.05
