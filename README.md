# eqtm
## file description:
./data/eqtmZscores/
2017-12-09-eQTLsFDR-et0.0-flipped.txt contains the eQTMs with FDR==0.0
2017-12-09-eQTLsFDR-gt0.0-flipped.txt contains the eQTMs with 0.05>FDR>0.0
2017-12-09-eQTLsFDR-gt0.05-flipped.txt with absolutely insignificant eqtms with FDR>0.05
random20000-eQTLsFDR-gt0.05-flipped.txt with a random subset of 20000 insignificant eqtms


## Bonder file with 8 features

KNN results:
Sensitivity: 0.7944518535732488 0.03201209834688265
Specificity: 0.8205457207090685 0.013116070543817448
AUC: 0.8771763812306269 0.011426978696967242

Ranfor results:
Sensitivity: 0.806242914678354 0.06868651128741192
Specificity: 0.7967631653992886 0.1351834463827146
AUC: 0.8864217057946683 0.01872067736585379

PCA+ranfor results:
Sensitivity: 0.8015005380683347 0.022371978096972403
Specificity: 0.8218047997663446 0.014043291264409008
AUC: 0.8808003617275166 0.010454404648194448


## eQTM - et 0.0
PCA+ranfor results:
Sensitivity: 0.5813720756026459 0.01893350908461518
Specificity: 0.9617158881647309 0.005568636343231043
AUC: 0.9269014886443153 0.002455317288639995

## eQTM - gt 0.0
PCA+ranfor results:
Sensitivity: 0.5119424094973878 0.015669131423281392
Specificity: 0.8603116428673548 0.0024415519870774026
AUC: 0.8078494912694445 0.004256705611510144

# train on each other datasets and test
# removing features from blood cell types


## eQTM - gt 0.05
