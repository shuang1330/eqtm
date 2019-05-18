
features = ["TssDistance","expressionMean","expressionVar","methyMean",
"methyVar","H4K5ac","H2A.Z","H2BK120ac","H3K79me2","H3K27ac",
"H2BK20ac","H3K14ac","H3K9ac","H4K20me1","H3K4ac","H2AK5ac",
"H3K9me3","H3K36me3","H3K4me1","H3K18ac","H3K23ac","H2BK5ac",
"H3K4me3","H2BK12ac","H3K23me2","H4K12ac","DNase","H2BK15ac",
"H3K9me1","H3K4me2","H3K27me3","H3K79me1","H2AK9ac","H3T11ph",
"H4K8ac","H4K91ac","H3K56ac"]

with open("/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/extras/model1_2_features.txt", "w") as f:
    f.write("\n".join(features))
