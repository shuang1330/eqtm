import pandas as pd
import argparse


def read_eqtm(eqtm_filepath):
    eqtm = pd.read_csv(eqtm_filepath, sep="\t")

    def flip(item):
        return -1 if item == "T" else 1
    eqtm["sign"] = [flip(item) for item in eqtm['AlleleAssessed']]
    eqtm["flippedZscore"] = eqtm["OverallZScore"] * eqtm["sign"]
    print(eqtm.shape, eqtm.head())
    return eqtm


def add_promoter_overlapRatio(eqtm, promoter_overlapRatio_filepath):
    overlap = pd.read_csv(promoter_overlapRatio_filepath, sep="\t", index_col=0)
    overlap_dict = overlap.T.to_dict()
    print(eqtm.columns)
    for col in overlap.columns:
        col_name = col+'_promoter'
        eqtm[col_name] = [overlap_dict[row][col] for row in eqtm['ProbeName'].values]
    print(eqtm.shape, eqtm.head())
    return eqtm


def add_cpg_overlapRatio(eqtm, cpg_overlapRatio_filepath):
    overlap = pd.read_csv(cpg_overlapRatio_filepath, sep="\t", index_col=0)
    overlap_dict = overlap.T.to_dict()
    for col in overlap.columns:
        def find_overlap_vector(item):
            if item in overlap_dict:
                return overlap_dict[item][col]
            else:
                return 0
        col_name = col
        eqtm[col_name] = [find_overlap_vector(item) for item in eqtm['SNPName'].values]
    print(eqtm.shape, eqtm.head())
    return eqtm


def add_tss_distance(eQTMs, gene_site_filepath):
    groupbyTss = pd.read_csv(gene_site_filepath, sep=",", index_col=0)

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
    eQTMs['TssDistance'] = eQTMs[['SNPChrPos', 'TssSite']].apply(calculateDis, axis=1)
    eQTMs['checkChr'] = eQTMs[['chr', 'SNPChr']].apply(checkChr, axis=1)
    # check whether they are from the same chromosome
    assert len(eQTMs['checkChr'].unique()) == 1
    print(eQTMs.shape, eQTMs.head())
    return eQTMs


def add_gene_meanVar(eqtm, gene_meanVar_filepath):
    gene_meanVar = pd.read_csv(gene_meanVar_filepath, sep=",", index_col=0, header=None).T.to_dict()
    eqtm["expressionMean"] = [gene_meanVar[gene][1] for gene in eqtm["ProbeName"]]
    eqtm["expressionVar"] = [gene_meanVar[gene][2] for gene in eqtm["ProbeName"]]
    print(eqtm.shape, eqtm.head())
    return eqtm


def add_methy_meanVar(eqtm, methy_meanVar_filepath):
    methy_meanVar = pd.read_csv(methy_meanVar_filepath, sep=",", index_col=0, header=None).T.to_dict()
    eqtm["methyMean"] = [methy_meanVar[gene][1] for gene in eqtm["SNPName"]]
    eqtm["methyVar"] = [methy_meanVar[gene][2] for gene in eqtm["SNPName"]]
    print(eqtm.shape, eqtm.head())
    return eqtm


if __name__ == "__main__":
    # parser = argparse.ArgumentParser()
    # parser.add_argument("--name", dest="name", type=str)
    # parser.add_argument("--overlapRatio_filename", dest="overlapRatio_filename", type=str)
    # args = parser.parse_args()
    # name = args.name
    # overlapRatio_filename = args.overlapRatio_filename
    for name in ["eqtms_metaTCGA_sig", "eqtms_metaTCGA_between0.05_0.5", "eqtms_metaTCGA_insig"]:
        print(name)
        eqtm_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/ORIGIN/%s.txt" % name
        save_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores" \
                        "/withExpressionTSSMethyCpgOverlapPromoter/%s_all_features_added.txt" % name
        cpg_overlapRatio_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores" \
                                    "/allCpgs/cpg_metaTCGA_all_overlapRatio.txt"
        gene_startEnd_position_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data" \
                                          "/features/TSSDistance/gene_startEndSite.txt"
        gene_meanVar_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/features/meanVar/metaTCGA_rna_meanVar.txt"
        methy_meanVar_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/features/meanVar/metaTCGA_methy_meanVar.txt"

        eqtm = read_eqtm(eqtm_filepath)
        eqtm = add_cpg_overlapRatio(eqtm, cpg_overlapRatio_filepath)
        eqtm = add_tss_distance(eqtm, gene_startEnd_position_filepath)
        eqtm = add_gene_meanVar(eqtm, gene_meanVar_filepath)
        eqtm = add_methy_meanVar(eqtm, methy_meanVar_filepath)
        eqtm.to_csv(save_filepath)