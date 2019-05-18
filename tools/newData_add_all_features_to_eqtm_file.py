import pandas as pd


def read_eqtm(eqtm_filepath):
    def return_line_with_correct_type(line):
        stuff = line.strip().split("\t")
        return [stuff[0], stuff[1], int(stuff[2]), int(stuff[3]),
                int(stuff[4]), int(stuff[5]), int(stuff[6]), float(stuff[7]),
                float(stuff[8]), float(stuff[9])]
    with open(eqtm_filepath, "r") as f:
        info = f.readlines()
        columns = info[0].strip().split("\t")
        data = [return_line_with_correct_type(line) for line in info[1:]]
        f.close()
    eqtm = pd.DataFrame(data=data, columns=columns)
    print(eqtm.head())
    # eqtm = pd.read_csv(eqtm_filepath, sep="\t", index_col=0)
    # print(eqtm.shape, eqtm.columns)
    return eqtm


def add_promoter_overlapRatio(eqtm, promoter_overlapRatio_filepath):
    overlap = pd.read_csv(promoter_overlapRatio_filepath, sep="\t", index_col=0)
    overlap_dict = overlap.T.to_dict()
    print(overlap_dict['ENSG00000006025.6_45885904_45886026'])

    for col in overlap.columns:
        def find_overlap_vector(item):
            if item in overlap_dict:
                # print(col)
                return overlap_dict[item][col]
            else:
                return 0
        col_name = col+'_promoter'
        eqtm[col_name] = [find_overlap_vector(row) for row in eqtm['GENE_ID'].values]
    print(eqtm.head())
    return eqtm


def add_cpg_overlapRatio(eqtm, cpg_overlapRatio_filepath):
    overlap = pd.read_csv(cpg_overlapRatio_filepath, sep="\t", index_col=0)
    overlap_dict = overlap.T.to_dict()
    print(overlap_dict['cg00777063'])

    for col in overlap.columns:
        def find_overlap_vector(item):
            if item in overlap_dict:
                return overlap_dict[item][col]
            else:
                return 0
        eqtm[col] = [find_overlap_vector(item) for item in eqtm['METHYL_label'].values]
    print(eqtm.head())
    return eqtm


if __name__ == "__main__":
    eqtm_filepath = \
        "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/new_blood_eqtm/eqtm_fibroblast.txt"
    save_filepath = \
        "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/withExpressionTSSMethyCpgOverlapPromoter/eqtm_fibroblast_allFeaturesAdded.txt"
    promoter_overlapRatio_filepath = \
        "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/allCpgs/cpg_fibroblast_gene_overlapRatio.txt"
    cpg_overlapRatio_filepath = \
        "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/allCpgs/cpg_fibroblast_overlapRatio.txt"
    eqtm = read_eqtm(eqtm_filepath)
    eqtm = add_promoter_overlapRatio(eqtm, promoter_overlapRatio_filepath)
    eqtm = add_cpg_overlapRatio(eqtm, cpg_overlapRatio_filepath)
    eqtm.to_csv(save_filepath)
