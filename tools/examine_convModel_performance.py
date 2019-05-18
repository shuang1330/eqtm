import __init__path
from lib.model.finetune_model import examineModel
import os
import pandas as pd
from lib.read.read_data import load_data
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    # locally
    project_rootdir = '/home/shuang/projects/development_eqtm'
    # calculon project root dir
    # project_rootdir = '/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm'

    # parser = argparse.ArgumentParser()
    # parser.add_argument("--eqtm_filepath", dest="eqtm_filepath",
    #                     default=".")
    # args = parser.parse_args()
    # print("Processing file: ", args.eqtm_filepath)
    # # raise NotImplementedError
    #
    # eqtm_path = os.path.join(project_rootdir, args.eqtm_filepath)
    #
    # model_name = 'ranfor'
    # exclude = ['SNPName', 'ProbeName']
    # keep = ['OverallZScore']
    # train_path = eqtm_path
    # test_path = train_path
    # examineModel(model_name, train_path, test_path,
    #              keep=keep, exclude=exclude, display=True)
    eqtm_dirpath = os.path.join(project_rootdir,
                                "data",
                                "eqtmZscores",
                                "withExpressionTSSMethyCpgOverlapGenePromoter")
    eqtm_names = os.listdir(eqtm_dirpath)
    exclude = ['SNPChr', 'PValue', 'SNPChrPos', 'ProbeChr',
               'ProbeCenterChrPos', 'CisTrans', 'SNPType', 'AlleleAssessed',
               'DatasetsWhereSNPProbePairIsAvailableAndPassesQC',
               'DatasetsZScores', 'DatasetsNrSamples',
               'IncludedDatasetsMeanProbeExpression',
               'IncludedDatasetsProbeExpressionVariance', 'HGNCName',
               'IncludedDatasetsCorrelationCoefficient', 'Meta-Beta (SE)',
               'Beta (SE)', 'FoldChange', 'FDR', 'checkChr', 'SNPName_ProbeName',
               'TssSite', 'chr']
    keep = ['SNPName', 'OverallZScore', 'ProbeName']

    # examine unique effects or not
    et_filepath = os.path.join(project_rootdir, "data", "eqtmZscores",
                               "withExpressionTSSMethyCpgOverlapPromoter",
                               "2017-12-09-eQTLsFDR-et0_withExpressionTssMethyOverlap_PromoterOverlap.txt")
    gt_filepath = os.path.join(project_rootdir, "data", "eqtmZscores",
                               "withExpressionTSSMethyCpgOverlapPromoter",
                               "2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_PromoterOverlap.txt")
    data_all = pd.concat([load_data(et_filepath, exclude=exclude, keep=keep, test_size=0).train.values,
                          load_data(gt_filepath, exclude=exclude, keep=keep, test_size=0).train.values])
    print(data_all.head())
    print(data_all.columns)
    raise NotImplementedError

    def check_direction(zscore):
        if zscore > 0:
            return 1
        else:
            return -1
    data_all["direction"] = [check_direction(item) for item in data_all["OverallZScore"]]
    data_groupby_count = data_all.groupby("SNPName")["direction"].agg(["count", np.var])  # one value causes NaN
    dups = []
    nums = []
    same_directions = []
    diff_directions = []

    def check_non_zero(item):
        if item > 0:
            return 1
        else:
            return 0
    for i in range(1, data_groupby_count["count"].max()):
        dups_num = data_groupby_count["count"][data_groupby_count["count"] == i].count()
        diff_num = sum([check_non_zero(item) for item in data_groupby_count["var"][data_groupby_count["count"]==i]])
        if dups_num >= 0:
            dups.append(i)
            nums.append(dups_num)
            same_directions.append((dups_num-diff_num)/data_all.shape[0])
            diff_directions.append(diff_num/data_all.shape[0])
    print(dups)
    print(nums)
    print(same_directions)
    print(diff_directions)
    p1 = plt.bar(dups, diff_directions, alpha=0.5, color="orange", label="Different directions")
    p2 = plt.bar(dups, same_directions, bottom=diff_directions, color="cornflowerblue",  alpha=0.5, label="Same direction")
    plt.legend()
    plt.savefig(os.path.join(project_rootdir, "output_fig",
                             "duplicate_number_distribution.png"))
    plt.show()

    raise NotImplementedError
    # Assumption3, examine different cell lines
    # test
    train_filename = "2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_PromoterOverlap.txt"
    train_path = os.path.join(project_rootdir,
                              "data",
                              "eqtmZscores",
                              "withExpressionTSSMethyCpgOverlapPromoter",
                              train_filename)
    test_path = os.path.join(project_rootdir,
                             "data",
                             "eqtmZscores",
                             "withExpressionTSSMethyCpgOverlapPromoter",
                             "2017-12-09-eQTLsFDR-et0.0-flipped_withExpressionTssMethy_nonBlood_PromoterOverlap.txt")

    model_savepath = os.path.join(project_rootdir, "model",
                                  train_filename+".pkl")
    examineModel(model_name='ranfor', trainPath=train_path,
                 testPath=test_path, keep=keep, exclude=exclude,
                 display=True, savefilename=model_savepath)

    raise NotImplementedError
    for eqtm_name in eqtm_names:
        eqtm_filepath = os.path.join(eqtm_dirpath, eqtm_name)
        print("Examining file:", eqtm_name, flush=True)
        model_savepath = os.path.join(project_rootdir,
                                      "model",
                                      "{}.pkl".format(eqtm_name))
        examineModel(model_name='ranfor', trainPath=eqtm_filepath,
                     testPath=eqtm_filepath, keep=keep, exclude=exclude,
                     display=True, savefilename = model_savepath)

