import __init__path
import xgboost as xgb
import os
from lib.read.read_allCpgFile import cpgFile
from lib.read.read_data import load_data
from sklearn.metrics import roc_auc_score


def xgb_model_fit(xgb, train_data,
                  predictors, useTrainCV=True,
                  cv_folds=5, early_stopping_round=50):
    if useTrainCV:
        xgb_param = xgb.get_xgb_param()


if __name__ == "__main__":
    project_rootdir = "/groups/umcg-gcc/tmp04/umcg-sli/development_eqtm"
    input_file = os.path.join(project_rootdir,
                              "data",
                              "eqtmZscores",
                              "withExpressionTSSMethyCpgOverlapGene",
                              "2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_withGeneOverlap.txt")
    exclude, keep = cpgFile().cpgFile_DefaultExcludeAndKeep()
    for iter in range(1):
        input_data = load_data(input_file, exclude=exclude, keep=keep, normalization=False)
        features = [col for col in input_data.train.values.columns if col not in ["OverallZScore"]]
        # model
        for max_depth in [18]:
            for learning_rate in [10e-2]:

                gbm = xgb.XGBClassifier(max_depth=max_depth,
                                        n_jobs=4,
                                        learning_rate=learning_rate).\
                    fit(input_data.train.values[features], input_data.train.labels)
                predictions = gbm.predict(input_data.test.values[features])
                print(max_depth, learning_rate, roc_auc_score(input_data.test.labels, predictions))
