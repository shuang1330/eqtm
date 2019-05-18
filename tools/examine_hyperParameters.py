import os
from lib.read.read_data import load_data
from lib.read.read_allCpgFile import cpgFile
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt

if __name__ == "__main__":

    project_rootdir = "/home/shuang/projects/development_eqtm"

    eqtm_filepath = os.path.join(
        project_rootdir,
        "data", "eqtmZscores", "withExpressionTSSMethyCpgOverlapGene",
        "2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_withGeneOverlap.txt"
    )
    exclude, keep = cpgFile().cpgFile_DefaultExcludeAndKeep()
    data = load_data(eqtm_filepath, keep=keep, exclude=exclude, test_size=0)
    print(data.train.values.columns)
    print(data.train.labels.head())
    features = [col for col in data.train.values.columns if col not in ["OverallZScore"]]
    data_train = data.train.values[features]
    data_test = data.test.values[features]

    assert "OverallZScore" not in data_train

    all_param = {"n_estimators": np.arange(2, 302, 20),
                 "max_depth": np.arange(1, 70, 3),
                 "min_samples_split": np.arange(0.1, 1.1, 0.1),
                 "min_samples_leaf": np.arange(1, 31, 2),
                 "max_leaf_nodes": np.arange(2, 32, 2),
                 "min_weight_fraction_leaf": np.arange(0.1, 0.5, 0.1),
                 "class_weight": ["balanced", {1: 2, 0: 1},
                                  {1: 4, 0: 1}, {1: 1, 0: 2},
                                  {1: 1, 0: 4}]}
    # all_param = {"n_estimators": np.arange(200, 301, 50),
    #               "max_leaf_nodes": np.arange(60, 121, 20)}
    num_cores = multiprocessing.cpu_count()

    for param in all_param.keys():
        model = GridSearchCV(estimator=RandomForestClassifier(),
                             param_grid={param: all_param[param]},
                             n_jobs=num_cores)
        model.fit(data_train, data.train.labels)
        f = plt.figure()
        ax = f.add_subplot(111)
        ax.scatter(model.cv_results_["param_{}".format(param)].data,
                   model.cv_results_["mean_test_score"])
        ax.set_title(param)
        plt.savefig("./development_eqtm/output_fig/param_effect/et/{}.png".
                    format(param))
        plt.show()
