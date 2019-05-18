import __init__path
import os
import lime
import lime.lime_tabular
from lib.read.read_data import load_data
import pickle
import copy
import collections
import numpy as np


def submodular_fn(explanations,
                  feature_value):
    z_words = set()
    for key in explanations.keys():
        exp = explanations[key]
        z_words = z_words.union([x[0] for x in exp])
    normalizer = sum([feature_value[w] for w in z_words])

    def fnz(x):
        all_words = set()
        for doc in x:
            key = str(doc)
            all_words = all_words.union([x[0] for x in explanations[key]])
        return sum([feature_value[w] for w in all_words]) / normalizer

    fnz.num_items = len(explanations)
    return fnz


def greedy(submodular_fn, k, chosen=[]):
    chosen = copy.deepcopy(chosen)
    all_items = [num for num in range(submodular_fn.num_items)]
    current_value = 0
    while len(chosen) != k:
        best_gain = 0
        best_item = all_items[0]
        for i in all_items:
            gain = submodular_fn(chosen + [i]) - current_value
            if gain > best_gain:
                best_gain = gain
                best_item = i
        chosen.append(best_item)
        all_items.remove(best_item)
        current_value += best_gain
    return chosen


def submodular_pick(exps1, B):
    def get_function(exps):
        feature_value = collections.defaultdict(float)
        for instance_id in exps.keys():
            exp = exps[instance_id]
            print(len(exp))
            for f, v in exp:
                feature_value[f] += np.abs(v)
        print(feature_value)
        for f in feature_value:
            feature_value[f] = np.sqrt(feature_value[f])
            print(feature_value[f])
            submodular = submodular_fn(exps, feature_value)
        return submodular

    fn1 = get_function(exps1)
    print(fn1)
    return greedy(fn1, B)


if __name__ == "__main__":
    # project_rootdir = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
    project_rootdir = "/home/shuang/projects/development_eqtm"
    gt_data_filepath = os.path.join(project_rootdir,
                                    "data",
                                    "eqtmZscores",
                                    "withExpressionTSSMethyCpgOverlapGene",
                                    "2017-12-09-eQTLsFDR-gt0_withExpressionTssMethyOverlap_withGeneOverlap.txt")
    et_data_filepath = os.path.join(project_rootdir,
                                    "data",
                                    "eqtmZscores",
                                    "withExpressionTSSMethyCpgOverlapGene",
                                    "2017-12-09-eQTLsFDR-et0_withExpressionTssMethyOverlap_withGeneOverlap.txt")

    exclude = ['SNPName', 'ProbeName', 'OverallZScore',
               'SNPChr', 'PValue', 'SNPChrPos', 'ProbeChr',
               'ProbeCenterChrPos', 'CisTrans', 'SNPType', 'AlleleAssessed',
               'DatasetsWhereSNPProbePairIsAvailableAndPassesQC',
               'DatasetsZScores', 'DatasetsNrSamples',
               'IncludedDatasetsMeanProbeExpression',
               'IncludedDatasetsProbeExpressionVariance', 'HGNCName',
               'IncludedDatasetsCorrelationCoefficient', 'Meta-Beta (SE)',
               'Beta (SE)', 'FoldChange', 'FDR', 'checkChr',
               'TssSite', 'chr', 'SNPName_ProbeName']
    keep = []
    gt_data = load_data(gt_data_filepath, exclude=exclude, keep=keep)
    et_data = load_data(et_data_filepath, exclude=exclude, keep=keep, test_size=0)

    model_filepath = os.path.join(project_rootdir,
                                  "model",
                                  "train_gt.pkl")
    model = pickle.load(open(model_filepath, "rb"))

    # create the explainer
    explainer = lime.lime_tabular.LimeTabularExplainer(training_data=gt_data.train.values.values,
                                                       training_labels=gt_data.train.labels.values,
                                                       feature_names=gt_data.train.values.columns.values,
                                                       class_names=np.array(["Negative", "Positive"]))

    # explain an instance
    example_map = explainer.explain_instance(et_data.train.values.loc[3],
                                     model.predict_proba,
                                     num_features=len(gt_data.train.values.columns)).as_map()
    print(example_map)

    # get the tp with confidence score == 1
    et_pred = model.predict_proba(et_data.train.values)
    et_conf_negative = et_pred[:, 0]
    et_conf_positive = et_pred[:, 1]
    negative_num = 0
    positive_num = 0
    negative_selected = []
    positive_selected = []
    for ind in range(et_data.train.values.shape[0]):
        if et_conf_negative[ind] > 0.9 and et_data.train.labels[ind] == 0\
                and negative_num <= 5:
            negative_num += 1
            exp = explainer.explain_instance(et_data.train.values.loc[ind],
                                             model.predict_proba,
                                             num_features=len(et_data.train.values.columns)).as_map()
            negative_selected.append([ind, exp])
            print(ind, et_conf_negative[ind], et_data.train.labels[ind])
        elif et_conf_positive[ind] > 0.9 and et_data.train.labels[ind] == 1 \
                and positive_num <= 5:
            positive_num += 1
            exp = explainer.explain_instance(et_data.train.values.loc[ind],
                                             model.predict_proba,
                                             num_features=len(et_data.train.values.columns)).as_map()
            positive_selected.append([ind, exp])
            print(ind, et_conf_positive[ind], et_data.train.labels[ind])


    # # explain all instances and pickle store the maps
    # all_instances_maps = {}
    # for i in range(et_data.train.values.shape[0]):
    #     print("Processing {}th instance...".format(i))
    #     instance = et_data.train.values.loc[i]
    #     map = explainer.explain_instance(instance,
    #                                      model.predict_proba,
    #                                      num_features=len(gt_data.train.values.columns)).as_map()
    #     print(map.keys())
    #     all_instances_maps.update({"%d"%i: map[1]})
    # mapfile_save_path = os.path.join(project_rootdir, "model_explanations", "et_all_maps.pkl")
    # pickle.dump(all_instances_maps, open(mapfile_save_path, "wb"))
    #
    # test_pickle = pickle.load(open(mapfile_save_path, "rb"))
    # for key in test_pickle.keys():
    #     print(key)
    #
    # exps_filepath = os.path.join(project_rootdir,
    #                              "model_explanations",
    #                              "et_all_maps.pkl")
    # explanations = pickle.load(open(exps_filepath, "rb"))
    #
    # # picked_instances = submodular_pick(explanations, 20)
    # feature_value = {}
    # for instance_id in explanations.keys():
    #     exp = explanations[instance_id]
    #     print(len(exp))
    #     for f, v in exp:
    #         if f in feature_value:
    #             feature_value[f] += np.abs(v)
    #         else:
    #             feature_value[f] = v
    # print(feature_value)
    # for f in feature_value.keys():
    #     feature_value[f] = np.sqrt(feature_value[f])
    #     print(feature_value[f])
    #     submodular = submodular_fn(explanations, feature_value)
    # items = greedy(submodular, 10)
    # save_item_filepath = os.path.join(project_rootdir,
    #                                   "model_explanations",
    #                                   "picked_items.txt")
    # with open(save_item_filepath, "w") as f:
    #     f.write('\t'.join([str(item) for item in items]))
    #     f.close()

