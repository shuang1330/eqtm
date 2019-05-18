import os
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
    project_rootdir = "/home/shuang/projects/development_eqtm"
    exps_filepath = os.path.join(project_rootdir,
                                 "model_explanations",
                                 "et_all_maps.pkl")
    explanations = pickle.load(open(exps_filepath, "rb"))
    print(len(explanations.keys()))

    # picked_instances = submodular_pick(explanations, 20)
    feature_value = {}
    for instance_id in explanations.keys():
        exp = explanations[instance_id]
      #  print(len(exp))
        for f, v in exp:
            if f in feature_value:
                feature_value[f] += np.abs(v)
            else:
                feature_value[f] = v
   # print(feature_value)
    for f in feature_value.keys():
        feature_value[f] = np.sqrt(feature_value[f])
    #    print(feature_value[f])
        submodular = submodular_fn(explanations, feature_value)
    items = greedy(submodular, 50)
    save_item_filepath = os.path.join(project_rootdir,
                                      "model_explanations",
                                      "picked_items.txt")
    with open(save_item_filepath, "w") as f:
        f.write('\t'.join([str(item) for item in items]))
        f.close()
