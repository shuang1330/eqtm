import __init__path
import pandas as pd
import os
# import tensorflow as tf
# from tensorflow.contrib.tensorboard.plugins import projector
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix
from sklearn import metrics
from sklearn.model_selection import train_test_split
import numpy as np
from scipy.stats import chisquare
from seaborn import clustermap
from scipy import interp
from sklearn.metrics import auc
from scipy.cluster.hierarchy import dendrogram, linkage
from  scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
import functools
print = functools.partial(print, flush=True)


def save_metadata(meta_filepath, num_significant, num_insignificant):
    with open(meta_filepath, "w") as f:
        for i in range(num_significant):
            f.write("1\n")
        for i in range(num_insignificant):
            f.write("0\n")
        f.close()


def generate_embedding(datapoints,
                       metadata_path,
                       save_dirpath,
                       save_filepath):
    sess = tf.Session()
    with tf.device("/cpu:0"):
        embedding = tf.Variable(tf.stack(datapoints, axis=0),
                                trainable=False,
                                name='embedding')
    tf.global_variables_initializer().run(session=sess)
    saver = tf.train.Saver()
    writer = tf.summary.FileWriter(save_dirpath, sess.graph)

    config = projector.ProjectorConfig()
    embed = config.embeddings.add()
    embed.tensor_name = 'embedding:0'
    embed.metadata_path = metadata_path
    print(metadata_path)

    projector.visualize_embeddings(writer, config)
    saver.save(sess,
               save_filepath,
               global_step=datapoints.shape[0])


def simple_classifier(train_data, train_label,
                      test_data, test_label):
    model = RandomForestClassifier(n_estimators=100)
    model.fit(train_data, train_label)
    pred = model.predict(test_data)
    tn, fp, fn, tp = confusion_matrix(test_label, pred).ravel()
    sensitivity = tp / (fn + tp)
    specificity = tn / (fp + tn)
    prods = model.predict_proba(test_data)[:, 1]
    fpr, tpr, _ = metrics.roc_curve(test_label, prods)
    score = metrics.auc(fpr, tpr)  # auc score
    print(round(sensitivity, 2), round(specificity, 2), round(score, 2))
    return sensitivity, specificity, score


def select_most_insignificant_eqtms(eqtm_filepath,
                                    most_insignificant_eqtm_filepath):
    eqtm = pd.read_csv(eqtm_filepath, sep="\t", index_col=0)
    if os.path.isfile(most_insignificant_eqtm_filepath[:-4]+"_names.txt"):
        with open(most_insignificant_eqtm_filepath[:-4]+"_names.txt", "r") as f:
            most_insig_names = [item.strip() for item in f.readlines()]
            f.close()
    else:
        with open(most_insignificant_eqtm_filepath, "r") as f:
            most_insig_names = [line.strip().split("\t")[3] for line in f.readlines()]
            print(most_insig_names[0])
            # raise NotImplementedError
            f.close()
        print(eqtm.loc[most_insig_names[0]])
        with open(most_insignificant_eqtm_filepath[:-4]+"_names.txt", "w") as f:
            f.write("\n".join([str(item)for item in most_insig_names]))
            f.close()
        print("Wrote names to file, ", most_insignificant_eqtm_filepath[:-4]+"_names.txt")
    return eqtm.loc[most_insig_names]


def measure_chisquare(dataset1, dataset2):
    measurement = {}
    for feature in dataset1.columns:
        hist1, _ = np.histogram(dataset1[feature], bins=100)
        hist2, _ = np.histogram(dataset2[feature], bins=100)
        measurement[feature] = chisquare(hist1, hist2)
    return measurement


def draw_feature_dff(dataset1, label1,
                     dataset2, label2,
                     savefig_path):
    for feature in dataset1.columns:
        plt.hist(dataset1[feature], bins=50, alpha=0.2,
                 label=label1, normed=True, stacked=True)
        plt.hist(dataset2[feature], bins=50, alpha=0.2,
                 label=label2, normed=True, stacked=True)
        plt.savefig(os.path.join(savefig_path, feature+".png"))
        plt.legend(loc="upper right")
        plt.title(feature)
        plt.show()
    return None


def measure_auc(classifier, test_data_values, test_data_labels):
    prods = classifier.predict_proba(test_data_values)[:, 1]
    fpr, tpr, _ = metrics.roc_curve(test_data_labels, prods)
    score = metrics.auc(fpr, tpr)
    return score


if __name__ == "__main__":
    # project_rootdir = "/home/shuang/projects/development_eqtm"
    # cluster
    project_rootdir = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm"
    data_dirpath = os.path.join(project_rootdir, "data",
                                "eqtmZscores", "allCpgs")
    data_significant = pd.read_csv(os.path.join(data_dirpath,
                                                "eqtm_FDR_smaller_than_0.05_bedtoolsFormat_overlapRatio.txt"),
                                   sep="\t", index_col=0)
    data_significant["label"] = 1
    features = [col for col in data_significant.columns if col not in ["label"]]
    data_insignificant_filepath = os.path.join(data_dirpath,
                                               "eqtm_FDR_larger_than_0.05_bedtoolsFormat_overlapRatio.txt")
    most_insignificant_filepath = os.path.join(data_dirpath,
                                               "eqtm_FDR_larger_than_0.5_bedtoolsFormat.txt")
    data_insignificant = select_most_insignificant_eqtms(data_insignificant_filepath,
                                                         most_insignificant_filepath).sample(
        n=data_significant.shape[0], random_state=0
    )
    data_insignificant["label"] = 0
    num_significant = data_significant.shape[0]
    num_insignificant = data_insignificant.shape[0]

    data_all = pd.concat([data_significant, data_insignificant], axis=0)
    print(data_all.columns)
    train_data, test_data_originally = train_test_split(data_all, test_size=0.25)
    test_data_new = pd.read_csv("/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/allCpgs/cpg_metaTCGA_overlapRatio.txt",
                                sep="\t")
    test_data_new["label"] = 1
    test_cpgs = data_all.index
    print(test_cpgs)
    isSeen = lambda x:True if x in test_cpgs else False
    test_data_new["isSeen"] = [isSeen(item) for item in test_data_new["Unnamed: 0"]]
    print("How many cpgs were not seen", test_data_new[test_data_new["isSeen"] == False].shape)
    # test_data = pd.concat([test_data_originally[features], test_data_new[features]])
    # test_label = pd.concat([test_data_originally["label"], test_data_new["label"]])
    # print(test_data.head())
    # print(test_label.head())
    print("Dataset with positive samples {}, negative samples {}".
          format(data_significant.shape[0], data_insignificant.shape[0]))
    print("Train data shape", train_data.shape)
    # print("Test data shape", test_data.shape)

    # predict whether significant or not, drawing AUROC
    print("Start to collect model results...", flush=True)
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    # for i in range(220):
    for i in range(1):
        auc_score = 0
        data_all = pd.concat([data_significant, data_insignificant.sample(n=data_significant.shape[0],
                                                                          random_state=i)],
                             axis=0)
        # train_data, test_data = train_test_split(data_all)
        model = RandomForestClassifier(n_estimators=100)
        model.fit(train_data[features], train_data["label"])
        probas_ = model.predict_proba(test_data_new[features])
        pred = model.predict(test_data_new[features][test_data_new["isSeen"] == False])
        print(pred)
        print(model.score(test_data_originally[features], test_data_originally["label"]))
        print(model.score(test_data_new[features][test_data_new["isSeen"] == False],
                          test_data_new["label"][test_data_new["isSeen"] == False]))
        print(len(probas_[:, 1][probas_[:, 1] > 0.4]))
        print(len(pred))
        # Compute ROC curve and area the curve
        # fpr, tpr, thresholds = metrics.roc_curve(test_label, probas_[:, 1])
        # tprs.append(interp(mean_fpr, fpr, tpr))
        # tprs[-1][0] = 0.0
        # roc_auc = auc(fpr, tpr)
        # print(i, roc_auc, flush=True)
        # aucs.append(roc_auc)
        # plt.plot(fpr, tpr, lw=1, alpha=0.3,
        #          label='ROC fold %d (AUC = %0.2f)' % (i, roc_auc))

    et_gt_filepath = "/groups/umcg-gcc/tmp03/umcg-sli/development_eqtm/data/eqtmZscores/withExpressionTSSMethyCpgOverlapPromoter/et_gt_all_celltypes.txt"
    et_gt = pd.read_csv(et_gt_filepath, sep=",")
    getDirection = lambda x:0 if x < 0 else 1
    et_gt["label"] = [getDirection(item) for item in et_gt["OverallZScore"]]

    test_subset = test_data_new[test_data_new["isSeen"] == False].copy()
    # test_subset["eqtm_combi"] = ["{}_{}".format(item[0], item[1]) for
    #                              item in zip(test_subset["SNPName"], test_subset["ProbeName"])]
    correct_cpgs = test_subset["Unnamed: 0"][test_subset["label"] == pred].values
    print("correct predicted significance: ", len(correct_cpgs))
    print(correct_cpgs)
    print(test_subset["Unnamed: 0"][test_subset["label"] != pred].values)

    meta_tcga_eqtms = pd.read_csv("/groups/umcg-gcc/tmp03/umcg-sli/temp_data/added_overlapRatio_metaTCGA.txt")

    isSig = lambda x:True if x in correct_cpgs else False
    meta_tcga_eqtms["isSig"] = [isSig(item) for item in meta_tcga_eqtms["SNPName"]]
    flipZscore = lambda zscore, allele: zscore if allele=="C" else -zscore
    meta_tcga_eqtms["flippedZscore"] = [flipZscore(item[0], item[1]) for item in zip(meta_tcga_eqtms["OverallZScore"],
                                                                                       meta_tcga_eqtms["AlleleAssessed"])]
    meta_tcga_eqtms["direction"] =[getDirection(zscore) for zscore in meta_tcga_eqtms["flippedZscore"]]

    model2 = RandomForestClassifier(n_estimators=100)
    model2.fit(et_gt[features], et_gt["label"])
    auc_sig = measure_auc(model2, meta_tcga_eqtms[features][meta_tcga_eqtms["isSig"]==True], meta_tcga_eqtms["direction"][meta_tcga_eqtms["isSig"]==True])

    print(auc_sig)
    raise NotImplementedError
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
             label='Luck', alpha=.8)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    plt.plot(mean_fpr, mean_tpr, color='b',
             label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
             lw=2, alpha=.8)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                     label=r'$\pm$ 1 std. dev.')

    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
    plt.savefig(os.path.join(project_rootdir, "output_fig", "sig_vs_insig_auc.png"))
    print("Fin.", flush=True)

    # model error analysis
    preds = model.predict(test_data[features])
    errors = test_data.index[test_data["label"] != preds]
    tn, fp, fn, tp = confusion_matrix(test_data["label"], preds).ravel()
    print(test_data.loc[errors])

    new_test_set = test_data.copy()
    new_test_set["isAllZero"] = new_test_set[features].apply(np.sum, axis=1)
    # all_zero_subset = new_test_set[new_test_set["isAllZero"]==0]

    # find all zeros
    data_all["isAllZero"] = data_all[features].apply(np.sum, axis=1)
    all_zero_subset = data_all[data_all["isAllZero"] == 0]
    print(all_zero_subset.shape)
    print(all_zero_subset[all_zero_subset["label"] > 0].shape)


    raise NotImplementedError

    # calculate zero ratios
    zero_ratios_sig = []
    zero_ratios_insig = []
    num_sig = data_significant.shape[0]
    num_insig = data_insignificant.shape[0]
    for feature in features:
        zero_ratios_sig.append(data_significant[feature][data_significant[feature]==0].count()/num_sig)
        zero_ratios_insig.append(data_insignificant[feature][data_insignificant[feature]==0].count()/num_insig)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(range(len(zero_ratios_sig)), zero_ratios_sig, alpha=0.2, label="Significant eQTMs")
    ax.bar(range(len(zero_ratios_insig)), zero_ratios_insig, alpha=0.2, label="Insignificant eQTMs")
    plt.xticks(range(len(zero_ratios_sig)), features, rotation='vertical')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(project_rootdir, "output_fig", "zeroRatios_diff_sigAndInsig.png"))
    plt.show()

    raise NotImplementedError

    # draw clustermaps
    # examine the zeros in the insignificant dataset
    for feature in features:
        null = data_insignificant.shape[0] - np.count_nonzero(data_insignificant[feature].values)
        if null > 0:
            print(feature, null)
    data_insignificant_nozero = data_insignificant
    data_insignificant_nozero[features] = data_insignificant_nozero[features] + 1e-6
    null = data_insignificant_nozero.size - np.count_nonzero(data_insignificant_nozero.values)

    # seaborn clustermap
    # all_data = pd.concat([data_significant, data_insignificant_nonzero.sample(data_significant.shape[0])])
    all_data = pd.concat([data_significant, data_insignificant_nozero.sample(data_significant.shape[0])])
    color_dic = {1: "c", 0: "azure"}
    all_data["is_significant"] = [color_dic[item] for item in all_data["label"].values]

    all_data_corr = all_data[features].T.corr()
    clustermap(all_data_corr, row_colors=all_data["is_significant"])
    plt.savefig(os.path.join(project_rootdir, "output_fig", "sig_insig_clustermap.png"))
    print("Saved clustermap to path:",
          os.path.join(project_rootdir, "output_fig", "sig_insig_clustermap.png"))
    raise NotImplementedError
    # draw the clustermap for correlation matrix


    lut = dict(zip(all_data["label"].unique(), "rbg"))
    all_data["Row_color"] = [lut[item] for item in all_data["label"]]
    g = clustermap(all_data[features],
                   row_colors=all_data["Row_color"])
    plt.savefig(os.path.join(project_rootdir,
                             "output_fig",
                             "clustermap_sig.png"))
    raise NotImplementedError
    # hierarchical clustering
    Z = linkage(all_data[features].values, 'ward')
    c, coph_dists = cophenet(Z, pdist(all_data[features].values))
    print("Cophenet coefficient: ", c)
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.set_title('Hierarchical Clustering Dendrogram (truncated)')
    ax.set_xlabel('sample index')
    ax.set_ylabel('distance')
    dendrogram(
        Z,
        truncate_mode='lastp',  # show only the last p merged clusters
        p=12,  # show only the last p merged clusters
        show_leaf_counts=False,  # otherwise numbers in brackets are counts
        leaf_rotation=90,
        leaf_font_size=12,
        show_contracted=True,  # to get a distribution impression in truncated branches
        labels=all_data.index.values
    )
    my_palette = plt.cm.get_cmap("Accent", 2)
    xlbls = ax.get_ymajorticklabels()
    num = -1
    # for lbl in xlbls:
    #     num += 1
    #     val = all_data["label"].values[num]
    #     print(val)
    #     lbl.set_color(my_palette(val))

    plt.savefig(os.path.join(project_rootdir, "output_fig", "dendrogram_sig.png"))

    raise NotImplementedError
    # chisquare differences in histograms
    chisquare_scores = measure_chisquare(data_significant,
                                         data_insignificant)
    raise NotImplementedError
    # draw feature differences
    savefig_path = os.path.join(project_rootdir,
                                "output_fig",
                                "sig_vs_insig")
    draw_feature_dff(data_significant, "Significant",
                     data_insignificant, "Insignificant",
                     savefig_path)
    raise NotImplementedError
    # tensorboard
    writer_dirpath = os.path.join(project_rootdir,
                                  "data",
                                  "tensorboard_writer")
    metadata_filepath = os.path.join(writer_dirpath,
                                     "all_eqtm_ratio_label_significance.tsv")
    # save_metadata(metadata_filepath, num_significant, num_insignificant)
    save_filepath = os.path.join(writer_dirpath, "all_eqtm_ratio_label_significance.ckpt")
    # generate_embedding(data_all.values,
    #                    "all_eqtm_ratio_label_significance.tsv",
    #                    writer_dirpath,
    #                    save_filepath)
    raise NotImplementedError